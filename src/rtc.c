#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>

#include <tc.h>

#define MAX_ERRSTR_LEN 512

#define MIN(a, b) ((a) <= (b) ? (a) : (b))
#define MAX(a, b) ((a) >= (b) ? (a) : (b))

struct cb_data {
    size_t n;
    const void **ds;
    size_t N;
    size_t K;
    SEXP result;
};

static tc_clustering_cb cb;


/* Adopted from "Writing R Extensions". */
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
       elmt = VECTOR_ELT(list, i);
       break;
    }
    return elmt;
}

SEXP
rtc_tree_clustering(
    SEXP Rds,
    SEXP RN,
    SEXP RK,
    SEXP limits,
    SEXP Rfragment_size,
    SEXP Ropts
) {
    SEXP rdim = getAttrib(Rds, R_DimSymbol);
    SEXP Rvec;
    SEXP Rrange;
    SEXP result;
    SEXP result_new;
    SEXP elt;
    size_t i = 0;
    size_t N = INTEGER(RN)[0];
    size_t K = INTEGER(RK)[0];
    size_t nprotected = 0;
    int res = 0;
    const void **ds = NULL;
    struct tc_param_def *param_def = NULL;
    struct cb_data cb_data;
    struct tc_opts opts = tc_default_opts;
    size_t burnin = 0;
    size_t nsamples = 0;
    char errstr[MAX_ERRSTR_LEN];

    errstr[0] = '\0';

    nsamples = INTEGER(getListElement(Ropts, "nsamples"))[0];
    burnin = INTEGER(getListElement(Ropts, "burnin"))[0];

    opts.nsamples = nsamples + burnin;
    opts.max_segments = INTEGER(getListElement(Ropts, "max.segments"))[0];
    opts.maxiter = INTEGER(getListElement(Ropts, "maxiter"))[0];

    result = PROTECT(allocVector(VECSXP, opts.nsamples));
    nprotected++;

    ds = calloc(K, sizeof(void *));
    if (ds == NULL) {
        snprintf(errstr, sizeof(errstr), "Allocation failed");
        goto cleanup;
    }

    for (size_t k = 0; k < K; k++) {
        Rvec = PROTECT(coerceVector(VECTOR_ELT(Rds, k), REALSXP));
        nprotected++;
        ds[k] = REAL(Rvec);
    }

    param_def = calloc(K, sizeof(struct tc_param_def));
    if (param_def == NULL) {
        snprintf(errstr, sizeof(errstr), "Allocation failed");
        goto cleanup;
    }

    for (size_t k = 0; k < K; k++) {
        param_def[k].type = TC_METRIC;
        param_def[k].size = TC_FLOAT64;
        if (!isNull(Rfragment_size)) {
            param_def[k].fragment_size =
                REAL(coerceVector(Rfragment_size, REALSXP))[k];
        }
        tc_param_def_init(&param_def[k], ds[k], N);
        if (!isNull(limits) && length(limits) >= k) {
            Rrange = PROTECT(coerceVector(VECTOR_ELT(limits, k), REALSXP));
            nprotected++;
            param_def[k].min.float64 = REAL(Rrange)[0];
            param_def[k].max.float64 = REAL(Rrange)[1];
            UNPROTECT(1);
            nprotected--;
        }
    }

    cb_data.ds = ds;
    cb_data.n = 0;
    cb_data.N = N;
    cb_data.K = K;
    cb_data.result = result;
    res = tc_clustering(ds, N, param_def, K, cb, &cb_data, &opts);
    if (res != 0) {
        snprintf(errstr, sizeof(errstr), "Clustering failed: %s", strerror(errno));
        goto cleanup;
    }

    nsamples = MIN(cb_data.n, nsamples);
    result_new = PROTECT(allocVector(VECSXP, nsamples));
    nprotected++;
    for (i = 0; i < nsamples; i++) {
        elt = VECTOR_ELT(result, i + (cb_data.n - nsamples));
        SET_VECTOR_ELT(result_new, i, elt);
        SET_VECTOR_ELT(result, i, R_NilValue);
    }
    result = result_new;

cleanup:
    UNPROTECT(nprotected);
    nprotected = 0;
    free(ds);
    ds = NULL;
    free(param_def);
    param_def = NULL;
    if (*errstr) error("%s", errstr);
    return result;
}

static bool
cb(const struct tc_tree *tree, double l, const void **ds, size_t N, void *data_)
{
    SEXP Rsegments;
    SEXP Rsegment;
    SEXP Rranges;
    SEXP Rrange;
    SEXP names;
    SEXP Rsegmentation;
    size_t S = 0;
    struct tc_segment *segments = NULL;
    struct tc_segment *segment = NULL;
    struct cb_data *data = NULL;
    size_t s = 0;
    size_t k = 0;
    size_t nprotected = 0;
    char errstr[MAX_ERRSTR_LEN];

    errstr[0] = '\0';
    data = data_;

    Rsegmentation = PROTECT(allocVector(VECSXP, 2));
    nprotected++;
    names = PROTECT(allocVector(STRSXP, 2));
    nprotected++;
    SET_STRING_ELT(names, 0, mkChar("segments"));
    SET_STRING_ELT(names, 1, mkChar("likelihood"));
    setAttrib(Rsegmentation, R_NamesSymbol, names);
    UNPROTECT(1);
    nprotected--;

    segments = tc_segments(tree, data->ds, data->N, &S);
    if (segments == NULL) {
        snprintf(errstr, sizeof(errstr), "Clustering failed: %s",
            strerror(errno));
        goto cleanup;
    }

    Rsegments = PROTECT(allocVector(VECSXP, S));
    nprotected++;
    for (s = 0; s < S; s++) {
        segment = &segments[s];

        Rsegment = PROTECT(allocVector(VECSXP, 3));
        nprotected++;
        names = PROTECT(allocVector(STRSXP, 3));
        nprotected++;
        SET_STRING_ELT(names, 0, mkChar("NX"));
        SET_STRING_ELT(names, 1, mkChar("V"));
        SET_STRING_ELT(names, 2, mkChar("ranges"));
        setAttrib(Rsegment, R_NamesSymbol, names);
        UNPROTECT(1);
        nprotected--;

        Rranges = PROTECT(allocVector(VECSXP, data->K));
        nprotected++;
        for (k = 0; k < data->K; k++) {
            Rrange = PROTECT(allocVector(REALSXP, 2));
            nprotected++;
            REAL(Rrange)[0] = segment->ranges[k].min;
            REAL(Rrange)[1] = segment->ranges[k].max;
            SET_VECTOR_ELT(Rranges, k, Rrange);
            UNPROTECT(1);
            nprotected--;
        }

        SET_VECTOR_ELT(Rsegment, 0, ScalarInteger(segment->NX));
        SET_VECTOR_ELT(Rsegment, 1, ScalarReal(segment->V));
        SET_VECTOR_ELT(Rsegment, 2, Rranges);
        SET_VECTOR_ELT(Rsegments, s, Rsegment);
        UNPROTECT(2);
        nprotected -= 2;
    }
    SET_VECTOR_ELT(Rsegmentation, 0, Rsegments);
    SET_VECTOR_ELT(Rsegmentation, 1, ScalarReal(l));
    SET_VECTOR_ELT(data->result, data->n++, Rsegmentation);

cleanup:
    if (segments != NULL) {
        tc_free_segments(segments, S);
        free(segments);
        segments = NULL;
    }
    UNPROTECT(nprotected);
    nprotected = 0;
    if (*errstr) error("%s", errstr);
    return true;
}

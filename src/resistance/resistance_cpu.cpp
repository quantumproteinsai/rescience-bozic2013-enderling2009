// Resistance model (resistance_core.h): replicates in parallel with OpenMP, C interface
// for Python (ctypes, python/resistance_model.py).
// Build: make build (g++ -O2 -fopenmp -ffp-contract=off -shared).
// Replicate r uses its own random stream (seed, r), so the results do not depend on
// the number of threads.

#include "resistance_core.h"
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {

int resistance_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// Simulates replicates first_rep .. first_rep + n_reps - 1. Outputs are arrays of n_reps
// (at_detect: n_reps x 4, types 00 10 01 11). Returns 0, or -1 for invalid parameters.
int resistance_run(const resistance::Params *p, double pre_growth, double sens_growth,
                   uint64_t seed, int64_t first_rep, int64_t n_reps,
                   int32_t *outcome, int32_t *restarts, double *t_detect, double *t_end,
                   double *at_detect, double *t_switch) {
    if (!(p->M >= 1.0) || !(p->dt > 0.0) || !(p->b > p->d) || p->k_exact < 1.0 || n_reps < 0 ||
        !(p->cost >= 0.0 && p->cost < 1.0))
        return -1;
    resistance::Rates r = {pre_growth, sens_growth, {0.0, 0.0, 0.0, 0.0}};
    for (int k = 0; k < resistance::NTYPES; k++) {
        double birth, death, growth;
        resistance::own_rates(*p, r, k, birth, death, growth);
        r.type_growth[k] = exp((birth - death) * p->dt);
    }
#pragma omp parallel for schedule(dynamic, 16)
    for (int64_t i = 0; i < n_reps; i++) {
        resistance::Result res;
        resistance::simulate_lesion(*p, r, seed, (uint64_t)(first_rep + i), res);
        outcome[i] = res.outcome;
        restarts[i] = res.restarts;
        t_detect[i] = res.t_detect;
        t_end[i] = res.t_end;
        t_switch[i] = res.t_switch;
        for (int k = 0; k < resistance::NOUT; k++) at_detect[4 * i + k] = res.at_detect[k];
    }
    return 0;
}

}  // extern "C"

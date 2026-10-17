// Evolution of resistance to mono or dual targeted therapy (Bozic et al., eLife 2013;
// 2:e00747, doi:10.7554/eLife.00747): one stochastic replicate = one lesion.
//
// Arithmetic is restricted to IEEE + - * / and the functions below (no libm, no fused
// multiply-add), so results are bit-reproducible across compilers and platforms (the
// same source also builds as a hardware kernel in a separate project). Compile with
// -ffp-contract=off and without -ffast-math.
//
// Model (paper, Materials and methods): continuous-time multitype birth-death process.
// Types 00, 10, 01, 11 (1 = resistant to drug 1 / drug 2). Pre-treatment: every cell
// divides at rate b and dies at rate d, from one sensitive cell until the lesion has M
// cells. Offspring of 00 is 10 with probability u n1, 01 with u n2, 11 with u n12;
// offspring of 10 is 11 with u (n2 + n12); of 01, 11 with u (n1 + n12). Treatment: cells
// sensitive to a drug in use divide at b' and die at d'; the others keep b and d.
//
// Simulation (as in the paper, hybrid): a type with fewer than k_exact cells evolves
// stochastically, one time step dt at a time (tau-leaping: Poisson numbers of births and
// deaths, Poisson mutants among the births); a larger type follows its deterministic
// exponential growth, and its mutants are Poisson draws with the expected number per step.
// The paper simulates small populations event by event; tau-leaping with b dt ~ 0.01 is
// a standard approximation, checked against the paper's analytic results (tests).
//
// A lesion whose founding cell lineage dies out before reaching M is restarted (the paper
// considers lesions that reach detection size). Outcome after treatment starts:
// eradicated (no cell left) or failed (a resistant type reaches fail_size cells, after
// which extinction is negligible), or undecided after t_max_treat days.

#ifndef RESISTANCE_CORE_H
#define RESISTANCE_CORE_H

#include <stdint.h>
#include <math.h>
#include <vector>

namespace resistance {

// T11X: fully resistant through ONE cross-resistance mutation (n12). Used only with a
// cost of resistance, which is paid per mutation (Table 2); without cost these cells
// are T11 and T11X stays empty (no random numbers drawn: results unchanged).
enum { T00 = 0, T10 = 1, T01 = 2, T11 = 3, T11X = 4, NTYPES = 5, NOUT = 4 };
enum { OUTCOME_ERADICATED = 0, OUTCOME_FAILED = 1, OUTCOME_UNDECIDED = 2 };

struct Params {
    double M;            // cells at the start of treatment (detection size)
    double u;            // point mutation rate per division
    double n1, n2, n12;  // resistance mutations: drug 1 only, drug 2 only, both
    double b, d;         // pre-treatment (and resistant) birth / death rates per day
    double bp, dp;       // treatment rates of cells sensitive to a drug in use
    double dt;           // time step (days)
    double k_exact;      // types below this size are stochastic
    double fail_size;    // a resistant type this large means treatment failure
    double t_max_treat;  // days of treatment before "undecided"
    int dual;            // 1: drugs 1 and 2; 0: drug 1 alone
    int exact_events;    // 1: small types event by event (Gillespie), as in the paper; 0: tau-leaping
    double cost;         // cost of resistance: a type with k resistances has birth rate
                         // (b - d)(1 - cost)^k + d (paper, Table 2); its growth factor per
                         // step is passed in Rates.type_growth
    int sequential;      // 1 (with dual): drug 1 until relapse (lesion back at M cells after
                         // shrinking below M), then drug 2 (paper, Figure 5); 0: simultaneous
};

struct Result {
    int outcome;              // OUTCOME_*
    int restarts;             // founding lineages that died out before reaching M
    double t_detect;          // days from the first cell to detection
    double t_end;             // days of treatment until eradication / failure
    double at_detect[NOUT];   // cells of types 00 10 01 11 when treatment starts (11 incl. T11X)
    double t_switch;          // sequential: days on drug 1 before the switch (-1: none)
};

// ----------------------------------------------------------------------------- RNG
// xoshiro256** (Blackman & Vigna), seeded with splitmix64 from (seed, replicate):
// replicate r always gets the same stream, whatever the thread or lane running it.
struct Rng {
    uint64_t s[4];
};

static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t splitmix64(uint64_t &x) {
    x += 0x9E3779B97F4A7C15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline void rng_seed(Rng &g, uint64_t seed, uint64_t replicate) {
    uint64_t x = seed * 0xD1B54A32D192ED03ULL + replicate;
    for (int i = 0; i < 4; i++) g.s[i] = splitmix64(x);
}

static inline uint64_t rng_next(Rng &g) {
    uint64_t result = rotl64(g.s[1] * 5, 7) * 9;
    uint64_t t = g.s[1] << 17;
    g.s[2] ^= g.s[0];
    g.s[3] ^= g.s[1];
    g.s[1] ^= g.s[2];
    g.s[0] ^= g.s[3];
    g.s[2] ^= t;
    g.s[3] = rotl64(g.s[3], 45);
    return result;
}

// Uniform in [0, 1) with 53 random bits (exact in double)
static inline double rng_uniform(Rng &g) {
    return (double)(rng_next(g) >> 11) * (1.0 / 9007199254740992.0);
}

// ----------------------------------------------------------------------------- exp(-x)
// exp(-x) for x >= 0 with + - * only: x = k ln2 + r, exp(-x) = 2^-k exp(-r), exp(-r) by
// a degree-13 Taylor series (r < ln2: truncation < 1e-16 relative). Below 2^-1000 -> 0.
static inline double exp_neg(double x) {
    const double LN2 = 0.69314718055994530942;
    if (x <= 0.0) return 1.0;
    if (x > 690.0) return 0.0;
    int k = 0;
    double r = x;
    while (r >= LN2) {            // at most 996 iterations; k ln2 subtracted exactly enough
        r -= LN2;
        k++;
    }
    double term = 1.0, sum = 1.0;
    for (int i = 1; i <= 13; i++) {
        term = term * (-r) / (double)i;
        sum += term;
    }
    double scale = 1.0;
    for (int i = 0; i < k; i++) scale *= 0.5;   // exact powers of two
    return sum * scale;
}

// ----------------------------------------------------------------------------- Poisson
// Poisson(mean) by inversion (exact up to rounding; at most 1000 terms). Means above 30
// use a normal approximation by the sum of 12 uniforms (rounded, >= 0): only reached by
// deterministic types, where the exact distribution of the mutant count is irrelevant.
static inline double poisson(Rng &g, double mean) {
    if (mean <= 0.0) return 0.0;
    if (mean > 30.0) {
        double s = 0.0;
        for (int i = 0; i < 12; i++) s += rng_uniform(g);
        double sd = 1.0;
        // sqrt(mean) by Newton iterations (no libm)
        for (int i = 0; i < 30; i++) sd = 0.5 * (sd + mean / sd);
        double v = mean + sd * (s - 6.0) + 0.5;
        if (v < 0.0) return 0.0;
        return (double)(int64_t)v;
    }
    double U = rng_uniform(g);
    double p = exp_neg(mean);
    double F = p;
    double k = 0.0;
    while (U > F && k < 1000.0) {
        k += 1.0;
        p = p * mean / k;
        F += p;
    }
    return k;
}

// ----------------------------------------------------------------------------- one lesion
// Type index bits = resistances (bit 0: drug 1, bit 1: drug 2; T11X resists both);
// drugs = bit mask of the drugs in use. A cell is sensitive if one drug in use is not
// in its resistances.
static inline int resistance_bits(int type) { return type == T11X ? 3 : type; }

static inline bool sensitive_to(int type, int drugs) { return (drugs & ~resistance_bits(type)) != 0; }

// Number of resistance mutations carried (for the cost of resistance)
static inline int resistances(int type) { return type == T11X ? 1 : (type & 1) + ((type >> 1) & 1); }

// Mutation probabilities per division, from type i to type j
static inline void mutation_table(const Params &p, double mu[NTYPES][NTYPES]) {
    for (int i = 0; i < NTYPES; i++)
        for (int j = 0; j < NTYPES; j++) mu[i][j] = 0.0;
    mu[T00][T10] = p.u * p.n1;
    mu[T00][T01] = p.u * p.n2;
    mu[T00][(p.cost > 0.0) ? T11X : T11] = p.u * p.n12;
    mu[T10][T11] = p.u * (p.n2 + p.n12);
    mu[T01][T11] = p.u * (p.n1 + p.n12);
}

// One time step for all types; returns the new sizes in x. growth[i] = exp((b_i - d_i) dt)
// for deterministic types.
static inline void step(Rng &g, const Params &p, double x[NTYPES], const double birth[NTYPES],
                        const double death[NTYPES], const double growth[NTYPES],
                        const double mu[NTYPES][NTYPES]) {
    double add[NTYPES] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double nx[NTYPES];
    for (int i = 0; i < NTYPES; i++) {
        double n = x[i];
        if (n <= 0.0) {
            nx[i] = 0.0;
            continue;
        }
        double births;
        if (n < p.k_exact) {
            // A deterministic type that has shrunk below k_exact has a fractional size:
            // round it stochastically (keeps the mean), or it could never reach 0
            double whole = (double)(int64_t)n;
            if (n != whole) n = whole + (rng_uniform(g) < n - whole ? 1.0 : 0.0);
            if (n <= 0.0) {
                nx[i] = 0.0;
                continue;
            }
            births = poisson(g, n * birth[i] * p.dt);
            double deaths = poisson(g, n * death[i] * p.dt);
            double v = n + births - deaths;
            nx[i] = v > 0.0 ? v : 0.0;
        } else {
            births = n * birth[i] * p.dt;        // expected births (for the mutants)
            nx[i] = n * growth[i];
        }
        for (int j = 0; j < NTYPES; j++) {
            if (mu[i][j] > 0.0) {
                double m = poisson(g, births * mu[i][j]);
                if (m > 0.0) {
                    // a mutant offspring replaces a same-type offspring
                    add[j] += m;
                    nx[i] = nx[i] > m ? nx[i] - m : 0.0;
                }
            }
        }
    }
    for (int i = 0; i < NTYPES; i++) x[i] = nx[i] + add[i];
}

// Exact alternative to step() for the small types (17 October 2026; CPU only, so it may
// use libm): within the step, types below k_exact evolve event by event (Gillespie);
// mutants of the deterministic types arrive as Poisson numbers of events at uniform
// times in the step (exact for an inhomogeneous Poisson process whose rate varies little
// within dt); a waiting time that passes an arrival is redrawn from it (memorylessness).
static inline void step_exact(Rng &g, const Params &p, double x[NTYPES], const double birth[NTYPES],
                              const double death[NTYPES], const double growth[NTYPES],
                              const double mu[NTYPES][NTYPES]) {
    bool small[NTYPES];
    std::vector<double> arrival_t;
    std::vector<int> arrival_type;
    for (int i = 0; i < NTYPES; i++) {
        double n = x[i];
        small[i] = n < p.k_exact;
        if (n <= 0.0) continue;
        if (small[i]) {
            double whole = (double)(int64_t)n;
            if (n != whole) x[i] = whole + (rng_uniform(g) < n - whole ? 1.0 : 0.0);
            continue;
        }
        double births = n * birth[i] * p.dt;
        x[i] = n * growth[i];
        for (int j = 0; j < NTYPES; j++) {
            if (mu[i][j] <= 0.0) continue;
            double m = poisson(g, births * mu[i][j]);
            x[i] = x[i] > m ? x[i] - m : 0.0;
            for (int k = 0; k < (int)m; k++) {
                arrival_t.push_back(rng_uniform(g) * p.dt);
                arrival_type.push_back(j);
            }
        }
    }
    // arrivals in time order (insertion sort: few per step)
    for (size_t a = 1; a < arrival_t.size(); a++)
        for (size_t c = a; c > 0 && arrival_t[c] < arrival_t[c - 1]; c--) {
            double tt = arrival_t[c]; arrival_t[c] = arrival_t[c - 1]; arrival_t[c - 1] = tt;
            int ty = arrival_type[c]; arrival_type[c] = arrival_type[c - 1]; arrival_type[c - 1] = ty;
        }
    double t = 0.0;
    size_t next = 0;
    for (;;) {
        double rate = 0.0;
        for (int i = 0; i < NTYPES; i++)
            if (small[i] && x[i] > 0.0) rate += (birth[i] + death[i]) * x[i];
        double t_arr = next < arrival_t.size() ? arrival_t[next] : p.dt;
        double tau = rate > 0.0 ? -log(1.0 - rng_uniform(g)) / rate : 2.0 * p.dt;
        if (t + tau >= t_arr) {                     // next arrival (or the end of the step) first
            if (next >= arrival_t.size()) break;
            t = t_arr;
            x[arrival_type[next]] += 1.0;
            next++;
            continue;
        }
        t += tau;
        // which small type, then birth or death
        double pick = rng_uniform(g) * rate, acc = 0.0;
        int i = 0;
        for (int k = 0; k < NTYPES; k++) {
            if (!(small[k] && x[k] > 0.0)) continue;
            i = k;
            acc += (birth[k] + death[k]) * x[k];
            if (pick < acc) break;
        }
        if (rng_uniform(g) * (birth[i] + death[i]) < birth[i]) {
            int child = i;
            double v = rng_uniform(g), cum = 0.0;
            for (int j = 0; j < NTYPES; j++) {
                cum += mu[i][j];
                if (mu[i][j] > 0.0 && v < cum) {
                    child = j;
                    break;
                }
            }
            x[child] += 1.0;
        } else {
            x[i] -= 1.0;
        }
    }
}

static inline double total_cells(const double x[NTYPES]) {
    double t = x[0] + x[1] + x[2] + x[3];   // same order as before T11X existed
    return t + x[T11X];
}

// growth factors exp((b - d) dt) are passed in (computed by the caller) so that the
// kernel needs no exp of a positive argument
struct Rates {
    double pre_growth;     // exp((b - d) dt)
    double sens_growth;    // exp((b' - d') dt), sensitive cells during treatment
    double type_growth[NTYPES];  // exp((b_k - d) dt), b_k = (b - d)(1 - cost)^k + d
};

// Rates of a type that is not sensitive to the drugs in use (or before treatment)
static inline void own_rates(const Params &p, const Rates &r, int type, double &birth, double &death,
                             double &growth) {
    int k = resistances(type);
    if (p.cost == 0.0 || k == 0) {
        birth = p.b;   // exactly b (no rounding through b - d)
    } else {
        double f = 1.0;
        for (; k > 0; k--) f *= 1.0 - p.cost;
        birth = (p.b - p.d) * f + p.d;
    }
    death = p.d;
    growth = p.cost == 0.0 ? r.pre_growth : r.type_growth[type];
}

static inline void set_treatment(const Params &p, const Rates &r, int drugs, double birth[NTYPES],
                                 double death[NTYPES], double growth[NTYPES]) {
    for (int i = 0; i < NTYPES; i++) {
        if (sensitive_to(i, drugs)) {
            birth[i] = p.bp;
            death[i] = p.dp;
            growth[i] = r.sens_growth;
        } else {
            own_rates(p, r, i, birth[i], death[i], growth[i]);
        }
    }
}

static inline void simulate_lesion(const Params &p, const Rates &r, uint64_t seed,
                                   uint64_t replicate, Result &out) {
    Rng g;
    rng_seed(g, seed, replicate);
    double mu[NTYPES][NTYPES];
    mutation_table(p, mu);

    double birth[NTYPES], death[NTYPES], growth[NTYPES];
    double x[NTYPES];
    out.restarts = 0;

    // Pre-treatment, restarted while the founding lineage dies out
    double t = 0.0;
    for (;;) {
        for (int i = 0; i < NTYPES; i++) {
            own_rates(p, r, i, birth[i], death[i], growth[i]);
            x[i] = 0.0;
        }
        x[T00] = 1.0;
        t = 0.0;
        double total = 1.0;
        while (total > 0.0 && total < p.M) {
            if (p.exact_events) step_exact(g, p, x, birth, death, growth, mu);
            else step(g, p, x, birth, death, growth, mu);
            t += p.dt;
            total = total_cells(x);
        }
        if (total >= p.M) break;
        out.restarts++;
    }
    out.t_detect = t;
    for (int i = 0; i < NOUT; i++) out.at_detect[i] = x[i];
    out.at_detect[T11] += x[T11X];

    // Treatment
    int drugs = p.dual ? (p.sequential ? 1 : 3) : 1;
    set_treatment(p, r, drugs, birth, death, growth);
    double tt = 0.0;
    out.t_switch = -1.0;
    bool shrunk = false;   // sequential: the lesion fell below M on drug 1
    out.outcome = OUTCOME_UNDECIDED;
    while (tt < p.t_max_treat) {
        double total = total_cells(x);
        if (total <= 0.0) {
            out.outcome = OUTCOME_ERADICATED;
            break;
        }
        if (p.dual && p.sequential && drugs == 1) {
            if (total < p.M) shrunk = true;
            else if (shrunk) {   // relapse on drug 1 (back at M cells): switch to drug 2
                drugs = 2;
                set_treatment(p, r, drugs, birth, death, growth);
                out.t_switch = tt;
            }
        } else {
            bool failed = false;
            for (int i = 0; i < NTYPES; i++)
                if (!sensitive_to(i, drugs) && x[i] >= p.fail_size) failed = true;
            if (failed) {
                out.outcome = OUTCOME_FAILED;
                break;
            }
        }
        if (p.exact_events) step_exact(g, p, x, birth, death, growth, mu);
        else step(g, p, x, birth, death, growth, mu);
        tt += p.dt;
    }
    out.t_end = tt;
}

}  // namespace resistance

#endif

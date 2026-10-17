// Reference implementation of the cancer stem cell cellular automaton of
// Poleszczuk & Enderling (Applied Mathematics 2014;5:144-152, arXiv:1309.6015), as
// published: cells updated one at a time in a random order each step (asynchronous).
// (Also the reference for a parallel sub-lattice version developed in a separate project.)
//
// Rules (2014 paper, "Tumor Growth Model"; death at division attempts as in Poleszczuk,
// Macklin & Enderling, Methods Mol Biol 2016, PMC6587968):
//   2-D square lattice, Moore neighbourhood (8); a cell with no free neighbour is quiescent.
//   Each step (1 h), each cell present at the start of the step, in random order:
//     if a neighbour is free:
//       with probability pd (= 1/24: 24-h cell cycle) attempt division:
//         stem cell: daughter is a stem cell with probability ps, else a non-stem cell
//                    with proliferation capacity rho_max; stem cells never die
//         non-stem:  dies with probability alpha, or if its capacity is 0; otherwise
//                    capacity - 1, and the daughter gets the same capacity
//       otherwise, with probability pm (= mu dt in cell widths) migrate to a free neighbour
//   The free neighbour is chosen uniformly among the free ones.
// Parameters of the paper's simulations: cct 24 h, rho_max 10, mu 100 um/day (10 cell
// widths/day), alpha 1 %, ps 0.1, one stem cell in the centre, 180 days: ~140,000 cells.
//
// Build: g++ -O2 -fopenmp -o ca_reference ca_reference.cpp
// Usage: ca_reference [tumours] [days] [seed] [alpha] [ps] [rho_max] [lattice]

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

struct Cell {
    int place;   // linear index in the lattice
    int rho;     // remaining divisions (non-stem); -1 for stem cells
};

struct Result {
    long cells, stem, quiescent;
    int radius_max;
};

static int VARIANT = 0;
static uint64_t DUMP_SEED = 0;   // CA_DUMP writes the lattice of the tumour with this seed
static double PM = 10.0 / 24.0;   // migration probability per step   // 0 as described; 1 newborn cells act in the step they are born

static Result simulate(uint64_t seed, int days, double alpha, double ps, int rho_max, int L) {
    const double pd = 1.0 / 24.0, pm = PM;
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<uint8_t> occupied(L * L, 0);
    std::vector<Cell> cells;
    std::vector<int> index(L * L, -1);   // cell index at each place (-1: empty)
    const int centre = (L / 2) * L + L / 2;
    cells.push_back({centre, -1});
    occupied[centre] = 1;
    index[centre] = 0;
    const int dr[8] = {-1, -1, -1, 0, 0, 1, 1, 1}, dc[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

    auto free_neighbour = [&](int place) -> int {
        int r = place / L, c = place % L, free_list[8], n = 0;
        for (int k = 0; k < 8; k++) {
            int rr = r + dr[k], cc = c + dc[k];
            if (rr < 0 || rr >= L || cc < 0 || cc >= L) continue;   // lattice edge: no space
            int q = rr * L + cc;
            if (!occupied[q]) free_list[n++] = q;
        }
        if (n == 0) return -1;
        return free_list[std::min(n - 1, (int)(U(gen) * n))];
    };

    std::vector<int> order;
    for (int step = 0; step < days * 24; step++) {
        int n0 = (int)cells.size();
        order.resize(n0);
        for (int i = 0; i < n0; i++) order[i] = i;
        std::shuffle(order.begin(), order.end(), gen);
        std::vector<char> dead(n0, 0);
        for (int oi = 0; oi < (int)order.size(); oi++) {
            int i = order[oi];
            if (i >= (int)dead.size()) dead.resize(i + 1, 0);
            if (dead[i]) continue;
            Cell &cell = cells[i];
            int target = free_neighbour(cell.place);
            if (target < 0) continue;                       // quiescent
            if (U(gen) < pd) {                               // division attempt
                if (cell.rho < 0) {                          // stem cell
                    Cell daughter = {target, U(gen) < ps ? -1 : rho_max};
                    occupied[target] = 1;
                    index[target] = (int)cells.size();
                    if (VARIANT == 1) order.push_back((int)cells.size());
                    cells.push_back(daughter);
                } else if (cell.rho > 0 && U(gen) >= alpha) {
                    cell.rho--;
                    Cell daughter = {target, cell.rho};
                    occupied[target] = 1;
                    index[target] = (int)cells.size();
                    if (VARIANT == 1) order.push_back((int)cells.size());
                    cells.push_back(daughter);
                } else {                                     // death
                    occupied[cell.place] = 0;
                    index[cell.place] = -1;
                    dead[i] = 1;
                }
            } else if (U(gen) < pm) {                        // migration
                occupied[cell.place] = 0;
                index[cell.place] = -1;
                cell.place = target;
                occupied[target] = 1;
                index[target] = i;
            }
        }
        // remove dead cells (keep the others in order)
        std::vector<Cell> alive;
        alive.reserve(cells.size());
        for (int i = 0; i < (int)cells.size(); i++)
            if (i >= (int)dead.size() || !dead[i]) alive.push_back(cells[i]);
        cells.swap(alive);
        for (int i = 0; i < (int)cells.size(); i++) index[cells[i].place] = i;
    }
    const char *dump = getenv("CA_DUMP");
    if (dump && seed == DUMP_SEED) {
        std::vector<uint8_t> img(L * L, 0);
        for (const Cell &cell : cells) img[cell.place] = cell.rho < 0 ? 255 : 128;
        FILE *f = fopen(dump, "wb");
        if (f) {
            fprintf(f, "P5\n%d %d\n255\n", L, L);
            fwrite(img.data(), 1, img.size(), f);
            fclose(f);
        }
    }
    Result res = {(long)cells.size(), 0, 0, 0};
    for (const Cell &cell : cells) {
        if (cell.rho < 0) res.stem++;
        int r = cell.place / L, c = cell.place % L, n = 0;
        for (int k = 0; k < 8; k++) {
            int rr = r + dr[k], cc = c + dc[k];
            if (rr >= 0 && rr < L && cc >= 0 && cc < L && !occupied[rr * L + cc]) n++;
        }
        if (n == 0) res.quiescent++;
        int dist = std::max(std::abs(r - L / 2), std::abs(c - L / 2));
        res.radius_max = std::max(res.radius_max, dist);
    }
    return res;
}

int main(int argc, char **argv) {
    int tumours = argc > 1 ? atoi(argv[1]) : 8;
    int days = argc > 2 ? atoi(argv[2]) : 180;
    uint64_t seed = argc > 3 ? strtoull(argv[3], nullptr, 10) : 1;
    double alpha = argc > 4 ? atof(argv[4]) : 0.01;
    double ps = argc > 5 ? atof(argv[5]) : 0.1;
    int rho_max = argc > 6 ? atoi(argv[6]) : 10;
    int L = argc > 7 ? atoi(argv[7]) : 1000;
    VARIANT = argc > 8 ? atoi(argv[8]) : 0;
    PM = argc > 9 ? atof(argv[9]) : PM;
    DUMP_SEED = seed * 1000003ULL;   // tumour 0
    std::vector<Result> res(tumours);
#pragma omp parallel for schedule(dynamic, 1)
    for (int t = 0; t < tumours; t++) res[t] = simulate(seed * 1000003ULL + t, days, alpha, ps, rho_max, L);
    double mean = 0, mean_stem = 0, mean_q = 0;
    int edge = 0;
    for (auto &r : res) {
        mean += r.cells;
        mean_stem += r.stem;
        mean_q += (double)r.quiescent / std::max(1L, r.cells);
        if (r.radius_max >= L / 2 - 1) edge++;
    }
    mean /= tumours;
    mean_stem /= tumours;
    mean_q /= tumours;
    double sd = 0;
    for (auto &r : res) sd += (r.cells - mean) * (r.cells - mean);
    sd = tumours > 1 ? std::sqrt(sd / (tumours - 1)) : 0;
    printf("tumours %d, days %d: cells mean %.0f (sd %.0f), stem mean %.1f, quiescent %.1f %%, touching edge %d\n",
           tumours, days, mean, sd, mean_stem, 100 * mean_q, edge);
    for (auto &r : res) printf("%ld ", r.cells);
    printf("\n");
    // machine-readable lines for scripts: one per tumour
    for (int t = 0; t < tumours; t++)
        printf("TUMOUR %d cells %ld stem %ld quiescent %ld edge %d\n", t, res[t].cells, res[t].stem,
               res[t].quiescent, res[t].radius_max >= L / 2 - 1 ? 1 : 0);
    return 0;
}

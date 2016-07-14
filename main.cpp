#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <ctime>

#include "circular_stat.h"


/*void generate_random_angles(size_t n, complex* angles) {
    for (size_t i=0; i<n; ++i) {
        double a = 2 * M_PI * rand() / double(RAND_MAX);
        angles[i] = complex(cos(a), sin(a));
    }
}*/


int main() {
    srand(time(NULL));

    const size_t n_perms = 12;
    const size_t n_f_stats = 8;

    double* shuffled_f_stat_mat = new double[n_perms*n_f_stats];
    double* zscore_mat = new double[n_perms*n_f_stats];

    for (size_t i=0; i<n_perms*n_f_stats; ++i) shuffled_f_stat_mat[i] = rand() / double(RAND_MAX);

    compute_zscores(shuffled_f_stat_mat, n_perms*n_f_stats, zscore_mat, n_perms*n_f_stats, n_perms);

    fputs("Shuffled f-stats:\n", stdout);
    for (size_t i=0; i<n_perms; ++i) {
        for (size_t j=0; j<n_f_stats; ++j) printf("%8.4lg", shuffled_f_stat_mat[i*n_f_stats+j]);
        fputc('\n', stdout);
    }
    fputs("z-scores:\n", stdout);
    for (size_t i=0; i<n_perms; ++i) {
        for (size_t j=0; j<n_f_stats; ++j) printf("%8.4lg", zscore_mat[i*n_f_stats+j]);
        fputc('\n', stdout);
    }

    delete[] zscore_mat;
    delete[] shuffled_f_stat_mat;

    return 0;
}

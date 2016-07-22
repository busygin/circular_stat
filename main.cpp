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

    constexpr size_t n_freqs{3};
    constexpr size_t n_bps{8};
    constexpr size_t n_events{10};
    constexpr size_t t_size{16};
    constexpr size_t n_wavelets{n_freqs*n_bps*n_events*t_size};
    std::complex<double>* wavelets = new std::complex<double>[n_wavelets];
    bool* recalls = new bool[n_events];
    recalls[0]=true; recalls[2]=true; recalls[3]=true; recalls[4]=true; recalls[7]=true; recalls[8]=true; recalls[9]=true;
    recalls[1]=false; recalls[5]=false; recalls[6]=false;

    for (size_t freq=0; freq<n_freqs; ++freq)
        for (size_t bp=0; bp<n_bps; ++bp)
            for (size_t event=0; event<n_events; ++event)
                for (size_t t=0; t<t_size; ++t)
                    wavelets[t+t_size*(event+n_events*(bp+n_bps*freq))] =
                            recalls[event] ? std::complex<double>{1.0,0.0} : std::complex<double>{-1.0,0.0};

    constexpr size_t n_bp_pairs{n_bps*(n_bps-1)/2};
    constexpr size_t n_features{n_freqs*n_bp_pairs};
    constexpr size_t n_ppc_output{n_features*n_events};
    double* ppc_output = new double[n_ppc_output];

    single_trial_ppc_all_features(recalls, n_events, wavelets, n_wavelets, ppc_output, n_ppc_output, n_freqs, n_bps, 10);

    size_t feature_idx{0};
    for (size_t freq=0; freq<n_freqs; ++freq)
        for (size_t bp1=1; bp1<n_bps; ++bp1)
            for (size_t bp2=0; bp2<bp1; ++bp2)
                for (size_t event=0; event<n_events; ++event) {
                    printf("%s feature: %lg\n", recalls[event]?"Recall":"Nonrecall", ppc_output[feature_idx++]);
                }

    delete[] ppc_output;
    delete[] recalls;
    delete[] wavelets;

//    const size_t n_perms = 12;
//    const size_t n_f_stats = 8;
//
//    double* shuffled_f_stat_mat = new double[n_perms*n_f_stats];
//    double* zscore_mat = new double[n_perms*n_f_stats];
//
//    for (size_t i=0; i<n_perms*n_f_stats; ++i) shuffled_f_stat_mat[i] = rand() / double(RAND_MAX);
//
//    compute_zscores(shuffled_f_stat_mat, n_perms*n_f_stats, zscore_mat, n_perms*n_f_stats, n_perms);
//
//    fputs("Shuffled f-stats:\n", stdout);
//    for (size_t i=0; i<n_perms; ++i) {
//        for (size_t j=0; j<n_f_stats; ++j) printf("%8.4lg", shuffled_f_stat_mat[i*n_f_stats+j]);
//        fputc('\n', stdout);
//    }
//    fputs("z-scores:\n", stdout);
//    for (size_t i=0; i<n_perms; ++i) {
//        for (size_t j=0; j<n_f_stats; ++j) printf("%8.4lg", zscore_mat[i*n_f_stats+j]);
//        fputc('\n', stdout);
//    }
//
//    delete[] zscore_mat;
//    delete[] shuffled_f_stat_mat;

    return 0;
}

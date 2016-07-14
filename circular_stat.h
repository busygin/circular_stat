#ifndef CIRCULAR_STAT_CIRCULAR_STAT_H
#define CIRCULAR_STAT_CIRCULAR_STAT_H

#include <cstddef>
#include <complex>


// void circ_diff(size_t n, complex* c1, complex* c2, complex* cdiff);
void circ_diff(std::complex<double>* c1, size_t n1, std::complex<double>* c2, size_t n2,
               std::complex<double>* cdiff, size_t n3);

std::complex<double> resultant_vector(std::complex<double>* c, size_t n);

double resultant_vector_length(std::complex<double>* c, size_t n);

// double circ_mean(size_t n, complex* cdiff);
std::complex<double> circ_mean(std::complex<double>* c, size_t n);

// void circ_diff_time_bins(size_t n, size_t bin_len, complex* c1, complex* c2, complex* cdiff, double* cdiff_means);
void circ_diff_time_bins(std::complex<double>* c1, size_t n1, std::complex<double>* c2, size_t n2,
                         std::complex<double>* cdiff, size_t n3, std::complex<double>* cdiff_means, size_t n_bins);

void compute_f_stat(std::complex<double>* phase_diff_mat, size_t n_phase_diffs, bool* recalls, size_t n_events, double* f_stat_mat, size_t n_f_stats);

void compute_zscores(double* mat, size_t n_mat, size_t n_perms);

void single_trial_ppc(
        std::complex<double>* wavelet1, size_t n_phases1,
        std::complex<double>* wavelet2, size_t n_phases2,
        double* ppcs, size_t n_ppcs, size_t n_events);

void single_trial_ppc_with_classes(
        bool* recalls, size_t n_events,
        std::complex<double>* wavelet1, size_t n_phases1,
        std::complex<double>* wavelet2, size_t n_phases2,
        double* ppcs, size_t n_ppcs);

void single_trial_ppc_all_features(
        bool* recalls, size_t n_events,
        std::complex<double>* wavelets, size_t n_wavelets,
        double* ppc_output, size_t n_ppc_output,
        size_t n_freqs, size_t n_bps);

#endif  // CIRCULAR_STAT_CIRCULAR_STAT_H

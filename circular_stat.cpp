#include <cmath>
#include <cstring>
#include <thread>
#include "circular_stat.h"


void circ_diff(std::complex<double>* c1, size_t n1, std::complex<double>* c2, size_t n2, std::complex<double>* cdiff, size_t n3) {
    for (size_t i=0; i<n1; ++i) {
        std::complex<double> d{c1[i].real()*c2[i].real()+c1[i].imag()*c2[i].imag(), c1[i].imag()*c2[i].real()-c1[i].real()*c2[i].imag()};
        cdiff[i] = d / std::abs(d);
    }
}

std::complex<double> resultant_vector(std::complex<double>* c, size_t n) {
    std::complex<double> s{0.0,0.0};
    for (size_t i=0; i<n; ++i) s += c[i];
    return s;
}

double resultant_vector_length(std::complex<double>* c, size_t n) {
    std::complex<double> s{0.0,0.0};
    for (size_t i=0; i<n; ++i) s += c[i];
    return std::abs(s);
}

std::complex<double> circ_mean(std::complex<double>* c, size_t n) {
    std::complex<double> s{0.0,0.0};
    for (size_t i=0; i<n; ++i) s += c[i];
    // return atan2(s.imag(),s.real());
    return s / std::abs(s);
}

void circ_diff_time_bins(std::complex<double>* c1, size_t n1, std::complex<double>* c2, size_t n2,
                         std::complex<double>* cdiff, size_t n3, std::complex<double>* cdiff_means, size_t n_bins) {
    circ_diff(c1, n1, c2, n1, cdiff, n1);
    size_t bin_len = n1 / n_bins;
    for (size_t i=0; i<n_bins; ++i, cdiff+=bin_len) {
        cdiff_means[i] = circ_mean(cdiff, bin_len);
    }
}

void compute_f_stat(std::complex<double>* phase_diff_mat, size_t n_phase_diffs, bool* recalls, size_t n_events, double* f_stat_mat, size_t n_f_stats) {
    size_t n_recalls{0};
    for (size_t i=0; i<n_events; ++i) {
        if (recalls[i]) ++n_recalls;
    }
    size_t n_non_recalls = n_events - n_recalls;

    size_t n_comps = n_phase_diffs / n_events;
    for (size_t i=0; i<n_comps; ++i) {
        std::complex<double>* pdm_i = phase_diff_mat + (i*n_events);
        std::complex<double> s_rec{0.0,0.0}, s_nrec{0.0,0.0};
        for (size_t j=0; j<n_events; ++j) {
            if (recalls[j]) {
                s_rec += pdm_i[j];
                // printf("Recall angle: %lg\n", atan2(pdm_i[j].imag(),pdm_i[j].real()));
            } else {
                s_nrec += pdm_i[j];
                // printf("Nonrecall angle: %lg\n", atan2(pdm_i[j].imag(),pdm_i[j].real()));
            }
            // if (j<10) printf("%lg ", atan2(pdm_i[j].imag(),pdm_i[j].real()));
        }
        double r_recalls = std::abs(s_rec);
        double r_non_recalls = std::abs(s_nrec);
        //f_stat_mat[i] = ((n_recalls-1)*(n_non_recalls-r_non_recalls)) / ((n_non_recalls-1)*(n_recalls-r_recalls));
        f_stat_mat[i] = ((n_non_recalls-1)*(n_recalls-r_recalls)) / ((n_recalls-1)*(n_non_recalls-r_non_recalls));
        // printf("\nn_recalls %ld, n_non_recalls %ld, r_recalls %lg r_non_recalls %lg fstat %lg\n", n_recalls, n_non_recalls, r_recalls, r_non_recalls, f_stat_mat[i]);
        // getc(stdin);
    }
}

inline double sqr(double x) { return x*x; }

void compute_zscores(double* mat, size_t n_mat, size_t n_perms) {
    size_t n_stats = n_mat / n_perms;
    double* ss = new double[n_stats];
    double* sos = new double[n_stats];

    memset(ss, 0, n_stats*sizeof(double));
    memset(sos, 0, n_stats*sizeof(double));
    for (size_t i=0; i<n_perms; ++i) {
        double* mat_i = mat + (i*n_stats);
        for (size_t j=0; j<n_stats; ++j) {
            double x = mat_i[j];
            ss[j] += x;
            sos[j] += sqr(x);
        }
    }

    for (size_t i=0; i<n_perms; ++i) {
        double* mat_i = mat + (i*n_stats);
        for (size_t j=0; j<n_stats; ++j) {
            double x = mat_i[j];
            double ss1 = ss[j] - x;
            double sos1 = sos[j] - sqr(x);
            double mu = ss1 / (n_perms-1);
            double sigma = sqrt((sos1-sqr(ss1)/(n_perms-1)) / (n_perms-2));
            mat_i[j] = (x-mu) / sigma;
        }
    }

    delete[] sos;
    delete[] ss;
}

// wavelet1 and wavele2 are n_events X tsize, and n_phases1=n_phases2=n_events*tsize
void single_trial_ppc(
        std::complex<double>* wavelet1, size_t n_phases1,
        std::complex<double>* wavelet2, size_t n_phases2,
        double* ppcs, size_t n_ppcs, size_t n_events)
{
    std::complex<double>* phase_diff = new std::complex<double>[n_phases1];
    circ_diff(wavelet1, n_phases1, wavelet2, n_phases1, phase_diff, n_phases1);
    size_t t_size = n_phases1 / n_events;
    std::complex<double>* phase_diff_diff = new std::complex<double>[t_size];
    memset(ppcs, 0, n_phases1*sizeof(double));
    for (size_t i=1; i<n_events; ++i) {
        std::complex<double>* phase_diff_i = phase_diff + (i*t_size);
        double* ppcs_i = ppcs + (i*t_size);
        for (size_t j=0; j<i; ++j) {
            std::complex<double>* phase_diff_j = phase_diff + (j*t_size);
            double* ppcs_j = ppcs + (j*t_size);
            circ_diff(phase_diff_i, t_size, phase_diff_j, t_size, phase_diff_diff, t_size);
            for (size_t k=0; k<t_size; ++k) {
                double cos_ij = phase_diff_diff[k].real();
                ppcs_i[k] += cos_ij;
                ppcs_j[k] += cos_ij;
            }
        }
    }
    for (size_t i=0; i<n_phases1; ++i) ppcs[i] /= n_events-1;
    delete[] phase_diff_diff;
    delete[] phase_diff;
}

void single_trial_ppc_with_classes(
        bool* recalls, size_t n_events,
        std::complex<double>* wavelet1, size_t n_phases1,
        std::complex<double>* wavelet2, size_t n_phases2,
        double* ppcs, size_t n_ppcs)
{
    std::complex<double>* phase_diff = new std::complex<double>[n_phases1];
    circ_diff(wavelet1, n_phases1, wavelet2, n_phases1, phase_diff, n_phases1);
    size_t t_size = n_phases1 / n_events;
    std::complex<double>* phase_diff_diff = new std::complex<double>[t_size];
    memset(ppcs, 0, n_phases1*sizeof(double));
    size_t n_recalls{0};
    for (size_t i=0; i<n_events; ++i) {
        bool recall_i = recalls[i];
        if (recall_i) ++n_recalls;
        std::complex<double>* phase_diff_i = phase_diff + (i*t_size);
        double* ppcs_i = ppcs + (i*t_size);
        for (size_t j=0; j<i; ++j) {
            bool recall_j = recalls[j];
            if (recall_i == recall_j) {
                std::complex<double> *phase_diff_j = phase_diff + (j * t_size);
                double *ppcs_j = ppcs + (j * t_size);
                circ_diff(phase_diff_i, t_size, phase_diff_j, t_size, phase_diff_diff, t_size);
                for (size_t k = 0; k < t_size; ++k) {
                    double cos_ij = phase_diff_diff[k].real();
                    ppcs_i[k] += cos_ij;
                    ppcs_j[k] += cos_ij;
                }
            }
        }
    }
    size_t n_non_recalls = n_events - n_recalls;
    for (size_t i=0; i<n_events; ++i) {
        bool recall = recalls[i];
        double* ppc_i = ppcs + (i*t_size);
        for (size_t t=0; t<t_size; ++t) {
            ppc_i[t] /= (recall ? n_recalls-1 : n_non_recalls-1);
        }
    }
    delete[] phase_diff_diff;
    delete[] phase_diff;
}

void compute_feature(bool* recalls, size_t n_events, size_t t_size,
                     std::complex<double>* wavelets_bp1, std::complex<double>* wavelets_bp2,
                     double* ppcs, double* ppc_output, size_t feature_idx, size_t thread_idx) {
    size_t m_size{n_events*t_size};
    ppcs += thread_idx*m_size;
    single_trial_ppc_with_classes(recalls, n_events, wavelets_bp1, m_size, wavelets_bp2, m_size, ppcs, m_size);
    double* ppc_output_f_bp1_bp2 = ppc_output + (feature_idx*n_events);
    for (size_t e=0; e<n_events; ++e) {
        double s{0.0};
        double* ppcs_e = ppcs + (e*t_size);
        for (size_t t=0; t<t_size; ++t) {
            s += ppcs_e[t];
        }
        ppc_output_f_bp1_bp2[e] = s / t_size;
    }
}

// wavelets: n_freqs X n_bps X n_events X t_size array (flattened)
// ppc_output: (n_freqs * n_bps*(n_bps-1)/2) X n_events (flattened)
void single_trial_ppc_all_features(
        bool* recalls, size_t n_events,
        std::complex<double>* wavelets, size_t n_wavelets,
        double* ppc_output, size_t n_ppc_output,
        size_t n_freqs, size_t n_bps, size_t n_threads)
{
    std::thread threads[n_threads];

    size_t t_size = n_wavelets / (n_freqs*n_bps*n_events);
    double* ppcs = new double[n_threads*n_events*t_size];
    size_t feature_idx{0};
    size_t thread_idx{0};
    for (size_t f=0; f<n_freqs; ++f) {
        printf("Frequency %ld\n", f);
        std::complex<double>* wavelets_f = wavelets + (f*n_bps*n_events*t_size);
        for (size_t bp1=1; bp1<n_bps; ++bp1) {
            printf("Bipolar pair %ld\n", bp1);
            std::complex<double>* wavelets_bp1 = wavelets_f + (bp1*n_events*t_size);
            for (size_t bp2=0; bp2<bp1; ++bp2) {
                std::complex<double>* wavelets_bp2 = wavelets_f + (bp2*n_events*t_size);
                threads[thread_idx] = std::thread(compute_feature, recalls, n_events, t_size, wavelets_bp1, wavelets_bp2, ppcs, ppc_output, feature_idx, thread_idx);
                ++thread_idx;
                ++feature_idx;
                if (thread_idx == n_threads) {
                    for (auto& tr : threads) tr.join();
                    thread_idx = 0;
                }
            }
        }
    }

    for (size_t i=0; i<thread_idx; ++i) threads[i].join();

    delete[] ppcs;
}

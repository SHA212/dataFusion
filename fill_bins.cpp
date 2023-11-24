#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> fill_bins(Rcpp::DataFrame spectra,
                              Rcpp::DataFrame bins,
                              double rt_bin_size,
                              double mz_bin_size,
                              double overlap = 0.1,
                              Rcpp::String summaryFunction = "max") {

  // add check for spectra col names and thorw error when not correct

  std::vector<double> ints = spectra["intensity"];

  std::vector<double> rts  = spectra["rt"];

  std::vector<double> mzs = spectra["mz"];

  std::vector<double> bin_rts = bins["rt"];

  std::vector<double> bin_mzs = bins["mz"];

  const int n_rt_bins = bin_rts.size();

  const int n_mz_bins = bin_mzs.size();

  // Rcpp::Rcout << n_rt_bins << std::endl;

  const int n_rts = rts.size();

  const int n_mzs = mzs.size();

  const int n_ints = ints.size();

  if (!(n_rts == n_ints && n_rts == n_mzs)) {
    throw std::range_error("rts, mzs, and ints must be the same size!");
  }

  double half_rt_bin = rt_bin_size / 2;

  double half_mz_bin = mz_bin_size / 2;

  std::vector<double> output(n_rt_bins, 0.0);

  const double* ints_ptr = ints.data();

  // loops over each bin
  for(int i = 0; i < n_rt_bins; i++){

    // Rcpp::Rcout << "Bin: " << i << std::endl;

    std::vector<double> ints_temp;

    double bin_rt_low = bin_rts[i] - half_rt_bin;
    double bin_rt_high = bin_rts[i] + half_rt_bin;

    double bin_mz_low = bin_mzs[i] - half_mz_bin;
    double bin_mz_high = bin_mzs[i] + half_mz_bin;

    bin_rt_low = bin_rt_low - (rt_bin_size * overlap);
    bin_rt_high = bin_rt_high + (rt_bin_size * overlap);

    bin_mz_low = bin_mz_low - (mz_bin_size * overlap);
    bin_mz_high = bin_mz_high + (mz_bin_size * overlap);

    // Make case when traces are outside last mz and rt bins

    std::vector<int> which_ints;

    for (int z = 0; z < n_ints; ++z) {
      if (rts[z] >= bin_rt_low && rts[z] <= bin_rt_high && mzs[z] >= bin_mz_low && mzs[z] <= bin_mz_high) {
        which_ints.push_back(z);
      }
    }

    int n_sel_ints = which_ints.size();

    // Rcpp::Rcout << n_sel_ints << std::endl;

    if (n_sel_ints > 0) {
      std::vector<double> ints_temp(n_sel_ints);
      double* ints_temp_ptr = ints_temp.data();

      for (const int& x : which_ints) {
        *(ints_temp_ptr++) = *(ints_ptr + x);
      }

      double summary_int;

      if (summaryFunction == "max") {
        auto maxElementIterator = std::max_element(ints_temp.begin(), ints_temp.end());
        summary_int = *maxElementIterator;

      } else {
        double sum = std::accumulate(ints_temp.begin(), ints_temp.end(), 0.0);
        summary_int = sum / static_cast<double>(ints_temp.size());
      }

      output[i] = summary_int;
    }
  } // end loop over each bin

  return output;
}

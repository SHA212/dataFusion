#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> fillMat(std::vector<double> ints,
                            std::vector<double> rts,
                            std::vector<double> mzs,
                            std::vector<double> bin_rts,
                            std::vector<double> bin_mzs,
                            double half_rt_bin,
                            double half_mz_bin ) {
  
  const int n_rt_bins = bin_rts.size();
  
  const int n_mz_bins = bin_mzs.size();
  
  const int n_bins = n_rt_bins * n_mz_bins;
  
  const int n_rts = rts.size();
  
  const int n_mzs = mzs.size();
  
  const int n_ints = ints.size();
  
  
  if (!(n_rts == n_ints && n_rts == n_mzs)) {
    throw std::range_error("rts, mzs, and ints must be the same size!");
  }
  
  std::vector<std::tuple<double, int>> output(n_bins);
  
  int x = 0; // go through sample
  
  for(int i = 0; i < n_rt_bins; i++){ //check each bin's associated rt
    
    bool in_rt_bin = (rts[x] >= bin_rts[i] - half_rt_bin);
    if (!in_rt_bin) --i;
    
    for(int j = 0; j < n_mz_bins; j++){ //check each rt's associated mz
      
      int output_bin = (i * n_mz_bins) + j;
      
      std::vector<double> ints_temp;
      
      bool in_the_bin = true;
      
      while(in_the_bin){
        
        if (x >= n_rts) break;
        in_the_bin = rts[x] <= (bin_rts[i] + half_rt_bin) && mzs[x] <= (bin_mzs[j] + half_mz_bin);
        
        if (!in_the_bin || mzs[x] < bin_mzs[j] - half_mz_bin){
          break;
        }
        
        if (rts[x] >= bin_rts[i] - half_rt_bin && rts[x] <=  bin_rts[i] + half_rt_bin){
          if (mzs[x] >= bin_mzs[j] - half_mz_bin && mzs[x] <=  bin_mzs[j] + half_mz_bin){
            ints_temp.push_back(ints[x]);
          }
        }
        
        x++;
        if (x >= n_rts) break;
      }
      
      int n_ints_temp = ints_temp.size();
      
      if (n_ints_temp > 0) {
        
        for (int a = 0; a < n_ints_temp; ++a){
          std::get<0>(output[output_bin]) = std::get<0>(output[output_bin]) + ints_temp[a];
        }
        
      } else {
        
        std::get<0>(output[output_bin]) = std::get<0>(output[output_bin]) + 0;
      }
      
      std::get<1>(output[output_bin]) = std::get<1>(output[output_bin]) + n_ints_temp;
      
      if (x >= n_rts) break;
      
    }
    if (x >= n_rts) break;
  }
  int output_size = output.size();
  
  std::vector<double> final_output(output_size, 0);
  
  for (int t = 0; t < output_size; t++) {
    
    if (std::get<0>(output[t]) > 0){
      final_output[t] = std::get<0>(output[t])/std::get<1>(output[t]);
    }
    else final_output[t] = 0;
  }
  return final_output;
}

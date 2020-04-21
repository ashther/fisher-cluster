#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List minLossSplitCpp(NumericMatrix dMtx, NumericMatrix min_loss_mtx, int n, int k) {
  List result;
  NumericVector::iterator it;
  NumericVector value(n - k + 1);
  
  // we construct dMtx(n, k) which mean D(k, n)
  // so be careful about row and column change between equation and code
  
  if (k == 2) {
    // when k is 2, loss value is the minimum value of 
    // D(1, j - 1) + D(j, n), which j is from 2 to n
    
    for (int j = k; j <= n; j++) {
      if (j == 2) {
        // when j is 2, loss value is D(1, 1) + D(2, n) = D(2, n)
        // 
        // because the reason 1: Cpp use 0-index and 
        // reason 2: dMtx rowname starts from 2, use idr = n - 2
        // 
        // because the reason 1: Cpp use 0-index, use idc = 2 - 1 = 1
        value[j - k] = dMtx(n - 2, 1);
        
      } else if (j == n) {
        // when j is n, loss value is D(1, n - 1) + D(n, n) = D(1, n - 1)
        // 
        // because the reason 1: Cpp use 0-index and 
        // reason 2: dMtx rowname starts from 2, use idr = (n - 1) - 2 = j - 3
        // 
        // because the reason 1: Cpp use 0-index, use idc = 1 - 1 = 0
        value[j - k] = dMtx(j - 3, 0);
        
      } else {
        // normally, loss value is D(1, j - 1) + D(j, n)
        // 
        // because the reason 1: Cpp use 0-index and 
        // reason 2: dMtx rowname starts from 2, 
        // for former, use idr = (j - 1) - 2 = j - 3
        // for latter, use idr = n - 2
        // 
        // because the reason 1: Cpp use 0-index
        // for former, use idc = 1 - 1 = 0
        // for latter, use idc = j - 1
        value[j - k] = dMtx(j - 3, 0) + dMtx(n - 2, j - 1);
      }
    }
  } else {
    // normally, loss value is the minimum value of 
    // L(j - 1, k - 1) + D(j, n), which j is from k, n
    
    for (int j = k; j <= n; j++) {
      if (j == n) {
        // when j is n, loss value is L(j - 1, k - 1) + D(n, n) = L(n - 1, k - 1)
        // 
        // because the reason 1: Cpp use 0-index and 
        // reason 2: L rowname starts from 3, use idr = (n - 1) - 1 - 2 = j - 4
        // 
        // because the reason 1: Cpp use 0-index and
        // reason 2: L colname starts from 2, use use idc = (k - 1) - 1 - 1 = k - 3
        if ((j - 4) < 0 || (k - 3) < 0) {
          value[j - k] = 0;
        } else {
          value[j - k] = min_loss_mtx(j - 4, k - 3);
        }
      } else {
        // normally, loss value is L(j - 1, k - 1) + D(j, n)
        // 
        // because the reason 1: Cpp use 0-index and 
        // reason 2: L rowname starts from 3, and dMtx rowname starts from 2
        // for former use idr = (j - 1) - 1 - 2 = j - 4
        // for latter use idr = n - 2
        // 
        // because the reason 1: Cpp use 0-index and
        // reason 2: L colname starts from 2, and dMtx colname starts from 1
        // for former use idc = (k - 1) - 1 - 1 = k - 3
        // for latter use idc = j - 1
        if ((j - 4) < 0 || (k - 3) < 0) {
          value[j - k] = dMtx(n - 2, j - 1);
        } else {
          value[j - k] = min_loss_mtx(j - 4, k - 3) + dMtx(n - 2, j - 1);
        }
      }
    }
  }
  
  double min_value = min(value);
  int idx = which_min(value) + k;
  result["value"] = min_value;
  result["idx"] = idx;
  return result;
}



/*** R
# minLossSplitCpp(dMtx, min_loss_mtx, 11, 3)
*/
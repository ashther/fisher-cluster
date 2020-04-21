// #include <RcppParallel.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
// using namespace RcppParallel;

double diameterMtx(NumericMatrix x) {
  int n_row = x.nrow();
  int n_col = x.ncol();
  NumericVector dist(n_row);
  NumericVector col_means(n_col);
  double col_means_square = 0;

  for (int j = 0; j < n_col; j++) {
    double temp = 0;
    for (int i = 0; i < n_row; i++) {
      temp += x(i, j);
    }
    col_means(j) = temp / n_row;
    col_means_square += pow(col_means(j), 2);
  }
  col_means_square = sqrt(col_means_square);

  // TODO this process can be parallel
  for (int i = 0; i < n_row; i++) {
    double x_square = 0;
    double xy = 0;
    for (int j = 0; j < n_col; j++) {
      x_square += pow(x(i, j), 2);
      xy += x(i, j) * col_means(j);
    }
    dist(i) = 1 - xy / (sqrt(x_square) * col_means_square);
  }

  return sum(pow(dist, 2));
}

// [[Rcpp::export]]
NumericMatrix dMtxCreateCpp(NumericMatrix dtm) {
  int n_row = dtm.nrow() - 1;
  int n_col = dtm.ncol();
  NumericMatrix dMtx(n_row, n_row);
  
  for (int k = 0; k < n_col; k++) {
    for (int n = k; n < n_row; n++) {
      dMtx(n, k) = diameterMtx(
        dtm(Range(k, n + 1), Range(0, n_col - 1))
      );
    }
  }
  
  return dMtx;
}

// use suger function, simple but slow ----------------------------------------

double diameterMtxSuger(NumericMatrix x) {
  int n_row = x.nrow();
  int n_col = x.ncol();
  NumericVector dist(n_row);
  NumericVector col_means(n_col);
  double col_means_square = 0;
  
  for (int j = 0; j < n_col; j++) {
    col_means(j) = mean(x(_, j));
  }
  col_means_square = sqrt(sum(pow(col_means, 2)));
  
  for (int i = 0; i < n_row; i++) {
    double x_square = sqrt(sum(pow(x(i, _), 2)));
    double xy = 0;
    for (int j = 0; j < n_col; j++) {
      xy += x(i, j) * col_means(j);
    }
    dist(i) = 1 - xy / (x_square * col_means_square);
  }
  
  return sum(pow(dist, 2));
}

NumericMatrix dMtxCreateSuger(NumericMatrix dtm) {
  int n_row = dtm.nrow() - 1;
  int n_col = dtm.ncol();
  NumericMatrix dMtx(n_row, n_row);
  
  for (int k = 0; k < n_col; k++) {
    for (int n = k; n < n_row; n++) {
      dMtx(n, k) = diameterMtxSuger(
        dtm(Range(k, n + 1), Range(0, n_col - 1))
      );
    }
  }
  
  return dMtx;
}

// try parallel ---------------------------------------------------------------
// can't resolve the problems about submatrix of RMatrix ----------------------

// struct DMTX: public Worker {
//   const int n_row;
//   const int n_col;
//   const RMatrix<double> dtm;
//   RMatrix<double> result;
//   
//   DMTX(const int n_row, const int n_col, const NumericMatrix dtm, NumericMatrix result) : 
//     n_row(n_row), n_col(n_col), dtm(dtm), result(result) {}
//   
//   void operator()(std::size_t begin, std::size_t end) {
//     
//     for (std::size_t k = begin; k < end; k++) {
//       for (int n = k; n < n_row; n++) {
//         // can't debug this error
//         result(n, k) = diameterMtx(dtm(Range(k, n + 1), Range(0, n_col - 1)));
//       }
//     }
//     
//   }
// };
// 
// NumericMatrix dMtxCreateParCpp(NumericMatrix dtm) {
//   int n_row = dtm.nrow() - 1;
//   int n_col = dtm.ncol();
//   NumericMatrix result(n_row, n_row);
//   
//   DMTX dMtx(n_row, n_col, dtm, result);
//   
//   parallelFor(0, n_col, dMtx);
//   
//   return result;
// };


/*** R
# x <- matrix(1:16, 4)
# sum(text2vec::dist2(rbind(x, colMeans(x)))[5, ] ^ 2)
# diameterMtx(x)

dtm <- matrix(rnorm(25), 5)
# dMtxCreate(dtm)
dMtxCreateCpp(dtm)
# dMtxCreateParCpp(dtm)
*/

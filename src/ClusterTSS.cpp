#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector clusterTSS(DataFrame df, DataFrame clusters) {
  int nRows = df.nrows();
  int prev;
  int nxt = -1;
  int i = 0;
  
  NumericVector rowNo(0);

  NumericVector iter = df["iterations.x"];
  NumericVector clust = clusters["clusters"];

  while (nxt < nRows - 1) {
    prev = nxt + 1;
    nxt = nxt + clust[i];
    
    int maxIter = 0;
    int maxIterIndex = 0;
    
    for(int j = prev; j <= nxt; j++){
      if(iter[j] > maxIter) {
        maxIter = iter[j];
        maxIterIndex = j - prev + 1;
      }
    }
    
    rowNo.push_back(prev + maxIterIndex);
    i++;
  }
    return rowNo;
}

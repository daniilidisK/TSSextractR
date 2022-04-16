#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector firstSTrand(DataFrame df) {
  int nRows = df.nrows();
  
  CharacterVector chr = df["V1"];
  NumericVector start = df["V2"];
  NumericVector end = df["V3"];
  CharacterVector strand = df["V6"];
  
  NumericVector newCord(0);
  
  for (int i = 0; i < nRows; i++) {
    if (strand[i] == "+") {
      newCord.push_back(start[i]);
    } else if (strand[i] == "-") {
      newCord.push_back(end[i]);
    }
  }
  return newCord;
}

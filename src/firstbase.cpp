#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List firstBase(DataFrame df) {
  int nRows = df.nrows();
  
  NumericVector start = df["V2"];
  NumericVector end = df["V3"];
  CharacterVector strand = df["V6"];
  CharacterVector cigar = df["V7"];
  
    for (int i = 0; i < nRows; i++) {
        if (strand[i] == "+") {
          end[i] = start[i] + 1;
          cigar[i] = "1M";
        } else if (strand[i] == "-") {
          start[i] = end[i] - 1;
          cigar[i] = "1M";
        }
    }
    return {start, end, cigar};
}

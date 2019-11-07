#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix transform_pts( Rcpp::NumericMatrix pts, Rcpp::NumericMatrix tmat ) {

	// Declare variables
	int i, j;
	int nrow = pts.nrow();
	int ncol = pts.ncol();

	// Create intermediate point matrix
	Rcpp::NumericMatrix pts_copy(nrow, ncol);

	// Copy matrix
	std::copy(pts.begin(), pts.end(), pts_copy.begin());

	// Matrix multiply
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < 3; j++) {
			pts_copy(i,j) = pts(i,0)*tmat(j,0) + pts(i,1)*tmat(j,1) + pts(i,2)*tmat(j,2) + tmat(j,3);
		}
	}

	return pts_copy;
}
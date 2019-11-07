#include <Rcpp.h>

// [[Rcpp::export]]
double min_pt_cloud_dist( Rcpp::NumericMatrix mat1, Rcpp::NumericMatrix mat2 ) {

	// Declare variables
	int i, j;
	int nrow1 = mat1.nrow();
	int nrow2 = mat2.nrow();
	double d_sum = 0;

	// For each row in mat1
	for (i = 0; i < nrow1; i++) {

		// Create vector to hold distances
		std::vector<double> pt_to_pt(nrow2);

		// Measure "distance" to each point, sum of abs difference in each dimension
		for (j = 0; j < nrow2; j++) {
			pt_to_pt[j] = std::abs(mat1(i,0) - mat2(j,0)) + std::abs(mat1(i,1) - mat2(j,1)) + std::abs(mat1(i,2) - mat2(j,2));
		}

		// Save minimum
		d_sum += *(pt_to_pt.begin() + std::distance(pt_to_pt.begin(), min_element(pt_to_pt.begin(), pt_to_pt.end())));
	}
	
	return d_sum;
}
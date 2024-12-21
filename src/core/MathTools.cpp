#include "MathTools.h"
#include "blas.h"

namespace vtls {

	void addArraysImag(int len, std::complex<double>* arr1, double* arr2targ) {
		for (int i = 0; i < len; i++)
			arr2targ[i] += std::imag(arr1[i]);
	}

	void scaMulArray(int len, std::complex<double> scalar, std::complex<double> *__restrict arr) {
		/*for (int i = 0; i < len; i++)
			arr[i] *= scalar;*/
		cblas_zscal(len, &scalar, arr, 1);
	}

	void scaMulArray(int len, double scalar, double* __restrict arr) {
		cblas_dscal(len, scalar, arr, 1);
	}

	void scaMulArray(int len, double scalar, std::complex<double> *__restrict arr) {
		/*for (int i = 0; i < len; i++)
			arr[i] *= scalar;*/
		cblas_zdscal(len, scalar, arr, 1);
	}

	void copyArray(int len, double *__restrict arr1, double *__restrict arr2) {
		/*for (int i = 0; i < len; i++)
			arr2[i] = arr1[i];*/
		cblas_dcopy(len, arr1, 1, arr2, 1);
	}

	void copyArray(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2) {
		/*for (int i = 0; i < len; i++)
			arr2[i] = arr1[i];*/
		cblas_zcopy(len, arr1, 1, arr2, 1);
	}

	void copyArray(int len, double* __restrict arr1, std::complex<double>* __restrict arr2) {
		for (int i = 0; i < len; i++)
			arr2[i] = arr1[i];
	}

	template <typename T>
	void copyArrayRe(int len, T* __restrict arr1, double* __restrict arr2){
		for (int i = 0; i < len; i++)
			arr2[i] = std::real(arr1[i]);
	}

	int findValue(int len, double *__restrict arr, double val) {
		for (int i = 0; i < len; i++)
			if (arr[i] >= val)
				return i;
		if (val > arr[len - 1])
			return len - 1;
		else
			return -1;
	}

	//idxs should be initialized by user
	void insertSort_idxs(int len, double *__restrict arr, int *__restrict idxs) {

		int i = 1, j;
		double tempA, tempI;

		while (i < len) {

			j = i;

			while (j > 0 && arr[j - 1] > arr[j]) {

				tempA = arr[j];
				tempI = idxs[j];

				arr[j] = arr[j - 1];
				idxs[j] = idxs[j - 1];

				arr[j - 1] = tempA;
				idxs[j - 1] = tempI;

				j--;

			}

			i++;

		}
	}
}

namespace vtlsPrnt {
	void printGraph(int n, double *__restrict arr) {
		double minVal = *std::min_element(arr, arr + n);
		double maxVal = *std::max_element(arr, arr + n);
		int * nArr = (int*) sq_malloc(n * sizeof(int));
		int k = n / 50 + 1;
		for (int i = 0; i < n; i += k) nArr[i] = (int)(((arr[i] - minVal) / (maxVal - minVal)) * 100.0);
		//printArray(n, nArr);
		for (int i = 0; i < n; i += k) {
			for (int j = 0; j < nArr[i]; j++)
				std::cout << "#";
			std::cout << std::endl;
		}

		sq_free(nArr);
	}

	void printGraph(int n, std::complex<double>* __restrict arr0) {
		double* arr = (double*) sq_malloc(n * sizeof(double));
		for (int i = 0; i < n; i++)
			arr[i] = std::real(arr0[i]);
		double minVal = *std::min_element(arr, arr + n);
		double maxVal = *std::max_element(arr, arr + n);
		int* nArr = (int*) sq_malloc(n * sizeof(int));
		int k = n / 50 + 1;
		for (int i = 0; i < n; i += k) nArr[i] = (int)(((arr[i] - minVal) / (maxVal - minVal)) * 100.0);
		//printArray(n, nArr);
		for (int i = 0; i < n; i += k) {
			for (int j = 0; j < nArr[i]; j++)
				std::cout << "#";
			std::cout << std::endl;
		}

		sq_free(arr);
		sq_free(nArr);
	}
}

namespace plotting {

	void GNUPlotter::update(int nPts, int nLines, double* x, double* y) {
		update(nPts, nLines, x, y, x[0], x[nPts - 1], vtls::min(nPts*nLines, y), vtls::max(nPts*nLines, y));
	}

	void GNUPlotter::update(int nPts, int nLines, double* y){
		double* x = (double*)sq_malloc(nPts * sizeof(double));
		for(int i = 0; i < nPts; i++)
			x[i] = i;

		update(nPts, nLines, x, y);

		sq_free(x);
	}

	void GNUPlotter::update(int nPts, int nLines, double* y, double ymin, double ymax){
		double* x = (double*)sq_malloc(nPts * sizeof(double));
		for(int i = 0; i < nPts; i++)
			x[i] = i;

		update(nPts, nLines, x, y, x[0], x[nPts - 1], ymin, ymax);

		sq_free(x);
	}

	void GNUPlotter::update(int nPts, int nLines, double* x, double* y, double xmin, double xmax, double ymin, double ymax){
		std::vector<double> xv(x, x+nPts), yv;

		gp << "set xrange [" << xmin << ":" << xmax << "]\n";
		gp << "set yrange [" << ymin << ":" << ymax << "]\n";

		gp << "plot ";
		for (int i = 0; i < nLines; i++)
			gp << "'-' with lines title '" << i << (i < nLines - 1 ? "'," : "'\n");
		gp.flush();
		for (int i = 0; i < nLines; i++){
			yv = std::vector<double>(y + i*nPts, y + (i + 1)*nPts);
			gp.send1d(boost::make_tuple(xv, yv));
		}
	}

}
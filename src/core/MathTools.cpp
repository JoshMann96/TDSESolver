#include "MathTools.h"
#include "blas.h"

namespace vtls {

	void addArraysImag(size_t len, const std::complex<double>* arr1, double* arr2targ) {
		for (size_t i = 0; i < len; i++)
			arr2targ[i] += std::imag(arr1[i]);
	}

	void scaMulArray(size_t len, std::complex<double> scalar, std::complex<double> *__restrict arr) {
		/*for (size_t i = 0; i < len; i++)
			arr[i] *= scalar;*/
		cblas_zscal(len, &scalar, arr, 1);
	}

	void scaMulArray(size_t len, double scalar, double* __restrict arr) {
		cblas_dscal(len, scalar, arr, 1);
	}

	void scaMulArray(size_t len, double scalar, std::complex<double> *__restrict arr) {
		/*for (size_t i = 0; i < len; i++)
			arr[i] *= scalar;*/
		cblas_zdscal(len, scalar, arr, 1);
	}

	void copyArray(size_t len, const double *__restrict arr1, double *__restrict arr2) {
		/*for (size_t i = 0; i < len; i++)
			arr2[i] = arr1[i];*/
		cblas_dcopy(len, arr1, 1, arr2, 1);
	}

	void copyArray(size_t len, const std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2) {
		/*for (size_t i = 0; i < len; i++)
			arr2[i] = arr1[i];*/
		cblas_zcopy(len, arr1, 1, arr2, 1);
	}

	void copyArray(size_t len, const double* __restrict arr1, std::complex<double>* __restrict arr2) {
		for (size_t i = 0; i < len; i++)
			arr2[i] = arr1[i];
	}

	template <typename T>
	void copyArrayRe(size_t len, const T* __restrict arr1, double* __restrict arr2){
		for (size_t i = 0; i < len; i++)
			arr2[i] = std::real(arr1[i]);
	}

	size_t findValue(size_t len, const double *__restrict arr, double val) {
		for (size_t i = 0; i < len; i++)
			if (arr[i] >= val)
				return i;
		if (val > arr[len - 1])
			return len - 1;
		else
			return 0;
	}

	//idxs should be initialized by user
	void insertSort_idxs(size_t len, double *__restrict arr, size_t *__restrict idxs) {

		size_t i = 1, j;
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

	std::unique_ptr<double[]> getPolynomialSmoothBoundary(size_t len, size_t inner, size_t outer, double rate) {
		std::unique_ptr<double[]> mask(new double[len], std::default_delete<double[]>());
		std::fill_n(mask.get(), len, 1.0);
		size_t size;
		double k;
		if (outer < inner) {
			size = inner - outer;
			for (size_t i = 0; i < size; i++) {
				k = (double)(i + 1) / size;
				mask[i + outer] = std::pow((924.0 * std::pow(k, 13) -
					6006.0 * std::pow(k, 12) +
					16380.0 * std::pow(k, 11) -
					24024.0 * std::pow(k, 10) +
					20020.0 * std::pow(k, 9) -
					9009.0 * std::pow(k, 8) +
					1716.0 * std::pow(k, 7)) , (rate));
			}
		}
		else {
			size = outer - inner;
			for (size_t i = 0; i < size; i++) {
				k = (double)(size - i) / size;
				mask[i + inner] = std::pow((924.0 * std::pow(k, 13) -
					6006.0 * std::pow(k, 12) +
					16380.0 * std::pow(k, 11) -
					24024.0 * std::pow(k, 10) +
					20020.0 * std::pow(k, 9) -
					9009.0 * std::pow(k, 8) +
					1716.0 * std::pow(k, 7)) , (rate));
			}
		}
		return mask;
	}
}

namespace vtlsPrnt {
	void printGraph(size_t n, const double *__restrict arr) {
		double minVal = *std::min_element(arr, arr + n);
		double maxVal = *std::max_element(arr, arr + n);
		size_t * nArr = (size_t*) sq_malloc(n * sizeof(size_t));
		size_t k = n / 50 + 1;
		for (size_t i = 0; i < n; i += k) nArr[i] = (size_t)(((arr[i] - minVal) / (maxVal - minVal)) * 100.0);
		//printArray(n, nArr);
		for (size_t i = 0; i < n; i += k) {
			for (size_t j = 0; j < nArr[i]; j++)
				std::cout << "#";
			std::cout << std::endl;
		}

		sq_free(nArr);
	}

	void printGraph(size_t n, const std::complex<double>* __restrict arr0) {
		double* arr = (double*) sq_malloc(n * sizeof(double));
		for (size_t i = 0; i < n; i++)
			arr[i] = std::real(arr0[i]);
		double minVal = *std::min_element(arr, arr + n);
		double maxVal = *std::max_element(arr, arr + n);
		size_t* nArr = (size_t*) sq_malloc(n * sizeof(size_t));
		size_t k = n / 50 + 1;
		for (size_t i = 0; i < n; i += k) nArr[i] = (size_t)(((arr[i] - minVal) / (maxVal - minVal)) * 100.0);
		//printArray(n, nArr);
		for (size_t i = 0; i < n; i += k) {
			for (size_t j = 0; j < nArr[i]; j++)
				std::cout << "#";
			std::cout << std::endl;
		}

		sq_free(arr);
		sq_free(nArr);
	}
}

namespace plotting {

	void GNUPlotter::update(size_t nPts, size_t nLines, const double* x, const double* y) {
		update(nPts, nLines, x, y, x[0], x[nPts - 1], vtls::min(nPts*nLines, y), vtls::max(nPts*nLines, y));
	}

	void GNUPlotter::update(size_t nPts, size_t nLines, const double* y){
		double* x = (double*)sq_malloc(nPts * sizeof(double));
		for(size_t i = 0; i < nPts; i++)
			x[i] = i;

		update(nPts, nLines, x, y);

		sq_free(x);
	}

	void GNUPlotter::update(size_t nPts, size_t nLines, const double* y, double ymin, double ymax){
		double* x = (double*)sq_malloc(nPts * sizeof(double));
		for(size_t i = 0; i < nPts; i++)
			x[i] = i;

		update(nPts, nLines, x, y, x[0], x[nPts - 1], ymin, ymax);

		sq_free(x);
	}

	void GNUPlotter::update(size_t nPts, size_t nLines, const double* x, const double* y, double xmin, double xmax, double ymin, double ymax){
		std::vector<double> xv(x, x+nPts), yv;

		gp << "set xrange [" << xmin << ":" << xmax << "]\n";
		gp << "set yrange [" << ymin << ":" << ymax << "]\n";

		gp << "plot ";
		for (size_t i = 0; i < nLines; i++)
			gp << "'-' with lines title '" << i << (i < nLines - 1 ? "'," : "'\n");
		gp.flush();
		for (size_t i = 0; i < nLines; i++){
			yv = std::vector<double>(y + i*nPts, y + (i + 1)*nPts);
			gp.send1d(boost::make_tuple(xv, yv));
		}
	}

}
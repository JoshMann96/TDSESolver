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



	Orthonormalizer::Orthonormalizer(size_t nPts, size_t nVecs) : nPts(static_cast<lapack_int>(nPts)), nVecs(static_cast<lapack_int>(nVecs)) {
		assert(nPts * nVecs < LAPACK_INT_MAX); // ensure we don't overflow lapack_int
		tau = (double*) sq_malloc(nPts * sizeof(double));
		getlwork(this->nPts, this->nVecs, &lwork);
		work = (double*) sq_malloc(lwork * sizeof(double));
	}

	Orthonormalizer::~Orthonormalizer() {
		if (work) sq_free(work);
		if (tau) sq_free(tau);
	}

	void Orthonormalizer::getlwork(lapack_int nPts, lapack_int nVecs, lapack_int* lwork) {
		double workSize;
		lapack_int info;
		*lwork = -1;
		LAPACK_dgeqrf(&nPts, &nVecs, nullptr, &nPts, nullptr, &workSize, lwork, &info);
		if (info != 0) {
			throw std::runtime_error("Error in LAPACK_dgeqrf (query): " + std::to_string(info));
		}
		if (workSize < 1) workSize = 1;
		*lwork = static_cast<lapack_int>(workSize);
	}

	void Orthonormalizer::orthonormalize(size_t nPts, size_t nVecs, double* __restrict vecs) {
		assert(nPts * nVecs < LAPACK_INT_MAX); // ensure we don't overflow lapack_int

		lapack_int nPts_int = static_cast<lapack_int>(nPts);
		lapack_int nVecs_int = static_cast<lapack_int>(nVecs);
		lapack_int lwork;
		getlwork(nPts_int, nVecs_int, &lwork);
		
		double* work = (double*) sq_malloc(lwork * sizeof(double));
		double* tau = (double*) sq_malloc(nPts * sizeof(double));

		orthonormalize(nPts_int, nVecs_int, vecs, tau, work, lwork);

		sq_free(work);
		sq_free(tau);
	}

	void Orthonormalizer::orthonormalize(lapack_int nPts, lapack_int nVecs, double* __restrict vecs, double* __restrict tau, double* __restrict work, lapack_int lwork){
		lapack_int info;
		LAPACK_dgeqrf(&nPts, &nVecs, vecs, &nPts, tau, work, &lwork, &info);
		if (info != 0) {
			sq_free(work);
			throw std::runtime_error("Error in LAPACK_dgeqrf: " + std::to_string(info));
		}
		LAPACK_dorgqr(&nPts, &nVecs, &nVecs, vecs, &nPts, tau, work, &lwork, &info);
		if (info != 0) {
			sq_free(work);
			throw std::runtime_error("Error in LAPACK_dorgqr: " + std::to_string(info));
		}
	}


	PolynomialExtrapolator::PolynomialExtrapolator(size_t nPts, size_t order, double stepFraction, const double* __restrict initialVector) : nPts(nPts), order(order), historyIndex(0, order) {
		history = (double*) sq_malloc(sizeof(double) * nPts * order);
		extrapStenc = (double*) sq_malloc(sizeof(double) * order * order);

		// initialize history
		if (initialVector)
			for (size_t i = 0; i < order; i++)
				cblas_dcopy(nPts, initialVector, 1, history + i, order);
		else
			std::fill_n(history, nPts * order, 0.0);

		// calculate extrapolation vector
		// fill with original matrix
		for (size_t m = 0; m < order; m++) {
			for (size_t q = 0; q < order; q++) {
				extrapStenc[m*order + q] = std::pow((m+stepFraction), q);
			}
		}

		// invert
		lapack_int info;
		lapack_int n = static_cast<lapack_int>(order);
		lapack_int* ipiv = (lapack_int*) sq_malloc(sizeof(lapack_int) * order);
		LAPACK_dgetrf(&n, &n, extrapStenc, &n, ipiv, &info);
		if(info != 0) {
			sq_free(ipiv);
			sq_free(extrapStenc);
			sq_free(history);
			throw std::runtime_error("PolynomialExtrapolator: dgetrf failed with info = " + std::to_string(info));
		}
		double* work = (double*) sq_malloc(sizeof(double) * order);
		LAPACK_dgetri(&n, extrapStenc, &n, ipiv, work, &n, &info);
		if(info != 0) {
			sq_free(work);
			sq_free(ipiv);
			sq_free(extrapStenc);
			sq_free(history);
			throw std::runtime_error("PolynomialExtrapolator: dgetri failed with info = " + std::to_string(info));
		}
		sq_free(work);
		sq_free(ipiv);

		// fill the matrix with permutations of the first row
		for (size_t m = 1; m < order; m++) {
			for (size_t lm = 0; lm < order; lm++) {
				extrapStenc[m*order + ((lm + order - m) % order)] = extrapStenc[lm];
			}
		}
	}

	void PolynomialExtrapolator::pushHistory(const double* __restrict vec) {
		cblas_dcopy(nPts, vec, 1, history + order - 1 - historyIndex, order);
		historyIndex++;
	}

	void PolynomialExtrapolator::extrapolate(double* __restrict targ) {
		for (size_t i = 0; i < nPts; i++)
			targ[i] = cblas_ddot(order, extrapStenc + historyIndex * order, 1, history + i * order, 1);
	}

	void PolynomialExtrapolator::printExtrapStenc() const {
		for (size_t m = 0; m < order; m++)
			vtlsPrnt::printArray(order, extrapStenc + m * order);
	}

	void PolynomialExtrapolator::printHistory() const {
		for (size_t m = 0; m < nPts; m++)
			vtlsPrnt::printArray(order, history + m * order);
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
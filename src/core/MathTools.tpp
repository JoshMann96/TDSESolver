#include <stddef.h>
namespace vtls{

	template <typename T>
	void copyArrayRe(size_t len, const T* __restrict arr1, double* __restrict arr2){
		for (size_t i = 0; i < len; i++)
			arr2[i] = std::real(arr1[i]);
	}
	

	template <class T>
	Orthonormalizer<T>::Orthonormalizer(size_t nPts, size_t nVecs) : nPts(static_cast<lapack_int>(nPts)), nVecs(static_cast<lapack_int>(nVecs)) {
		static_assert(std::is_same_v<T, double> || std::is_same_v<T, lapack_complex_double>, "Orthonormalizer only supports double or lapack_complex_double types");
		assert(nPts * nVecs < LAPACK_INT_MAX); // ensure we don't overflow lapack_int
		tau = (T*) sq_malloc(nPts * sizeof(T));
		getlwork(this->nPts, this->nVecs, &lwork);
		work = (T*) sq_malloc(lwork * sizeof(T));
	}

	template <class T>
	Orthonormalizer<T>::~Orthonormalizer() {
		if (work) sq_free(work);
		if (tau) sq_free(tau);
	}

	template <class T>
	void Orthonormalizer<T>::getlwork(lapack_int nPts, lapack_int nVecs, lapack_int* lwork) {
		double workSize;
		lapack_int info;
		*lwork = -1;

		// select appropriate routine based on type T
		if constexpr (std::is_same_v<T, double>){
			LAPACK_dgeqrf(&nPts, &nVecs, nullptr, &nPts, nullptr, &workSize, lwork, &info);
			if (workSize < 1) workSize = 1;
		}
		else if constexpr (std::is_same_v<T, lapack_complex_double>){
			std::complex<double> workSize_c;
			LAPACK_zgeqrf(&nPts, &nVecs, nullptr, &nPts, nullptr, reinterpret_cast<lapack_complex_double*>(&workSize_c), lwork, &info);
			if (workSize_c.real() < 1) workSize_c = 1;
			workSize = workSize_c.real();
		}
		else
			throw std::runtime_error("Unsupported type for Orthonormalizer");

		if (info != 0) {
			throw std::runtime_error("Error in LAPACK_?geqrf (query): " + std::to_string(info));
		}
		*lwork = static_cast<lapack_int>(workSize);
	}

	template <class T>
	void Orthonormalizer<T>::orthonormalize(size_t nPts, size_t nVecs, T* __restrict vecs) {
		assert(nPts * nVecs < LAPACK_INT_MAX); // ensure we don't overflow lapack_int

		lapack_int nPts_int = static_cast<lapack_int>(nPts);
		lapack_int nVecs_int = static_cast<lapack_int>(nVecs);
		lapack_int lwork;
		getlwork(nPts_int, nVecs_int, &lwork);
		
		T* work = (T*) sq_malloc(lwork * sizeof(T));
		T* tau = (T*) sq_malloc(nPts * sizeof(T));

		orthonormalize(nPts_int, nVecs_int, vecs, tau, work, lwork);

		sq_free(work);
		sq_free(tau);
	}

	template <class T>
	void Orthonormalizer<T>::orthonormalize(lapack_int nPts, lapack_int nVecs, T* __restrict vecs, T* __restrict tau, T* __restrict work, lapack_int lwork){
		lapack_int info;
		
		if constexpr (std::is_same_v<T, double>)
			LAPACK_dgeqrf(&nPts, &nVecs, vecs, &nPts, tau, work, &lwork, &info);
		else if constexpr (std::is_same_v<T, lapack_complex_double>)
			LAPACK_zgeqrf(&nPts, &nVecs, vecs, &nPts, tau, work, &lwork, &info);
		else
			throw std::runtime_error("Unsupported type for Orthonormalizer");
		if (info != 0) {
			sq_free(work);
			throw std::runtime_error("Error in LAPACK_?geqrf: " + std::to_string(info));
		}

		if constexpr (std::is_same_v<T, double>)
			LAPACK_dorgqr(&nPts, &nVecs, &nVecs, vecs, &nPts, tau, work, &lwork, &info);
		else if constexpr (std::is_same_v<T, lapack_complex_double>)
			LAPACK_zungqr(&nPts, &nVecs, &nVecs, vecs, &nPts, tau, work, &lwork, &info);
		else
			throw std::runtime_error("Unsupported type for Orthonormalizer");
		if (info != 0) {
			sq_free(work);
			throw std::runtime_error("Error in LAPACK_?orgqr: " + std::to_string(info));
		}
	}

    template<class T>
	Convolver<T>::Convolver(size_t len) : len(len){
		temp1 = reinterpret_cast<std::complex<double>*>(sq_malloc(sizeof(fftw_complex)*len));
		temp2 = reinterpret_cast<std::complex<double>*>(sq_malloc(sizeof(fftw_complex)*len));

		mtx.lock();

		fftw_plan_with_nthreads(omp_get_max_threads());
		//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

		assert(len <= INT_MAX);
		int lenInt = static_cast<int>(len);
		fp = fftw_plan_dft(1, &lenInt, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1), FFTW_FORWARD, FFTW_ESTIMATE);
		bp = fftw_plan_dft(1, &lenInt, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2), FFTW_BACKWARD, FFTW_ESTIMATE);

		mtx.unlock();
	}

	template<class T>
	Convolver<T>::~Convolver(){
		sq_free(temp1);
		sq_free(temp2);
		mtx.lock();
		fftw_destroy_plan(fp);
		fftw_destroy_plan(bp);
		mtx.unlock();
	}

	template<class T>
	void Convolver<T>::compute(const T* arr1, const T* arr2, T* targ){
		if (len < 1)
			throw std::runtime_error("Convolver: length not set.");

		vtls::copyArray(len, arr1, temp1);
		vtls::copyArray(len, arr2, temp2);

		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1));
		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2));

		vtls::seqMulArrays(len, temp1, temp2);

		fftw_execute_dft(bp, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2));

		vtls::scaMulArrayRe(len, 1.0/len, temp2, targ);
	}


	template<class T>
	MaskConvolver<T>::MaskConvolver(size_t len, const T* maskIn) : len(len){
		mask = reinterpret_cast<std::complex<double>*>(sq_malloc(sizeof(fftw_complex)*len));
		temp = reinterpret_cast<std::complex<double>*>(sq_malloc(sizeof(fftw_complex)*len));

		mtx.lock();

		fftw_plan_with_nthreads(omp_get_max_threads());
		//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

		assert(len <= INT_MAX);
		int lenInt = static_cast<int>(len);
		fp = fftw_plan_dft(1, &lenInt, reinterpret_cast<fftw_complex*>(mask), reinterpret_cast<fftw_complex*>(mask), FFTW_FORWARD, FFTW_ESTIMATE);
		bp = fftw_plan_dft(1, &lenInt, reinterpret_cast<fftw_complex*>(temp), reinterpret_cast<fftw_complex*>(temp), FFTW_BACKWARD, FFTW_ESTIMATE);

		mtx.unlock();

		vtls::copyArray(len, maskIn, mask);
		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(mask), reinterpret_cast<fftw_complex*>(mask));
	}

	template<class T>
	MaskConvolver<T>::~MaskConvolver(){
		sq_free(mask);
		sq_free(temp);
		mtx.lock();
		fftw_destroy_plan(fp);
		fftw_destroy_plan(bp);
		mtx.unlock();
	}

	template<class T>
	void MaskConvolver<T>::compute(const T* arr, T* targ){
		vtls::copyArray(len, arr, temp);

		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp), reinterpret_cast<fftw_complex*>(temp));

		vtls::seqMulArrays(len, mask, temp);

		fftw_execute_dft(bp, reinterpret_cast<fftw_complex*>(temp), reinterpret_cast<fftw_complex*>(temp));

		vtls::scaMulArrayRe(len, 1.0/len, temp, targ);
	}

	template<class T>
	void MaskConvolver<T>::compute(T* arr){
		vtls::copyArray(len, arr, temp);

		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp), reinterpret_cast<fftw_complex*>(temp));

		vtls::seqMulArrays(len, mask, temp);

		fftw_execute_dft(bp, reinterpret_cast<fftw_complex*>(temp), reinterpret_cast<fftw_complex*>(temp));

		vtls::scaMulArrayRe(len, 1.0/len, temp, arr);
	}
}
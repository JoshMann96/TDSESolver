#include <stddef.h>
namespace vtls{

	template <typename T>
	void copyArrayRe(size_t len, const T* __restrict arr1, double* __restrict arr2){
		for (size_t i = 0; i < len; i++)
			arr2[i] = std::real(arr1[i]);
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
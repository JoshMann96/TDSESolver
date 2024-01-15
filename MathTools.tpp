namespace vtls{
    template<class T>
	Convolver<T>::Convolver(int len) : len(len){
		temp1 = reinterpret_cast<std::complex<double>*>(fftw_malloc(sizeof(fftw_complex)*len));
		temp2 = reinterpret_cast<std::complex<double>*>(fftw_malloc(sizeof(fftw_complex)*len));

		fftw_plan_with_nthreads(omp_get_max_threads());
		//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

		fp = fftw_plan_dft(1, &len, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1), FFTW_FORWARD, FFTW_PATIENT);
		bp = fftw_plan_dft(1, &len, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2), FFTW_BACKWARD, FFTW_PATIENT);
	}

	template<class T>
	Convolver<T>::~Convolver(){
		fftw_free(temp1);
		fftw_free(temp2);
		fftw_destroy_plan(fp);
		fftw_destroy_plan(bp);
	}

	template<class T>
	void Convolver<T>::compute(T* arr1, T* arr2, T* targ){
		vtls::copyArray(len, arr1, temp1);
		vtls::copyArray(len, arr2, temp2);

		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1));
		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2));

		vtls::seqMulArrays(len, temp1, temp2);

		fftw_execute_dft(bp, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2));

		vtls::scaMulArrayRe(len, 1.0/len, temp2, targ);
	}


	template<class T>
	MaskConvolver<T>::MaskConvolver(int len, T* constArr) : len(len){
		temp1 = reinterpret_cast<std::complex<double>*>(fftw_malloc(sizeof(fftw_complex)*len));
		temp2 = reinterpret_cast<std::complex<double>*>(fftw_malloc(sizeof(fftw_complex)*len));

		fftw_plan_with_nthreads(omp_get_max_threads());
		//std::cout << "Assigned FFTW threads: " << fftw_planner_nthreads() << std:: endl;

		fp = fftw_plan_dft(1, &len, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1), FFTW_FORWARD, FFTW_PATIENT);
		bp = fftw_plan_dft(1, &len, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2), FFTW_BACKWARD, FFTW_PATIENT);

		vtls::copyArray(len, constArr, temp1);
		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp1), reinterpret_cast<fftw_complex*>(temp1));
	}

	template<class T>
	MaskConvolver<T>::~MaskConvolver(){
		fftw_free(temp1);
		fftw_free(temp2);
		fftw_destroy_plan(fp);
		fftw_destroy_plan(bp);
	}

	template<class T>
	void MaskConvolver<T>::compute(T* arr, T* targ){
		vtls::copyArray(len, arr, temp2);

		fftw_execute_dft(fp, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2));

		vtls::seqMulArrays(len, temp1, temp2);

		fftw_execute_dft(bp, reinterpret_cast<fftw_complex*>(temp2), reinterpret_cast<fftw_complex*>(temp2));

		vtls::scaMulArrayRe(len, 1.0/len, temp2, targ);
	}
}
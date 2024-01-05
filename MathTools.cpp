#include "stdafx.h"


namespace vtls {

	template<typename T>
	Convolver<T>::Convolver(int len) : len(len){
		temp1 = fftw_malloc(sizeof(fftw_complex)*len);
		temp2 = fftw_malloc(sizeof(fftw_complex)*len);
		fp = fftw_plan_dft(1, &len, temp1, temp1, FFTW_FORWARD, FFTW_PATIENT);
		bp = fftw_plan_dft(1, &len, temp2, temp2, FFTW_BACKWARD, FFTW_PATIENT);
	}

	template<typename T>
	Convolver<T>::~Convolver(){
		fftw_free(temp1);
		fftw_free(temp2);
		fftw_destroy_plan(fp);
		fftw_destroy_plan(bp);
	}

	template<typename T>
	void Convolver<T>::compute(T* arr1, T* arr2, T* targ){
		vtls::copyArray(len, arr1, temp1);
		vtls::copyArray(len, arr2, temp2);

		fftw_execute_dft(fp, temp1, temp1);
		fftw_execute_dft(fp, temp2, temp2);

		vtls::seqMulArrays(len, temp1, temp2);

		fftw_execute_dft(bp, temp2, temp2);

		vtls::scaMulArrayRe(len, 1.0/len, temp2, targ);
	}



	template<typename T>
	MaskConvolver<T>::MaskConvolver(int len, T* constArr) : len(len){
		temp1 = fftw_malloc(sizeof(fftw_complex)*len);
		temp2 = fftw_malloc(sizeof(fftw_complex)*len);
		fp = fftw_plan_dft(1, &len, temp1, temp1, FFTW_FORWARD, FFTW_PATIENT);
		bp = fftw_plan_dft(1, &len, temp2, temp2, FFTW_BACKWARD, FFTW_PATIENT);

		vtls::copyArray(len, constArr, temp1);
		fftw_execute_dft(fp, temp1, temp1);
	}

	template<typename T>
	MaskConvolver<T>::~MaskConvolver(){
		fftw_free(temp1);
		fftw_free(temp2);
		fftw_destroy_plan(fp);
		fftw_destroy_plan(bp);
	}

	template<typename T>
	void MaskConvolver<T>::compute(T* arr, T* targ){
		vtls::copyArray(len, arr, temp2);

		fftw_execute_dft(fp, temp2, temp2);

		vtls::seqMulArrays(len, temp1, temp2);

		fftw_execute_dft(bp, temp2, temp2);

		vtls::scaMulArrayRe(len, 1.0/len, temp2, targ);
	}



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
		int * nArr = new int[n];
		int k = n / 50 + 1;
		for (int i = 0; i < n; i += k) nArr[i] = (int)(((arr[i] - minVal) / (maxVal - minVal)) * 100.0);
		//printArray(n, nArr);
		for (int i = 0; i < n; i += k) {
			for (int j = 0; j < nArr[i]; j++)
				std::cout << "#";
			std::cout << std::endl;
		}
	}

	void printGraph(int n, std::complex<double>* __restrict arr0) {
		double* arr = new double[n];
		for (int i = 0; i < n; i++)
			arr[i] = std::real(arr0[i]);
		double minVal = *std::min_element(arr, arr + n);
		double maxVal = *std::max_element(arr, arr + n);
		int* nArr = new int[n];
		int k = n / 50 + 1;
		for (int i = 0; i < n; i += k) nArr[i] = (int)(((arr[i] - minVal) / (maxVal - minVal)) * 100.0);
		//printArray(n, nArr);
		for (int i = 0; i < n; i += k) {
			for (int j = 0; j < nArr[i]; j++)
				std::cout << "#";
			std::cout << std::endl;
		}
	}
}
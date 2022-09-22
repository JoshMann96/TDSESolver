#include "stdafx.h"


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

	int findValue(int len, double *__restrict arr, double val) {
		for (int i = 0; i < len; i++)
			if (arr[i] >= val)
				return i;
		if (val > arr[len - 1])
			return len - 1;
		else
			return -1;
	}

	void mkl_ddcon(double h[], int inch, double x[], int incx, double y[], int incy, int nh, int nx, int iy0, int ny, int id) {
		int status = VSL_STATUS_OK, error;
		VSLConvTaskPtr task, *task_ptr = &task;

		vsldConvNewTask1D(task_ptr, VSL_CONV_MODE_DIRECT, nh, nx, ny);
		vslConvSetStart(task, &iy0);
		vslConvSetDecimation(task, &id);
		status = vsldConvExec1D(task, h, inch, x, incx, y, incy);

		error = vslConvDeleteTask(task_ptr);

		if (status != VSL_STATUS_OK) {
			printf("ERROR: ddcon(): bad status=%d\n", status);
			exit(1);
		}
		if (error != 0) {
			printf("ERROR: ddcon(): failed to destroy the task descriptor\n");
			exit(1);
		}
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
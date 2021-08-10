#include "stdafx.h"


namespace vtls {

	/*void addArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2, std::complex<double> *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr1[i] + arr2[i];
	}*/

	/*
	void addArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2targ) {
		for (int i = 0; i < len; i++)
			arr2targ[i] += arr1[i];
	}

	void addArrays(int len, double *__restrict arr1, double *__restrict arr2targ) {
		for (int i = 0; i < len; i++)
			arr2targ[i] += arr1[i];
	}*/

	/*void seqMulArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2, std::complex<double> *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr1[i] * arr2[i];
	}

	void seqMulArrays(int len, double *__restrict arr1, std::complex<double> *__restrict arr2, std::complex<double> *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr1[i] * arr2[i];
	}

	void seqMulArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2targ) {
		for (int i = 0; i < len; i++)
			arr2targ[i] *= arr1[i];
	}

	void seqMulArrays(int len, double *__restrict arr1, double *__restrict arr2targ) {
		for (int i = 0; i < len; i++)
			arr2targ[i] *= arr1[i];
	}

	void seqMulArrays(int len, double* __restrict arr1, double* __restrict arr2, double*__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr1[i] * arr2[i];
	}*/

	/*void scaMulArray(int len, std::complex<double> scalar, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr[i] * scalar;
	}*/

	void scaMulArray(int len, std::complex<double> scalar, std::complex<double> *__restrict arr) {
		/*for (int i = 0; i < len; i++)
			arr[i] *= scalar;*/
		cblas_zscal(len, &scalar, arr, 1);
	}

	/*void scaMulArray(int len, double scalar, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr[i] * scalar;
	}

	void scaMulArray(int len, double scalar, double *__restrict arr, double *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr[i] * scalar;
	}*/

	void scaMulArray(int len, double scalar, double* __restrict arr) {
		cblas_dscal(len, scalar, arr, 1);
	}

	void scaMulArray(int len, double scalar, std::complex<double> *__restrict arr) {
		/*for (int i = 0; i < len; i++)
			arr[i] *= scalar;*/
		cblas_zdscal(len, scalar, arr, 1);
	}

	/*void scaMulArrayRe(int len, std::complex<double> scalar, std::complex<double> *__restrict arr, double *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = std::real(arr[i] * scalar);
	}*/

	/*void normSqr(int len, std::complex<double> *__restrict arr, double *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = std::norm(arr[i]);
	}*/

	/*void abs(int len, std::complex<double> *__restrict arr, double *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = std::abs(arr[i]);
	}*/

	/*void normalizeSqrNorm(int len, std::complex<double> *__restrict arr, double dx) {
		double *  tarr = new double[len];
		normSqr(len, arr, tarr);
		scaMulArray(len, 1.0 / std::sqrt(vtlsInt::simps(len, tarr, dx)), arr);
		delete tarr;
	}*/

	/*void downSampleLinearInterpolateNoEdge(int len, double *__restrict arr, int newLen, double *__restrict targ) {
		double step = (double)(len - 1) / newLen;
		double curPos = step / 2.0;
		for (int i = 0; i < newLen; i++) {
			targ[i] = (1.0 - std::fmod(curPos, 1))*arr[(int)curPos] + (fmod(curPos, 1))*arr[(int)curPos + 1 * ((int)curPos != (len - 1))];
			curPos += step;
		}
	}

	void downSampleLinearInterpolateEdge(int len, double *__restrict arr, int newLen, double *__restrict targ) {
		double step = (double)(len - 1) / (newLen - 1);
		double curPos = 0.0;
		for (int i = 0; i < newLen; i++) {
			targ[i] = (1.0 - std::fmod(curPos, 1))*arr[(int)curPos] + (fmod(curPos, 1))*arr[(int)curPos + 1 * ((int)curPos != (len - 1))];
			curPos += step;
		}
	}

	void downSampleLinearInterpolateEdge(int len, std::complex<double> *__restrict arr, int newLen, std::complex<double> *__restrict targ) {
		double step = (double)(len - 1) / (newLen - 1);
		double curPos = 0.0;
		for (int i = 0; i < newLen; i++) {
			targ[i] = (1.0 - std::fmod(curPos, 1))*arr[(int)curPos] + (fmod(curPos, 1))*arr[(int)curPos + 1 * ((int)curPos != (len - 1))];
			curPos += step;
		}
	}*/

	/*void linearInterpolate(int l1, double *__restrict x1, double *__restrict y1, int l2, double *__restrict x2, double *__restrict y2) {
		int curPos = 0;
		for (int i = 0; i < l2; i++) {
			while (x2[i] > x1[curPos] && curPos < l1)
				curPos++;
			if (curPos > 0 && curPos < l1)
				y2[i] = y1[curPos - 1] + (y1[curPos] - y1[curPos - 1])*(x2[i] - x1[curPos - 1]) / (x1[curPos] - x1[curPos - 1]);
			else if (curPos < 1)
				y2[i] = y1[0];
			else
				y2[i] = y1[l1 - 1];
		}
	}*/

	/*void linspace(int len, double min, double max, double *__restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = (max - min) * i / (double)(len - 1) + min;
	}

	double* linspace(int len, double min, double max) {
		double * ret = new double[len];
		for (int i = 0; i < len; i++)
			ret[i] = (max - min) * i / (double)(len - 1) + min;
		return ret;
	}*/

	/*void scaAddArray(int len, double scalar, double *__restrict arr) {
		for (int i = 0; i < len; i++)
			arr[i] += scalar;
	}

	void scaAddArray(int len, std::complex<double> scalar, std::complex<double> *__restrict arr) {
		for (int i = 0; i < len; i++)
			arr[i] += scalar;
	}*/

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

	/*std::complex<double> firstDerivative(int len, std::complex<double> *__restrict arr, int pos, double dx) {
		if (pos == 0)
			return (arr[1] - arr[0]) / dx;
		else if (pos == len - 1)
			return (arr[pos] - arr[pos-1]) / dx;
		else
			return (arr[pos + 1] - arr[pos - 1]) / (2.0*dx);
	}*/

	/*void firstDerivative(int len, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ, double dx) {
		targ[0] = (arr[1] - arr[0]) / dx;
		for (int i = 1; i < len - 1; i++)
			targ[i] = (arr[i + 1] - arr[i - 1]) / (2.0*dx);
		targ[len - 1] = (arr[len - 1] - arr[len - 2]) / dx;
	}

	void firstDerivative(int len, double *__restrict arr, double *__restrict targ, double dx) {
		targ[0] = (arr[1] - arr[0]) / dx;
		for (int i = 1; i < len - 1; i++)
			targ[i] = (arr[i + 1] - arr[i - 1]) / (2.0*dx);
		targ[len - 1] = (arr[len - 1] - arr[len - 2]) / dx;
	}*/

	/*void secondDerivative(int len, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ, double dx) {
		double dx2 = dx * dx;
		targ[0] = (arr[1] - arr[0])*2.0 / dx2;
		for (int i = 1; i < len - 1; i++)
			targ[i] = (arr[i - 1] - 2.0*arr[i] + arr[i + 1]) / dx2;
		targ[len - 1] = (arr[len - 2] - arr[len - 1])*2.0 / dx2;
	}*/

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



namespace vtlsInt {
	/*double trapz(int len, double *__restrict arr, double dx) {
		if (len <= 1)
			return 0;

		double sum = (arr[0] + arr[len - 1]) / 2.0;
		for (int i = 1; i < len - 1; i++) {
			sum += arr[i];
		}
		return sum * dx;
	}

	std::complex<double> trapz(int len, std::complex<double> *__restrict arr, double dx) {
		if (len <= 1)
			return 0;

		std::complex<double> sum = (arr[0] + arr[len - 1]) / 2.0;
		for (int i = 1; i < len - 1; i++)
			sum += arr[i];
		return sum * dx;
	}*/

	/*std::complex<double> trapzMul(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2, double dx) {
		if (len <= 1)
			return 0;

		std::complex<double> sum = (arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) / 2.0;
		for (int i = 1; i < len - 1; i++)
			sum += arr1[i] * arr2[i];
		return sum * dx;
	}

	std::complex<double> trapzMul(int len, std::complex<double> *__restrict arr1, double *__restrict arr2, double dx) {
		if (len <= 1)
			return 0;

		std::complex<double> sum = (arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) / 2.0;
		for (int i = 1; i < len - 1; i++)
			sum += arr1[i] * arr2[i];
		return sum * dx;
	}*/

	/*double simps(int len, double *__restrict arr, double dx) {
		if (len % 2) {
			double sum = (arr[len - 1] + arr[0]) / 4.0;
			for (int i = 1; i < len - 1; i += 2)
				sum += arr[i];
			sum *= 2.0;
			for (int i = 2; i < len - 1; i += 2)
				sum += arr[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			double sum = (5.0*(arr[0] + arr[len - 1]) + 13.0*(arr[1] + arr[len - 2])) / 12.0;
			for (int i = 2; i < len - 2; i++)
				sum += arr[i];
			return sum * dx;
		}
	}

	std::complex<double> simps(int len, std::complex<double> *__restrict arr, double dx) {
		if (len % 2) {
			std::complex<double> sum = (arr[len - 1] + arr[0]) / 4.0;
			for (int i = 1; i < len - 1; i += 2)
				sum += arr[i];
			sum *= 2.0;
			for (int i = 2; i < len - 1; i += 2)
				sum += arr[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			std::complex<double> sum =
				(
					5.0*(arr[0] + arr[len - 1]) +
					13.0*(arr[1] + arr[len - 2])
					) / 12.0;
			for (int i = 2; i < len - 2; i++)
				sum += arr[i];
			return sum * dx;
		}
	}*/

	/*std::complex<double> simpsMul(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2, double dx) {
		if (len % 2) {
			std::complex<double> sum = (arr1[len - 1] * arr2[len - 1] + arr1[0] * arr2[0]) / 4.0;
			for (int i = 1; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			sum *= 2.0;
			for (int i = 2; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			std::complex<double> sum =
				(
					5.0*(arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) +
					13.0*(arr1[1] * arr2[1] + arr1[len - 2] * arr2[len - 2])
					) / 12.0;
			for (int i = 2; i < len - 2; i++)
				sum += arr1[i] * arr2[i];
			return sum * dx;
		}
	}

	double simpsMul(int len, double *__restrict arr1, double *__restrict arr2, double dx) {
		if (len % 2) {
			double sum = (arr1[len - 1] * arr2[len - 1] + arr1[0] * arr2[0]) / 4.0;
			for (int i = 1; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			sum *= 2.0;
			for (int i = 2; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			double sum =
				(
					5.0*(arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) +
					13.0*(arr1[1] * arr2[1] + arr1[len - 2] * arr2[len - 2])
					) / 12.0;
			for (int i = 2; i < len - 2; i++)
				sum += arr1[i] * arr2[i];
			return sum * dx;
		}
	}*/

	/*void cumIntRectLeft(int len, double *__restrict arr, double dx, double *__restrict targ) {
		double sum = 0.0;
		for (int i = 0; i < len; i++) {
			sum += arr[i] * dx;
			targ[i] = sum;
		}
	}

	void cumIntRectRight(int len, double *__restrict arr, double dx, double *__restrict targ) {
		double sum = 0.0;
		for (int i = 0; i < len; i++) {
			targ[i] = sum;
			sum += arr[i] * dx;
		}
	}*/

	/*void cumIntTrapz(int len, double *__restrict arr, double dx, double *__restrict targ) {
		double sum = 0.0;
		for (int i = 0; i < len - 1; i++) {
			targ[i] = sum;
			sum += (arr[i] + arr[i + 1])*(dx / 2.0);
		}
		targ[len - 1] = sum;
	}

	void cumIntTrapz(int len, std::complex<double> *__restrict arr, double dx, std::complex<double> *__restrict targ) {
		std::complex<double> sum = 0.0;
		for (int i = 0; i < len - 1; i++) {
			targ[i] = sum;
			sum += (arr[i] + arr[i + 1])*(dx / 2.0);
		}
		targ[len - 1] = sum;
	}*/
}



namespace vtlsPrnt {
	/*void printArray(int n, double *__restrict arr) {
		std::cout << "[";
		for (int i = 0; i < n; i++) {
			std::cout << arr[i];
			if (i != n - 1) {
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
	}

	void printArray(int n, int *__restrict arr) {
		std::cout << "[";
		for (int i = 0; i < n; i++) {
			std::cout << arr[i];
			if (i != n - 1) {
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
	}

	void printArray(int n, std::complex<double> *__restrict arr) {
		std::cout << "[";
		for (int i = 0; i < n; i++) {
			std::cout << arr[i];
			if (i != n - 1) {
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
	}*/

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
}
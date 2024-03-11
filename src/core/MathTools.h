#pragma once
#include "CORECommonHeader.h"
#include <cblas.h>
#include <omp.h>
#include <fftw3.h>
#include "exprtk.hpp"

namespace vtlsInt {
	// Riemann integration.

	template <typename T, typename U>
	decltype(std::declval<T&>()* std::declval<U&>()) rSum(int len, T* __restrict arr, U dx) {
		if (len <= 1)
			return 0;

		T sum = 0;
		for (int i = 0; i < len; i++) {
			sum += arr[i];
		}
		return sum * dx;
	}

	template <typename T, typename U, typename V>
	decltype(std::declval<T&>()* std::declval<U&>()* std::declval<V&>()) rSumMul(int len, T* __restrict arr1, U* __restrict arr2, V dx) {
		if (len <= 1)
			return 0;

		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0;
		for (int i = 0; i < len; i++)
			sum += arr1[i] * arr2[i];
		return sum * dx;
	}

	// Trapezoidal integration.
	template <typename T, typename U>
	decltype(std::declval<T&>() * std::declval<U&>()) trapz(int len, T* __restrict arr, U dx) {
		if (len <= 1)
			return 0;

		T sum = (arr[0] + arr[len - 1]) / 2.0;
		for (int i = 1; i < len - 1; i++) {
			sum += arr[i];
		}
		return sum * dx;
	}

	// Multiplies the elements of two arrays and then trapezoidally integrates the result.
	template <typename T, typename U, typename V>
	decltype(std::declval<T&>()* std::declval<U&>()* std::declval<V&>()) trapzMul(int len, T* __restrict arr1, U* __restrict arr2, V dx) {
		if (len <= 1)
			return 0;

		decltype(std::declval<T&>() * std::declval<U&>()) sum = (arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) / 2.0;
		for (int i = 1; i < len - 1; i++)
			sum += arr1[i] * arr2[i];
		return sum * dx;
	}

	// Simpson rule integration.
	template <typename T, typename U>
	decltype(std::declval<T&>()* std::declval<U&>()) simps(int len, T* __restrict arr, U dx) {
		if (len % 2) {
			T sum = (arr[len - 1] + arr[0]) / 4.0;
			for (int i = 1; i < len - 1; i += 2)
				sum += arr[i];
			sum *= 2.0;
			for (int i = 2; i < len - 1; i += 2)
				sum += arr[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			T sum = (5.0 * (arr[0] + arr[len - 1]) + 13.0 * (arr[1] + arr[len - 2])) / 12.0;
			for (int i = 2; i < len - 2; i++)
				sum += arr[i];
			return sum * dx;
		}
	}

	template <typename T, typename U, typename V>
	decltype(std::declval<T&>()* std::declval<U&>()* std::declval<V&>()) simpsMul(int len, T* __restrict arr1, U* __restrict arr2, V dx) {
		if (len % 2) {
			decltype(std::declval<T&>() * std::declval<U&>()) sum = (arr1[len - 1] * arr2[len - 1] + arr1[0] * arr2[0]) / 4.0;
			for (int i = 1; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			sum *= 2.0;
			for (int i = 2; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			decltype(std::declval<T&>() * std::declval<U&>()) sum =
				(
					5.0 * (arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) +
					13.0 * (arr1[1] * arr2[1] + arr1[len - 2] * arr2[len - 2])
					) / 12.0;
			for (int i = 2; i < len - 2; i++)
				sum += arr1[i] * arr2[i];
			return sum * dx;
		}
	}
	
	// Cumulative integration using left points as values.
	template <typename T, typename U>
	void cumIntRectLeft(int len, T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0.0;
		for (int i = 0; i < len; i++) {
			sum += arr[i] * dx;
			targ[i] = sum;
		}
	}

	// Cumulative integration using right points as values.
	template <typename T, typename U>
	void cumIntRectRight(int len, T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0.0;
		for (int i = 0; i < len; i++) {
			targ[i] = sum;
			sum += arr[i] * dx;
		}
	}

	// Cumulative integration using trapezoidal rule.
	template <typename T, typename U>
	void cumIntTrapz(int len, T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0.0;
		for (int i = 0; i < len - 1; i++) {
			targ[i] = sum;
			sum += (arr[i] + arr[i + 1]) * (dx / 2.0);
		}
		targ[len - 1] = sum;
	}

	// Trapezoidal integration.
	//double trapz(int len, double* __restrict arr, double dx);
	//std::complex<double> trapz(int len, std::complex<double>* __restrict arr, double dx);
	// Multiplies the elements of two arrays and then trapezoidally integrates the result.
	//std::complex<double> trapzMul(int len, std::complex<double>* __restrict arr1, std::complex<double>* __restrict arr2, double dx);
	//std::complex<double> trapzMul(int len, std::complex<double>* __restrict arr1, double* __restrict arr2, double dx);
	// Simpson rule integration.
	//double simps(int len, double* __restrict arr, double dx);
	//std::complex<double> simps(int len, std::complex<double>* __restrict arr, double dx);
	// Same as trapzMul, but for simps.
	//std::complex<double> simpsMul(int len, std::complex<double>* __restrict arr1, std::complex<double>* __restrict arr2, double dx);
	//double simpsMul(int len, double* __restrict arr1, double* __restrict arr2, double dx);
	// Cumulative integration using left points as values.
	//void cumIntRectLeft(int len, double* __restrict arr, double dx, double* __restrict targ);
	// Cumulative integration using right points as values.
	//void cumIntRectRight(int len, double* __restrict arr, double dx, double* __restrict targ);
	// Cumulative integration using trapezoidal rule.
	//void cumIntTrapz(int len, double* __restrict arr, double dx, double* __restrict targ);
	//void cumIntTrapz(int len, std::complex<double>* __restrict arr, double dx, std::complex<double>* __restrict targ);
};

namespace vtls {
	template <class T>
	class Convolver{
	private:
		int len;
		fftw_plan fp, bp;
		std::complex<double> *temp1, *temp2;
	public:
		Convolver<T>(){}
		Convolver<T>(int len);
		~Convolver();
		void compute(T* arr1, T* arr2, T* targ);
	};

	template <class T>
	class MaskConvolver :
		public Convolver<T>{
	private:
		int len;
		fftw_plan fp, bp;
		std::complex<double> *temp1, *temp2;
	public:
		MaskConvolver<T>(int len, T* constArr);
		~MaskConvolver();
		void compute(T* arr, T* targ);
	};

	// Multiplies triag*diag*triag, with triag being a Hermitian matrix (upper triangular rep, column major) and diag being diagonal
	template <typename T, typename U>
	void mulTriagDiagTriag(int len, T* triag, U* diag, decltype(std::declval<T&>()* std::declval<U&>())* targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) csum = 0.0;
		for (int i = 0; i < len; i++) {
			for (int j = i; j < len; j++) {
				//targ[i,j] = sum_k triag[i,k]diag[k]triag[k,l]
				csum = 0.0;
				for (int k = 0; k < i; k++)
					csum += std::conj(triag[k + (i * (i + 1)) / 2]) * diag[k] * triag[k + (j * (j + 1)) / 2];
				for(int k = i; k < j; k++)
					csum += triag[i + (k * (k + 1)) / 2] * diag[k] * triag[k + (j * (j + 1)) / 2];
				for (int k = j; k < len; k++)
					csum += triag[i + (k * (k + 1)) / 2] * diag[k] * std::conj(triag[j + (k * (k + 1)) / 2]);
				targ[i + (j * (j + 1)) / 2] = csum;
			}
		}
	}

	// Adds two arrays, taking only the imaginary component of the first
	void addArraysImag(int len, std::complex<double>* arr1, double* arr2targ);

	template <typename T, typename U>
	void evalMathExpr(int len, const char* var, T* vals, std::string expr, U* res) {
		typedef exprtk::symbol_table<T> symbol_table_t;
		typedef exprtk::expression<T>   expression_t;
		typedef exprtk::parser<T>       parser_t;

		T cval;
		symbol_table_t symbol_table;
		symbol_table.add_variable(var, cval);
		symbol_table.add_constants();

		expression_t expression;
		expression.register_symbol_table(symbol_table);

		parser_t parser;
		parser.compile(expr, expression);

		for (int i = 0; i < len; i++) {
			cval = vals[i];
			res[i] = (U)(expression.value());
		}
	}

	template <typename T, typename U>
	void polyEval(int len, T* x, int nPoly, U* __restrict polyCoeffs, decltype(std::declval<T&>()* std::declval<U&>())* y) {
		for (int i = 0; i < len; i++)
			y[i] = boost::math::tools::evaluate_polynomial(polyCoeffs, x[i], nPoly);
	}

	// Adds two arrays into a third array
	template <typename T, typename U, typename V>
	void addArrays(int len, T* __restrict arr1, U* __restrict arr2, V* __restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr1[i] + arr2[i];
	}

	// Adds two arrays, storing the result in the second array
	template <typename T, typename U>
	void addArrays(int len, T* __restrict arr1, U* __restrict arr2targ) {
		for (int i = 0; i < len; i++)
			arr2targ[i] += arr1[i];
	}

	// Sequentially multiplies two arrays into a third array
	template <typename T, typename U, typename V>
	void seqMulArrays(int len, const T* __restrict arr1, const U* __restrict arr2, V* __restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr1[i] * arr2[i];
	}
	// Sequentially multiplies two arrays, storing the result in the second array
	template <typename T, typename U>
	void seqMulArrays(int len, T* __restrict arr1, U* __restrict arr2targ) {
		for (int i = 0; i < len; i++)
			arr2targ[i] *= arr1[i];
	}

	// Multiplies an array by a scalar multiple into a second array
	template <typename T, typename U, typename V>
	void scaMulArray(int len, T scalar, U* __restrict arr, V* __restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = arr[i] * scalar;
	}
	// Multiplies an array by a scalar multiple
	void scaMulArray(int len, std::complex<double> scalar, std::complex<double>* __restrict arr);
	void scaMulArray(int len, double scalar, std::complex<double>* __restrict arr);
	void scaMulArray(int len, double scalar, double* __restrict arr);

	// Multiplies an array by a scalar multiple and stores the real component of the result in the second array 
	template <typename T, typename U, typename V>
	void scaMulArrayRe(int len, T scalar, U* __restrict arr, V* __restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = std::real(arr[i] * scalar);
	}

	// Gets the norm squared of an array point-by-point
	template <typename T>
	void normSqr(int len, T* __restrict arr, double* __restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = std::norm(arr[i]);
	}

	// Gets the absolute value of an array
	template <typename T>
	void abs(int len, T* __restrict arr, double* __restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = std::abs(arr[i]);
	}

	// Normalizes an array such that its total norm is 1
	template <typename T>
	void normalizeSqrNorm(int len, T* __restrict arr, double dx) {
		double* tarr = new double[len];
		normSqr(len, arr, tarr);
		scaMulArray(len, 1.0 / std::sqrt(vtlsInt::simps(len, tarr, dx)), arr);
		delete[] tarr;
	}

	// Gets square norm of an array
	template <typename T>
	double getNorm(int len, T* __restrict arr, double dx) {
		double sm = 0.0;
		for (int i = 0; i < len; i++)
			sm += std::pow(std::abs(arr[i]), 2);
		return sm *= dx;
	}

	// Sets square norm of an array
	template <typename T>
	void setNorm(int len, T* __restrict arr, double dx, double norm) {
		scaMulArray(len, std::sqrt(norm / getNorm(len, arr, dx)), arr);
	}

	// Linear interpolation
	template <typename T>
	void downSampleLinearInterpolateEdge(int len, T* __restrict arr, int newLen, T* __restrict targ) {
		double step = (double)(len - 1) / (newLen - 1);
		double curPos = 0.0;
		for (int i = 0; i < newLen; i++) {
			targ[i] = (1.0 - std::fmod(curPos, 1)) * arr[(int)curPos] + (fmod(curPos, 1)) * arr[(int)curPos + 1 * ((int)curPos != (len - 1))];
			curPos += step;
		}
	}

	template <typename T>
	void downSampleLinearInterpolateNoEdge(int len, T* __restrict arr, int newLen, T* __restrict targ) {
		double step = (double)(len - 1) / newLen;
		double curPos = step / 2.0;
		for (int i = 0; i < newLen; i++) {
			targ[i] = (1.0 - std::fmod(curPos, 1)) * arr[(int)curPos] + (fmod(curPos, 1)) * arr[(int)curPos + 1 * ((int)curPos != (len - 1))];
			curPos += step;
		}
	}

	template <typename T>
	void linearInterpolate(int l1, double* __restrict x1, T* __restrict y1, int l2, double* __restrict x2, T* __restrict y2) {
		int curPos = 0;
		for (int i = 0; i < l2; i++) {
			while (x2[i] > x1[curPos] && curPos < l1)
				curPos++;
			if (curPos > 0 && curPos < l1)
				y2[i] = y1[curPos - 1] + (y1[curPos] - y1[curPos - 1]) * (x2[i] - x1[curPos - 1]) / (x1[curPos] - x1[curPos - 1]);
			else if (curPos < 1)
				y2[i] = y1[0];
			else
				y2[i] = y1[l1 - 1];
		}
	}

	template <typename T, typename U>
	U linearInterpolate(int len, U* __restrict arr, T xStart, T dx, T samp) {
		int ix = (int)((samp - xStart) / dx);
		T t = (samp - xStart) / dx - ix;
		if (ix < 0) {
			t += ix;
			ix = 0;
		}
		else if (ix >= len - 1) {
			t += ix - len + 2;
			ix = len - 2;
		}
		return arr[ix] + t * (arr[ix + 1] - arr[ix]);
	}

	template <typename T>
	void linspace(int len, T min, T max, T* __restrict targ) {
		for (int i = 0; i < len; i++)
			targ[i] = (max - min) * i / (T)(len - 1) + min;
	}

	template <typename T>
	T* linspace(int len, T min, T max) {
		T* ret = new T[len];
		for (int i = 0; i < len; i++)
			ret[i] = (max - min) * i / (T)(len - 1) + min;
		return ret;
	}

	// Adds a scalar value to each component of the array
	template <typename T, typename U>
	void scaAddArray(int len, T scalar, U* __restrict arr) {
		for (int i = 0; i < len; i++)
			arr[i] += scalar;
	}

	// Copies the array
	void copyArray(int len, double* __restrict arr1, double* __restrict arr2);
	void copyArray(int len, std::complex<double>* __restrict arr1, std::complex<double>* __restrict arr2);
	void copyArray(int len, double* __restrict arr1, std::complex<double>* __restrict arr2);

	template <typename T>
	void copyArrayRe(int len, T* __restrict arr1, double* __restrict arr2);

	// Takes first derivative at one point
	template <typename T>
	T firstDerivative(int len, T* __restrict arr, int pos, double dx) {
		if (pos == 0)
			return (arr[1] - arr[0]) / dx;
		else if (pos == len - 1)
			return (arr[pos] - arr[pos - 1]) / dx;
		else
			return (arr[pos + 1] - arr[pos - 1]) / (2.0 * dx);
	}

	// Takes first derivative across entire array
	template <typename T>
	void firstDerivative(int len, T* __restrict arr, T* __restrict targ, double dx) {
		targ[0] = (arr[1] - arr[0]) / dx;
		for (int i = 1; i < len - 1; i++)
			targ[i] = (arr[i + 1] - arr[i - 1]) / (2.0 * dx);
		targ[len - 1] = (arr[len - 1] - arr[len - 2]) / dx;
	}

	// Takes second derivative across entire array
	template <typename T>
	void secondDerivative(int len, T* __restrict arr, T* __restrict targ, double dx) {
		double dx2 = dx * dx;
		targ[0] = (arr[1] - arr[0]) * 2.0 / dx2;
		for (int i = 1; i < len - 1; i++)
			targ[i] = (arr[i - 1] - 2.0 * arr[i] + arr[i + 1]) / dx2;
		targ[len - 1] = (arr[len - 2] - arr[len - 1]) * 2.0 / dx2;
	}

	// Finds the location of a value
	int findValue(int len, double* __restrict arr, double val);

	void insertSort_idxs(int len, double* __restrict arr, int* __restrict idxs);

	//void addArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2, std::complex<double> *__restrict targ);
	//void addArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2targ);
	//void addArrays(int len, double *__restrict arr1, double *__restrict arr2targ);
	// Multiplies two arrays together sequentially.
	//void seqMulArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2, std::complex<double> *__restrict targ);
	//void seqMulArrays(int len, double *__restrict arr1, std::complex<double> *__restrict arr2, std::complex<double> *__restrict targ);
	//void seqMulArrays(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2targ);
	//void seqMulArrays(int len, double *__restrict arr1, double *__restrict arr2targ);
	//void seqMulArrays(int len, double* __restrict arr1, double* __restrict arr2, double *__restrict targ);
	// Multiplies an array with a scalar.
	//void scaMulArray(int len, std::complex<double> scalar, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ);
	//void scaMulArray(int len, double scalar, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ);
	//void scaMulArray(int len, double scalar, double *__restrict arr, double *__restrict targ);
	//void scaMulArrayRe(int len, std::complex<double> scalar, std::complex<double> *__restrict arr, double *__restrict targ);
	// Gets the magnitude squared of each element of an array.
	//void normSqr(int len, std::complex<double> *__restrict arr, double *__restrict targ);
	// Gets the absolute value of each element in an array.
	//void abs(int len, std::complex<double> *__restrict arr, double *__restrict targ);
	// Normalizes array (wave function).
	//void normalizeSqrNorm(int len, std::complex<double> *__restrict arr, double dx);
	// Linear interpolation.
	//void downSampleLinearInterpolateEdge(int len, double *__restrict arr, int newLen, double *__restrict targ);
	//void downSampleLinearInterpolateEdge(int len, std::complex<double> *__restrict arr, int newLen, std::complex<double> *__restrict targ);
	//void downSampleLinearInterpolateNoEdge(int len, double *__restrict arr, int newLen, double *__restrict targ);
	// Linear interpolation using new grid (y1,x1 sampled at x2 written to y2). Assumes ordered in increasing order in x1, x2
	//void linearInterpolate(int l1, double *__restrict x1, double *__restrict y1, int l2, double *__restrict x2, double * y2);
	// Linspace
	//void linspace(int len, double min, double max, double *__restrict targ);
	//double* linspace(int len, double min, double max);
	// Adds a scalar to each element of an array.
	//void scaAddArray(int len, double scalar, double *__restrict arr);
	//void scaAddArray(int len, std::complex<double> scalar, std::complex<double> *__restrict arr);
	// Copies the elements of two arrays.
	//void copyArray(int len, double *__restrict arr1, double *__restrict arr2);
	//void copyArray(int len, std::complex<double> *__restrict arr1, std::complex<double> *__restrict arr2);
	//void copyArray(int len, double* __restrict arr1, std::complex<double>* __restrict arr2);
	// Calculates first derivative (midpoint method except at edges).
	//std::complex<double> firstDerivative(int len, std::complex<double> *__restrict arr, int pos, double dx);
	//void firstDerivative(int len, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ, double dx);
	//void firstDerivative(int len, double *__restrict arr, double *__restrict targ, double dx);
	// Calculates second derivative.
	//void secondDerivative(int len, std::complex<double> *__restrict arr, std::complex<double> *__restrict targ, double dx);
	//int findValue(int len, double *__restrict arr, double val);
	// Performs convolution. Recommended to use the contents of this function as a template, as many optimizations can potentially be done for special cases.
	//void mkl_ddcon(double h[], int inch, double x[], int incx, double y[], int incy, int nh, int nx, int iy0, int ny, int id);
	//void insertSort_idxs(int len, double* __restrict arr, int* __restrict idxs);
};

namespace vtlsPrnt {
	// Prints contents of array.
	template <typename T>
	void printArray(int n, T* __restrict arr) {
		std::cout << "[";
		for (int i = 0; i < n; i++) {
			std::cout << arr[i];
			if (i != n - 1) {
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
	}

	// Uses text to graph an array.
	void printGraph(int n, double* __restrict arr);
	void printGraph(int n, std::complex<double>* __restrict arr);

	// Prints contents of array.
	/*void printArray(int n, double *__restrict arr);
	void printArray(int n, int *__restrict arr);
	void printArray(int n, std::complex<double> *__restrict arr);*/
	// Uses text to graph an array.
	//void printGraph(int n, double *__restrict arr);

	template <typename T>
	void saveArray(int n, const char* fil, T* data) {
		std::fstream fid(fil, std::ios::out | std::ios::binary);
		fid.write(reinterpret_cast<char*>(n), sizeof(int));
		fid.write(reinterpret_cast<char*>(data), sizeof(T) * n);
		fid.close();
	}
}

#include "MathTools.tpp"
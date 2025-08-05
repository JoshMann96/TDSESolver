/**
 * @file MathTools.h
 * @brief Mathematical tools for numerical calculations, including integration, convolution, and array manipulation.
 */
#pragma once
#include "CORECommonHeader.h"
#include "blas.h"
#include <omp.h>
#include <mutex>
#include "CyclicInt.h"

/**
 * @namespace vtlsInt
 * @brief Numerical integration methods.
 */
namespace vtlsInt {

	/**
	 * Riemannian integration of an array.
	 * @tparam T The type of the array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to sum.
	 * @param dx The increment (step size).
	 */
	template <typename T, typename U>
	decltype(std::declval<T&>()* std::declval<U&>()) sum(size_t len, const T* __restrict arr, U dx) {
		T sum = 0;
		for (size_t i = 0; i < len; i++) {
			sum += arr[i];
		}
		return sum * dx;
	}

	/**
	 * Riemannian integration of the elementwise product of two arrays (scaled dot product).
	 * @tparam T The type of the first array elements.
	 * @tparam U The type of the second array elements.
	 * @tparam V The type of the increment (dx).
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to sum.
	 * @param arr2 (in) The second array to sum.
	 * @param dx The increment (step size).
	 */
	template <typename T, typename U, typename V>
	decltype(std::declval<T&>()* std::declval<U&>()* std::declval<V&>()) innerProduct(size_t len, const T* __restrict arr1, const U* __restrict arr2, V dx) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0;
		for (size_t i = 0; i < len; i++)
			sum += arr1[i] * arr2[i];
		return sum * dx;
	}

	inline double innerProduct(size_t len, const double* __restrict arr1, const double* __restrict arr2, double dx) {
		return cblas_ddot(len, arr1, 1, arr2, 1) * dx;
	}

	inline std::complex<double> innerProduct(size_t len, const std::complex<double>* __restrict arr1, const std::complex<double>* __restrict arr2, double dx) {
		return cblas_zdotu(len, arr1, 1, arr2, 1) * dx;
	}

	/**
	 * Riemannian integration of the elementwise product of the conjugate of the first array and the second array (dot product).
	 * @tparam T The type of the first array elements.
	 * @tparam U The type of the second array elements.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to sum (to be conjugated).
	 * @param arr2 (in) The second array to sum.
	 * @param dx The increment (step size).
	 */
	template <typename T, typename U, typename V>
	decltype(std::declval<T&>()* std::declval<U&>()* std::declval<V&>()) conjugateInnerProduct(size_t len, const T* __restrict arr1, const U* __restrict arr2, V dx) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0;
		for (size_t i = 0; i < len; i++)
			sum += std::conj(arr1[i]) * arr2[i];
		return sum * dx;
	}

	inline double conjugateInnerProduct(size_t len, const double* __restrict arr1, const double* __restrict arr2, double dx) {
		return cblas_ddot(len, arr1, 1, arr2, 1) * dx;
	}

	inline std::complex<double> conjugateInnerProduct(size_t len, const std::complex<double>* __restrict arr1, const std::complex<double>* __restrict arr2, double dx) {
		return cblas_zdotc(len, arr1, 1, arr2, 1) * dx;
	}

	/**
	 * Trapezoidal integration of an array.
	 * @tparam T The type of the array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to sum.
	 * @param dx The increment (step size).
	 */
	template <typename T, typename U>
	decltype(std::declval<T&>() * std::declval<U&>()) trapz(size_t len, const T* __restrict arr, U dx) {
		T sum = (arr[0] + arr[len - 1]) / 2.0;
		for (size_t i = 1; i < len - 1; i++) {
			sum += arr[i];
		}
		return sum * dx;
	}

	// Multiplies the elements of two arrays and then trapezoidally integrates the result.
	/**
	 * Trapezoidal integration of the elementwise product of two arrays.
	 * @tparam T The type of the first array elements.
	 * @tparam U The type of the second array elements.
	 * @tparam V The type of the increment (dx).
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to sum.
	 * @param arr2 (in) The second array to sum.
	 * @param dx The increment (step size).
	 */
	template <typename T, typename U, typename V>
	decltype(std::declval<T&>()* std::declval<U&>()* std::declval<V&>()) trapzMul(size_t len, const T* __restrict arr1, const U* __restrict arr2, V dx) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = (arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) / 2.0;
		for (size_t i = 1; i < len - 1; i++)
			sum += arr1[i] * arr2[i];
		return sum * dx;
	}

	/**
	 * Simpson's rule integration of an array.
	 * @tparam T The type of the array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to sum.
	 * @param dx The increment (step size).
	 */
	template <typename T, typename U>
	decltype(std::declval<T&>()* std::declval<U&>()) simps(size_t len, const T* __restrict arr, U dx) {
		if (len % 2) {
			T sum = (arr[len - 1] + arr[0]) / 4.0;
			for (size_t i = 1; i < len - 1; i += 2)
				sum += arr[i];
			sum *= 2.0;
			for (size_t i = 2; i < len - 1; i += 2)
				sum += arr[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			T sum = (5.0 * (arr[0] + arr[len - 1]) + 13.0 * (arr[1] + arr[len - 2])) / 12.0;
			for (size_t i = 2; i < len - 2; i++)
				sum += arr[i];
			return sum * dx;
		}
	}

	/**
	 * Simpson's rule integration of the elementwise product of two arrays.
	 * @tparam T The type of the first array elements.
	 * @tparam U The type of the second array elements.
	 * @tparam V The type of the increment (dx).
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to sum.
	 * @param arr2 (in) The second array to sum.
	 * @param dx The increment (step size).
	 */
	template <typename T, typename U, typename V>
	decltype(std::declval<T&>()* std::declval<U&>()* std::declval<V&>()) simpsMul(size_t len, const T* __restrict arr1, const U* __restrict arr2, V dx) {
		if (len % 2) {
			decltype(std::declval<T&>() * std::declval<U&>()) sum = (arr1[len - 1] * arr2[len - 1] + arr1[0] * arr2[0]) / 4.0;
			for (size_t i = 1; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			sum *= 2.0;
			for (size_t i = 2; i < len - 1; i += 2)
				sum += arr1[i] * arr2[i];
			return sum * dx * 2.0 / 3.0;
		}
		else {
			decltype(std::declval<T&>() * std::declval<U&>()) sum =
				(
					5.0 * (arr1[0] * arr2[0] + arr1[len - 1] * arr2[len - 1]) +
					13.0 * (arr1[1] * arr2[1] + arr1[len - 2] * arr2[len - 2])
					) / 12.0;
			for (size_t i = 2; i < len - 2; i++)
				sum += arr1[i] * arr2[i];
			return sum * dx;
		}
	}
	/**
	 * Cumulative Riemannian integration where the input array defines the value of the intervals to the left of each point (the first element of the result is arr[0]*dx).
	 * 
	 * \f$s_i = s_{i-1} + arr_i dx\f$
	 * @tparam T The type of the input array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to integrate.
	 * @param dx The increment (step size).
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U>
	void cumIntRectLeft(size_t len, const T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0.0;
		for (size_t i = 0; i < len; i++) {
			sum += arr[i] * dx;
			targ[i] = sum;
		}
	}

	/**
	 * Cumulative Riemannian integration where the input array defines the value of the intervals to the right of each point (the first element of the result is 0).
	 * 
	 * \f$s_i = s_{i-1} + arr_{i-1} dx\f$
	 * @tparam T The type of the input array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to integrate.
	 * @param dx The increment (step size).
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U>
	void cumIntRectRight(size_t len, const T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0.0;
		for (size_t i = 0; i < len; i++) {
			targ[i] = sum;
			sum += arr[i] * dx;
		}
	}

	/**
	 * Cumulative trapezoidal integration where the input array is on the same grid as the output. The first element of the result is 0.
	 * 
	 * \f$s_i = s_{i-1} + \frac{dx}{2} (arr_{i-1} + arr_i)\f$
	 * @tparam T The type of the input array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to integrate.
	 * @param dx The increment (step size).
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U>
	void cumIntTrapz(size_t len, const T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = 0.0;
		for (size_t i = 0; i < len - 1; i++) {
			targ[i] = sum;
			sum += (arr[i] + arr[i + 1]) * (dx / 2.0);
		}
		targ[len - 1] = sum;
	}

	/**
	 * Cumulative trapezoidal integration where the output array is one half step to the right of the input array.
	 * The first element of the result is \f$s_0 = \frac{dx}{4} (arr_0 + arr_1)\f$.
	 * 
	 * \f$s_i = s_{i-1} + \frac{dx}{4} (arr_{i-1} + 2 arr_{i} + arr_{i+1})\f$
	 * @tparam T The type of the input array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to integrate.
	 * @param dx The increment (step size).
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U>
	void cumIntTrapzToRight(size_t len, const T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) sum = (arr[0]+arr[1]) * dx / 4.0;
		for (size_t i = 0; i < len - 2; i++) {
			targ[i] = sum;
			sum += (arr[i] + 2.0*arr[i+1] + arr[i+2]) * (dx / 4.0);
		}
		targ[len - 2] = sum;
		targ[len - 1] = targ[len-2] + (arr[len-2] + arr[len-1]) * (dx / 4.0);
	}

	/**
	 * Cumulative trapezoidal integration where the output array is one half step to the left of the input array.
	 * The first element of the result is 0.
	 * 
	 * \f$s_i = s_{i-1} + \frac{dx}{4} (arr_{i-2} + 2 arr_{i-1} + arr_{i})\f$
	 * @tparam T The type of the input array elements.
	 * @tparam U The type of the increment (dx).
	 * @param len The length of the array.
	 * @param arr (in) The array to integrate.
	 * @param dx The increment (step size).
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U>
	void cumIntTrapzToLeft(size_t len, const T* __restrict arr, U dx, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		targ[0] = 0.0;
		decltype(std::declval<T&>() * std::declval<U&>()) sum = (arr[0]+arr[1]) * dx / 4.0;
		for (size_t i = 1; i < len - 1; i++) {
			targ[i] = sum;
			sum += (arr[i-1] + 2.0*arr[i] + arr[i+1]) * (dx / 4.0);
		}
		targ[len - 1] = sum;
	}
};

/**
 * @namespace vtls
 * @brief General mathematical tools for numerical calculations, including convolution, polynomial evaluation, and array manipulation.
 */
namespace vtls {

	/**
	 * @brief A class for performing convolution using FFTW.
	 * @tparam T The type of the input arrays (e.g., double, std::complex<double>).
	 * @details This class uses FFTW to perform the convolution between two arrays.
	 */
	template <class T>
	class Convolver{
	private:
		size_t len = -1;
		fftw_plan fp, bp;
		std::complex<double> *temp1, *temp2;
	public:
		/// Default constructor for Convolver. Do not call this constructor directly.
		Convolver<T>();

		/**
		 * Constructor for Convolver.
		 * @param len The length of the arrays to be convolved.
		 */
		Convolver<T>(size_t len);

		~Convolver();

		/**
		 * Performs the convolution of two arrays and stores the result in a target array.
		 * @param arr1 (in) The first input array to convolve.
		 * @param arr2 (in) The second input array to convolve.
		 * @param targ (out) The target array to store the result.
		 */
		void compute(const T* arr1, const T* arr2, T* targ);
	};

	/**
	 * @brief A class for efficiently performing convolutions of arrays with a constant mask using FFTW.
	 * @tparam T The type of the input arrays (e.g., double, std::complex<double>).
	 * @details This class uses FFTW to perform the convolution of an array with a constant mask, which is defined at construction.
	 * The mask is only transformed into reciprocal space once at construction.
	 */
	template <class T>
	class MaskConvolver{
	private:
		size_t len;
		fftw_plan fp, bp;
		std::complex<double> *mask, *temp;
	public:
		/**
		 * Constructor for MaskConvolver.
		 * @param len The length of the arrays to be convolved.
		 * @param maskIn (in) The mask to be used for the convolution.
		 */
		MaskConvolver<T>(size_t len, const T* maskIn);

		~MaskConvolver();

		/**
		 * Performs the out-of-place convolution of an array with the constant mask and stores the result in a target array.
		 * @param arr (in) The input array to convolve.
		 * @param targ (out) The target array to store the result of the convolution.
		 */
		void compute(const T* arr, T* targ);

		/**
		 * Performs the in-place convolution of an array with the constant mask, overwriting the input array with the result.
		 * @param arr (in/out) The input array to convolve, which will be overwritten with the result.
		 */
		void compute(T* arr);
	};

	/**
	 * Performs the matrix product \f$H D H\f$, where \f$H\f$ is a Hermitian matrix represented in upper triangular form and \f$D\f$ is a diagonal matrix represented as a 1D array.
	 * @tparam T The type of the elements in the Hermitian matrix.
	 * @tparam U The type of the elements in the diagonal matrix.
	 * @param len The length of the Hermitian matrix and diagonal matrix (number of rows or columns).
	 * @param hermitTriag (in) The Hermitian matrix represented in upper triangular form (column-major order, len*(len-1)/2 elements.
	 * @param diag (in) The diagonal matrix represented as a 1D array, len elements.
	 * @param targ (out) The target array to store the result of the matrix product, which will be a Hermitian matrix in upper triangular form (len*(len+1)/2 elements).
	 */
	template <typename T, typename U>
	void mulHermitDiagHermit(size_t len, const T* hermitTriag, const U* diag, decltype(std::declval<T&>()* std::declval<U&>())* targ) {
		decltype(std::declval<T&>() * std::declval<U&>()) csum = 0.0;
		for (size_t i = 0; i < len; i++) {
			for (size_t j = i; j < len; j++) {
				//targ[i,j] = sum_k triag[i,k]diag[k]triag[k,l]
				csum = 0.0;
				for (size_t k = 0; k < i; k++)
					csum += std::conj(hermitTriag[k + (i * (i + 1)) / 2]) * diag[k] * hermitTriag[k + (j * (j + 1)) / 2];
				for(size_t k = i; k < j; k++)
					csum += hermitTriag[i + (k * (k + 1)) / 2] * diag[k] * hermitTriag[k + (j * (j + 1)) / 2];
				for (size_t k = j; k < len; k++)
					csum += hermitTriag[i + (k * (k + 1)) / 2] * diag[k] * std::conj(hermitTriag[j + (k * (k + 1)) / 2]);
				targ[i + (j * (j + 1)) / 2] = csum;
			}
		}
	}

	// Adds two arrays, taking only the imaginary component of the first
	/**
	 * Adds the imaginary component of one array to the real component of the second array, storing the result in the second array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array, whose imaginary component will be added.
	 * @param arr2targ (in/out) The second array, which will be modified to store the result of the addition.
	 */
	void addArraysImag(size_t len, const std::complex<double>* arr1, double* arr2targ);

	/// @deprecated This function has been removed due to exprtk functionality being removed.
	template <typename T, typename U>
	void evalMathExpr(size_t len, const char* var, T* vals, std::string expr, U* res) {
		throw std::runtime_error("exprtk functionality has been removed for 'vtls::evalMathExpr'. Future versions should implements lambdas from pybind11.");
		/*typedef exprtk::symbol_table<T> symbol_table_t;
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

		for (size_t i = 0; i < len; i++) {
			cval = vals[i];
			res[i] = (U)(expression.value());
		}*/
	}

	/**
	 * Evaluates a polynomial at each point in an array using the given polynomial coefficients.
	 * @tparam T The type of the sample points.
	 * @tparam U The type of the polynomial coefficients.
	 * @param len The length of the array of sample points and the target array.
	 * @param x (in) The array of sample points at which to evaluate the polynomial.
	 * @param nPoly The number of polynomial coefficients (should be 2 for linear, 3 for quadratic, etc.).
	 * @param polyCoeffs (in) The array of polynomial coefficients, of length nPoly.
	 * @param y (out) The target array to store the results of the polynomial evaluation, which will be of length len.
	 */
	template <typename T, typename U>
	void polyEval(size_t len, const T* x, size_t nPoly, const U* __restrict polyCoeffs, decltype(std::declval<T&>()* std::declval<U&>())* y) {
		for (size_t i = 0; i < len; i++)
			y[i] = boost::math::tools::evaluate_polynomial(polyCoeffs, x[i], nPoly);
	}

	/**
	 * Adds a scalar multiple of one array to a second array, storing the result in the second array.
	 * @tparam T The type of the scalar.
	 * @tparam U The type of the first array.
	 * @tparam V The type of the second array.
	 * @param len The length of the arrays.
	 * @param scalar The scalar value to multiply the first array by.
	 * @param arr1 (in) The first array, which will be multiplied by the scalar.
	 * @param arr2targ (in/out) The second array, which will contain the result.
	 */
	template <typename T, typename U, typename V>
	void scaMulAddArrays(size_t len, T scalar, const U* __restrict arr1, V* __restrict arr2targ) {
		for (size_t i = 0; i < len; i++)
			arr2targ[i] += arr1[i] * scalar;
	}

	/**
	 * Adds a scalar multiple of one array to another array, storing the result in a target array.
	 * @tparam T The type of the scalar.
	 * @tparam U The type of the first array.
	 * @tparam V The type of the second array.
	 * @tparam W The type of the target array.
	 * @param len The length of the arrays.
	 * @param scalar The scalar value to multiply the first array by.
	 * @param arr1 (in) The first array, which will be multiplied by the scalar.
	 * @param arr2 (in) The second array, which will be added to the scaled first array.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U, typename V, typename W>
	void scaMulAddArrays(size_t len, T scalar, const U* __restrict arr1, const V* __restrict arr2, W* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = arr1[i] * scalar + arr2[i];
	}

	/**
	 * Adds the real component of a scalar multiple of one array to another array.
	 * @tparam T The type of the scalar.
	 * @tparam U The type of the first array.
	 * @param len The length of the arrays.
	 * @param scalar The scalar value to multiply the first array by.
	 * @param arr (in) The first array, which will be multiplied by the scalar.
	 * @param targ (out) The second array, which will be added to with the result.
	 */
	template <typename T, typename U>
	void scaMulAddArraysRe(size_t len, T scalar, const U* __restrict arr, double* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] += std::real(arr[i] * scalar);
	}

	/**
	 * Adds the imaginary component of a scalar multiple of one array to another array.
	 * @tparam T The type of the scalar.
	 * @tparam U The type of the first array.
	 * @param len The length of the arrays.
	 * @param scalar The scalar value to multiply the first array by.
	 * @param arr (in) The first array, which will be multiplied by the scalar.
	 * @param targ (out) The second array, which will be added to with the result.
	 */
	template <typename T, typename U>
	void scaMulAddArraysIm(size_t len, T scalar, const U* __restrict arr, double* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] += std::imag(arr[i] * scalar);
	}

	/**
	 * Adds two arrays together, storing the result in a third array.
	 * @tparam T The type of the first array.
	 * @tparam U The type of the second array.
	 * @tparam V The type of the target array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to add.
	 * @param arr2 (in) The second array to add.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U, typename V>
	void addArrays(size_t len, const T* __restrict arr1, const U* __restrict arr2, V* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = arr1[i] + arr2[i];
	}

	/**
	 * Averages two arrays, storing the result in a third array.
	 * @tparam T The type of the first array.
	 * @tparam U The type of the second array.
	 * @tparam V The type of the target array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to average.
	 * @param arr2 (in) The second array to average.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U, typename V>
	void averageArrays(size_t len, const T* __restrict arr1, const U* __restrict arr2, V* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = (arr1[i] + arr2[i])*0.5;
	}

	/**
	 * Averages an array with a target array, storing the result in the target array.
	 * @tparam T The type of the first array.
	 * @tparam U The type of the target array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to average.
	 * @param arr2targ (in/out) The target array, which will contain the result of the averaging.
	 */
	template <typename T, typename U>
	void averageArrays(size_t len, const T* __restrict arr1, U* __restrict arr2targ) {
		for (size_t i = 0; i < len; i++)
			arr2targ[i] = (arr1[i] + arr2targ[i]) * 0.5;
	}

	/**
	 * Adds the elements of one array to another, storing the result in the second array.
	 * @tparam T The type of the first array.
	 * @tparam U The type of the second array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to add.
	 * @param arr2targ (in/out) The second array, which will contain the result of the addition.
	 */
	template <typename T, typename U>
	void addArrays(size_t len, const T* __restrict arr1, U* __restrict arr2targ) {
		for (size_t i = 0; i < len; i++)
			arr2targ[i] += arr1[i];
	}

	/**
	 * Multiplies the elements of two arrays together, storing the result in a third array.
	 * @tparam T The type of the first array.
	 * @tparam U The type of the second array.
	 * @tparam V The type of the target array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to multiply.
	 * @param arr2 (in) The second array to multiply.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U, typename V>
	void seqMulArrays(size_t len, const T* __restrict arr1, const U* __restrict arr2, V* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = arr1[i] * arr2[i];
	}
	
	/**
	 * Multiplies the elements of two arrays together, storing the result in the second array.
	 * @tparam T The type of the first array.
	 * @tparam U The type of the second array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to multiply.
	 * @param arr2targ (in/out) The second array, which will contain the result of the multiplication.
	 */
	template <typename T, typename U>
	void seqMulArrays(size_t len, const T* __restrict arr1, U* __restrict arr2targ) {
		for (size_t i = 0; i < len; i++)
			arr2targ[i] *= arr1[i];
	}

	/**
	 * Multiplies the elements of two arrays together and adds the result to a target array.
	 * @tparam T The type of the first array.
	 * @tparam U The type of the second array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to multiply.
	 * @param arr2 (in) The second array to multiply.
	 * @param targ (out) The target array to add the result of the multiplication.
	 */
	template<typename T, typename U>
	void seqMulAddArrays(size_t len, const T* __restrict arr1, const U* __restrict arr2, decltype(std::declval<T&>()* std::declval<U&>())* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] += arr1[i] * arr2[i];
	}

	/**
	 * Multiplies an array by a scalar multiple and stores the result in a target array.
	 * @tparam T The type of the scalar.
	 * @tparam U The type of the input array.
	 * @tparam V The type of the target array.
	 * @param len The length of the arrays.
	 * @param scalar The scalar value to multiply the input array by.
	 * @param arr (in) The input array to multiply.
	 * @param targ (out) The target array to store the result of the multiplication.
	 * @note There are type-specific in-place overloads which, in implementation, use cblas for further optimization.
	 */
	template <typename T, typename U, typename V>
	void scaMulArray(size_t len, T scalar, const U* __restrict arr, V* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = arr[i] * scalar;
	}

	/**
	 * Multiplies an array by a scalar multiple and stores the result in the same array.
	 * @param len The length of the array.
	 * @param scalar The scalar value to multiply the array by.
	 * @param arr (in/out) The array to multiply, which will be modified to store the result of the multiplication.
	 */
	void scaMulArray(size_t len, double scalar, double* __restrict arr);

	/// @copydoc vtls::scaMulArray(size_t,double,double*)
	void scaMulArray(size_t len, double scalar, std::complex<double>* __restrict arr);

	/// @copydoc vtls::scaMulArray(size_t,double,double*)
	void scaMulArray(size_t len, std::complex<double> scalar, std::complex<double>* __restrict arr);

	/**
	 * Multiplies an array by a scalar multiple and stores the real part of the result in a target array.
	 * @tparam T The type of the scalar.
	 * @tparam U The type of the input array.
	 * @tparam V The type of the target array.
	 * @param len The length of the arrays.
	 * @param scalar The scalar value to multiply the input array by.
	 * @param arr (in) The input array to multiply.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T, typename U, typename V>
	void scaMulArrayRe(size_t len, T scalar, const U* __restrict arr, V* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = std::real(arr[i] * scalar);
	}

	/**
	 * Calculates the square norm of each element in an array and stores the result in a target array.
	 * 
	 * \f$ \mathrm{normSqr}(arr_i) = |arr_i|^2 \f$
	 * @tparam T The type of the input array elements.
	 * @param len The length of the array.
	 * @param arr (in) The input array to be processed.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T>
	void normSqr(size_t len, const T* __restrict arr, double* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = std::norm(arr[i]);
	}

	/**
	 * Calculates the absolute value of each element in an array and stores the result in a target array.
	 * 
	 * \f$ \mathrm{abs}(arr_i) = |arr_i| \f$
	 * @tparam T The type of the input array elements.
	 * @param len The length of the array.
	 * @param arr (in) The input array to be processed.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T>
	void abs(size_t len, const T* __restrict arr, double* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = std::abs(arr[i]);
	}

	/**
	 * Normalizes an array by dividing each element by the the L2 norm of the array assuming uniform grid spacing \f$dx\f$.
	 * 
	 * \f$ arr_i = \frac{arr_i}{\sqrt{\sum_j |arr_j|^2 \cdot dx}} \f$
	 * @tparam T The type of the input array elements.
	 * @param len The length of the array.
	 * @param arr (in/out) The input array to be normalized, which will be modified to store the normalized values.
	 * @param dx The grid spacing.
	 * @note This function allocates and subsequently frees a temporary array to store the squared norms. Consider implementing a version which accepts a user-allocated workspace if this will be used often.
	 */
	template <typename T>
	void normalizeSqrNorm(size_t len, T* __restrict arr, double dx) {
		double* temp = (double*) sq_malloc(sizeof(double)*len);
		normSqr(len, arr, temp);
		scaMulArray(len, 1.0 / std::sqrt(vtlsInt::simps(len, temp, dx)), arr);
		sq_free(temp);
	}

	/**
	 * Calculates the L2 norm squared of an array assuming uniform grid spacing \f$dx\f$.
	 * 
	 * \f$ \mathrm{norm} = \sum_j |arr_j|^2 \cdot dx \f$
	 * @tparam T The type of the input array elements.
	 * @param len The length of the array.
	 * @param arr (in) The input array to be processed.
	 * @param dx The grid spacing.
	 * @return The L2 norm of the array.
	 */
	template <typename T>
	double getNorm(size_t len, const T* __restrict arr, double dx) {
		double sm = 0.0;
		for (size_t i = 0; i < len; i++)
			sm += std::norm(arr[i]);
			//sm += std::pow(std::abs(arr[i]), 2);
		return sm *= dx;
	}

	inline double getNorm(size_t len, const double* __restrict arr, double dx) {
		double val = cblas_dnrm2(len, arr, 1);
		return val*val * dx;
	}
	
	inline double getNorm(size_t len, const std::complex<double>* __restrict arr, double dx) {
		double val = cblas_dznrm2(len, arr, 1);
		return val*val * dx;
	}

	/**
	 * Sets the L2 norm of an array to a specified value assuming uniform grid spacing \f$dx\f$.
	 * 
	 * @tparam T The type of the input array elements.
	 * @param len The length of the array.
	 * @param arr (in/out) The input array to be normalized, which will be modified to store the normalized values.
	 * @param dx The grid spacing.
	 * @param norm The target L2 norm (squared).
	 */
	template <typename T>
	void setNorm(size_t len, T* __restrict arr, double dx, double norm) {
		scaMulArray(len, std::sqrt(norm / getNorm(len, arr, dx)), arr);
	}

	/**
	 * Linearly interpolates an array to a new length, including the edge points. The first and last points are the same.
	 * @tparam T The type of the array elements.
	 * @param len The length of the input array.
	 * @param arr (in) The input array to be interpolated.
	 * @param newLen The length of the target array.
	 * @param targ (out) The target array to store the result of the interpolation.
	 */
	template <typename T>
	void linearInterpolateEdge(size_t len, const T* __restrict arr, size_t newLen, T* __restrict targ) {
		double step = (double)(len - 1) / (newLen - 1);
		double curPos = 0.0;
		for (size_t i = 0; i < newLen; i++) {
			targ[i] = (1.0 - std::fmod(curPos, 1)) * arr[(size_t)curPos] + (fmod(curPos, 1)) * arr[(size_t)curPos + 1 * ((size_t)curPos != (len - 1))];
			curPos += step;
		}
	}

	/**
	 * Linearly interpolates an array to a new length, excluding the edge points. The first and last points are half a (new) step away from the edge.
	 * @tparam T The type of the array elements.
	 * @param len The length of the input array.
	 * @param arr (in) The input array to be interpolated.
	 * @param newLen The length of the target array.
	 * @param targ (out) The target array to store the result of the interpolation.
	 */
	template <typename T>
	void linearInterpolateNoEdge(size_t len, const T* __restrict arr, size_t newLen, T* __restrict targ) {
		double step = (double)(len - 1) / newLen;
		double curPos = step / 2.0;
		for (size_t i = 0; i < newLen; i++) {
			targ[i] = (1.0 - std::fmod(curPos, 1)) * arr[(size_t)curPos] + (fmod(curPos, 1)) * arr[(size_t)curPos + 1 * ((size_t)curPos != (len - 1))];
			curPos += step;
		}
	}

	/**
	 * Linearly interpolates an array \f$y_1(x_1)\f$ to a new array \f$y_2(x_2)\f$.
	 * @tparam T The type of the array elements.
	 * @param l1 The length of the first array.
	 * @param x1 (in) The x-values of the first array.
	 * @param y1 (in) The y-values of the first array.
	 * @param l2 The length of the second array.
	 * @param x2 (in) The x-values of the second array.
	 * @param y2 (out) The target array to store the result of the interpolation.
	 */
	template <typename T>
	void linearInterpolate(size_t l1, const double* __restrict x1, const T* __restrict y1, size_t l2, const double* __restrict x2, T* __restrict y2) {
		size_t curPos = 0;
		for (size_t i = 0; i < l2; i++) {
			while (curPos < l1 && x2[i] > x1[curPos])
				curPos++;
			if (curPos > 0 && curPos < l1)
				y2[i] = y1[curPos - 1] + (y1[curPos] - y1[curPos - 1]) * (x2[i] - x1[curPos - 1]) / (x1[curPos] - x1[curPos - 1]);
			else if (curPos < 1)
				y2[i] = y1[0];
			else
				y2[i] = y1[l1 - 1];
		}
	}

	/**
	 * Interpolates a single point from an array using linear interpolation.
	 * @tparam T The type of sample points.
	 * @tparam U The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to interpolate.
	 * @param xStart The starting x-value of the array.
	 * @param dx The increment (step size).
	 * @param samp The sample point at which to interpolate.
	 * @return The interpolated value at the sample point.
	 */
	template <typename T, typename U>
	U linearInterpolate(size_t len, const U* __restrict arr, T xStart, T dx, T samp) {
		size_t ix = (size_t)((samp - xStart) / dx);
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

	/**
	 * Generates a linearly spaced array. The edge points are included.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param min The minimum (first) value of the array.
	 * @param max The maximum (last) value of the array.
	 * @param targ (out) The target array to store the result.
	 */
	template <typename T>
	void linspace(size_t len, T min, T max, T* __restrict targ) {
		for (size_t i = 0; i < len; i++)
			targ[i] = (max - min) * i / (T)(len - 1) + min;
	}

	/**
	 * Generates a linearly spaced array. The edge points are included.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param min The minimum (first) value of the array.
	 * @param max The maximum (last) value of the array.
	 * @return The generated array.
	 */
	template <typename T>
	std::vector<T> linspace(size_t len, T min, T max) {
		std::vector<T> ret = std::vector<T>(len);
		for (size_t i = 0; i < len; i++)
			ret[i] = (max - min) * i / (T)(len - 1) + min;
		return ret;
	}

	/**
	 * Adds a scalar to each element of an array.
	 * @tparam T The type of the scalar.
	 * @tparam U The type of the array elements.
	 * @param len The length of the array.
	 * @param scalar The scalar value to add to each element.
	 * @param arr (in/out) The array to which the scalar will be added.
	 */
	template <typename T, typename U>
	void scaAddArray(size_t len, T scalar, U* __restrict arr) {
		for (size_t i = 0; i < len; i++)
			arr[i] += scalar;
	}

	/**
	 * Copies the elements of one array to another.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The source array to copy from.
	 * @param arr2 (out) The target array to copy to.
	 */
	void copyArray(size_t len, const double* __restrict arr1, double* __restrict arr2);

	/// @copydoc copyArray(size_t len, const double* __restrict arr1, double* __restrict arr2)
	void copyArray(size_t len, const std::complex<double>* __restrict arr1, std::complex<double>* __restrict arr2);

	/// @copydoc copyArray(size_t len, const double* __restrict arr1, double* __restrict arr2)
	void copyArray(size_t len, const double* __restrict arr1, std::complex<double>* __restrict arr2);

	/**
	 * Copies the elements of one array to another, taking the real part of the source array.
	 * @tparam T The type of the source array elements.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The source array to copy from.
	 * @param arr2 (out) The target array to copy to.
	 */
	template <typename T>
	void copyArrayRe(size_t len, const T* __restrict arr1, double* __restrict arr2);

	/**
	 * Copies the elements of one array to another, taking the complex conjugate of the source array.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The source array to copy from.
	 * @param arr2 (out) The target array to copy to.
	 */
	inline void copyArrayConj(size_t len, const std::complex<double>* __restrict arr1, std::complex<double>* __restrict arr2) {
		for (size_t i = 0; i < len; i++)
			arr2[i] = std::conj(arr1[i]);
	}

	/**
	 * Copies the elements of one array to another, taking the complex conjugate of the source array and multiplying by a scalar.
	 * @param len The length of the arrays.
	 * @param s The scalar to multiply the conjugated elements by.
	 * @param arr1 (in) The source array to copy from.
	 * @param arr2 (out) The target array to copy to.
	 */
	inline void copyArrayConj(size_t len, std::complex<double> s, const std::complex<double>* __restrict arr1, std::complex<double>* __restrict arr2) {
		for (size_t i = 0; i < len; i++)
			arr2[i] = s * std::conj(arr1[i]);
	}

	/**
	 * Evaluates the first derivative of an array at a specified position.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to differentiate.
	 * @param pos The position at which to evaluate the derivative.
	 * @param dx The increment (step size).
	 * @return The value of the first derivative at the specified position.
	 */
	template <typename T>
	T firstDerivative(size_t len, const T* __restrict arr, size_t pos, double dx) {
		if (pos == 0)
			return (arr[1] - arr[0]) / dx;
		else if (pos == len - 1)
			return (arr[pos] - arr[pos - 1]) / dx;
		else
			return (arr[pos + 1] - arr[pos - 1]) / (2.0 * dx);
	}

	/**
	 * Evaluates the first derivative of an array at each point.
	 * On the edges, the forward or backward difference is used. Otherwise a central difference is used.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to differentiate.
	 * @param targ (out) The target array to store the result.
	 * @param dx The increment (step size).
	 */
	template <typename T>
	void firstDerivative(size_t len, const T* __restrict arr, T* __restrict targ, double dx) {
		targ[0] = (arr[1] - arr[0]) / dx;
		for (size_t i = 1; i < len - 1; i++)
			targ[i] = (arr[i + 1] - arr[i - 1]) / (2.0 * dx);
		targ[len - 1] = (arr[len - 1] - arr[len - 2]) / dx;
	}

	/**
	 * Evaluates the second derivative of an array at each point.
	 * On the edges, the forward or backward difference is used. Otherwise a central difference is used.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to differentiate.
	 * @param targ (out) The target array to store the result.
	 * @param dx The increment (step size).
	 */
	template <typename T>
	void secondDerivative(size_t len, const T* __restrict arr, T* __restrict targ, double dx) {
		double dx2 = dx * dx;
		targ[0] = (arr[1] - arr[0]) * 2.0 / dx2;
		for (size_t i = 1; i < len - 1; i++)
			targ[i] = (arr[i - 1] - 2.0 * arr[i] + arr[i + 1]) / dx2;
		targ[len - 1] = (arr[len - 2] - arr[len - 1]) * 2.0 / dx2;
	}

	/**
	 * Finds the index of the first element in an array that is greater than or equal to a specified value.
	 * The array must be sorted in ascending order.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to search.
	 * @param val The value to search for.
	 * @return The index of the first element greater than or equal to the specified value, len-1 if the value is greater than all elements, or 0 if the value is less than all elements.
	 */
	size_t findValue(size_t len, const double* __restrict arr, double val);

	/**
	 * Sorts the input array in ascending order using insertion sort.
	 * @param len The length of the array.
	 * @param arr (in/out) The array to sort, which will be modified to store the sorted values.
	 * @param idxs (out) The target array to store the indices of the sorted values.
	 */
	void insertSort_idxs(size_t len, double* __restrict arr, size_t* __restrict idxs);

	/**
	 * Finds the maximum value in an array.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to search.
	 * @return The maximum value in the array.
	 */
	template <typename T>
	T max(size_t len, const T* __restrict arr) {
		T mx = arr[0];
		for (size_t i = 1; i < len; i++)
			if (arr[i] > mx)
				mx = arr[i];
		return mx;
	}

	/**
	 * Finds the minimum value in an array.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to search.
	 * @return The minimum value in the array.
	 */
	template <typename T>
	T min(size_t len, const T* __restrict arr) {
		T mn = arr[0];
		for (size_t i = 1; i < len; i++)
			if (arr[i] < mn)
				mn = arr[i];
		return mn;
	}

	/**
	 * Tests if the maximum absolute difference between two arrays exceeds a specified threshold.
	 * @tparam T The type of the array elements.
	 * @param len The length of the arrays.
	 * @param arr1 (in) The first array to compare.
	 * @param arr2 (in) The second array to compare.
	 * @param maxDiff The maximum allowed absolute difference.
	 * @return True if the maximum absolute difference exceeds maxDiff, false otherwise.
	 */
	template <typename T>
	bool testMaxAbsDiffExceedsThresh(size_t len, const T* __restrict arr1, const T* __restrict arr2, decltype(std::abs(std::declval<T&>())) maxDiff) {
		for (size_t i = 0; i < len; i++)
			if (std::abs(arr1[i] - arr2[i]) > maxDiff)
				return true;
		return false;
	}

	/**
	 * Checks if an array contains any NaN (Not a Number) values.
	 * @tparam T The type of the array elements.
	 * @param len The length of the array.
	 * @param arr (in) The array to check for NaN values.
	 * @return True if the array contains any NaN values, false otherwise.
	 */
	template <typename T>
	bool hasNaN(size_t len, const T* __restrict arr) {
		for (size_t i = 0; i < len; i++)
			if (std::isnan(arr[i]))
				return true;
		return false;
	}

	/// @copydoc hasNaN(size_t len, const T* __restrict arr)
	inline bool hasNaN(size_t len, const std::complex<double>* __restrict arr) {
		for (size_t i = 0; i < len; i++)
			if (std::isnan(std::real(arr[i])) || std::isnan(std::imag(arr[i])))
				return true;
		return false;
	}

	template <typename T>
	bool hasNegative(size_t len, const T* __restrict arr) {
		for (size_t i = 0; i < len; i++)
			if (arr[i] < 0)
				return true;
		return false;
	}

	/**
	 * Returns an array of 1's except for within the external region where it is a 13th order
	 * polynomial that is 1 at the inner boundary and 0 at the outer boundary. The polynomial is
	 * then raised to the power of \a rate. The result is continuous up to 7th order at the inner and outer boundaries.
	 * The side of the boundary is determined by the order of inner and outer.
	 * 
	 * @warning This function allocates memory for the array. It is the responsibility of the caller to free this memory.
	 * @param len The number of elements in the vector.
	 * @param inner The inner boundary position.
	 * @param outer The outer boundary position.
	 * @param rate The exponent of the polynomial. Larger values make for a stronger decay.
	 * @return A unique_ptr to an array of doubles representing the smoothed decay mask.
	*/
	std::unique_ptr<double[]> getPolynomialSmoothBoundary(size_t len, size_t inner, size_t outer, double rate);

	template <class T>
	struct Orthonormalizer {
	private:
		T *tau = nullptr, *work = nullptr;
		lapack_int nPts, nVecs, lwork;

		/**
		 * Orthonormalizes a set of vectors using QR factorization via Lapack's ?geqrf and ?(or/un)gqr routines.
		 * @param nPts The length of each vector.
		 * @param nVecs The number of vectors to orthonormalize.
		 * @param vecs (in/out) The array of vectors to orthonormalize. The vectors are stored in a 1D array in row-major order.
		 * @param tau (out) The array to store the scalar factors for the Householder reflections.
		 * @param work (out) The workspace array for the QR factorization.
		 * @param lwork The length of the workspace array.
		 */
		static void orthonormalize(lapack_int nPts, lapack_int nVecs, T* __restrict vecs, T* __restrict tau, T* __restrict work, lapack_int lwork);
		
		/**
		 * Calculates the optimal size of the workspace array for the QR factorization.
		 * @param nPts The length of each vector.
		 * @param nVecs The number of vectors to orthonormalize.
		 * @param lwork (out) The length of the workspace array.
		 * @note This function uses Lapack's ?geqrf to determine the optimal size of the workspace array.
		 */
		static void getlwork(lapack_int nPts, lapack_int nVecs, lapack_int* lwork);
	public:
		/**
		 * Constructor. Use this as opposed to the static orthonormalize function to manage working memory.
		 * @param nPts The length of each vector.
		 * @param nVecs The number of vectors to orthonormalize.
		 */
		Orthonormalizer(size_t nPts, size_t nVecs);
		~Orthonormalizer();

		/**
		 * Orthonormalizes a set of vectors using QR factorization via Lapack's ?geqrf and ?(or/un)gqr routines.
		 * @param nPts The length of each vector.
		 * @param nVecs The number of vectors to orthonormalize.
		 * @param vecs (in/out) The array of vectors to orthonormalize. The vectors are stored in a 1D array in row-major order.
		 * @note This function allocates then frees working memory for the QR factorization.
		 */
		static void orthonormalize(size_t nPts, size_t nVecs, T* __restrict vecs);

		/**
		 * Orthonormalizes a set of vectors using QR factorization via Lapack's ?geqrf and ?(or/un)gqr routines.
		 * @param vecs (in/out) The array of vectors to orthonormalize. The vectors are stored in a 1D array in row-major order.
		 */
		void orthonormalize (T* __restrict vecs) {
			orthonormalize(nPts, nVecs, vecs, tau, work, lwork);
		};
	};

	/// A class used for fitting an array's history to a polynomial to project (or extrapolate) to the future time.
	/// I took this idea from Octopus, which apparently uses this method for time propagation without SCF.
	struct PolynomialExtrapolator {
	private:
		size_t nPts, order;
		double* extrapStenc = nullptr;
		double* history = nullptr;
		cyclic_int<size_t> historyIndex;
	public:
		/**
		 * Constructor.
		 * @param nPts The number of points in the vector to extrapolate.
		 * @param order The order of the polynomial (minus one) to use for extrapolation. order = length of history used.
		 * @param stepFraction The fraction of the step to extrapolate. For example, 1.0 extrapolates to the next step, 0.5 extrapolates to the middle of the next step. 0.0 would return the last point.
		 * @param initialVector (in) An optional initial vector to use for the first extrapolation. If not provided, the first extrapolation will be zero.
		 */
		PolynomialExtrapolator(size_t nPts, size_t order, double stepFraction, const double* __restrict initialVector = nullptr);

		~PolynomialExtrapolator(){
			if (extrapStenc) sq_free(extrapStenc);
			if (history) sq_free(history);
		};

		/**
		 * Pushes a new vector onto the history for extrapolation.
		 * @param vec (in) The vector to push onto the history.
		 */
		void pushHistory(const double* __restrict vec);

		/**
		 * Fills the history with a vector, overwriting all previous history.
		 * @param vec (in) The vector to fill the history with.
		 */
		void fillHistory(const double* __restrict vec);

		/**
		 * Extrapolates the next vector using the history and the polynomial coefficients.
		 * @param targ (out) The target array to store the extrapolated vector.
		 */
		void extrapolate(double* __restrict targ);

		/// Prints the history of vectors used for extrapolation.
		void printHistory() const;

		/// Prints the extrapolation stencil used for extrapolation. Each row is a permutation of the first to minimize reorganizing the vector at each step.
		void printExtrapStenc() const;
	};

	namespace masks{

		struct step{
			/**
			 * @param nPts The number of points in the mask.
			 * @param center The centroid of the sigmoid.
			 * @param stepLength The sigmoid width. 
			 * 	The sign determines the direction of the step: positive values create a step from 0 to 1, negative values create a step from 1 to 0.
			 * 	If zero, a simple Heaviside step function is used. Using -0.0 is distinguished as a step from 1 to 0.
			 * @param mask (out) The target array to store the mask values.
			 */
			virtual void operator()(size_t nPts, size_t center, double stepLength, double* __restrict mask) const = 0;
		};

		inline void heaviside(size_t nPts, size_t center, int direction, double* __restrict mask) {
			center = std::clamp(center, (size_t)0, nPts - 1);
			if (direction > 0) {
				std::fill_n(mask, center, 0.0);
				std::fill_n(mask + center, nPts - center, 1.0);
			}
			else {
				std::fill_n(mask, center, 1.0);
				std::fill_n(mask + center, nPts - center, 0.0);
			}
		}

		// SIGMOID
		/**
		 * Generates a sigmoid-smoothed step function:
		 * \f[
		 * mask_i = \frac{1}{1 + e^{2 \cdot \frac{i - center}{maskLength}}}
		 * \f]
		 * 
		 * @copydoc step::operator()(size_t nPts, size_t center, double stepLength, double* __restrict mask)
		 */
		inline void sigmoid(size_t nPts, size_t center, double stepLength, double* __restrict mask) {
			center = std::clamp(center, (size_t)0, nPts - 1);

			if (std::abs(stepLength) > 0.0)
				for (size_t i = 0; i < nPts; i++)
					mask[i] = 1.0 / (1.0 + std::exp(-2.0 * ((double)i - (double)center) / stepLength));
			else
				heaviside(nPts, center, std::signbit(stepLength) ? -1 : 1, mask);
		}
		struct sigmoid_c : step {
			void operator()(size_t nPts, size_t center, double stepLength, double* __restrict mask) const override {
				sigmoid(nPts, center, stepLength, mask);
			}
		};

		// POLYNOMIAL SMOOTHED
		/**
		 * Generates a 13th order polynomial-smoothed step function. Continuous up to 7th order at each boundary.
		 * 
		 * @copydoc step::operator()(size_t nPts, size_t center, double stepLength, double* __restrict mask)
		 */
		inline void poly13(size_t nPts, size_t center, double stepLength, double* __restrict mask) {
			int direction = std::signbit(stepLength) ? -1 : 1;
			size_t stepLengthA = (size_t)std::abs(stepLength);
			stepLengthA = std::clamp(stepLengthA, (size_t)0, nPts - 1);
			center = std::clamp(center, stepLengthA/2, nPts - 1 - stepLengthA/2);
			size_t left = center - stepLengthA / 2;
			size_t right = center + stepLengthA / 2;
			
			if (stepLengthA > 0){
				std::fill_n(mask, left, direction == 1 ? 0.0 : 1.0);
				std::fill_n(mask + right, nPts - right, direction == 1 ? 1.0 : 0.0);
				double k;
				for (size_t i = left; i < right; i++) {
					k = (double)(i - left) / (double)(right - left);
					k = direction == 1 ? k : 1.0 - k; // reverse k for negative stepLength
					mask[i] = (
							924.0*std::pow(k, 13) -
							6006.0*std::pow(k, 12) +
							16380.0*std::pow(k, 11) -
							24024.0*std::pow(k, 10) +
							20020.0*std::pow(k, 9) -
							9009.0*std::pow(k, 8) +
							1716.0*std::pow(k, 7)
						);
				}
			}
			else
				heaviside(nPts, center, direction, mask);
		}
		struct poly13_c : step {
			void operator()(size_t nPts, size_t center, double stepLength, double* __restrict mask) const override {
				poly13(nPts, center, stepLength, mask);
			}
		};

		// BIDIRECTIONAL MASK
		inline void biMask(size_t nPts, size_t leftCenter, size_t rightCenter, double maskLength, double* __restrict mask, step &stepFunc){
			assert(leftCenter <= rightCenter);
			
			double* tempMask = (double*) sq_malloc(sizeof(double) * nPts);
			
			stepFunc(nPts, leftCenter, maskLength, tempMask);
			stepFunc(nPts, rightCenter, -maskLength, mask);
			if (!std::signbit(maskLength)) // positive, 0-1-0, multiply masks
				vtls::seqMulArrays(nPts, tempMask, mask);
			else // negative, 1-0-1, add masks
				vtls::addArrays(nPts, tempMask, mask);

			sq_free(tempMask);
		} 

		/**
		 * Generates a sigmoid-smoothed mask:
		 * \f[
		 * mask_i = \frac{1}{(1 + e^{2 \cdot \frac{leftCenter - i}{maskLength}}) \cdot (1 + e^{2 \cdot \frac{i - rightCenter}{maskLength}})}
		 * \f]
		 * 
		 * @param nPts The number of points in the mask.
		 * @param leftCenter The centroid of the left-sided sigmoid.
		 * @param rightCenter The centroid of the right-sided sigmoid.
		 * @param maskLength The sigmoid width. 
		 * 	The sign determines the shape of the mask: positive creates a mask from 0 to 1 then 0, and negative creates a mask from 1 to 0 then 1.
		 * 	If zero, a simple Heaviside step function is used.
		 * @param mask (out) The target array to store the mask values.
		 */
		inline void biSigmoid(size_t nPts, size_t leftCenter, size_t rightCenter, double maskLength, double* __restrict mask) {
			sigmoid_c f;
			biMask(nPts, leftCenter, rightCenter, maskLength, mask, f);
		}

		/**
		 * Generates a 13th order polynomial-smoothed mask. Continuous up to 7th order at each boundary.
		 * @param nPts The number of points in the mask.
		 * @param leftCenter The centroid of the left-sided polynomial.
		 * @param rightCenter The centroid of the right-sided polynomial.
		 * @param maskLength The polynomial width.
		 * The sign determines the shape of the mask: positive creates a mask from 0 to 1 then 0, and negative creates a mask from 1 to 0 then 1.
		 * If zero, a simple Heaviside step function is used.
		 * @param mask (out) The target array to store the mask values.
		 */
		inline void biPoly13(size_t nPts, size_t leftCenter, size_t rightCenter, double maskLength, double* __restrict mask) {
			poly13_c f;
			biMask(nPts, leftCenter, rightCenter, maskLength, mask, f);
		}
	};
};

/**
 * @namespace vtlsPrnt
 * @brief A namespace for printing data to the console.
 */
namespace vtlsPrnt {
	/**
	 * Prints the contents of an array to stdout.
	 * @tparam T The type of the array elements.
	 * @param n The length of the array.
	 * @param arr (in) The array to print.
	 */
	template <typename T>
	void printArray(size_t n, const T* __restrict arr) {
		std::cout << "[";
		for (size_t i = 0; i < n; i++) {
			std::cout << arr[i];
			if (i != n - 1) {
				std::cout << ",";
			}
		}
		std::cout << "]" << std::endl;
	}

	/**
	 * Plots the contents of an array in the console using text.
	 * @param n The length of the array.
	 * @param arr (in) The array to plot.
	 */
	void printGraph(size_t n, const double* __restrict arr);

	/// @copydoc printGraph(size_t n, const double* __restrict arr)
	void printGraph(size_t n, const std::complex<double>* __restrict arr);

	/**
	 * Saves an array to a binary file.
	 * @tparam T The type of the array elements.
	 * @param n The length of the array.
	 * @param fil The filename to save the array to.
	 * @param data (in) The array to save.
	 */
	template <typename T>
	void saveArray(size_t n, const char* fil, const T* data) {
		std::fstream fid(fil, std::ios::out | std::ios::binary);
		fid.write(reinterpret_cast<char*>(n), sizeof(size_t));
		fid.write(reinterpret_cast<char*>(data), sizeof(T) * n);
		fid.close();
	}
};

#include "gnuplot-iostream.h"

/**
 * @namespace plotting
 * @brief A namespace for plotting data using GNUPlot.
 */
namespace plotting{
	/**
	 * A class for plotting data using GNUPlot.
	 * @details This class uses the gnuplot-iostream library to plot data using GNUPlot.
	 */
	class GNUPlotter {
	private:
		Gnuplot gp;
		const char* ylabel = nullptr;
		const char* xlabel = nullptr;
	public:
		/// Default constructor for initializing a GNUPlotter.
		GNUPlotter(){};

		GNUPlotter(const char* xlabel, const char* ylabel) : xlabel(xlabel), ylabel(ylabel) {
			gp << "set grid\n";
			if (ylabel)
				gp << "set ylabel '" << ylabel << "'\n";
			if (xlabel)
				gp << "set xlabel '" << xlabel << "'\n";
		};

		/**
		 * Constructor for initializing a GNUPlotter with data.
		 * @param nPts The number of points in the data.
		 * @param nLines The number of lines in the data.
		 * @param x (in) The x-values of the data.
		 * @param y (in) The y-values of the data.
		 */
		GNUPlotter(size_t nPts, size_t nLines, const double* x, const double* y){update(nPts, nLines, x, y);};

		/**
		 * Constructor for initializing a GNUPlotter with data. The x-values are assumed to be the index of the data.
		 * @param nPts The number of points in the data.
		 * @param nLines The number of lines in the data.
		 * @param y (in) The y-values of the data.
		 */
		GNUPlotter(size_t nPts, size_t nLines, const double* y){update(nPts, nLines, y);};

		/**
		 * Updates the plot with new data.
		 * @param nPts The number of points in the data.
		 * @param nLines The number of lines in the data.
		 * @param x (in) The x-values of the data.
		 * @param y (in) The y-values of the data.
		 * @param xmin The minimum value of the x-axis to display.
		 * @param xmax The maximum value of the x-axis to display.
		 * @param ymin The minimum value of the y-axis to display.
		 * @param ymax The maximum value of the y-axis to display.
		 */
		void update(size_t nPts, size_t nLines, const double* x, const double* y, double xmin, double xmax, double ymin, double ymax);

		/**
		 * Updates the plot with new data.
		 * @param nPts The number of points in the data.
		 * @param nLines The number of lines in the data.
		 * @param x (in) The x-values of the data.
		 * @param y (in) The y-values of the data.
		 */
		void update(size_t nPts, size_t nLines, const double* x, const double* y);

		/**
		 * Updates the plot with new data. The x-values are assumed to be the index of the data.
		 * @param nPts The number of points in the data.
		 * @param nLines The number of lines in the data.
		 * @param y (in) The y-values of the data.
		 * @param ymin The minimum value of the y-axis to display.
		 * @param ymax The maximum value of the y-axis to display.
		 */
		void update(size_t nPts, size_t nLines, const double* y, double ymin, double ymax);

		/**
		 * Updates the plot with new data. The x-values are assumed to be the index of the data.
		 * @param nPts The number of points in the data.
		 * @param nLines The number of lines in the data.
		 * @param y (in) The y-values of the data.
		 */
		void update(size_t nPts, size_t nLines, const double* y);
	};
};

#include "MathTools.tpp"
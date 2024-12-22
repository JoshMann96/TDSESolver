#pragma once
#include <complex>
#include "PhysCon.h"
#include "blas.h"

namespace FDBCs
{
    template <typename T>
    class CyclicArray
    {
		template <typename U>
		friend class CyclicArray; // allow access to other CyclicArray instances' private members when template types differ
    private:
        int size, idx;
        T* arr;
		int localIdx(int i) { return (idx + i) % size; };
    public:
        CyclicArray(int size) : size(size), idx(size-1), arr(static_cast<T *>(sq_malloc(sizeof(T) * size))) {};
		CyclicArray(int size, T val) : CyclicArray(size) { std::fill_n(arr, size, val); };
        ~CyclicArray() { sq_free(arr); };
        void push(T val) { step(); arr[idx] = val; };
		void pushBack(T val) { stepBack(); arr[idx] = val; };
        void set(int i, T val) { arr[localIdx(i)] = val; };
        void step() { idx = (idx + 1) % size; };
		void stepBack() { idx = (idx - 1 + size) % size; };
        T get(int i) { return arr[localIdx(i)]; };
        void print() { for (int i = 0; i < size; i++) std::cout << arr[i] << " "; std::cout << std::endl; };

		void mul(T val) { for (int i = 0; i < size; i++) arr[i] *= val; };

		template <typename U>
		decltype(std::declval<T&>()* std::declval<U&>()) inner(CyclicArray<U>* arr) {
			if (this->size != arr->size)
				throw std::invalid_argument("Arrays must be of the same size.");

			if (this->idx == arr->idx) // aligned arrays lead to faster code
				return this->aligned_inner(arr);
			else{
				decltype(std::declval<T&>()* std::declval<U&>()) sum = 0;
				for (int i = 0; i < this->size; i++)
					sum += this->get(i) * arr->get(i);
				return sum;
			}
		}

		template <typename U>
		decltype(std::declval<T&>()* std::declval<U&>()) aligned_inner(CyclicArray<U>* arr) {
			if (this->size != arr->size)
				throw std::invalid_argument("Arrays must be of the same size.");

			decltype(std::declval<T&>()* std::declval<U&>()) sum = 0;
			if (this->idx == arr->idx) // aligned arrays lead to faster code
				for (int i = 0; i < this->size; i++)
					sum += this->arr[i] * arr->arr[i];
			else
				throw std::invalid_argument("Arrays must be aligned.");
			return sum;
		}

		template <typename U>
		decltype(std::declval<T&>()* std::declval<U&>()) inner(U* arr) {
			decltype(std::declval<T&>()* std::declval<U&>()) sum = 0;
			if constexpr (std::is_same_v<U, T>) //types are the same
				if constexpr (std::is_same_v<T, std::complex<double>>){ // complex double inner product
					std::complex<double> temp;
					cblas_zdotu_sub(this->size - this->idx, this->arr + idx, 1, arr, 1, &temp);
					sum += temp;
					cblas_zdotu_sub(this->idx, this->arr, 1, arr + this->size - this->idx, 1, &temp);
					sum += temp;
				}
			else
				for (int i = 0; i < this->size; i++)
					sum += this->get(i) * arr[i];
			return sum;
		}
    };

	enum class BCSide {
		LEFT, RIGHT
	};

	class BoundaryCondition
	{
	public:
		virtual std::complex<double> getLHSEle() = 0; // LHS matrix element at the boundary position
		virtual std::complex<double> getLHSAdjEle() = 0; // LHS matrix element adjacent to the boundary position
		virtual void getRHS(std::complex<double>* psibd, std::complex<double>* psiad, double vb, std::complex<double>* res, int nElec) = 0; // RHS value for the condition
		virtual void finishStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb) = 0;
		virtual void prepareStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb) = 0;
		virtual void fillHistory(std::complex<double>* psibd, double* kin) = 0;
	};

	/*
	Boundary condition which is the same for all orbitals
	*/
	class CommonBC :
		public BoundaryCondition
	{
	public:
		void getRHS(std::complex<double>* psibd, std::complex<double>* psiad, double vb, std::complex<double>* res, int nElec) {
			for (int i = 0; i < nElec; i++)
				res[i] = getRHS(vb);
		}
		void fillHistory(std::complex<double>* psibd, double* kin) { return; };
		void prepareStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb) {};
		void finishStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb) {};
		virtual std::complex<double> getRHS(double vb) = 0; // RHS value for the condition
	};

	class DirichletBC :
		public CommonBC
	{
	private:
		std::complex<double> bdVal;
	public:
		DirichletBC(std::complex<double> bdVal) : bdVal(bdVal) {};
		std::complex<double> getLHSEle() { return 1.0; };
		std::complex<double> getLHSAdjEle() { return 0.0; };
		std::complex<double> getRHS(double vb) { return bdVal; };
	};

	class NeumannBC :
		public CommonBC
	{
	private:
		std::complex<double> bdDer;
		double dx;
		double direction;
	public:
		NeumannBC(std::complex<double> bdDer, double dx, BCSide side) : bdDer(bdDer), dx(dx), direction(side == BCSide::LEFT ? 1 : -1) {};
		std::complex<double> getLHSEle() { return -1.0 / dx * direction; };
		std::complex<double> getLHSAdjEle() { return 1.0 / dx * direction; };
		std::complex<double> getRHS(double vb) { return bdDer; };
	};

	// Homogeneous Discrete Transparent Boundary Condition
	class HDTransparentBC :
		public BoundaryCondition
	{
	private:
        int order, nElec;
        double dx, dt;
        CyclicArray<std::complex<double>> **psis;
        std::complex<double> kernel0, *kernel;

        void calcKernel(double vb);
	public:

		HDTransparentBC(int order, int nElec, double dx, double dt);
        ~HDTransparentBC();
		std::complex<double> getLHSEle() { return -kernel0; };
		std::complex<double> getLHSAdjEle() { return 1.0; };
		void getRHS(std::complex<double>* psibd, std::complex<double>* psiad, double vb, std::complex<double>* res, int nElec);
        void finishStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb){
			for (int i = 0; i < nElec; i++){
				psis[i]->stepBack();
				psis[i]->mul(std::exp(PhysCon::im/PhysCon::hbar*vb*dt));
			}
		};
		void prepareStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb){
			for (int i = 0; i < nElec; i++){
				psis[i]->set(0, psibd[i]);
			}
			calcKernel(vb);
		};
		void printKernel() { for (int i = 0; i < order; i++) std::cout << kernel[i] << " "; std::cout << std::endl; };
		void fillHistory(std::complex<double>* psibd, double* kin) { 
			for (int i = 0; i < nElec; i++)
				for (int j = 0; j < order; j++)
					psis[i]->set(j, psibd[i]*std::exp(-PhysCon::im/PhysCon::hbar*kin[i]*(j*dt)));
		};
	};

	// Inhomogeneous Discrete Transparent Boundary Condition (incoming current)
	class IDTransparentBC :
		public HDTransparentBC
	{
	};

}
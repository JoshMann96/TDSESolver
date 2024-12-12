#pragma once
#include <complex>
#include "PhysCon.h"

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
    public:
        CyclicArray(int size) : size(size), idx(0), arr((T*)sq_malloc(sizeof(T) * size)) {};
		CyclicArray(int size, T val) : size(size), idx(0), arr((T*)sq_malloc(sizeof(T) * size)) { std::fill_n(arr, size, val); };
        ~CyclicArray() { sq_free(arr); };
        void push(T val) { arr[idx] = val; step(); };
        void set(int i, T val) { arr[(idx + i) % size] = val; };
        void step() { idx = (idx + 1) % size; };
        T get(int i) { return arr[(idx + i) % size]; };
        void print() { for (int i = 0; i < size; i++) std::cout << arr[i] << " "; std::cout << std::endl; };

		void mul(T val) { for (int i = 0; i < size; i++) arr[i] *= val; };

		template <typename U>
		decltype(std::declval<T&>()* std::declval<U&>()) inner(CyclicArray<U>* arr) {
			if (this->size != arr->size)
				throw std::invalid_argument("Arrays must be of the same size.");

			decltype(std::declval<T&>()* std::declval<U&>()) sum = 0;
			if (this->idx == arr->idx) // aligned arrays lead to faster code
				for (int i = 0; i < this->size; i++)
					sum += this->arr[i] * arr->arr[i];
			else
				for (int i = 0; i < this->size; i++)
					sum += this->get(i) * arr->get(i);
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
        void iterate(std::complex<double> psibd, std::complex<double> psiad, double vb) {};
	};

	// Unique Boundary Condition which depends on the individual orbital
	class UniqueBC :
		public BoundaryCondition
	{
	public:
		virtual void finishStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb) = 0;
		virtual void prepareStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb) = 0;
		virtual void fillHistory(std::complex<double> psibd, double kin) = 0; // fill the history assuming an eigenstate with local kinetic energy kin
	};

	// Homogeneous Discrete Transparent Boundary Condition
	class HDTransparentBC :
		public UniqueBC
	{
	private:
        int order, nElec;
        double dx, dt;
        CyclicArray<std::complex<double>> **psis, *kernel;
        std::complex<double> kernel0;

        void calcKernel(double vb);
	public:

		HDTransparentBC(int order, int nElec, double dx, double dt);
        ~HDTransparentBC();
		std::complex<double> getLHSEle() { return -kernel0; };
		std::complex<double> getLHSAdjEle() { return 1.0; };
		void getRHS(std::complex<double>* psibd, std::complex<double>* psiad, double vb, std::complex<double>* res, int nElec);
        void finishStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb){
			for (int i = 0; i < nElec; i++){
				psis[i]->mul(std::exp(PhysCon::im/PhysCon::hbar*vb*dt));
				psis[i]->push(psibd[i]);
			}
		};
		void prepareStep(std::complex<double>* psibd, std::complex<double>* psiad, double vb){
			calcKernel(vb);
			kernel->step();
		};
		void printKernel() { kernel->print(); };
		void fillHistory(std::complex<double> psibd, double kin) { throw std::runtime_error("HDTransparentBC::fillHistory not implemented."); };
	};

	// Inhomogeneous Discrete Transparent Boundary Condition (incoming current)
	class IDTransparentBC :
		public HDTransparentBC
	{
	};

	inline bool isUnique(BoundaryCondition* bc) {return dynamic_cast<UniqueBC*>(bc) != nullptr;}
	inline bool isCommon(BoundaryCondition* bc) {return dynamic_cast<CommonBC*>(bc) != nullptr;}

}
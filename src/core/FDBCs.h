#pragma once
#include <complex>

namespace FD_BCs
{
    template <typename T>
    class CyclicArray
    {
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

		template <typename U, typename V>
		static decltype(std::declval<U&>()* std::declval<V&>()) inner(CyclicArray<U>* arr1, CyclicArray<V>* arr2) {
			if (arr1->size != arr2->size)
				throw std::invalid_argument("Arrays must be of the same size.");

			decltype(std::declval<U&>()* std::declval<V&>()) sum = 0;
			if (arr1->idx == arr2->idx) // aligned arrays lead to faster code
				for (int i = 0; i < arr1->size; i++)
					sum += arr1->arr[i] * arr2->arr[i];
			else
				for (int i = 0; i < arr1->size; i++)
					sum += arr1->get(i) * arr2->get(i);
			return sum;
		}
    };

	enum class BC_Side {
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
		NeumannBC(std::complex<double> bdDer, double dx, BC_Side side) : bdDer(bdDer), dx(dx), direction(side == BC_Side::LEFT ? 1 : -1) {};
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
	};

	// Inhomogeneous Discrete Transparent Boundary Condition (incoming current)
	class IDTransparentBC :
		public HDTransparentBC
	{
	};
}
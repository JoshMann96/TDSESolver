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
        ~CyclicArray() { sq_free(arr); };
        void push(T val) { arr[idx] = val; step(); };
        void set(int i, T val) { arr[(idx + i) % size] = val; };
        void step() { idx = (idx + 1) % size; };
        T get(int i) { return arr[(idx + i) % size]; };
        void print() { for (int i = 0; i < size; i++) std::cout << get(i) << " "; std::cout << std::endl; };
    };

	enum class BC_Side {
		LEFT, RIGHT
	};

	/*
	Calculates LHS matrix elements and RHS value for a boundary condition.
	*/
	class FD_BC 
	{
	public:
		virtual std::complex<double> getLHSEle(std::complex<double> psibd, std::complex<double> psiad, double vb) = 0; // LHS matrix element at the boundary position
		virtual std::complex<double> getLHSAdjEle(std::complex<double> psibd, std::complex<double> psiad, double vb) = 0; // LHS matrix element adjacent to the boundary position
		virtual std::complex<double> getRHS(std::complex<double> psibd, std::complex<double> psiad, double vb) = 0; // RHS value for the condition
        virtual void iterate(std::complex<double> psibd, std::complex<double> psiad, double vb) = 0;
	};

	class DirichletBC :
		public FD_BC
	{
	private:
		std::complex<double> bdVal;
	public:
		DirichletBC(std::complex<double> bdVal) : bdVal(bdVal) {};
		std::complex<double> getLHSEle(std::complex<double> psibd, std::complex<double> psiad, double vb) { return 1.0; };
		std::complex<double> getLHSAdjEle(std::complex<double> psibd, std::complex<double> psiad, double vb) { return 0.0; };
		std::complex<double> getRHS(std::complex<double> psibd, std::complex<double> psiad, double vb) { return bdVal; };
        void iterate(std::complex<double> psibd, std::complex<double> psiad, double vb) {};
	};

	class NeumannBC :
		public FD_BC
	{
	private:
		std::complex<double> bdDer;
		double dx;
		double direction;
	public:
		NeumannBC(std::complex<double> bdDer, double dx, BC_Side side) : bdDer(bdDer), dx(dx), direction(side == BC_Side::LEFT ? 1 : -1) {};
		std::complex<double> getLHSEle(std::complex<double> psibd, std::complex<double> psiad, double vb) { return -1.0 / dx * direction; };
		std::complex<double> getLHSAdjEle(std::complex<double> psibd, std::complex<double> psiad, double vb) { return 1.0 / dx * direction; };
		std::complex<double> getRHS(std::complex<double> psibd, std::complex<double> psiad, double vb) { return bdDer; };
        void iterate(std::complex<double> psibd, std::complex<double> psiad, double vb) {};
	};

	// Homogeneous Discrete Transparent Boundary Condition
	class HDTransparentBC :
		public FD_BC
	{
	private:
        int order;
        double dx, dt;
        CyclicArray<std::complex<double>> *psis, *kernel;
        std::complex<double> kernel0;

	public:

        void calcKernel(double vb);

		HDTransparentBC(int order, double dx, double dt);
        ~HDTransparentBC();
		std::complex<double> getLHSEle(std::complex<double> psibd, std::complex<double> psiad, double vb) { return 1.0; };
		std::complex<double> getLHSAdjEle(std::complex<double> psibd, std::complex<double> psiad, double vb) { return 1.0; };
		std::complex<double> getRHS(std::complex<double> psibd, std::complex<double> psiad, double vb) { return 0.0; };
        void iterate(std::complex<double> psibd, std::complex<double> psiad, double vb){ kernel->step(); };
        void printKernel(){kernel->print();}
	};

	// Inhomogeneous Discrete Transparent Boundary Condition (incoming current)
	class IDTransparentBC :
		public FD_BC
	{
	};
}
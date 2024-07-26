#include "Measurers.h"
#include "KineticOperator.h"
#include "PhysCon.h"
#include <iostream>
#include <stdexcept>
#include <system_error>
#include "blas_externs.h"


namespace Measurers {

	std::fstream openFile(const char* fil) {
		return std::fstream(fil, std::ios::out | std::ios::binary);
	}

	DoubleConst::DoubleConst(double c, const char* filName, const char* fol){
		DoubleConst::c = c;
		const char *ext = ".dat";
		int l1 = std::strlen(fol), l2 = std::strlen(filName), l3 = std::strlen(ext);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], filName, l2);
		strcpy(&nfil[l1 + l2], ext);
		/*std::stringstream ss;
		ss << fol << filName << ".dat";
		std::string str = ss.str();
		const char* nfil = str.c_str();*/

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&c), sizeof c);
		fil.close();
	}

	DoubleConst::~DoubleConst() {
		kill();
	}

	int DoubleConst::measure(int step, std::complex<double>* psi, double* v, double t) { return 1; }

	void DoubleConst::terminate() {}


	NElec::NElec(int* nElec, const char* fol)
		: nElec(nElec) {
			int l1 = std::strlen(fol), l2 = std::strlen(fname);
			nfil = new char[l1 + l2 + 1];
			strncpy(nfil, fol, l1);
			strcpy(&nfil[l1], fname);
		}

	NElec::~NElec() {
		kill();
		if(nfil)
			delete[] nfil; nfil = nullptr;
	}

	int NElec::measure(int step, std::complex<double>* psi, double* v, double t) { 
		if(first){
			first = 0;
			
			fil = openFile(nfil);
			delete[] nfil; nfil = nullptr;
			if (!fil) {
				std::cout << "Could not open file: " << nfil;
				std::cin.ignore();
			}

			fil.write(reinterpret_cast<char*>(&index), sizeof(int));
			fil.write(reinterpret_cast<char*>(nElec), sizeof(int));
		}
		
		return 1; 
	}

	void NElec::terminate() {fil.close();}

	//Header
	Header::Header(const char* title, const char* fol) {
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));

		fil.write(title, sizeof(title));

		fil.close();
	}

	Header::~Header() {
		kill();
	}

	int Header::measure(int step, std::complex<double>* psi, double* v, double t) { return 1; }

	void Header::terminate() {}

	//NPts
	NPts::NPts(int nPts, const char* fol) {
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nPts), sizeof(int));

		fil.close();
	}

	NPts::~NPts() {
		kill();
	}

	void NPts::terminate() {}

	int NPts::measure(int step, std::complex<double>* psi, double* v, double t) { return 1; }

	//NSteps
	NSteps::NSteps(int numSteps, const char* fol) {
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&numSteps), sizeof(int));
		fil.close();
	}

	NSteps::~NSteps() {kill();}

	void NSteps::terminate() {}

	int NSteps::measure(int step, std::complex<double> * psi, double * v, double t) {return 0;}

	//DX
	DX::DX(double dx, const char* fol) {
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&dx), sizeof(double));

		fil.close();
	}

	DX::~DX() {
		kill();
	}

	int DX::measure(int step, std::complex<double>* psi, double* v, double t) { return 1; }

	void DX::terminate() {}

	//DT
	DT::DT(double dt, const char* fol) {
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&dt), sizeof(double));

		fil.close();
	}

	DT::~DT() {
		kill();
	}

	int DT::measure(int step, std::complex<double>* psi, double* v, double t) { return 1; }

	void DT::terminate() {}

	//XS
	XS::XS(int len, double* xs, const char* fol) {
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*len);

		fil.close();
	}

	XS::~XS() {
		kill();
	}

	int XS::measure(int step, std::complex<double>* psi, double* v, double t) { return 1; }

	void XS::terminate() {}


	//TS
	TS::TS(const char* fol) {
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	TS::~TS() {
		kill();
	}

	int TS::measure(int step, std::complex<double> * psi, double * v, double t) {
		fil.write(reinterpret_cast<char*>(&t), sizeof(double));
		return 0;
	}

	void TS::terminate() {
		fil.close();
	}


	OrigPot::OrigPot(int n, const char* fol) {
		OrigPot::n = n;
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);
		//std::stringstream ss;
		//ss << fol << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&n), sizeof(int));
	}

	OrigPot::~OrigPot() {
		kill();
	}

	int OrigPot::measure(int step, std::complex<double> * psi, double * v, double t) {
		fil.write(reinterpret_cast<char*>(&v[0]), sizeof(double)*n);
		return 1;
	}

	void OrigPot::terminate() {fil.close();}

	//Psi2t
	Psi2t::Psi2t(int nPts, int nx, int nt, int numSteps, double maxT, double * x, int* nElec, const char* fol) :
		nPts(nPts), nx(nx), nt(nt), numSteps(numSteps), nElec(nElec), curIdx(0)
	{
		measSteps = (int*) sq_malloc(sizeof(int)*numSteps);
		vtls::linspace(nt, 0, numSteps - 1, measSteps);

		xs = (double*) sq_malloc(sizeof(double)*nx);
		vtls::downSampleLinearInterpolateEdge(nPts, x, nx, xs);
		ts = (double*) sq_malloc(sizeof(double)*nt);
		psi2b = (double*) sq_malloc(sizeof(double)*nPts);
		psi2s = (double*) sq_malloc(sizeof(double)*nx);

		vtls::linspace(nt, 0.0, maxT, ts);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nx), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nt), sizeof(int));
	}

	Psi2t::~Psi2t() {
		kill();

		sq_free(psi2b);
		sq_free(psi2s);
		sq_free(xs);
		sq_free(ts);
		sq_free(measSteps);
	}

	int Psi2t::measure(int step, std::complex<double> * psi, double * v, double t) {
		while(step == measSteps[curIdx]){
			for(int i = 0; i < *nElec; i++){
				vtls::normSqr(nPts, &psi[i*nPts], psi2b);
				vtls::downSampleLinearInterpolateEdge(nPts, psi2b, nx, psi2s);
				fil.write(reinterpret_cast<char*>(&psi2s[0]), sizeof(double)*nx);
			}

			curIdx++;
			if(curIdx >= nt)
				return 1;
		}
		return 0;
	}

	void Psi2t::terminate() {
		fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*nx);
		fil.write(reinterpret_cast<char*>(&ts[0]), sizeof(double)*nt);
		fil.close();
	}

	ExpectE::ExpectE(int nPts, double dx, int* nElec, const char* fol, KineticOperators::KineticOperator ** kin) : nPts(nPts), dx(dx), kin(kin), nElec(nElec) {
		rho = (double*) sq_malloc(sizeof(double)*nPts);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	ExpectE::~ExpectE() {
		kill();

		sq_free(rho);
	}

	int ExpectE::measure(int step, std::complex<double> * psi, double * v, double t) {
		double ex;
		for(int i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], rho);
			ex = vtlsInt::rSumMul(nPts, rho, v, dx) + (*kin)->evaluateKineticEnergy(&psi[i*nPts]);

			fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		}

		return 0;
	}

	void ExpectE::terminate() {
		fil.close();
	}


	ExpectX::ExpectX(int nPts, double* xs, double dx, int* nElec, const char* fol) :
		nPts(nPts), dx(dx), nElec(nElec), x(xs)
	 {
		scratch = (double*) sq_malloc(sizeof(double)*nPts);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	ExpectX::~ExpectX() {
		kill();

		sq_free(scratch);
	}

	int ExpectX::measure(int step, std::complex<double> * psi, double * v, double t) {
		for(int i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], scratch);
			double ex = vtlsInt::simpsMul(nPts, x, scratch, dx);
			fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		}
		return 0;
	}

	void ExpectX::terminate() {
		fil.close();
	}


	ExpectP::ExpectP(int len, double dx, int* nElec, const char* fol) :
		nPts(len), dx(dx), nElec(nElec)
	 {
		scratch1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*len);
		scratch2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*len);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	ExpectP::~ExpectP() {
		kill();

		sq_free(scratch1);
		sq_free(scratch2);
	}

	int ExpectP::measure(int step, std::complex<double> * psi, double * v, double t) {
		double ex;
		for(int i = 0; i < *nElec; i++){
			vtls::firstDerivative(nPts, &psi[i*nPts], scratch1, dx);
			for (int j = 0; j < nPts; j++)
				scratch2[j] = std::conj(psi[i*nPts + j]);
			ex = std::imag(vtlsInt::simpsMul(nPts, scratch2, scratch1, dx))*PhysCon::hbar;
			fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		}
		return 0;
	}

	void ExpectP::terminate() {
		fil.close();
	}


	ExpectA::ExpectA(int nPts, double dx, int* nElec, const char* fol) :
		nPts(nPts), dx(dx), nElec(nElec)
		 {
		scratch1 = (double*) sq_malloc(sizeof(double)*nPts);
		scratch2 = (double*) sq_malloc(sizeof(double)*nPts);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	ExpectA::~ExpectA() {
		kill();

		sq_free(scratch1);
		sq_free(scratch2);
	}

	int ExpectA::measure(int step, std::complex<double> * psi, double * v, double t) {
		double ex;
		vtls::firstDerivative(nPts, v, scratch1, dx);
		for(int i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], scratch2);
			ex = vtlsInt::simpsMul(nPts, scratch2, scratch1, dx)*(-1.0 / PhysCon::me);
			fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		}

		return 0;
	}

	void ExpectA::terminate() {
		fil.close();
	}


	//TotProb
	TotProb::TotProb(int nPts, double dx, int* nElec, const char* fol) :
		nPts(nPts), dx(dx), nElec(nElec)
	{
		psi2 = (double*) sq_malloc(sizeof(double)*nPts);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	TotProb::~TotProb() {
		kill();

		sq_free(psi2);
	}

	int TotProb::measure(int step, std::complex<double> * psi, double * v, double t) {
		double sum;
		for(int i = 0; i < *nElec; i++){	
			vtls::normSqr(nPts, &psi[i*nPts], psi2);
			sum = vtlsInt::simps(nPts, psi2, dx);
			fil.write(reinterpret_cast<char*>(&sum), sizeof(double));
		}
		return 0;
	}

	void TotProb::terminate() {
		fil.close();
	}


	VDProbCurrent::VDProbCurrent(int nPts, double dx, int* nElec, int vdPos, int vdNum, const char* name, const char* fol) :
		nPts(nPts), dx(dx), vdPos(vdPos), vdNum(vdNum), nElec(nElec)
	 {
		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1+l2], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&vdPos), sizeof(int));
	}

	VDProbCurrent::~VDProbCurrent() {
		kill();
	}

	int VDProbCurrent::measure(int step, std::complex<double> * psi, double * v, double t) {
		std::complex<double> der;
		double j;
		for(int i = 0; i < *nElec; i++){
			der = vtls::firstDerivative(nPts, &psi[i*nPts], vdPos, dx);
			j = std::imag(PhysCon::hbar / (2.0 * PhysCon::me)*(std::conj(psi[i*nPts + vdPos])*der - psi[i*nPts + vdPos] * std::conj(der)));
			fil.write(reinterpret_cast<char*>(&j), sizeof(double));
		}

		return 0;
	}

	void VDProbCurrent::terminate() {
		fil.close();
	}

	PsiT::PsiT(int nPts, double meaT, int *nElec, int vdNum, const char* name, const char* fol) :
		nElec(nElec), vdNum(vdNum), meaT(meaT), nPts(nPts)
	 {
		PsiT::meaT = meaT;
		PsiT::vdNum = vdNum;
		PsiT::nPts = nPts;

		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1 + l2], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&meaT), sizeof(double));
	}

	PsiT::~PsiT() {
		kill();
	}

	int PsiT::measure(int step, std::complex<double> * psi, double * v, double t) {
		if ((!done && t >= meaT)) {
			for(int i = 0; i < *nElec; i++)
				fil.write(reinterpret_cast<char*>(&psi[i*nPts]), sizeof(std::complex<double>)* nPts);
			done = 1;
			return 1;
		}
		
		return 0;

	}

	void PsiT::terminate() {
		fil.close();
	}


	PotT::PotT(int n, double meaT, int vdNum, const char* name, const char* fol) {
		PotT::meaT = meaT;
		PotT::vdNum = vdNum;
		PotT::n = n;

		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1 + l2], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&meaT), sizeof(double));
	}

	PotT::~PotT() {
		kill();
	}

	int PotT::measure(int step, std::complex<double> * psi, double * v, double t) {
		if (!done && t >= meaT) {
			fil.write(reinterpret_cast<char*>(v), sizeof(double)*n);
			done = 1;
			return 1;
		}
		return 0;
	}

	void PotT::terminate() {
		fil.close();
	}

	VDPsi::VDPsi(int* nElec, int vdPos, int vdNum, const char* name, const char* fol) : 
		nElec(nElec), vdPos(vdPos), vdNum(vdNum)
	{
		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1 + l2], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&vdPos), sizeof(int));
	}

	VDPsi::~VDPsi() {
		kill();
	}

	int VDPsi::measure(int step, std::complex<double> * psi, double * v, double t) {
		for(int i = 0; i < *nElec; i++)
			fil.write(reinterpret_cast<char*>(&psi[i*nPts + vdPos]), sizeof(std::complex<double>));
		return 0;
	}

	void VDPsi::terminate() {
		fil.close();
	}

	VDPot::VDPot(int vdPos, int vdNum, const char* name, const char* fol) :
		vdPos(vdPos), vdNum(vdNum)
	 {
		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1 + l2], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&vdPos), sizeof(int));
	}

	VDPot::~VDPot() {
		kill();
	}

	int VDPot::measure(int step, std::complex<double> * psi, double * v, double t) {
		fil.write(reinterpret_cast<char*>(&v[vdPos]), sizeof(double));
		return 0;
	}

	void VDPot::terminate() {
		fil.close();
	}


	VDFluxSpec::VDFluxSpec(int nPts, int vdPos, int vdNum, int* nElec, int nSamp, double emax, double tmax, const char* name, const char* fol) :
		vdPos(vdPos), vdNum(vdNum), nElec(nElec), nSamp(nSamp), tmax(tmax), nPts(nPts), dw(emax / PhysCon::hbar / nSamp)
	 {
		phss = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		temp = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);

		phaseCalcExpMul = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		for(int i = 0; i < nSamp; i++)
			phaseCalcExpMul[i] = PhysCon::im * dw * (double)i; //to be multiplied by t then exponentiated later
		
		cumPotPhs = 1;

		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1 + l2], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&vdPos), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nSamp), sizeof(int));
		fil.write(reinterpret_cast<char*>(&emax), sizeof(double));
	}

	VDFluxSpec::~VDFluxSpec() {
		kill();

		if(wfcs0)
			sq_free(wfcs0); wfcs0 = nullptr;
		if(wfcs1)
			sq_free(wfcs1); wfcs1 = nullptr;
		sq_free(phss);
		sq_free(phaseCalcExpMul);
		sq_free(temp);
	}

	int VDFluxSpec::measure(int step, std::complex<double>* psi, double* v, double t) {
		if (first) {
			if(wfcs0)
				sq_free(wfcs0);
			if(wfcs1)
				sq_free(wfcs1);
			wfcs0 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp * *nElec);
			wfcs1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp * *nElec);

			for (int i = 0; i < nSamp * *nElec; i++) {
				wfcs0[i] = 0;
				wfcs1[i] = 0;
			}
			
			first = 0;
			ct = t;
		}

		cumPotPhs *= std::exp(PhysCon::im * (v[vdPos] + v[vdPos + 1]) / (2.0 * PhysCon::hbar) * (t - ct));
		double winMul;
		if (t < tukeyAl / 2 * tmax)
			winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * t / (tukeyAl * tmax)));
		else if (t > (1.0 - tukeyAl / 2) * tmax)
			winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * (tmax - t) / (tukeyAl * tmax)));
		else
			winMul = 1.0;

		//exp(i t dw (idx))*cumPotPhs*winMul, expanded to hopefully vectorize better
		cblas_zcopy(nSamp, phaseCalcExpMul, 1, phss, 1); 	// phss = 		i dw (idx)
		cblas_zdscal(nSamp, t, phss, 1); 					// phss = 		i dw (idx) t
		for(int i = 0; i < nSamp; i++)
			phss[i] = std::exp(phss[i]);					// phss = exp(	i dw (idx) t)
		std::complex<double> cpwm = cumPotPhs * winMul;
		cblas_zscal(nSamp, &cpwm, phss, 1);
			
		for(int i = 0; i < *nElec; i++){
			//wfcs0[i0 + i] += psip0 * phss[i]
			//wfcs1[i0 + i] += psip1 * phss[i]
			cblas_zaxpy(nSamp, &psi[i*nPts + vdPos  ], phss, 1, &wfcs0[i*nSamp], 1);
			cblas_zaxpy(nSamp, &psi[i*nPts + vdPos+1], phss, 1, &wfcs1[i*nSamp], 1);
		}

		ct = t;

		return 0;
	}

	void VDFluxSpec::terminate() {
		fil.write(reinterpret_cast<char*>(&wfcs0[0]), *nElec * nSamp * sizeof(std::complex<double>));
		fil.write(reinterpret_cast<char*>(&wfcs1[0]), *nElec * nSamp * sizeof(std::complex<double>));
		fil.close();
	}



	Vfunct::Vfunct(int nPts, int nx, int nt, int numSteps, double maxT, double * x, const char* fol) :
		nPts(nPts), nx(nx), nt(nt), maxT(maxT), curIdx(0)
	{
		measSteps = (int*) sq_malloc(sizeof(int)*numSteps);
		vtls::linspace(nt, 0, numSteps - 1, measSteps);

		xs = (double*) sq_malloc(sizeof(double)*nx);
		vtls::downSampleLinearInterpolateEdge(nPts, x, nx, xs);
		ts = (double*) sq_malloc(sizeof(double)*nt);
		vs = (double*) sq_malloc(sizeof(double)*nx);

		vtls::linspace(nt, 0.0, maxT, ts);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nx), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nt), sizeof(int));
	}

	Vfunct::~Vfunct() {
		kill();

		sq_free(vs);
		sq_free(xs);
		sq_free(ts);
		sq_free(measSteps);
	}

	int Vfunct::measure(int step, std::complex<double> * psi, double * v, double t) {
		while(step == measSteps[curIdx]){
			vtls::downSampleLinearInterpolateEdge(nPts, v, nx, vs);
			fil.write(reinterpret_cast<char*>(&vs[0]), sizeof(double)*nx);
			
			curIdx++;
			if(curIdx >= nt)
				return 1;
		}
		return 0;
	}

	void Vfunct::terminate() {
		fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*nx);
		fil.write(reinterpret_cast<char*>(&ts[0]), sizeof(double)*nt);
		fil.close();
	}

	ExpectE0::ExpectE0(int nPts, double dx, int* nElec, const char* fol, KineticOperators::KineticOperator ** kin) : 
		nPts(nPts), dx(dx), kin(kin), nElec(nElec)
	{
		rho = (double*) sq_malloc(sizeof(double)*nPts);

		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	ExpectE0::~ExpectE0() {
		kill();

		sq_free(rho);
	}

	int ExpectE0::measure(int step, std::complex<double>* psi, double* v, double t) {
		if(first){
			first = 0;
			double ex;
			for(int i = 0; i < *nElec; i++){
				vtls::normSqr(nPts, &psi[i*nPts], rho);
				ex = vtlsInt::rSumMul(nPts, rho, v, dx) + (*kin)->evaluateKineticEnergy(&psi[i*nPts]);
				fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
			}
		}
		return 1;
	}

	void ExpectE0::terminate() {
		fil.close();
	}


	WfcRhoWeights::WfcRhoWeights(int* nElec, double** weights, const char* fol) : 
		nElec(nElec), weights(weights)
	{
		int l1 = std::strlen(fol), l2 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + 1];
		strncpy(nfil, fol, l1);
		strcpy(&nfil[l1], fname);

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
	}

	WfcRhoWeights::~WfcRhoWeights() {
		kill();
	}

	int WfcRhoWeights::measure(int step, std::complex<double>* psi, double* v, double t) {
		if (first) {
			if (*weights){
				first = 0;
				fil.write(reinterpret_cast<char*>(nElec), sizeof(int));
				fil.write(reinterpret_cast<char*>(*weights), sizeof(double)* *nElec);
			}
			else{
				throw std::runtime_error("Weights not set for WfcRhoWeights.");
			}
		}
		return 1;
	}

	void WfcRhoWeights::terminate() {
		fil.close();
	}


	BasicMeasurers::BasicMeasurers(int nPts, int numSteps, double dx, double dt, double * xs, const char* fol) {
		meas.push_back(new NPts(nPts, fol));
		//meas.push_back(new Header(title_4char, fol));
		meas.push_back(new NSteps(numSteps, fol));
		meas.push_back(new DX(dx, fol));
		meas.push_back(new DT(dt, fol));
		//meas.push_back(new XS(nPts, xs, fol));
		//meas.push_back(new TS(fol));
	}

	BasicMeasurers::~BasicMeasurers() { kill(); }

	int BasicMeasurers::measure(int step, std::complex<double> * psi, double * v, double t) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(step, psi, v, t) == 1) {
				(*it)->kill();
				it = meas.erase(it);
			}
			else
				++it;
		}
		return 0;
	}

	void BasicMeasurers::terminate() {
		for(Measurer* m : meas)
			delete m;
		meas.clear();
	}


	MeasurementManager::MeasurementManager(const char* fname) {
		MeasurementManager::fname = fname;
	}

	MeasurementManager::~MeasurementManager() {
		kill();
	}

	void MeasurementManager::addMeasurer(Measurer * m) {
		meas.push_back(m);
	}

	int MeasurementManager::measure(int step, std::complex<double> * psi, double * v, double t) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(step, psi, v, t) == 1) {
				(*it)->kill();
				it = meas.erase(it);
			}
			else
				++it;
		}
		//std::cout << std::flush;
		return 0;
	}

	void MeasurementManager::terminate() {
		for(Measurer* m : meas)
			m->kill();
		meas.clear();
	}
}
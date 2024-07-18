#include "Measurers.h"
#include <iostream>
#include <stdexcept>
#include <system_error>
#include "blas_externs.h"


namespace Measurers {

	std::fstream openFile(const char* fil) {
		return std::fstream(fil, std::ios::out | std::ios::binary);
	}

	DoubleConst::DoubleConst(double c, const char* filName, const char* fol) {
		DoubleConst::c = c;
		int l1 = std::strlen(fol), l2 = std::strlen(filName);
		char* nfil = new char[l1 + l2 + 4];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], filName, l2);
		strcpy(&nfil[l1 + l2], ".dat");
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

	int DoubleConst::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

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

	int NElec::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { 
		if(first){
			first = 0;
			
			fil = openFile(nfil);
			delete[] nfil; nfil = nullptr;
			if (!fil) {
				std::cout << "Could not open file: " << nfil;
				std::cin.ignore();
			}

			fil.write(reinterpret_cast<char*>(&index), sizeof(int));
			fil.write(reinterpret_cast<char*>(&(*nElec)), sizeof(int));
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

	int Header::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

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

	int NPts::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

	//NSteps
	NSteps::NSteps(const char* fol) {
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

	NSteps::~NSteps() {kill();}

	void NSteps::terminate() {
		fil.write(reinterpret_cast<char*>(&steps), sizeof(int));
		fil.close();
	}

	int NSteps::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		if (t != tmea) {
			steps++;
			tmea = t;
		}
		return 0;
	}

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

	int DX::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

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

	int DT::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

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

	int XS::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

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

	int TS::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		if (t != tmea) {
			fil.write(reinterpret_cast<char*>(&t), sizeof(double));
			tmea = t;
		}
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

	int OrigPot::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		if (!measured) {
			fil.write(reinterpret_cast<char*>(&v[0]), sizeof(double)*n);
			measured = 1;
		}
		return 1;
	}

	void OrigPot::terminate() {fil.close();}

	//Psi2t
	Psi2t::Psi2t(int n, int nx, int nt, double maxT, double dt, double * x, const char* fol) {
		Psi2t::n = n;
		Psi2t::nx = nx;
		Psi2t::nt = nt;
		Psi2t::maxT = maxT;
		Psi2t::dt = dt;

		interval = maxT / (double)nt;

		xs = new double[nx];
		vtls::downSampleLinearInterpolateEdge(n, x, nx, xs);
		ts = new double[nt];
		psi2b = new double[n];
		psi2s = new double[nx];

		mulT = interval;

		curIts = 0;

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
		fil.write(reinterpret_cast<char*>(&nx), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nt), sizeof(int));
	}

	Psi2t::~Psi2t() {
		kill();

		delete[] psi2b;
		delete[] psi2s;
		delete[] xs;
		delete[] ts;
	}

	int Psi2t::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		int k = (int)(t / interval) - curIts;
		if ((k == 0 || std::abs(t-mulT)<dt/8) && curIts<=nt) {
			if (std::abs(t - mulT) >= dt / 8 && curIts == nt)
				curIts++;
			else {
				vtls::normSqr(n, psi, psi2b);
				vtls::downSampleLinearInterpolateEdge(n, psi2b, nx, psi2s);

				fil.write(reinterpret_cast<char*>(&psi2s[0]), sizeof(double)*nx);
				if (std::abs(t - mulT) >= dt / 8) {
					ts[curIts] = t;
					curIts++;
					mulT = t;
				}
			}
		}
		else if (k > 1) {
			std::cout << "Downsampled measurers were not able to keep up with time step (there needs to be more time steps in the simulation than sampled)." << std::endl;
			throw -1;
		}
		return 0;
	}

	void Psi2t::terminate() {
		fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*nx);
		fil.write(reinterpret_cast<char*>(&ts[0]), sizeof(double)*nt);
		fil.close();
	}

	ExpectE::ExpectE(int len, double dx, const char* fol) {
		nPts = len;
		ExpectE::dx = dx;
		rho = new double[nPts];

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

	ExpectE::~ExpectE() {
		kill();

		delete[] rho;
	}

	int ExpectE::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		/*vtls::secondDerivative(nPts, psi, scratch1, dx);
		vtls::scaMulArray(nPts, -PhysCon::hbar*PhysCon::hbar / (2.0*PhysCon::me), scratch1);
		vtls::seqMulArrays(nPts, v, psi, scratch2);
		vtls::addArrays(nPts, scratch2, scratch1);
		for (int i = 0; i < nPts; i++)
			scratch2[i] = std::conj(psi[i]);
		double ex = std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));*/

		vtls::normSqr(nPts, psi, rho);
		double ex = vtlsInt::rSumMul(nPts, rho, v, dx) + kin->evaluateKineticEnergy(psi);

		fil.write(reinterpret_cast<char*>(&ex), sizeof(double));

		return 0;
	}

	void ExpectE::terminate() {
		fil.close();
	}


	ExpectX::ExpectX(int len, double* xs, double dx, const char* fol) {
		nPts = len;
		x = xs;
		ExpectX::dx = dx;
		scratch = new double[len];

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

	ExpectX::~ExpectX() {
		kill();

		delete[] scratch;
	}

	int ExpectX::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::normSqr(nPts, psi, scratch);
		double ex = vtlsInt::simpsMul(nPts, x, scratch, dx);
		fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		return 0;
	}

	void ExpectX::terminate() {
		fil.close();
	}


	ExpectP::ExpectP(int len, double dx, const char* fol) {
		nPts = len;
		ExpectP::dx = dx;
		scratch1 = new std::complex<double>[len];
		scratch2 = new std::complex<double>[len];

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

	ExpectP::~ExpectP() {
		kill();

		delete[] scratch1;
		delete[] scratch2;
	}

	int ExpectP::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::firstDerivative(nPts, psi, scratch1, dx);
		for (int i = 0; i < nPts; i++)
			scratch2[i] = std::conj(psi[i]);
		double ex = std::imag(vtlsInt::simpsMul(nPts, scratch2, scratch1, dx))*PhysCon::hbar;
		fil.write(reinterpret_cast<char*>(&ex), sizeof(double));

		return 0;
	}

	void ExpectP::terminate() {
		fil.close();
	}


	ExpectA::ExpectA(int len, double dx, const char* fol) {
		nPts = len;
		ExpectA::dx = dx;
		scratch1 = new double[len];
		scratch2 = new double[len];

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

	ExpectA::~ExpectA() {
		kill();

		delete[] scratch1;
		delete [] scratch2;
	}

	int ExpectA::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::firstDerivative(nPts, v, scratch1, dx);
		vtls::normSqr(nPts, psi, scratch2);
		double ex = vtlsInt::simpsMul(nPts, scratch2, scratch1, dx)*(-1.0 / PhysCon::me);
		fil.write(reinterpret_cast<char*>(&ex), sizeof(double));

		return 0;
	}

	void ExpectA::terminate() {
		fil.close();
	}


	//TotProb
	TotProb::TotProb(int n, double dx, const char* fol) {
		TotProb::n = n;
		TotProb::dx = dx;
		psi2 = new double[n];

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

	TotProb::~TotProb() {
		kill();

		delete[] psi2;
	}

	int TotProb::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::normSqr(n, psi, psi2);
		double sum = vtlsInt::simps(n, psi2, dx);
		fil.write(reinterpret_cast<char*>(&sum), sizeof(double));

		return 0;
	}

	void TotProb::terminate() {
		fil.close();
	}


	VDProbCurrent::VDProbCurrent(int n, double dx, int vdPos, int vdNum, const char* name, const char* fol) {
		VDProbCurrent::n = n;
		VDProbCurrent::dx = dx;
		VDProbCurrent::pos = vdPos;
		VDProbCurrent::vdNum = vdNum;

		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1+l2], fname);
		//std::stringstream ss;
		//ss << fol << vdNum << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&pos), sizeof(int));
	}

	VDProbCurrent::~VDProbCurrent() {
		kill();
	}

	int VDProbCurrent::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		std::complex<double> der = vtls::firstDerivative(n, psi, pos, dx);
		double j = std::imag(PhysCon::hbar / (2.0 * PhysCon::me)*(std::conj(psi[pos])*der - psi[pos] * std::conj(der)));
		fil.write(reinterpret_cast<char*>(&j), sizeof(double));

		return 0;
	}

	void VDProbCurrent::terminate() {
		fil.close();
	}

	PsiT::PsiT(int n, double meaT, int vdNum, const char* name, const char* fol) {
		PsiT::meaT = meaT;
		PsiT::vdNum = vdNum;
		PsiT::n = n;

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

	int PsiT::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		if ((!done && t >= meaT) || t == curTime) {
			fil.write(reinterpret_cast<char*>(psi), sizeof(std::complex<double>)* n);
			done = 1;
			curTime = t;
		}
		else if ((done && t >= meaT) && t != curTime)
			return 1;
		
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

	int PotT::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
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

	VDPsi::VDPsi(int vdPos, int vdNum, const char* name, const char* fol) {
		VDPsi::pos = vdPos;
		VDPsi::vdNum = vdNum;

		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1 + l2], fname);
		//std::stringstream ss;
		//ss << fol << vdNum << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&pos), sizeof(int));
	}

	VDPsi::~VDPsi() {
		kill();
	}

	int VDPsi::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		fil.write(reinterpret_cast<char*>(&psi[pos]), sizeof(std::complex<double>));
		return 0;
	}

	void VDPsi::terminate() {
		fil.close();
	}

	VDPot::VDPot(int vdPos, int vdNum, const char* name, const char* fol) {
		VDPot::pos = vdPos;
		VDPot::vdNum = vdNum;

		std::string tempstr = std::to_string(vdNum);
		const char* nm = tempstr.c_str();
		int l1 = std::strlen(fol), l2 = std::strlen(nm), l3 = std::strlen(fname);
		char* nfil = new char[l1 + l2 + l3 + 1];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], nm, l2);
		strcpy(&nfil[l1 + l2], fname);
		//std::stringstream ss;
		//ss << fol << vdNum << fname;
		//std::string str = ss.str();
		//const char* nfil = str.c_str();

		fil = openFile(nfil);
		delete[] nfil;
		if (!fil) {
			std::cout << "Could not open file: " << nfil;
			std::cin.ignore();
		}

		fil.write(reinterpret_cast<char*>(&index), sizeof(int));
		fil.write(reinterpret_cast<char*>(&vdNum), sizeof(int));
		fil.write(reinterpret_cast<char*>(&name), 4);
		fil.write(reinterpret_cast<char*>(&pos), sizeof(int));
	}

	VDPot::~VDPot() {
		kill();
	}

	int VDPot::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		fil.write(reinterpret_cast<char*>(&v[pos]), sizeof(double));
		return 0;
	}

	void VDPot::terminate() {
		fil.close();
	}


	VDFluxSpec::VDFluxSpec(int vdPos, int vdNum, int* nElec, int nsamp, double emax, double tmax, const char* name, const char* fol) {
		VDFluxSpec::pos = vdPos;
		VDFluxSpec::vdNum = vdNum;
		VDFluxSpec::dw = emax / PhysCon::hbar / nsamp;
		VDFluxSpec::nelecPtr = nElec;
		VDFluxSpec::nsamp = nsamp;
		VDFluxSpec::tmax = tmax;

		phss = new std::complex<double>[nsamp];
		temp = new std::complex<double>[nsamp];

		phaseCalcExpMul = new std::complex<double>[nsamp];
		for(int i = 0; i < nsamp; i++)
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
		fil.write(reinterpret_cast<char*>(&pos), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nsamp), sizeof(int));
		fil.write(reinterpret_cast<char*>(&emax), sizeof(double));
	}

	VDFluxSpec::~VDFluxSpec() {
		kill();

		if(wfcs0)
			delete[] wfcs0; wfcs0 = nullptr;
		if(wfcs1)
			delete[] wfcs1; wfcs1 = nullptr;
		delete[] phss;
		delete[] phaseCalcExpMul;
		delete[] temp;
	}

	int VDFluxSpec::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) {
		if (first) {
			nElec = *nelecPtr;

			if(wfcs0)
				delete[] wfcs0;
			if(wfcs1)
				delete[] wfcs1;
			wfcs0 = new std::complex<double>[nsamp * nElec];
			wfcs1 = new std::complex<double>[nsamp * nElec];

			for (int i = 0; i < nsamp * nElec; i++) {
				wfcs0[i] = 0;
				wfcs1[i] = 0;
			}

			celec = 0;
			first = 0;
			ct = 0;
		}
		if (celec == 0) {
			cumPotPhs *= std::exp(PhysCon::im * (v[pos] + v[pos + 1]) / (2.0 * PhysCon::hbar) * (t - ct));
			double winMul;
			if (t < tukeyAl / 2 * tmax)
				winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * t / (tukeyAl * tmax)));
			else if (t > (1.0 - tukeyAl / 2) * tmax)
				winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * (tmax - t) / (tukeyAl * tmax)));
			else
				winMul = 1;

			//exp(i t dw (idx))*cumPotPhs*winMul, expanded to hopefully vectorize better
			cblas_zcopy(nsamp, phaseCalcExpMul, 1, phss, 1); 	// phss = 		i dw (idx)
			cblas_zdscal(nsamp, t, phss, 1); 					// phss = 		i dw (idx) t
			for(int i = 0; i < nsamp; i++)
				phss[i] = std::exp(phss[i]);					// phss = exp(	i dw (idx) t)
			std::complex<double> cpwm = cumPotPhs * winMul;
			cblas_zscal(nsamp, &cpwm, phss, 1);					// phss = exp(	i dw (idx) t)*cumPotPhs*winMul
			//used to be this
			/*for (int i = 0; i < nsamp; i++)
				phss[i] = std::exp(PhysCon::im * dw * t * (double)i)*cumPotPhs*winMul;*/
			ct = t;
		}
		if (ct != t) {
			std::cout << "VDFluxSpec: Inconsistent timing." << std::endl;
			return 1;
		}
		//wfcs0[i0 + i] += psip0 * phss[i]
		//wfcs1[i0 + i] += psip1 * phss[i]
		cblas_zaxpy(nsamp, &psi[pos  ], phss, 1, &wfcs0[celec*nsamp], 1);
		cblas_zaxpy(nsamp, &psi[pos+1], phss, 1, &wfcs1[celec*nsamp], 1);
		//was this
		/*for (int i = 0; i < nsamp; i++) {
			wfcs0[celec * nsamp + i] += psi[pos] * phss[i];
			wfcs1[celec * nsamp + i] += psi[pos+1] * phss[i];
		}*/
		celec++;
		if (celec == nElec){
			celec = 0;
		}

		return 0;
	}

	void VDFluxSpec::terminate() {
		fil.write(reinterpret_cast<char*>(&wfcs0[0]), nElec * nsamp * sizeof(std::complex<double>));
		fil.write(reinterpret_cast<char*>(&wfcs1[0]), nElec * nsamp * sizeof(std::complex<double>));
		fil.close();
	}



	Vfunct::Vfunct(int n, int nx, int nt, double maxT, double * x, const char* fol) {
		Vfunct::n = n;
		Vfunct::nx = nx;
		Vfunct::nt = nt;
		Vfunct::maxT = maxT;

		interval = maxT / (double)nt;

		xs = new double[nx];
		vtls::downSampleLinearInterpolateEdge(n, x, nx, xs);
		ts = new double[nt];
		vs = new double[nx];

		curIts = 0;

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
		fil.write(reinterpret_cast<char*>(&nx), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nt), sizeof(int));
	}

	Vfunct::~Vfunct() {
		kill();

		delete[] vs;
		delete[] xs;
		delete[] ts;
	}

	int Vfunct::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		int k = (int)(t / interval) - curIts;
		if (k == 0 && curIts<nt) {
			vtls::downSampleLinearInterpolateEdge(n, v, nx, vs);

			fil.write(reinterpret_cast<char*>(&vs[0]), sizeof(double)*nx);
			ts[curIts] = t;
			curIts++;
		}
		else if (k > 1)
			while (k > 1) {
				vtls::downSampleLinearInterpolateEdge(n, v, nx, vs);

				fil.write(reinterpret_cast<char*>(&vs[0]), sizeof(double)*nx);
				ts[curIts] = t;
				curIts++;
				k = (int)(t / interval) - curIts;
			}

		return 0;
	}

	void Vfunct::terminate() {
		fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*nx);
		fil.write(reinterpret_cast<char*>(&ts[0]), sizeof(double)*nt);
		fil.close();
	}

	ExpectE0::ExpectE0(int len, double dx, const char* fol) {
		nPts = len;
		ExpectE0::dx = dx;
		rho = new double[nPts];

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

	ExpectE0::~ExpectE0() {
		kill();

		delete[] rho;
	}

	int ExpectE0::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) {
		if(first || t == tmea){
			first = 0;
			tmea = t;

			/*vtls::secondDerivative(nPts, psi, scratch1, dx);
			vtls::scaMulArray(nPts, -PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me), scratch1);
			vtls::seqMulArrays(nPts, v, psi, scratch2);
			vtls::addArrays(nPts, scratch2, scratch1);
			for (int i = 0; i < nPts; i++)
				scratch2[i] = std::conj(psi[i]);
			double ex = std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));*/

			vtls::normSqr(nPts, psi, rho);
			double ex = vtlsInt::rSumMul(nPts, rho, v, dx) + kin->evaluateKineticEnergy(psi);

			fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		}
		else if(!first && t!=tmea)
			return 1;
		return 0;
	}

	void ExpectE0::terminate() {
		fil.close();
	}


	WfcRhoWeights::WfcRhoWeights(int* nelecPtr, int nPts, double dx, WfcToRho::Weight* wght, const char* fol) : nelecPtr(nelecPtr), nPts(nPts), dx(dx), wght(wght) {
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

	WfcRhoWeights::~WfcRhoWeights() {
		kill();
	}

	int WfcRhoWeights::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator * kin) {
		if (first) {
			first = 0;
			int nElec = *nelecPtr;
			fil.write(reinterpret_cast<char*>(&nElec), sizeof(int));
			double* energies = new double[nElec];
			double* wghts = new double[nElec];
			WfcToRho::calcEnergies(nElec, nPts, dx, psi, v, kin, energies);
			wght->calcWeights(nElec, energies, wghts);

			fil.write(reinterpret_cast<char*>(&wghts[0]), sizeof(double)* nElec);

			delete[] energies;
			delete[] wghts;

			return 0;
		}
		return 1;
	}

	void WfcRhoWeights::terminate() {
		fil.close();
	}


	BasicMeasurers::BasicMeasurers(int nPts, double dx, double dt, double * xs, const char* fol) {
		meas.push_back(new NPts(nPts, fol));
		//meas.push_back(new Header(title_4char, fol));
		meas.push_back(new NSteps(fol));
		meas.push_back(new DX(dx, fol));
		meas.push_back(new DT(dt, fol));
		//meas.push_back(new XS(nPts, xs, fol));
		//meas.push_back(new TS(fol));
	}

	BasicMeasurers::~BasicMeasurers() { kill(); }

	int BasicMeasurers::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(psi, v, t, kin) == 1) {
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

	int MeasurementManager::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(psi, v, t, kin) == 1) {
				(*it)->kill();
				it = meas.erase(it);
			}
			else
				++it;
		}
		//std::cout << std::flush;
		return 0;
	}

	int MeasurementManager::measureMany(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin, int nElec, int nPts){
		for ( auto it = meas.begin(); it != meas.end(); ){
			for(int j = 0; j < nElec; j++){
				if( (*it)->measure(&psi[nPts*j], v, t, kin) == 1) {
					(*it)->kill();
					it = meas.erase(it);
					--it;
					break;
				}
			}
			++it;
		}
		/*for(int i = meas.size()-1; i>=0; i--){
			for(int j = 0; j < nElec; j++){
				if (meas[i]->measure(&psi[nPts*j], v, t, kin) == 1){
					meas[i]->kill();
					meas.erase(meas.begin() + i);
					break;
				}
			}
		}*/
		return 0;
	}

	void MeasurementManager::terminate() {
		meas.clear();
	}
}
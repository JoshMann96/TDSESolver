#include "stdafx.h"


/*Measurer::Measurer() {}
Measurer::~Measurer() {}
void Measurer::measure(std::complex<double> * psi, double * v, double t) {}
void Measurer::terminate() {}*/


//Constant

namespace Measurers {

	std::fstream openFile(const char* fil) {
		return std::fstream(fil, std::ios::out | std::ios::binary);
	}

	DoubleConst::DoubleConst(double c, const char* filName, const char* fol) {
		DoubleConst::c = c;
		int l1 = std::strlen(fol), l2 = std::strlen(filName);
		char* nfil = new char[l1 + l2 + 10];
		strncpy(nfil, fol, l1);
		strncpy(&nfil[l1], filName, l2);
		strcpy(&nfil[l1 + l2], ".tdsePART");
		/*std::stringstream ss;
		ss << fol << filName << ".tdsePART";
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
		if (fil) fil.close();
	}

	DoubleConst::~DoubleConst() {
		terminate();
	}

	int DoubleConst::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

	void DoubleConst::terminate() {
		if (fil) fil.close();
	}


	ElectronNumber::ElectronNumber(int nelec, const char* fol) {
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
		fil.write(reinterpret_cast<char*>(&nelec), sizeof(int));
		fil.close();
	}

	ElectronNumber::~ElectronNumber() {
		terminate();
	}

	int ElectronNumber::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

	void ElectronNumber::terminate() {
		if (fil) fil.close();
	}

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
	}

	Header::~Header() {
		terminate();
	}

	int Header::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

	void Header::terminate() {
		if (fil) fil.close();
	}

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
	}

	NPts::~NPts() {
		terminate();
	}

	void NPts::terminate() {
		if (fil) fil.close();
	}

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

	NSteps::~NSteps() { terminate(); }

	void NSteps::terminate() {
		fil.write(reinterpret_cast<char*>(&steps), sizeof(int));
		if (fil) fil.close();
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
	}

	DX::~DX() {
		terminate();
	}

	int DX::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

	void DX::terminate() {
		if (fil) fil.close();
	}

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
	}

	DT::~DT() {
		terminate();
	}

	int DT::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

	void DT::terminate() {
		if (fil) fil.close();
	}

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
	}

	XS::~XS() {
		terminate();
	}

	int XS::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) { return 1; }

	void XS::terminate() {
		if (fil) fil.close();
	}


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
		terminate();
	}

	int TS::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		if (t != tmea) {
			fil.write(reinterpret_cast<char*>(&t), sizeof(double));
			tmea = t;
		}
		return 0;
	}

	void TS::terminate() {
		if (fil) fil.close();
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
		terminate();
	}

	int OrigPot::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		if (!measured) {
			fil.write(reinterpret_cast<char*>(&v[0]), sizeof(double)*n);
			measured = 1;
			return 1;
		}
		return 0;
	}

	void OrigPot::terminate() {
		if (fil) fil.close();
	}

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
		terminate();
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
		else if (k > 1)
			throw "Downsampled measurers were not able to keep up with time step (there needs to be more time steps in the simulation than sampled).";
		return 0;
	}

	void Psi2t::terminate() {
		if (fil) {
			fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*nx);
			fil.write(reinterpret_cast<char*>(&ts[0]), sizeof(double)*nt);
			fil.close();
		}

	}

	
	WignerQPD::WignerQPD(int n, int nx, int np, int nt, double minP, double maxP, double maxT, double windowWidth, double * x, const char* fol) {
		WignerQPD::n = n;
		WignerQPD::nx = nx;
		WignerQPD::np = np;
		WignerQPD::nt = nt;
		WignerQPD::minP = minP;
		WignerQPD::maxP = maxP;
		WignerQPD::maxT = maxT;
		WignerQPD::x = x;

		WignerQPD::windowWid = (int)(windowWidth / (x[n - 1] - x[0])*n);

		interval = maxT / (double)nt;

		xs = new double[nx];
		vtls::downSampleLinearInterpolateEdge(n, x, nx, xs);
		ps = new double[np];
		vtls::linspace(np, minP, maxP, ps);
		ts = new double[nt];
		convEnt = new std::complex<double>[n];
		convOut = new std::complex<double>[n];
		output = new std::complex<double>[nx];

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
		fil.write(reinterpret_cast<char*>(&np), sizeof(int));
		fil.write(reinterpret_cast<char*>(&nt), sizeof(int));
	}

	WignerQPD::~WignerQPD() {
		terminate();
	}

	void WignerQPD::initializeConvolution(std::complex<double> * staticElement) {
		status = VSL_STATUS_OK;
		task_ptr = &task;

		vslzConvNewTaskX1D(task_ptr, VSL_CONV_MODE_FFT, n, n, n, staticElement, 1);
		int nPts2 = n / 2;
		vslConvSetStart(task, &nPts2);
		int one = 1;
		vslConvSetDecimation(task, &one);
	}

	int WignerQPD::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		double dx = (x[n - 1] - x[0]) / n;
		int tk = (int)(t / interval) - curIts;
		double sum;
		std::complex<double> entr;
		int pos;
		int maxk;
		if (tk == 0 && curIts < nt) {
			for (int i = 0; i < nx; i++) {
				pos = vtls::findValue(n, x, xs[i]);
				for (int j = 0; j < np; j++) {
					sum = 0;
					maxk = std::min(std::min(windowWid, pos), n - pos);
					entr = 2.0*PhysCon::im*ps[j]*dx / PhysCon::hbar;
					for (int k = 0; k < maxk; k++)
						sum += std::real(std::conj(psi[pos + k])*psi[pos - k] * std::exp(entr*(double)k));
					sum *= dx / (PhysCon::pi*PhysCon::hbar);
					fil.write(reinterpret_cast<char*>(&sum), sizeof(double));
				}
			}
			ts[curIts] = t;
			curIts++;
		}
		else if (tk > 1)
			throw "Downsampled mreasurers were not able to keep up with time step (there needs to be more time steps in the simulation than sampled).";
		return 0;
	}

	/* SPACE CONVOLUTION METHOD -- WORK ON DIRECT METHOD TO SEE IF IT MAKES SENSE
	void WignerQPD::measure(std::complex<double> * psi, double * v, double t) {
		int tk = (int)(t / interval) - curIts;
		if (tk == 0 && curIts < nt) {
			initializeConvolution(psi);
			for (int k = 0; k < np; k++) {
				for (int i = 0; i < n; i++)
					convEnt[i] = std::conj(psi[i])*std::exp(2.0*PhysCon::im*ps[k] * x[i] / PhysCon::hbar);

				status = vslzConvExecX1D(task, convEnt, 1, convOut, 1);
				if (status != VSL_STATUS_OK) {
					std::cout << "ERROR: WignerQPD convolution: bad status=" << status << std::endl;
					std::cin.ignore();
					exit(1);
				}
				//Assuming for now that the 2x parameter is just 1x
				vtls::downSampleLinearInterpolateEdge(n, convOut, nx, output);
				for (int i = 0; i < nx; i++)
					output[i] *= std::exp(-2.0*PhysCon::im*ps[k] * xs[i] / PhysCon::hbar) / (PhysCon::pi*PhysCon::hbar);
				fil.write(reinterpret_cast<char*>(&output[0]), sizeof(std::complex<double>)*nx);
			}
			status = vslConvDeleteTask(task_ptr);
			if (status != VSL_STATUS_OK) {
				std::cout << "ERROR: WignerQPD convolution: bad status=" << status << std::endl;
				std::cin.ignore();
				exit(1);
			}
			ts[curIts] = t;
			curIts++;
		}
		else if (tk > 1)
			throw "Downsampled mreasurers were not able to keep up with time step (there needs to be more time steps in the simulation than sampled).";
	}
	*/

	void WignerQPD::terminate() {
		if (fil) {
			fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*nx);
			fil.write(reinterpret_cast<char*>(&ps[0]), sizeof(double)*np);
			fil.write(reinterpret_cast<char*>(&ts[0]), sizeof(double)*nt);
			fil.close();
		}
	}

	ExpectE::ExpectE(int len, double dx, const char* fol) {
		nPts = len;
		ExpectE::dx = dx;
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

	ExpectE::~ExpectE() {
		terminate();
	}

	int ExpectE::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::secondDerivative(nPts, psi, scratch1, dx);
		vtls::scaMulArray(nPts, -PhysCon::hbar*PhysCon::hbar / (2.0*PhysCon::me), scratch1);
		vtls::seqMulArrays(nPts, v, psi, scratch2);
		vtls::addArrays(nPts, scratch2, scratch1);
		for (int i = 0; i < nPts; i++)
			scratch2[i] = std::conj(psi[i]);
		double ex = std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));

		fil.write(reinterpret_cast<char*>(&ex), sizeof(double));

		return 0;
	}

	void ExpectE::terminate() {
		if (fil) fil.close();
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
		terminate();
	}

	int ExpectX::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::normSqr(nPts, psi, scratch);
		double ex = vtlsInt::simpsMul(nPts, x, scratch, dx);
		fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		return 0;
	}

	void ExpectX::terminate() {
		if (fil) fil.close();
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
		terminate();
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
		if (fil) fil.close();
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
		terminate();
	}

	int ExpectA::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::firstDerivative(nPts, v, scratch1, dx);
		vtls::normSqr(nPts, psi, scratch2);
		double ex = vtlsInt::simpsMul(nPts, scratch2, scratch1, dx)*(-1.0 / PhysCon::me);
		fil.write(reinterpret_cast<char*>(&ex), sizeof(double));

		return 0;
	}

	void ExpectA::terminate() {
		if (fil) fil.close();
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
		terminate();
	}

	int TotProb::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		vtls::normSqr(n, psi, psi2);
		double sum = vtlsInt::simps(n, psi2, dx);
		fil.write(reinterpret_cast<char*>(&sum), sizeof(double));

		return 0;
	}

	void TotProb::terminate() {
		if (fil) fil.close();
	}


	VDProbCurrent::VDProbCurrent(int n, double dx, int vdPos, int vdNum, const char* name, const char* fol) {
		VDProbCurrent::n = n;
		VDProbCurrent::dx = dx;
		VDProbCurrent::pos = vdPos;
		VDProbCurrent::vdNum = vdNum;

		const char* nm = std::to_string(vdNum).c_str();
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
		terminate();
	}

	int VDProbCurrent::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		std::complex<double> der = vtls::firstDerivative(n, psi, pos, dx);
		double j = std::imag(PhysCon::hbar / (2.0 * PhysCon::me)*(std::conj(psi[pos])*der - psi[pos] * std::conj(der)));
		fil.write(reinterpret_cast<char*>(&j), sizeof(double));

		return 0;
	}

	void VDProbCurrent::terminate() {
		if (fil) fil.close();
	}

	PsiT::PsiT(int n, double meaT, int vdNum, const char* name, const char* fol) {
		PsiT::meaT = meaT;
		PsiT::vdNum = vdNum;
		PsiT::n = n;

		const char* nm = std::to_string(vdNum).c_str();
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
		terminate();
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
		if (fil) fil.close();
	}


	PotT::PotT(int n, double meaT, int vdNum, const char* name, const char* fol) {
		PotT::meaT = meaT;
		PotT::vdNum = vdNum;
		PotT::n = n;

		const char* nm = std::to_string(vdNum).c_str();
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
		terminate();
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
		if (fil) fil.close();
	}

	VDPsi::VDPsi(int vdPos, int vdNum, const char* name, const char* fol) {
		VDPsi::pos = vdPos;
		VDPsi::vdNum = vdNum;

		const char* nm = std::to_string(vdNum).c_str();
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
		terminate();
	}

	int VDPsi::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		fil.write(reinterpret_cast<char*>(&psi[pos]), sizeof(std::complex<double>));
		return 0;
	}

	void VDPsi::terminate() {
		if (fil) fil.close();
	}

	VDPot::VDPot(int vdPos, int vdNum, const char* name, const char* fol) {
		VDPot::pos = vdPos;
		VDPot::vdNum = vdNum;

		const char* nm = std::to_string(vdNum).c_str();
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
		terminate();
	}

	int VDPot::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		fil.write(reinterpret_cast<char*>(&v[pos]), sizeof(double));
		return 0;
	}

	void VDPot::terminate() {
		if (fil) fil.close();
	}


	VDFluxSpec::VDFluxSpec(int vdPos, int vdNum, int* nelec, int nsamp, double emax, double tmax, const char* name, const char* fol) {
		VDFluxSpec::pos = vdPos;
		VDFluxSpec::vdNum = vdNum;
		VDFluxSpec::dw = emax / PhysCon::hbar / nsamp;
		VDFluxSpec::nelecPtr = nelec;
		VDFluxSpec::nsamp = nsamp;
		VDFluxSpec::tmax = tmax;

		phss = new std::complex<double>[nsamp];
		cumPotPhs = 1;

		const char* nm = std::to_string(vdNum).c_str();
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
		terminate();
	}

	int VDFluxSpec::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) {
		if (first) {
			nelec = nelecPtr[0];
			wfcs0 = new std::complex<double>[nsamp * nelec];
			wfcs1 = new std::complex<double>[nsamp * nelec];
			for (int i = 0; i < nsamp * nelec; i++) {
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

			for (int i = 0; i < nsamp; i++)
				phss[i] = std::exp(PhysCon::im * dw * t * (double)i)*cumPotPhs*winMul;
			ct = t;
		}
		if (ct != t) {
			std::cout << "VDFluxSpec: Inconsistent timing." << std::endl;
			return 1;
		}

		for (int i = 0; i < nsamp; i++) {
			wfcs0[celec * nsamp + i] += psi[pos] * phss[i];
			wfcs1[celec * nsamp + i] += psi[pos+1] * phss[i];
		}

		celec++;
		if (celec == nelec)
			celec = 0;

		return 0;
	}

	void VDFluxSpec::terminate() {
		fil.write(reinterpret_cast<char*>(&wfcs0[0]), nelec * nsamp * sizeof(std::complex<double>));
		fil.write(reinterpret_cast<char*>(&wfcs1[0]), nelec * nsamp * sizeof(std::complex<double>));
		if (fil) fil.close();
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
		terminate();
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
		if (fil) {
			fil.write(reinterpret_cast<char*>(&xs[0]), sizeof(double)*nx);
			fil.write(reinterpret_cast<char*>(&ts[0]), sizeof(double)*nt);
			fil.close();
		}

	}

	ExpectE0::ExpectE0(int len, double dx, const char* fol) {
		nPts = len;
		ExpectE0::dx = dx;
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

	ExpectE0::~ExpectE0() {
		terminate();
	}

	int ExpectE0::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator* kin) {
		if(first || t == tmea){
			first = 0;
			tmea = t;

			vtls::secondDerivative(nPts, psi, scratch1, dx);
			vtls::scaMulArray(nPts, -PhysCon::hbar * PhysCon::hbar / (2.0 * PhysCon::me), scratch1);
			vtls::seqMulArrays(nPts, v, psi, scratch2);
			vtls::addArrays(nPts, scratch2, scratch1);
			for (int i = 0; i < nPts; i++)
				scratch2[i] = std::conj(psi[i]);
			double ex = std::real(vtlsInt::simpsMul(nPts, scratch1, scratch2, dx));

			fil.write(reinterpret_cast<char*>(&ex), sizeof(double));
		}
		else if(!first && t!=tmea)
			return 1;
		return 0;
	}

	void ExpectE0::terminate() {
		if (fil) fil.close();
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
		terminate();
	}

	int WfcRhoWeights::measure(std::complex<double>* psi, double* v, double t, KineticOperators::KineticOperator * kin) {
		if (first) {
			first = 0;
			int nelec = nelecPtr[0];
			fil.write(reinterpret_cast<char*>(&nelec), sizeof(int));
			double* energies = new double[nelec];
			double* wghts = new double[nelec];
			WfcToRho::calcEnergies(nelec, nPts, dx, psi, v, kin, energies);
			wght->calcWeights(nelec, energies, wghts);

			fil.write(reinterpret_cast<char*>(&wghts[0]), sizeof(double)* nelec);

			delete[] energies, wghts;

			return 0;
		}
		return 1;
	}

	void WfcRhoWeights::terminate() {
		if (fil) fil.close();
	}


	BasicMeasurers::BasicMeasurers(const char * title_4char, int nPts, double dx, double dt, double * xs, const char* fol) {
		meas.push_back(new NPts(nPts, fol));
		meas.push_back(new Header(title_4char, fol));
		meas.push_back(new NSteps(fol));
		meas.push_back(new DX(dx, fol));
		meas.push_back(new DT(dt, fol));
		meas.push_back(new XS(nPts, xs, fol));
		meas.push_back(new TS(fol));
	}

	BasicMeasurers::~BasicMeasurers() { terminate(); }

	int BasicMeasurers::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		for (int i = meas.size() - 1; i >= 0; i--)
			if (meas.at(i)->measure(psi, v, t, kin) == 1) {
				meas.at(i)->terminate();
				meas.erase(meas.begin() + i);
			}
		return 0;
				
	}

	void BasicMeasurers::terminate() {
		for (uint i = 0; i < meas.size(); i++)
			meas[i]->terminate();
	}


	MeasurementManager::MeasurementManager(const char* fname) {
		MeasurementManager::fname = fname;
	}

	MeasurementManager::~MeasurementManager() {
		terminate();
	}

	void MeasurementManager::addMeasurer(Measurer * m) {
		//meas.push_back(m);
		if (m->isHeavy())
			heavies.push_back(m);
		else
			lights.push_back(m);
	}

	int MeasurementManager::measure(std::complex<double> * psi, double * v, double t, KineticOperators::KineticOperator* kin) {
		/*for (int i = meas.size()-1; i >= 0; i--) {
			if (meas[i]->measure(psi, v, t) == 1)
				meas.erase(meas.begin() + i);
		}
		return 0;*/

		if (doParallel == -1) {
			// Only work in parallel if we have more than 1 heavy
			if (heavies.size() > 1)
				doParallel = 1;
			else
				doParallel = 0;
		}

		if (pool == NULL && doParallel)
			pool = new ThreadPool(std::fmin((int)heavies.size()+1, std::thread::hardware_concurrency()));

		if (doParallel) {
			//Do heavies in parallel with eachother
			for (int i = 0; i < (int)heavies.size(); i++)
				res.emplace_back(
					pool->enqueue([i, this, psi, v, t, kin] {
						return heavies.at(i)->measure(psi, v, t, kin);
						}
					)
				);
			//Do lights together, in parallel with heavies
			res.emplace_back(
				pool->enqueue([this, psi, v, t, kin] {
					for (int i = lights.size() - 1; i >= 0; i--) {
						if (lights[i]->measure(psi, v, t, kin) == 1) {
							lights[i]->terminate();
							lights.erase(lights.begin() + i);
						}
					}
					return 0;
				}
				)
			);

			for (int i = (int)heavies.size() - 1; i >= 0; i--) {
				if (res[i].get() == 1) {
					heavies[i]->terminate();
					heavies.erase(heavies.begin() + i);
				}
			}
			res.back().get();
			res.clear();

			return 0;
		}
		else {
			//Do light measurements
			for (int i = lights.size()-1; i >= 0; i--) {
				if (lights[i]->measure(psi, v, t, kin) == 1) {
					lights[i]->terminate();
					lights.erase(lights.begin() + i);
				}
			}
			//Do heavy measurement if it exists
			if (heavies.size())
				if (heavies[0]->measure(psi, v, t, kin) == 1) {
					heavies[0]->terminate();
					heavies.erase(heavies.begin());
				}
			return 0;
		}
	}

	void MeasurementManager::terminate() {
		/*const char* tfiln;
		std::fstream tfil;
		fil = std::fstream(fname, std::ios::out | std::ios::binary);
		for (int i = 0; i < (int)meas.size(); i++) {
		tfiln = meas[i]->terminate();
		tfil = std::fstream(tfiln, std::ios::in || std::ios::binary);
		int size = (int)(tfil.tellg()) + sizeof(int);
		fil.write(reinterpret_cast<char*>(&size), sizeof(int));
		fil << tfil.rdbuf();
		t//fill.close();
		tfil.clear();
		}
		//fill.close();
		//DELETE INDIVIDUAL MEASURER FILES
		*/
		for (int i = 0; i < (int)(lights.size()); i++)
			lights [i] ->terminate();
		for (int i = 0; i < (int)(heavies.size()); i++)
			heavies[i]->terminate();
	}
}
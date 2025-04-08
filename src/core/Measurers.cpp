#include "Measurers.h"
#include "KineticOperator.h"
#include "PhysCon.h"
#include <iostream>
#include <stdexcept>
#include <system_error>
#include "blas.h"


namespace Measurers {

	std::fstream openFile(const char* fil) {
		std::fstream fstrm = std::fstream(fil, std::ios::out | std::ios::binary);
		if(!fstrm)
			throw std::system_error(errno, std::system_category(), "Could not open file: " + std::string(fil));
		return fstrm;
	}

	std::fstream openFile(std::initializer_list<const char*> args) {
		std::string fil("");
		for( auto ele : args )
			fil += ele;
		std::fstream fstrm = openFile(fil.c_str());
		return fstrm;
	}

	std::fstream openFile(std::list<const char*> args) {
		std::string fil("");
		for( auto ele : args )
			fil += ele;
		std::fstream fstrm = openFile(fil.c_str());
		return fstrm;
	}

	DoubleConst::DoubleConst(double c, const char* name, const char* fol) : Measurer(-2, {fol, name}), c(c){
		write(&c, sizeof(double));
		close();
	}


	NElec::NElec(int* nElec, const char* fol) : nElec(nElec), fol(fol), Measurer(22) {}

	MeasurerStatus NElec::measure(int step, const std::complex<double> * psi, const double* v, double t) { 
		if(first && *nElec > 0){
			first = false;
			
			open({fol, fname});

			write(&index, sizeof(int));
			write(nElec, sizeof(int));

			close();
		}
		
		return MeasurerStatus::ALL_DONE; 
	}

	
	Header::Header(const char* title, const char* fol) : Measurer(-1, {fol, fname}) {
		write(title, sizeof(title));
		close();
	}

	
	NPts::NPts(int nPts, const char* fol) : Measurer(0, {fol, fname}) {
		write(&nPts, sizeof(int));
		close();
	}


	NSteps::NSteps(const char* fol) : Measurer(1, {fol, fname}), steps(0) {}
	
	NSteps::~NSteps() {
		write(&steps, sizeof(int));
	}

	MeasurerStatus NSteps::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		steps++;
		return MeasurerStatus::SUCCESS;
	}
	

	DX::DX(double dx, const char* fol) : Measurer(2, {fol, fname}) {
		write(&dx, sizeof(double));
		close();
	}

	
	DT::DT(double dt, const char* fol) : Measurer(3, {fol, fname}) {
		write(&dt, sizeof(double));
		close();
	}

	
	XS::XS(int len, const double* xs, const char* fol) : Measurer(4, {fol, fname}) {
		write(&xs[0], sizeof(double)*len);
		close();
	}

	
	TS::TS(const char* fol) : Measurer(5, {fol, fname}){}

	MeasurerStatus TS::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		write(&t, sizeof(double));
		return MeasurerStatus::SUCCESS;
	}


	OrigPot::OrigPot(int n, const char* fol) : Measurer(6, {fol, fname}), n(n){
		write(&n, sizeof(int));
	}

	MeasurerStatus OrigPot::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		write(v, sizeof(double)*n);
		return MeasurerStatus::ALL_DONE;
	}


	Psi2t::Psi2t(int nPts, int nx, int nt, int numSteps, const double * x, const int* nElec, const char* fol) :
		nPts(nPts), nx(nx), nt(nt), numSteps(numSteps), nElec(nElec), curIdx(0), Measurer(9, {fol, fname})
	{
		measSteps = (int*) sq_malloc(sizeof(int)*numSteps);
		vtls::linspace(nt, 0, numSteps - 1, measSteps);

		xs = (double*) sq_malloc(sizeof(double)*nx);
		vtls::linearInterpolateEdge(nPts, x, nx, xs);
		ts = (double*) sq_malloc(sizeof(double)*nt);
		psi2b = (double*) sq_malloc(sizeof(double)*nPts);
		psi2s = (double*) sq_malloc(sizeof(double)*nx);

		write(&nx, sizeof(int));
		write(&nt, sizeof(int));
	}

	Psi2t::~Psi2t() {
		write(xs, sizeof(double)*nx);
		write(ts, sizeof(double)*nt);

		sq_free(psi2b);
		sq_free(psi2s);
		sq_free(xs);
		sq_free(ts);
		sq_free(measSteps);
	}

	MeasurerStatus Psi2t::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		while(step == measSteps[curIdx]){
			for(int i = 0; i < *nElec; i++){
				vtls::normSqr(nPts, &psi[i*nPts], psi2b);
				vtls::linearInterpolateEdge(nPts, psi2b, nx, psi2s);
				write(psi2s, sizeof(double)*nx);
			}

			ts[curIdx] = t;

			curIdx++;
			if(curIdx >= nt)
				return MeasurerStatus::ALL_DONE;
		}
		return MeasurerStatus::SUCCESS;
	}


	ExpectE::ExpectE(int nPts, double dx, const int* nElec, const char* fol, KineticOperators::KineticOperator * const* kin) : 
		nPts(nPts), dx(dx), kin(kin), nElec(nElec), Measurer(10, {fol, fname})
	{
		rho = (double*) sq_malloc(sizeof(double)*nPts);
	}

	ExpectE::~ExpectE() {
		sq_free(rho);
	}

	MeasurerStatus ExpectE::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		double ex;
		for(int i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], rho);
			ex = vtlsInt::rSumMul(nPts, rho, v, dx) / vtlsInt::rSum(nPts, rho, dx) + (*kin)->evaluateKineticEnergy(&psi[i*nPts]);

			write(&ex, sizeof(double));
		}

		return MeasurerStatus::SUCCESS;
	}


	ExpectX::ExpectX(int nPts, const double* xs, double dx, const int* nElec, const char* fol) :
		nPts(nPts), dx(dx), nElec(nElec), x(xs), Measurer(11, {fol, fname})
	 {
		scratch = (double*) sq_malloc(sizeof(double)*nPts);
	}

	ExpectX::~ExpectX() {
		sq_free(scratch);
	}

	MeasurerStatus ExpectX::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		for(int i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], scratch);
			double ex = vtlsInt::simpsMul(nPts, x, scratch, dx);
			write(&ex, sizeof(double));
		}
		return MeasurerStatus::SUCCESS;
	}


	ExpectP::ExpectP(int len, double dx, const int* nElec, const char* fol) :
		nPts(len), dx(dx), nElec(nElec), Measurer(12, {fol, fname})
	 {
		scratch1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*len);
		scratch2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*len);
	}

	ExpectP::~ExpectP() {
		sq_free(scratch1);
		sq_free(scratch2);
	}

	MeasurerStatus ExpectP::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		double ex;
		for(int i = 0; i < *nElec; i++){
			vtls::firstDerivative(nPts, &psi[i*nPts], scratch1, dx);
			for (int j = 0; j < nPts; j++)
				scratch2[j] = std::conj(psi[i*nPts + j]);
			ex = std::imag(vtlsInt::simpsMul(nPts, scratch2, scratch1, dx))*PhysCon::hbar;
			write(&ex, sizeof(double));
		}
		return MeasurerStatus::SUCCESS;
	}


	ExpectA::ExpectA(int nPts, double dx, const int* nElec, const char* fol) :
		nPts(nPts), dx(dx), nElec(nElec), Measurer(13, {fol, fname})
		 {
		scratch1 = (double*) sq_malloc(sizeof(double)*nPts);
		scratch2 = (double*) sq_malloc(sizeof(double)*nPts);
	}

	ExpectA::~ExpectA() {
		sq_free(scratch1);
		sq_free(scratch2);
	}

	MeasurerStatus ExpectA::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		double ex;
		vtls::firstDerivative(nPts, v, scratch1, dx);
		for(int i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], scratch2);
			ex = vtlsInt::simpsMul(nPts, scratch2, scratch1, dx)*(-1.0 / PhysCon::me);
			write(&ex, sizeof(double));
		}

		return MeasurerStatus::SUCCESS;
	}


	TotProb::TotProb(int nPts, double dx, const int* nElec, const char* fol) :
		nPts(nPts), dx(dx), nElec(nElec), Measurer(16, {fol, fname})
	{
		psi2 = (double*) sq_malloc(sizeof(double)*nPts);
	}

	TotProb::~TotProb() {
		sq_free(psi2);
	}

	MeasurerStatus TotProb::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		double sum;
		for(int i = 0; i < *nElec; i++){	
			vtls::normSqr(nPts, &psi[i*nPts], psi2);
			sum = vtlsInt::simps(nPts, psi2, dx);
			write(&sum, sizeof(double));
		}
		return MeasurerStatus::SUCCESS;
	}


	VDProbCurrent::VDProbCurrent(int nPts, double dx, const int* nElec, int vdPos, int vdNum, const char* name, const char* fol) :
		nPts(nPts), dx(dx), vdPos(vdPos), nElec(nElec), Measurer(14, {fol, std::to_string(vdNum).c_str(), fname})
	 {
		write(&vdNum, sizeof(int));
		write(&name, 4);
		write(&vdPos, sizeof(int));
	}

	MeasurerStatus VDProbCurrent::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		std::complex<double> der;
		double j;
		for(int i = 0; i < *nElec; i++){
			der = vtls::firstDerivative(nPts, &psi[i*nPts], vdPos, dx);
			j = std::imag(PhysCon::hbar / (2.0 * PhysCon::me)*(std::conj(psi[i*nPts + vdPos])*der - psi[i*nPts + vdPos] * std::conj(der)));
			write(&j, sizeof(double));
		}

		return MeasurerStatus::SUCCESS;
	}


	PsiT::PsiT(int nPts, double meaT, const int *nElec, int vdNum, const char* name, const char* fol) :
		nElec(nElec), meaT(meaT), nPts(nPts), Measurer(19, {fol, std::to_string(vdNum).c_str(), fname})
	 {
		write(&vdNum, sizeof(int));
		write(&name, 4);
		write(&meaT, sizeof(double));
	}

	MeasurerStatus PsiT::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		if ((!done && t >= meaT)) {
			write(psi, sizeof(std::complex<double>)* nPts * *nElec);
			done = true;
			return MeasurerStatus::ALL_DONE;
		}
		
		return MeasurerStatus::SUCCESS;
	}


	PotT::PotT(int n, double meaT, int vdNum, const char* name, const char* fol) :
		n(n), meaT(meaT), Measurer(20, {fol, std::to_string(vdNum).c_str(), fname})
	{
		write(&vdNum, sizeof(int));
		write(&name, 4);
		write(&meaT, sizeof(double));
	}

	MeasurerStatus PotT::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		if (!done && t >= meaT) {
			write(v, sizeof(double)*n);
			done = true;
			return MeasurerStatus::ALL_DONE;
		}
		return MeasurerStatus::SUCCESS;
	}


	VDPsi::VDPsi(const int* nElec, int vdPos, int vdNum, const char* name, const char* fol) : 
		nElec(nElec), vdPos(vdPos), Measurer(15, {fol, std::to_string(vdNum).c_str(), fname})
	{
		write(&vdNum, sizeof(int));
		write(&name, 4);
		write(&vdPos, sizeof(int));
	}

	MeasurerStatus VDPsi::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		for(int i = 0; i < *nElec; i++)
			write(&psi[i*nPts + vdPos], sizeof(std::complex<double>));
		return MeasurerStatus::SUCCESS;
	}


	VDPot::VDPot(int vdPos, int vdNum, const char* name, const char* fol) :
		vdPos(vdPos), Measurer(21, {fol, std::to_string(vdNum).c_str(), fname}), vdNum(vdNum)
	 {
		write(&vdNum, sizeof(int));
		write(&name, 4);
		write(&vdPos, sizeof(int));
	}

	MeasurerStatus VDPot::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		write(&v[vdPos], sizeof(double));
		return MeasurerStatus::SUCCESS;
	}


	VDFluxSpec::VDFluxSpec(int nPts, int vdPos, int vdNum, const int* nElec, int nSamp, double emax, double tmax, const char* name, const char* fol) :
		vdPos(vdPos),  nElec(nElec), nSamp(nSamp), tmax(tmax), nPts(nPts), dw(emax / PhysCon::hbar / nSamp),
		Measurer(24, {fol, std::to_string(vdNum).c_str(), fname})
	 {
		phss = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		temp = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);

		phaseCalcExpMul = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		for(int i = 0; i < nSamp; i++)
			phaseCalcExpMul[i] = PhysCon::im * dw * (double)i; //to be multiplied by t then exponentiated later
		
		cumPotPhs = 1;

		write(&vdNum, sizeof(int));
		write(&name, 4);
		write(&vdPos, sizeof(int));
		write(&nSamp, sizeof(int));
		write(&emax, sizeof(double));
	}

	VDFluxSpec::~VDFluxSpec() {
		write(wfcs0, *nElec * nSamp * sizeof(std::complex<double>));
		write(wfcs1, *nElec * nSamp * sizeof(std::complex<double>));

		if(wfcs0)
			sq_free(wfcs0); wfcs0 = nullptr;
		if(wfcs1)
			sq_free(wfcs1); wfcs1 = nullptr;
		sq_free(phss);
		sq_free(phaseCalcExpMul);
		sq_free(temp);
	}

	MeasurerStatus VDFluxSpec::measure(int step, const std::complex<double> * psi, const double* v, double t) {
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
			
			first = false;
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

		return MeasurerStatus::SUCCESS;
	}


	Vfunct::Vfunct(int potNum, int nPts, int nx, int nt, int numSteps, double maxT, const double * x, const char* fol) :
		nPts(nPts), nx(nx), nt(nt), maxT(maxT), curIdx(0), Measurer(17, {fol, std::to_string(potNum).c_str(), fname})
	{
		measSteps = (int*) sq_malloc(sizeof(int)*numSteps);
		vtls::linspace(nt, 0, numSteps - 1, measSteps);

		xs = (double*) sq_malloc(sizeof(double)*nx);
		vtls::linearInterpolateEdge(nPts, x, nx, xs);
		ts = (double*) sq_malloc(sizeof(double)*nt);
		vs = (double*) sq_malloc(sizeof(double)*nx);

		vtls::linspace(nt, 0.0, maxT, ts);

		write(&nx, sizeof(int));
		write(&nt, sizeof(int));
	}

	Vfunct::~Vfunct() {
		write(xs, sizeof(double)*nx);
		write(ts, sizeof(double)*nt);

		sq_free(vs);
		sq_free(xs);
		sq_free(ts);
		sq_free(measSteps);
	}

	MeasurerStatus Vfunct::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		while(step == measSteps[curIdx]){
			vtls::linearInterpolateEdge(nPts, v, nx, vs);
			write(vs, sizeof(double)*nx);
			
			curIdx++;
			if(curIdx >= nt)
				return MeasurerStatus::ALL_DONE;
		}
		return MeasurerStatus::SUCCESS;
	}


	ExpectE0::ExpectE0(int nPts, double dx, const int* nElec, const char* fol, KineticOperators::KineticOperator * const* kin) : 
		nPts(nPts), dx(dx), kin(kin), nElec(nElec), Measurer(23, {fol, fname})
	{
		rho = (double*) sq_malloc(sizeof(double)*nPts);
	}

	ExpectE0::~ExpectE0() {
		sq_free(rho);
	}

	
	MeasurerStatus ExpectE0::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		if(first){
			first = false;
			double ex;
			for(int i = 0; i < *nElec; i++){
				vtls::normSqr(nPts, &psi[i*nPts], rho);
				ex = vtlsInt::rSumMul(nPts, rho, v, dx) + (*kin)->evaluateKineticEnergy(&psi[i*nPts]);
				write(&ex, sizeof(double));
			}
		}
		return MeasurerStatus::ALL_DONE;
	}


	WfcRhoWeights::WfcRhoWeights(const int* nElec, double * const * weights, const char* fol) : 
		nElec(nElec), weights(weights), Measurer(18, {fol, fname}){ needsDens=true; }

	MeasurerStatus WfcRhoWeights::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		if (first) {
			if (*weights){
				first = false;
				write(nElec, sizeof(int));
				write(*weights, sizeof(double)* *nElec);
			}
			else{
				throw std::runtime_error("Weights not set for WfcRhoWeights.");
			}
		}
		return MeasurerStatus::ALL_DONE;
	}


	BasicMeasurers::BasicMeasurers(int nPts, double dx, double dt, const char* fol)
	{
		meas.push_back(new NPts(nPts, fol));
		meas.push_back(new NSteps(fol));
		meas.push_back(new DX(dx, fol));
		meas.push_back(new DT(dt, fol));
	}

	BasicMeasurers::~BasicMeasurers() {
		for(Measurer* m : meas)
			delete m;
		meas.clear();
	}

	MeasurerStatus BasicMeasurers::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(step, psi, v, t) == MeasurerStatus::ALL_DONE) {
				delete (*it);
				it = meas.erase(it);
			}
			else
				++it;
		}
		return MeasurerStatus::SUCCESS;
	}


	DensityPlotter::DensityPlotter(int nPts, const int *nElec, double dx, const double *xs, WfcToRho::Density *const dens, double * const * wght, int stepsPerPlot, bool pause):
		nPts(nPts), nElec(nElec), dens(dens), wght(wght), dx(dx), xs(xs), pause(pause), stepsPerPlot(stepsPerPlot)
	{
		needsDens = true;
		plotter = new plotting::GNUPlotter();
		rho = (double*) sq_malloc(sizeof(double)*nPts);
	}

	DensityPlotter::~DensityPlotter(){
		delete plotter;
		sq_free(rho);
	}
	
	MeasurerStatus DensityPlotter::measure(int step, const std::complex<double> * psi, const double* v, double t){
		if(step%stepsPerPlot == 0){
			dens->calcRho(nPts, *nElec, dx, *wght, psi, rho);
			plotter->update(nPts, 1, xs, rho);

			if(pause)
				std::cin.get();
		}
		
		return MeasurerStatus::SUCCESS;
	}


	PotentialPlotter::PotentialPlotter(int nPts, const double *xs, int stepsPerPlot, bool pause):
		nPts(nPts), xs(xs), pause(pause), stepsPerPlot(stepsPerPlot)
	{
		plotter = new plotting::GNUPlotter();
	}

	PotentialPlotter::~PotentialPlotter(){
		delete plotter;
	}

	MeasurerStatus PotentialPlotter::measure(int step, const std::complex<double> * psi, const double* v, double t){
		if(step%stepsPerPlot == 0){
			plotter->update(nPts, 1, xs, v);

			if(pause)
				std::cin.get();
		}
		
		return MeasurerStatus::SUCCESS;
	}


	MeasurementManager::MeasurementManager(const char* fname) {
		MeasurementManager::fname = fname;
	}

	MeasurementManager::~MeasurementManager() {
		for(Measurer* m : meas)
			delete m;
		meas.clear();
	}

	void MeasurementManager::addMeasurer(Measurer * m) {
		meas.push_back(m);
	}

	MeasurerStatus MeasurementManager::measure(int step, const std::complex<double> * psi, const double* v, double t) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(step, psi, v, t) == 1) {
				delete (*it);
				it = meas.erase(it);
			}
			else
				++it;
		}
		//std::cout << std::flush;
		return MeasurerStatus::SUCCESS;
	}
}
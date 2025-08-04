#include "Measurers.h"
#include "KineticOperator.h"
#include "PhysCon.h"
#include <iostream>
#include <stdexcept>
#include <system_error>
#include "blas.h"
#include <filesystem>


namespace Measurers {

	std::fstream openFile(const std::string fil) {
		mtx.lock();

		std::filesystem::create_directories(std::filesystem::path(fil).parent_path());

		std::fstream fstrm = std::fstream(fil, std::ios::out | std::ios::binary);
		if(!fstrm)
			throw std::system_error(errno, std::system_category(), "Could not open file: " + fil);

		mtx.unlock();

		return fstrm;
	}

	DoubleConst::DoubleConst(double c, const std::string name, const std::string fol) : Measurer(-2, fol, name), c(c){
		write(&c, sizeof(double));
		close();
	}


	NElec::NElec(size_t* nElec, const std::string fol) : nElec(nElec), fol(fol), Measurer(22) {}

	MeasurerStatus NElec::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) { 
		if(first && *nElec > 0){
			first = false;
			
			open({makeDirectory(fol), fname});

			write(nElec, sizeof(size_t));

			close();
		}
		
		return MeasurerStatus::ALL_DONE; 
	}

	
	Header::Header(const std::string title, const std::string fol) : Measurer(-1, fol, fname) {
		assert(title.length() == 8);

		write(title.c_str(), 8);
		close();
	}

	
	NPts::NPts(size_t nPts, const std::string fol) : Measurer(0, fol, fname) {
		write(&nPts, sizeof(size_t));
		close();
	}


	NSteps::NSteps(const std::string fol) : Measurer(1, fol, fname), steps(0) {}
	
	NSteps::~NSteps() {
		write(&steps, sizeof(size_t));
	}

	MeasurerStatus NSteps::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		steps++;
		return MeasurerStatus::SUCCESS;
	}
	

	DX::DX(double dx, const std::string fol) : Measurer(2, fol, fname) {
		write(&dx, sizeof(double));
		close();
	}

	
	DT::DT(double dt, const std::string fol) : Measurer(3, fol, fname) {
		write(&dt, sizeof(double));
		close();
	}

	
	XS::XS(size_t len, const double* xs, const std::string fol) : Measurer(4, fol, fname) {
		write(&xs[0], sizeof(double)*len);
		close();
	}

	
	TS::TS(const std::string fol) : Measurer(5, fol, fname){}

	MeasurerStatus TS::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		write(&t, sizeof(double));
		return MeasurerStatus::SUCCESS;
	}


	OrigPot::OrigPot(size_t n, const std::string fol) : Measurer(6, fol, fname), n(n){
		write(&n, sizeof(size_t));
	}

	MeasurerStatus OrigPot::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		write(v, sizeof(double)*n);
		return MeasurerStatus::ALL_DONE;
	}


	Psi2t::Psi2t(size_t nPts, size_t nx, size_t nt, size_t numSteps, const double * x, const size_t* nElec, const std::string fol) :
		nPts(nPts), nx(nx), nt(nt), numSteps(numSteps), nElec(nElec), curIdx(0), Measurer(9, fol, fname)
	{
		measSteps = (size_t*) sq_malloc(sizeof(size_t)*nt);
		vtls::linspace(nt, (size_t)0, (size_t)numSteps, measSteps);

		xs = (double*) sq_malloc(sizeof(double)*nx);
		vtls::linearInterpolateEdge(nPts, x, nx, xs);
		ts = (double*) sq_malloc(sizeof(double)*nt);
		psi2b = (double*) sq_malloc(sizeof(double)*nPts);
		psi2s = (double*) sq_malloc(sizeof(double)*nx);

		write(&nx, sizeof(size_t));
		write(&nt, sizeof(size_t));
	}

	Psi2t::~Psi2t() {
		if(curIdx != nt){
			std::cerr << "Warning: Psi2t measurer terminated before all measurements were made. Expected " << nt << " measurements, but only " << curIdx << " were made." << std::endl;
			std::cerr << "\t Padding with zeros." << std::endl;

			std::fill_n(psi2s, nx, 0.0);
			while(curIdx < nt){
				for(size_t i = 0; i < *nElec; i++)
					write(psi2s, sizeof(double)*nx);
				ts[curIdx] = 0.0;
				curIdx++;
			}
		}

		write(xs, sizeof(double)*nx);
		write(ts, sizeof(double)*nt);

		sq_free(psi2b);
		sq_free(psi2s);
		sq_free(xs);
		sq_free(ts);
		sq_free(measSteps);
	}

	MeasurerStatus Psi2t::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		while(step >= measSteps[curIdx]){
			for(size_t i = 0; i < *nElec; i++){
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


	ExpectE::ExpectE(size_t nPts, double dx, const size_t* nElec, const std::string fol, KineticOperators::KineticOperator * const* kin) : 
		nPts(nPts), dx(dx), kin(kin), nElec(nElec), Measurer(10, fol, fname)
	{
		rho = (double*) sq_malloc(sizeof(double)*nPts);
	}

	ExpectE::~ExpectE() {
		sq_free(rho);
	}

	MeasurerStatus ExpectE::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		double ex;
		for(size_t i = 0; i < *nElec; i++){
			ex = (*kin)->evaluateEnergy(&psi[i*nPts], v);

			write(&ex, sizeof(double));
		}

		return MeasurerStatus::SUCCESS;
	}


	ExpectX::ExpectX(size_t nPts, const double* xs, double dx, const size_t* nElec, const std::string fol, size_t minPos, size_t maxPos, double maskLength) :
		nPts(nPts), dx(dx), nElec(nElec), x(xs), Measurer(11, fol, fname)
	 {
		scratch = (double*) sq_malloc(sizeof(double)*nPts);
		mask = (double*) sq_malloc(sizeof(double)*nPts);
		vtls::masks::biSigmoid(nPts, minPos, maxPos, maskLength/dx, mask);
	}

	ExpectX::~ExpectX() {
		sq_free(scratch);
	}

	MeasurerStatus ExpectX::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		for(size_t i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], scratch);
			vtls::seqMulArrays(nPts, mask, scratch);
			double ex = vtlsInt::simpsMul(nPts, x, scratch, dx) / vtlsInt::simps(nPts, scratch, dx);
			write(&ex, sizeof(double));
		}
		return MeasurerStatus::SUCCESS;
	}


	ExpectP::ExpectP(size_t len, double dx, const size_t* nElec, const std::string fol, size_t minPos, size_t maxPos, double maskLength) :
		nPts(len), dx(dx), nElec(nElec), Measurer(12, fol, fname)
	 {
		scratch1 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*len);
		scratch2 = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*len);
		mask = (double*) sq_malloc(sizeof(double)*nPts);
		vtls::masks::biSigmoid(nPts, minPos, maxPos, maskLength/dx, mask);
	}

	ExpectP::~ExpectP() {
		sq_free(scratch1);
		sq_free(scratch2);
	}

	MeasurerStatus ExpectP::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		double ex;
		for(size_t i = 0; i < *nElec; i++){
			vtls::firstDerivative(nPts, &psi[i*nPts], scratch1, dx);
			for (size_t j = 0; j < nPts; j++)
				scratch2[j] = std::conj(psi[i*nPts + j]) * mask[j];
			ex = std::imag(vtlsInt::simpsMul(nPts, scratch2, scratch1, dx))*PhysCon::hbar / vtls::getNorm(nPts, scratch2, dx);
			write(&ex, sizeof(double));
		}
		return MeasurerStatus::SUCCESS;
	}


	ExpectA::ExpectA(size_t nPts, double dx, const size_t* nElec, const std::string fol, size_t minPos, size_t maxPos, double maskLength) :
		nPts(nPts), dx(dx), nElec(nElec), Measurer(13, fol, fname)
		 {
		scratch1 = (double*) sq_malloc(sizeof(double)*nPts);
		scratch2 = (double*) sq_malloc(sizeof(double)*nPts);
		mask = (double*) sq_malloc(sizeof(double)*nPts);
		vtls::masks::biSigmoid(nPts, minPos, maxPos, maskLength/dx, mask);
	}

	ExpectA::~ExpectA() {
		sq_free(scratch1);
		sq_free(scratch2);
	}

	MeasurerStatus ExpectA::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		double ex;
		vtls::firstDerivative(nPts, v, scratch1, dx);
		for(size_t i = 0; i < *nElec; i++){
			vtls::normSqr(nPts, &psi[i*nPts], scratch2);
			vtls::seqMulArrays(nPts, mask, scratch2);
			ex = vtlsInt::simpsMul(nPts, scratch2, scratch1, dx)*(-1.0 / PhysCon::me)/vtlsInt::simps(nPts, scratch2, dx);
			write(&ex, sizeof(double));
		}

		return MeasurerStatus::SUCCESS;
	}


	TotProb::TotProb(size_t nPts, double dx, const size_t* nElec, const std::string fol) :
		nPts(nPts), dx(dx), nElec(nElec), Measurer(16, fol, fname)
	{
		psi2 = (double*) sq_malloc(sizeof(double)*nPts);
	}

	TotProb::~TotProb() {
		sq_free(psi2);
	}

	MeasurerStatus TotProb::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		double sum;
		for(size_t i = 0; i < *nElec; i++){	
			vtls::normSqr(nPts, &psi[i*nPts], psi2);
			sum = vtlsInt::sum(nPts, psi2, dx);
			write(&sum, sizeof(double));
		}
		return MeasurerStatus::SUCCESS;
	}


	VDProbCurrent::VDProbCurrent(size_t nPts, double dx, const size_t* nElec, size_t vdPos, int vdNum, const std::string name, const std::string fol) :
		nPts(nPts), dx(dx), vdPos(vdPos), nElec(nElec), Measurer(14, fol, std::to_string(vdNum) + fname)
	 {
		assert(name.length() == 4);
		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&vdPos, sizeof(size_t));
	}

	MeasurerStatus VDProbCurrent::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		std::complex<double> der;
		double j;
		for(size_t i = 0; i < *nElec; i++){
			der = vtls::firstDerivative(nPts, &psi[i*nPts], vdPos, dx);
			j = std::imag(PhysCon::hbar / (2.0 * PhysCon::me)*(std::conj(psi[i*nPts + vdPos])*der - psi[i*nPts + vdPos] * std::conj(der)));
			write(&j, sizeof(double));
		}

		return MeasurerStatus::SUCCESS;
	}


	PsiT::PsiT(size_t nPts, double meaT, const size_t *nElec, int vdNum, const std::string name, const std::string fol) :
		nElec(nElec), meaT(meaT), nPts(nPts), Measurer(19, fol, std::to_string(vdNum) + fname)
	 {
		assert(name.length() == 4);
		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&meaT, sizeof(double));
	}

	MeasurerStatus PsiT::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		if ((!done && t >= meaT)) {
			write(psi, sizeof(std::complex<double>)* nPts * *nElec);
			done = true;
			return MeasurerStatus::ALL_DONE;
		}
		
		return MeasurerStatus::SUCCESS;
	}


	PotT::PotT(size_t n, double meaT, int vdNum, const std::string name, const std::string fol) :
		n(n), meaT(meaT), Measurer(20, fol, std::to_string(vdNum) + fname)
	{
		assert(name.length() == 4);
		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&meaT, sizeof(double));
	}

	MeasurerStatus PotT::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		if (!done && t >= meaT) {
			write(v, sizeof(double)*n);
			done = true;
			return MeasurerStatus::ALL_DONE;
		}
		return MeasurerStatus::SUCCESS;
	}


	VDPsi::VDPsi(const size_t* nElec, size_t vdPos, int vdNum, const std::string name, const std::string fol) : 
		nElec(nElec), vdPos(vdPos), Measurer(15, fol, std::to_string(vdNum) + fname)
	{
		assert(name.length() == 4);
		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&vdPos, sizeof(size_t));
	}

	MeasurerStatus VDPsi::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		for(size_t i = 0; i < *nElec; i++)
			write(&psi[i*nPts + vdPos], sizeof(std::complex<double>));
		return MeasurerStatus::SUCCESS;
	}


	VDPot::VDPot(size_t vdPos, int vdNum, const std::string name, const std::string fol) :
		vdPos(vdPos), Measurer(21, fol, std::to_string(vdNum) + fname), vdNum(vdNum)
	 {
		assert(name.length() == 4);
		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&vdPos, sizeof(size_t));
	}

	MeasurerStatus VDPot::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		write(&v[vdPos], sizeof(double));
		return MeasurerStatus::SUCCESS;
	}


	VDFluxSpec::VDFluxSpec(size_t nPts, double dx, double dt, size_t vdPos, int vdNum, const size_t* nElec, size_t nSamp, double emax, KineticOperators::KineticOperator** kinOp, double tmax, const std::string name, const std::string fol) :
		nElec(nElec), dx(dx), dt(dt), emax(emax), nSamp(nSamp), tmax(tmax), nPts(nPts), kinOp(kinOp),
		Measurer(24, fol, std::to_string(vdNum) + fname)
	 {
		assert(name.length() == 4);

		// set right-sided derivative by default, left-sided if on right boundary
		if(vdPos == nPts-1){
			vdpL = nPts-2;
			vdpR = nPts-1;
		}
		else{
			vdpL = vdPos;
			vdpR = vdPos+1;
		}

		phsL = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		phsR = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		scaledPhsL = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		scaledPhsR = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		kineticEnergies = (double*) sq_malloc(sizeof(double)*nSamp);
		sqrtGroupVelocities = (double*) sq_malloc(sizeof(double)*nSamp);

		std::fill_n(phsL, nSamp, 1.0);
		std::fill_n(phsR, nSamp, 1.0);
		vtls::linspace(nSamp, 0.0, emax, kineticEnergies);

		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&vdPos, sizeof(size_t));
		write(&nSamp, sizeof(size_t));
	}

	VDFluxSpec::~VDFluxSpec() {
		// combine left and right vds for directional flux
		// calculate wavenumbers
		double* ks = (double*) sq_malloc(sizeof(double) * nSamp);
		for(int i = 0; i < nSamp; i++)
			ks[i] = getWavenumber(kineticEnergies[i]);

		// momentum space wavefunctions
		std::complex<double>* psik = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>) * (2*nSamp-1) * *nElec);
		for(size_t i = 0; i < *nElec; i++){
			// negative wavenumbers
			for(size_t j = 1; j < nSamp; j++)
				psik[i*(2*nSamp-1) + (nSamp-1)-j] = -std::complex<double>(0,0.5) / std::sin(ks[j]*dx) / std::sqrt(2.0*PhysCon::pi * PhysCon::hbar) * ( // is ks supposed to be in the numerator?
					std::exp(0.5*PhysCon::im*ks[j]*dx) * wfcsL[i*nSamp + j] - std::exp(-0.5*PhysCon::im*ks[j]*dx) * wfcsR[i*nSamp + j]);
			// zero wavenumber
			psik[i*(2*nSamp-1) + nSamp-1] = 0.5 / std::sqrt(2.0*PhysCon::pi * PhysCon::hbar) * (wfcsL[i*nSamp] + wfcsR[i*nSamp]); // zero wavenumber is average of left and right wavefunctions... gets overwritten by velocity anyway
			// positive wavenumbers (flip sign in exponent)
			for(size_t j = 1; j < nSamp; j++)
				psik[i*(2*nSamp-1) + (nSamp-1) + j] = std::complex<double>(0,0.5) / std::sin(ks[j]*dx) / std::sqrt(2.0*PhysCon::pi * PhysCon::hbar) * ( 
					std::exp(-0.5*PhysCon::im*ks[j]*dx) * wfcsL[i*nSamp + j] - std::exp(0.5*PhysCon::im*ks[j]*dx) * wfcsR[i*nSamp + j]);
		}

		// array of signed momenta
		double* signedMomenta = (double*) sq_malloc(sizeof(double) * (2*nSamp-1));
		for(size_t i = 1; i < nSamp; i++){
			signedMomenta[nSamp-1-i] = -ks[i];
			signedMomenta[nSamp-1+i] = ks[i];
		}
		signedMomenta[nSamp-1] = 0.0;

		// array of signed kinetic energies
		double* signedEnergies = (double*) sq_malloc(sizeof(double) * (2*nSamp-1));
		for(size_t i = 1; i < nSamp; i++){
			signedEnergies[nSamp-1-i] = -getEnergy(ks[i]);
			signedEnergies[nSamp-1+i] = getEnergy(ks[i]);
		}
		signedEnergies[nSamp-1] = 0.0;

		// write to file
		write(signedEnergies, sizeof(double) * (2*nSamp-1));
		write(signedMomenta, sizeof(double) * (2*nSamp-1));
		write(psik, sizeof(std::complex<double>) * (2*nSamp-1) * *nElec);

		if(wfcsL)
			sq_free(wfcsL); wfcsL = nullptr;
		if(wfcsR)
			sq_free(wfcsR); wfcsR = nullptr;
		sq_free(phsL);
		sq_free(scaledPhsL);
		sq_free(phsR);
		sq_free(scaledPhsR);
		sq_free(kineticEnergies);
		sq_free(ks);
		sq_free(psik);
		sq_free(signedMomenta);
		sq_free(signedEnergies);
		sq_free(sqrtGroupVelocities);
	}

	MeasurerStatus VDFluxSpec::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		if (first) {
			if(wfcsL)
				sq_free(wfcsL);
			if(wfcsR)
				sq_free(wfcsR);
			wfcsL = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp * *nElec);
			wfcsR = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp * *nElec);

			std::fill_n(wfcsL, nSamp * *nElec, 0.0);
			std::fill_n(wfcsR, nSamp * *nElec, 0.0);
			
			first = false;
			ct = t;
			tstart = t;

			timeEvolutionType = (*kinOp)->getTimeEvolutionType();

			// ensure max kinetic energy does not exceed grid density
			try{
				if(getWavenumber(emax)*dx > PhysCon::pi/2.0)
					emax = getEnergy(PhysCon::pi/(2.0*dx)*0.9999);
			}
			catch(const std::runtime_error& e){
				emax = std::clamp(emax, 0.0, getEnergy(PhysCon::pi/(2.0*dx)*0.9999));
			}

			std::fill_n(phsL, nSamp, 1.0);
			std::fill_n(phsR, nSamp, 1.0);
			vtls::linspace(nSamp, 0.0, emax, kineticEnergies);
		}

		// calculate Tukey window value
		double winMul;
		if (t < tukeyAl / 2 * tmax)
			winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * t / (tukeyAl * tmax)));
		else if (t > (1.0 - tukeyAl / 2) * tmax)
			winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * (tmax - t) / (tukeyAl * tmax)));
		else
			winMul = 1.0;

		// accumulate phase for each kinetic energy
		switch(timeEvolutionType){
			case KineticOperators::TimeEvolutionType::PSEUDOSPECTRAL:
				advancePhaseOS(nSamp, t - ct, v[vdpL], kineticEnergies, phsL);
				advancePhaseOS(nSamp, t - ct, v[vdpR], kineticEnergies, phsR);
				break;
			case KineticOperators::TimeEvolutionType::CRANK_NICOLSON:
				advancePhaseCN(nSamp, t - ct, v[vdpL], kineticEnergies, phsL);
				advancePhaseCN(nSamp, t - ct, v[vdpR], kineticEnergies, phsR);
				break;
			default:
				throw std::runtime_error("Unknown time evolution type in VDFluxSpec measurer.");
		}
		// apply the window function and time step to the phase
		vtls::scaMulArray(nSamp, winMul*(t-ct), phsL, scaledPhsL);
		vtls::scaMulArray(nSamp, winMul*(t-ct), phsR, scaledPhsR);
		// apply the square-root group velocity to the phase
		fillSqrtGroupVelocities((v[vdpL] + v[vdpR])/2.0);
		vtls::seqMulArrays(nSamp, sqrtGroupVelocities, scaledPhsL);
		vtls::seqMulArrays(nSamp, sqrtGroupVelocities, scaledPhsR);

		for(size_t i = 0; i < *nElec; i++){
			//wfcs0[i0 + i] += psip0 * phss[i]
			//wfcs1[i0 + i] += psip1 * phss[i]
			cblas_zaxpy(nSamp, &psi[i*nPts + vdpL], scaledPhsL, 1, &wfcsL[i*nSamp], 1);
			cblas_zaxpy(nSamp, &psi[i*nPts + vdpR], scaledPhsR, 1, &wfcsR[i*nSamp], 1);
		}

		ct = t;

		return MeasurerStatus::SUCCESS;
	}


	VDUnidirectionalFluxSpec::VDUnidirectionalFluxSpec(size_t nPts, double dx, double dt, size_t vdPos, int vdNum, const size_t* nElec, size_t nSamp, double emax, KineticOperators::KineticOperator** kinOp, double tmax, std::string name, const std::string fol) :
		nElec(nElec), dx(dx), dt(dt), emax(emax), nSamp(nSamp), nPts(nPts), tmax(tmax), kinOp(kinOp),
		Measurer(24, fol, std::to_string(vdNum) + fname)
	 {
		assert(name.length() == 4);
		assert(vdPos <= nPts-1 && vdPos >= 0);
		vdp = vdPos;

		kineticEnergies = (double*) sq_malloc(sizeof(double)*nSamp);
		phs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		scaledPhs = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp);
		sqrtGroupVelocities = (double*) sq_malloc(sizeof(double)*nSamp);

		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&vdPos, sizeof(size_t));
		write(&nSamp, sizeof(size_t));
	}

	VDUnidirectionalFluxSpec::~VDUnidirectionalFluxSpec(){
		// write energies
		write(kineticEnergies, sizeof(double)*nSamp);

		// write wavenumbers
		double* ks = (double*) sq_malloc(sizeof(double) * nSamp);
		for(int i = 0; i < nSamp; i++)
			ks[i] = getWavenumber(kineticEnergies[i]);
		write(ks, sizeof(double)*nSamp);
		sq_free(ks);

		// write scaled Fourier transform
		vtls::scaMulArray(nSamp, 1.0/std::sqrt(2.0*PhysCon::pi * PhysCon::hbar), psift);
		write(psift, sizeof(std::complex<double>)*nSamp * *nElec);

		if(psift)
			sq_free(psift);
		sq_free(kineticEnergies);
		sq_free(phs);
		sq_free(scaledPhs);
		sq_free(sqrtGroupVelocities);
	}

	MeasurerStatus VDUnidirectionalFluxSpec::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		if (first) {
			if(psift)
				sq_free(psift);
			psift = (std::complex<double>*) sq_malloc(sizeof(std::complex<double>)*nSamp * *nElec);
			std::fill_n(psift, nSamp * *nElec, 0.0);

			first = false;
			ct = t;
			tstart = t;

			timeEvolutionType = (*kinOp)->getTimeEvolutionType();

			try{
				if(getWavenumber(emax)*dx > PhysCon::pi/2.0)
					emax = getEnergy(PhysCon::pi/(2.0*dx)*0.9999);
			}
			catch(const std::runtime_error& e){
				emax = std::clamp(emax, 0.0, getEnergy(PhysCon::pi/(2.0*dx)*0.9999));
			}

			std::fill_n(phs, nSamp, 1.0);
			vtls::linspace(nSamp, 0.0, emax, kineticEnergies);
		}

		// calculate Tukey window value
		double winMul;
		if (t < tukeyAl / 2 * tmax)
			winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * t / (tukeyAl * tmax)));
		else if (t > (1.0 - tukeyAl / 2) * tmax)
			winMul = 0.5 * (1 - std::cos(2.0 * PhysCon::pi * (tmax - t) / (tukeyAl * tmax)));
		else
			winMul = 1.0;

		// accumulate phase for each kinetic energy
		switch(timeEvolutionType){
			case KineticOperators::TimeEvolutionType::PSEUDOSPECTRAL:
				VDFluxSpec::advancePhaseOS(nSamp, t - ct, v[vdp], kineticEnergies, phs);
				break;
			case KineticOperators::TimeEvolutionType::CRANK_NICOLSON:
				VDFluxSpec::advancePhaseCN(nSamp, t - ct, v[vdp], kineticEnergies, phs);
				break;
			default:
				throw std::runtime_error("Unknown time evolution type in VDFluxSpec measurer.");
		}
		// apply the window function and time step to the phase
		vtls::scaMulArray(nSamp, winMul*(t-ct), phs, scaledPhs);
		// apply the square-root group velocity to the phase
		fillSqrtGroupVelocities(v[vdp]);
		vtls::seqMulArrays(nSamp, sqrtGroupVelocities, scaledPhs);

		for(size_t i = 0; i < *nElec; i++)
			cblas_zaxpy(nSamp, &psi[i*nPts + vdp], scaledPhs, 1, &psift[i*nSamp], 1);

		ct = t;

		return MeasurerStatus::SUCCESS;
	}


	VDClassicalFluxSpec::VDClassicalFluxSpec(size_t nPts, double dx, double dt, size_t vdPos, int vdNum, const size_t* nElec, size_t nSamp, double emax, const std::string name, const std::string fol) :
		nElec(nElec), dx(dx), dt(dt), emax(emax), nSamp(nSamp), nPts(nPts),
		Measurer(24, fol, std::to_string(vdNum) + fname)
	 {
		assert(name.length() == 4);

		// set right-sided derivative by default, left-sided if on right boundary
		if(vdPos == nPts-1){
			vdpL = nPts-2;
			vdpR = nPts-1;
		}
		else{
			vdpL = vdPos;
			vdpR = vdPos+1;
		}
		momenta = (double*) sq_malloc(sizeof(double)*(2*nSamp-1));

		mommax = std::sqrt(2.0*PhysCon::me*emax)/PhysCon::hbar;
		vtls::linspace(2*nSamp-1, -mommax, mommax, momenta);
		dk = momenta[1] - momenta[0];

		write(&vdNum, sizeof(int));
		write(name.c_str(), 4);
		write(&vdPos, sizeof(size_t));
		write(&nSamp, sizeof(size_t));
	}

	VDClassicalFluxSpec::~VDClassicalFluxSpec() {
		// write momenta
		write(momenta, sizeof(double)*(2*nSamp-1));
		
		// write yields
		if(!yields){
			std::cerr << "Warning: VDClassicalFluxSpec measurer terminated before measurements were made. Zeros will be written." << std::endl;
			yields = (double*) sq_malloc(sizeof(double)*(2*nSamp-1)*(*nElec));
			std::fill_n(yields, (2*nSamp-1)*(*nElec), 0.0);
		}
		write(yields, sizeof(double)*(2*nSamp-1)*(*nElec));

		if(momenta)
			sq_free(momenta); momenta = nullptr;
		if(yields)
			sq_free(yields); yields = nullptr;
	}

	MeasurerStatus VDClassicalFluxSpec::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t){
		if(first){
			if(yields)
				sq_free(yields);
			yields = (double*) sq_malloc(sizeof(double)*(2*nSamp-1)*(*nElec));
			std::fill_n(yields, (2*nSamp-1)* (*nElec), 0.0);
			ct = t;
			first = false;
		}
		double dt = t - ct;
		double ccur, crho, cmom;
		for(size_t i = 0; i < *nElec; i++){
			// calculate probability current
			ccur = 0.5 * PhysCon::hbar / PhysCon::me * std::imag(
				std::conj(psi[i*nPts+vdpL] + psi[i*nPts+vdpR])*
				(psi[i*nPts+vdpR] - psi[i*nPts+vdpL]) / (dx));
			// density
			crho = 0.25*std::norm(psi[i*nPts+vdpL] + psi[i*nPts+vdpR]);
			// momentum
			cmom = PhysCon::me / PhysCon::hbar * ccur / crho;

			if(std::abs(cmom) <= mommax) // map momentum to index, integrate
				yields[i*(2*nSamp-1) + (size_t)(((cmom + mommax) * (2*nSamp-2)) / (2*mommax) + 0.5)] += ccur * dt / dk;
		}
		ct = t;

		return MeasurerStatus::SUCCESS;
	}

	
	DownsampleMeasurer::DownsampleMeasurer(size_t nPts, size_t nx, size_t nt, size_t numSteps, const double* x, int measIndex, const std::string fname, const std::string fol) : 
		nPts(nPts), nx(nx), nt(nt), numSteps(numSteps), curIdx(0), Measurer(measIndex, fol, fname)
	{
		measSteps = (size_t*) sq_malloc(sizeof(size_t)*nt);
		vtls::linspace(nt, (size_t)0, (size_t)numSteps, measSteps);

		xs = (double*) sq_malloc(sizeof(double)*nx);
		vtls::linearInterpolateEdge(nPts, x, nx, xs);
		ts = (double*) sq_malloc(sizeof(double)*nt);
		data = (double*) sq_malloc(sizeof(double)*nx);

		write(&nx, sizeof(size_t));
		write(&nt, sizeof(size_t));
	}

	DownsampleMeasurer::~DownsampleMeasurer() {
		if(curIdx < nt){
			std::cerr << "Warning: DownsampleMeasurer terminated before all measurements were made. Expected " << nt << " measurements, but only " << curIdx << " were made." << std::endl;
			std::cerr << "\t Padding with zeros." << std::endl;

			std::fill_n(data, nx, 0.0);
			while(curIdx < nt){
				write(data, sizeof(double)*nx);
				ts[curIdx] = 0.0;
				curIdx++;
			}
		}

		write(xs, sizeof(double)*nx);
		write(ts, sizeof(double)*nt);

		sq_free(data);
		sq_free(xs);
		sq_free(ts);
		sq_free(measSteps);
	}

	MeasurerStatus DownsampleMeasurer::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		const double* source = dataSource(rho, cur, v);
		while(step >= measSteps[curIdx]){
			vtls::linearInterpolateEdge(nPts, source, nx, data);
			write(data, sizeof(double)*nx);

			ts[curIdx] = t;
			curIdx++;
			if(curIdx >= nt)
				return MeasurerStatus::ALL_DONE;
		}
		return MeasurerStatus::SUCCESS;
	}


	ExpectE0::ExpectE0(size_t nPts, double dx, const size_t* nElec, const std::string fol, KineticOperators::KineticOperator * const* kin) : 
		nPts(nPts), dx(dx), kin(kin), nElec(nElec), Measurer(23, fol, fname)
	{
		rho = (double*) sq_malloc(sizeof(double)*nPts);
	}

	ExpectE0::~ExpectE0() {
		sq_free(rho);
	}

	
	MeasurerStatus ExpectE0::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		if(first){
			first = false;
			double ex;
			for(size_t i = 0; i < *nElec; i++){
				ex = (*kin)->evaluateEnergy(&psi[i*nPts], v);
				write(&ex, sizeof(double));
			}
		}
		return MeasurerStatus::ALL_DONE;
	}


	WfcRhoWeights::WfcRhoWeights(const size_t* nElec, double * const * weights, const std::string fol) : 
		nElec(nElec), weights(weights), Measurer(18, fol, fname){}

	MeasurerStatus WfcRhoWeights::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		if (first) {
			if (*weights){
				first = false;
				write(nElec, sizeof(size_t));
				write(*weights, sizeof(double)* *nElec);
			}
			else{
				throw std::runtime_error("Weights not set for WfcRhoWeights.");
			}
		}
		return MeasurerStatus::ALL_DONE;
	}


	BasicMeasurers::BasicMeasurers(size_t nPts, double dx, double dt, const std::string fol)
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

	MeasurerStatus BasicMeasurers::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(step, psi, rho, cur, v, t) == MeasurerStatus::ALL_DONE) {
				delete (*it);
				it = meas.erase(it);
			}
			else
				++it;
		}
		return MeasurerStatus::SUCCESS;
	}


	QuantityPlotter::QuantityPlotter(size_t nPts, const double *xs, const char* ylabel, bool plotChange, size_t stepsPerPlot, bool pause):
		nPts(nPts), xs(xs), pause(pause), stepsPerPlot(stepsPerPlot)
	{
		if(plotChange){
			q0 = (double*) sq_malloc(sizeof(double)*nPts);
			tq = (double*) sq_malloc(sizeof(double)*nPts);
		}
		plotter = new plotting::GNUPlotter("x (m)", ylabel);
	}

	QuantityPlotter::~QuantityPlotter(){
		if(q0)
			sq_free(q0);
		if(tq)
			sq_free(tq);
		delete plotter;
	}

	void QuantityPlotter::plotQuantity(size_t step, const double* quant){
		if(step%stepsPerPlot == 0){
			if(q0){ // plot difference
				if(first){
					vtls::copyArray(nPts, quant, q0);
					first = false;
				}
				vtls::scaMulAddArrays(nPts, -1.0, q0, quant, tq);
				plotter->update(nPts, 1, xs, tq);
			}
			else
				plotter->update(nPts, 1, xs, quant);

			if(pause)
				std::cin.get();
		}
	}


	MeasurementManager::MeasurementManager(const std::string fname) : fname(fname) {}

	MeasurementManager::~MeasurementManager() {
		for(Measurer* m : meas)
			delete m;
		meas.clear();
	}

	void MeasurementManager::addMeasurer(Measurer * m) {
		meas.push_back(m);
	}

	MeasurerStatus MeasurementManager::measure(size_t step, const std::complex<double> * psi, const double * rho, const double* cur, const double* v, double t) {
		for ( auto it = meas.begin(); it != meas.end(); ){
			if( (*it)->measure(step, psi, rho, cur, v, t) == MeasurerStatus::ALL_DONE) {
				delete (*it);
				it = meas.erase(it);
			}
			else
				++it;
		}
		return MeasurerStatus::SUCCESS;
	}
}
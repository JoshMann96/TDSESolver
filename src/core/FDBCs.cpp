#include "CORECommonHeader.h"
#include "MathTools.h"
#include "KineticOperator.h"
#include "FDBCs.h"

namespace FDBCs{
    UniformHDTransparentBC::UniformHDTransparentBC(size_t order, size_t nElec, double dx, double dt) : order(order), nElec(nElec), dx(dx / PhysCon::a0), dt(dt / PhysCon::hbar * PhysCon::auE_ha) {
        if (order < 1)
            throw std::invalid_argument("Order of HDTransparentBC must be greater than 0.");

        psis = new CyclicArray<std::complex<double>>*[nElec];
        for (size_t i = 0; i < nElec; i++){
            psis[i] = new CyclicArray<std::complex<double>>(order, 0.0);
        }

        kernel = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * order);
    }

    UniformHDTransparentBC::~UniformHDTransparentBC() {
        for (size_t i = 0; i < nElec; i++){
            delete psis[i];
        }
        delete[] psis;
        sq_free(kernel);
    }

    void UniformHDTransparentBC::calcKernel(double vb, double dt){
        if(!kernelCalculated && dt == 0.0){ // standard calculation, use internal dt
            kernelCalculated = true;
            kernelVb = vb;
            dt = this->dt;
        }
        else if (dt != 0.0){ // particular calculation with specified dt, kernel will need to be re-calculated later
            kernelCalculated = false;
        }
        else if (abs(kernelVb-vb) > 1e-5/PhysCon::auE_ha)
            throw std::runtime_error("Potential at UniformHDTransparentBC is not constant. Consider using a different boundary condition.");
        else
            return;
        
        double rr = 4.0*dx*dx/dt;
        double sig = 2.0*dx*dx*vb;
        double phi = std::atan(2.0*rr*(sig+2.0)/(rr*rr-4.0*sig-sig*sig));
        double mu = (rr*rr+4*sig+sig*sig)/std::sqrt((rr*rr+sig*sig)*(rr*rr+(sig+4.0)*(sig+4.0)));
        std::complex<double> lam = std::exp(PhysCon::im * phi);
        std::complex<double> al = 0.5*PhysCon::im * std::exp(0.5*PhysCon::im*phi) * std::pow((rr*rr+sig*sig)*(rr*rr+(sig+4.0)*(sig+4.0)), 0.25);
        
        kernel0 = (1.0-PhysCon::im*rr/2.0 + sig/2.0) - al;
        kernel[0] = (1.0+PhysCon::im*rr/2.0 + sig/2.0) + al*std::exp(-PhysCon::im*phi) * mu;
        kernel[1] =  al * std::exp(-PhysCon::im*2.0*phi) * 0.5 * (mu*mu - 1.0);
        for (size_t i = 2; i < order; i++)
            kernel[i] = (2.0*i-1.0)/(i+1.0) * mu / lam * kernel[i-1] - (i-2.0)/(i+1.0) / (lam*lam) * kernel[i-2];
    }

    void UniformHDTransparentBC::getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, size_t nElec){
        if (this->nElec != nElec)
            throw std::invalid_argument("Number of electrons in HDTransparentBC does not match the number of electrons in the system.");

        for (size_t i = 0; i < nElec; i++)
            res[i] = psis[i]->inner(kernel) - psiad[i];
    }

    void UniformHDTransparentBC::prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb){
        for (size_t i = 0; i < nElec; i++)
            psis[i]->set(0, psibd[i]);

        calcKernel(vb);
    }

    void UniformHDTransparentBC::finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb){
        //std::complex<double> phs = 1.0/phasePerStep(vb);
        for (size_t i = 0; i < nElec; i++)
            psis[i]->stepBack();
	}

    void UniformHDTransparentBC::fillHistory(const std::complex<double>* psibd, const std::complex<double>* historicalPhaseAdvance, double vb) {
        std::complex<double> phs;
        for (size_t i = 0; i < nElec; i++){
            phs = 1.0;
            for (size_t j = 0; j < order; j++){
                psis[i]->set(j, psibd[i]*phs);
                phs /= historicalPhaseAdvance[i];
            }
        }
    };

    std::complex<double> UniformHDTransparentBC::getSteadyRHS(double kin){
        return 0.0; // homogeneous, no source
    }

    std::complex<double> UniformHDTransparentBC::getSteadyLHSEle(double kin){
        if(kin >= 0){
            double k = KineticOperators::CrankNicolson::wavenumberFromKineticEnergy(kin, dx*PhysCon::a0, 1.0);
            return 1.0 - 0.5*(std::exp(PhysCon::im*k*dx*PhysCon::a0)    + kin / (PhysCon::hbar*PhysCon::hbar/(2.0*PhysCon::me*dx*dx * PhysCon::a0 * PhysCon::a0)));
            //return 0.5*std::exp(-PhysCon::im*k*dx*PhysCon::a0);
        }
        else{
            double k = KineticOperators::CrankNicolson::wavenumberFromKineticEnergy(-kin, dx*PhysCon::a0, 1.0);
            return 1.0 - 0.5*(std::exp(-k*dx*PhysCon::a0)               + kin / (PhysCon::hbar*PhysCon::hbar/(2.0*PhysCon::me*dx*dx * PhysCon::a0 * PhysCon::a0)));
        }
    }

    std::complex<double> UniformHDTransparentBC::getSteadyLHSAdjEle(double kin){
        return -0.5;
    }

    std::complex<double> UniformHDTransparentBC::getSteadyRHS_PA(std::complex<double> phaseAdvance, double k0, double vb, double dt){
        return 0.0; // homogeneous, no source
    }

    std::complex<double> UniformHDTransparentBC::getSteadyLHSEle_PA(std::complex<double> phaseAdvance, double k0, double vb, double dt){
        calcKernel(vb, dt / PhysCon::hbar * PhysCon::auE_ha);

        std::complex<double> phs = phaseAdvance;
        std::complex<double> sum = -phs * kernel0;
        for (size_t i = 0; i < order; i++){
            phs /= phaseAdvance;
            sum -= phs * kernel[i];
        }
        return sum;
    }

    std::complex<double> UniformHDTransparentBC::getSteadyLHSAdjEle_PA(std::complex<double> phaseAdvance, double k0, double vb, double dt){
        return 1.0 + phaseAdvance;
    }

    UniformIDTransparentBC::UniformIDTransparentBC(size_t order, size_t nElec, double dx, double dt, const double* k0, double vb) : UniformHDTransparentBC(order, nElec, dx, dt) {
        phaseAdvance = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        phs = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        adjphs = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        ihpsi = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);

        hompsi = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);

        for (size_t i = 0; i < nElec; i++){
            // TODO: pull from input or from KineticOperator instead of CN result
            phaseAdvance[i] = 1.0-0.5*PhysCon::im*this->dt*(1.0/this->dx/this->dx*(1.0-std::cos(k0[i]*PhysCon::a0*this->dx)) + vb/PhysCon::auE_ha);
            phaseAdvance[i] /= std::conj(phaseAdvance[i]);

            adjphs[i] = std::exp(PhysCon::im*k0[i]*PhysCon::a0*this->dx);

            ihpsi[i] = 1.0;
        }
        std::fill_n(phs, nElec, 1.0);
    }

    void UniformIDTransparentBC::prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb){
        for (size_t i = 0; i < nElec; i++)
            psis[i]->set(0, psibd[i] - ihpsi[i]*phs[i]); // subtract off inhomogeneous component

        if(!kernelCalculated){
            calcKernel(vb);

            kernelCalculated = true;
            kernelVb = vb;
        }
        else if (abs(kernelVb-vb) > 1e-5*PhysCon::auE_ha)
            throw std::runtime_error("Potential at UniformHDTransparentBC is not constant. Consider using a different boundary condition.");
        
    }

    void UniformIDTransparentBC::getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, size_t nElec){
        if (this->nElec != nElec)
            throw std::invalid_argument("Number of electrons in IDTransparentBC does not match the number of electrons in the system.");

        for (size_t i = 0; i < nElec; i++)
            res[i] = psis[i]->inner(kernel) - psiad[i] + ihpsi[i]*phs[i]*(adjphs[i]+phaseAdvance[i]*(adjphs[i] - kernel0)); // add inhomogeneous component of present step
    }

    void UniformIDTransparentBC::finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb) {
        UniformHDTransparentBC::finishStep(psibd, psiad, vb);

        for (size_t i = 0; i < nElec; i++)
            phs[i] *= phaseAdvance[i];
    }

    void UniformIDTransparentBC::fillHistory(const std::complex<double>* psibd, const std::complex<double>* historicalPhaseAdvance, double vb) {
        std::complex<double>* dpsibd = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        for (size_t i = 0; i < nElec; i++)
            dpsibd[i] = psibd[i] - ihpsi[i];

        UniformHDTransparentBC::fillHistory(dpsibd, historicalPhaseAdvance, vb);
        
        sq_free(dpsibd);
    }

    std::complex<double> UniformIDTransparentBC::getSteadyRHS(double kin){
        if(kin >= 0){
            double k = KineticOperators::CrankNicolson::wavenumberFromKineticEnergy(kin, dx*PhysCon::a0, 1.0);
            return -PhysCon::im * std::sin(k * dx * PhysCon::a0);
        }
        else{
            double k = KineticOperators::CrankNicolson::wavenumberFromKineticEnergy(-kin, dx*PhysCon::a0, 1.0);
            return std::sinh(k * dx * PhysCon::a0);
        }
    }

    std::complex<double> UniformIDTransparentBC::getSteadyRHS_PA(std::complex<double> phaseAdvance, double k0, double vb, double dt){
        std::complex<double> sum = getSteadyLHSEle_PA(phaseAdvance, k0, vb, dt);
                                // Interpret negative k0 as exponential growth (thereby decaying in the external domain) instead of wave
        sum += (phaseAdvance+1.0) * std::exp((k0>=0.0 ? PhysCon::im : 1.0)*k0*dx*PhysCon::a0);

        return sum;
    }
}

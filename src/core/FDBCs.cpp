#include "CORECommonHeader.h"
#include "MathTools.h"
#include "FDBCs.h"

namespace FDBCs{
    UniformHDTransparentBC::UniformHDTransparentBC(int order, int nElec, double dx, double dt) : order(order), nElec(nElec), dx(dx / PhysCon::a0), dt(dt / PhysCon::hbar * PhysCon::auE_ha) {
        if (order < 1)
            throw std::invalid_argument("Order of HDTransparentBC must be greater than 0.");

        psis = new CyclicArray<std::complex<double>>*[nElec];
        for (int i = 0; i < nElec; i++){
            psis[i] = new CyclicArray<std::complex<double>>(order, 0.0);
        }

        kernel = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * order);
    }

    UniformHDTransparentBC::~UniformHDTransparentBC() {
        for (int i = 0; i < nElec; i++){
            delete psis[i];
        }
        delete[] psis;
        sq_free(kernel);
    }

    void UniformHDTransparentBC::calcKernel(double vb){
        if(!kernelCalculated){
            kernelCalculated = true;
            kernelVb = vb;
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
        for (int i = 2; i < order; i++)
            kernel[i] = (2.0*i-1.0)/(i+1.0) * mu / lam * kernel[i-1] - (i-2.0)/(i+1.0) / (lam*lam) * kernel[i-2];
    }

    /*void UniformHDTransparentBC::calcCVPsis(double vb){
        for (int i = 0; i < nElec; i++){
            std::complex<double> phs = 1.0;
            for (int j = 0; j < order-1; j++){
                cvpsis[i]->set(j, phs * psis[i]->get(j));
                phs *= (2.0 + PhysCon::im*dt*(vb - vbs->get(j+1))*(1.0+psis[i]->get(j)/psis[i]->get(j+1))) /\
                    (2.0 - PhysCon::im*dt*(vb - vbs->get(j+1))*(1.0+psis[i]->get(j+1)/psis[i]->get(j)));
            }
            cvpsis[i]->set(order-1, phs * psis[i]->get(order-1));
        }
    }*/

    void UniformHDTransparentBC::getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, int nElec){
        if (this->nElec != nElec)
            throw std::invalid_argument("Number of electrons in HDTransparentBC does not match the number of electrons in the system.");

        for (int i = 0; i < nElec; i++)
            res[i] = psis[i]->inner(kernel) - psiad[i];
    }

    void UniformHDTransparentBC::prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb){
        for (int i = 0; i < nElec; i++)
            psis[i]->set(0, psibd[i]);

        calcKernel(vb);
    }

    void UniformHDTransparentBC::finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb){
        //std::complex<double> phs = 1.0/phasePerStep(vb);
        for (int i = 0; i < nElec; i++)
            psis[i]->stepBack();
	}

    void UniformHDTransparentBC::fillHistory(const std::complex<double>* psibd, const std::complex<double>* historicalPhaseAdvance, double vb) {
        std::complex<double> phs;
        for (int i = 0; i < nElec; i++){
            phs = 1.0;
            for (int j = 0; j < order; j++){
                psis[i]->set(j, psibd[i]*phs);
                phs /= historicalPhaseAdvance[i];
            }
        }
    };

    std::complex<double> UniformHDTransparentBC::getSteadyRHS(std::complex<double> phaseAdvance, double k0, double vb){
        return 0.0; // homogeneous, no source
    }

    std::complex<double> UniformHDTransparentBC::getSteadyLHSEle(std::complex<double> phaseAdvance, double k0, double vb){
        calcKernel(vb);

        std::complex<double> phs = phaseAdvance;
        std::complex<double> sum = -phs * kernel0;
        for (int i = 0; i < order; i++){
            phs /= phaseAdvance;
            sum -= phs * kernel[i];
        }
        return sum;
    }

    std::complex<double> UniformHDTransparentBC::getSteadyLHSAdjEle(std::complex<double> phaseAdvance, double k0, double vb){
        return 1.0 + phaseAdvance;
    }

    UniformIDTransparentBC::UniformIDTransparentBC(int order, int nElec, double dx, double dt, std::complex<double>* psibd, double* k0, double vb) : UniformHDTransparentBC(order, nElec, dx, dt) {
        phaseAdvance = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        phs = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        adjphs = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        ihpsi = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);

        hompsi = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);

        for (int i = 0; i < nElec; i++){
            // TODO: pull from input or from KineticOperator instead of CN result
            phaseAdvance[i] = 1.0-0.5*PhysCon::im*this->dt*(1.0/this->dx/this->dx*(1.0-std::cos(k0[i]*PhysCon::a0*this->dx)) + vb/PhysCon::auE_ha);
            phaseAdvance[i] /= std::conj(phaseAdvance[i]);

            adjphs[i] = std::exp(PhysCon::im*k0[i]*PhysCon::a0*this->dx);

            ihpsi[i] = psibd[i];
        }
        std::fill_n(phs, nElec, 1.0);
    }

    void UniformIDTransparentBC::prepareStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb){
        for (int i = 0; i < nElec; i++)
            psis[i]->set(0, psibd[i] - ihpsi[i]*phs[i]); // subtract off inhomogeneous component

        if(!kernelCalculated){
            calcKernel(vb);

            kernelCalculated = true;
            kernelVb = vb;
        }
        else if (abs(kernelVb-vb) > 1e-5*PhysCon::auE_ha)
            throw std::runtime_error("Potential at UniformHDTransparentBC is not constant. Consider using a different boundary condition.");
        
    }

    void UniformIDTransparentBC::getRHS(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb, std::complex<double>* res, int nElec){
        if (this->nElec != nElec)
            throw std::invalid_argument("Number of electrons in IDTransparentBC does not match the number of electrons in the system.");

        for (int i = 0; i < nElec; i++)
            res[i] = psis[i]->inner(kernel) - psiad[i] + ihpsi[i]*phs[i]*(adjphs[i]+phaseAdvance[i]*(adjphs[i] - kernel0)); // add inhomogeneous component of present step
    }

    void UniformIDTransparentBC::finishStep(const std::complex<double>* psibd, const std::complex<double>* psiad, double vb) {
        UniformHDTransparentBC::finishStep(psibd, psiad, vb);

        for (int i = 0; i < nElec; i++)
            phs[i] *= phaseAdvance[i];
    }

    void UniformIDTransparentBC::fillHistory(const std::complex<double>* psibd, const std::complex<double>* historicalPhaseAdvance, double vb) {
        std::complex<double>* dpsibd = (std::complex<double>*)sq_malloc(sizeof(std::complex<double>) * nElec);
        for (int i = 0; i < nElec; i++)
            dpsibd[i] = psibd[i] - ihpsi[i];

        UniformHDTransparentBC::fillHistory(dpsibd, historicalPhaseAdvance, vb);
        
        sq_free(dpsibd);
    }

    std::complex<double> UniformIDTransparentBC::getSteadyRHS(std::complex<double> phaseAdvance, double k0, double vb){
        calcKernel(vb);

        std::complex<double> ihPhs = phaseAdvance;
        std::complex<double> sum = 0.0;
        for (int i = 0; i < order; i++){
            ihPhs /= phaseAdvance;
            sum -= ihPhs * kernel[i];
        }
                                        // Interpret negative k0 as exponential growth (thereby decaying in the external domain) instead of wave
        sum += (1.0+phaseAdvance) * std::exp((k0>=0.0 ? PhysCon::im : 1.0)*k0*dx*PhysCon::a0) - phaseAdvance*kernel0;

        return sum;
    }
}

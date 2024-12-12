#include "CORECommonHeader.h"
#include "MathTools.h"
#include "FDBCs.h"

namespace FDBCs{
    HDTransparentBC::HDTransparentBC(int order, int nElec, double dx, double dt) : order(order), nElec(nElec), dx(dx / PhysCon::a0), dt(dt / PhysCon::hbar * PhysCon::auE_ha) {
        if (order < 1)
            throw std::invalid_argument("Order of HDTransparentBC must be greater than 0.");

        psis = new CyclicArray<std::complex<double>>*[nElec];
        for (int i = 0; i < nElec; i++)
            psis[i] = new CyclicArray<std::complex<double>>(order, 0.0);
        kernel = new CyclicArray<std::complex<double>>(order, 0.0);
    }

    HDTransparentBC::~HDTransparentBC() {
        for (int i = 0; i < nElec; i++)
            delete psis[i];
        delete[] psis;
        delete kernel;
    }

    void HDTransparentBC::calcKernel(double vb){
        double rr = 4.0*dx*dx/dt;
        double sig = 2*dx*dx*vb;
        double phi = std::atan(2*rr*(sig+2.0)/(rr*rr-4.0*sig-sig*sig));
        double mu = (rr*rr+4*sig+sig*sig)/std::sqrt((rr*rr+sig*sig)*(rr*rr+(sig+4.0)*(sig+4.0)));
        std::complex<double> lam = std::exp(PhysCon::im * phi);
        std::complex<double> al = 0.5*PhysCon::im * std::exp(0.5*PhysCon::im*phi) * std::pow((rr*rr+sig*sig)*(rr*rr+(sig+4.0)*(sig+4.0)), 0.25);
        
        kernel0 = (1.0-PhysCon::im*rr/2.0 + sig/2.0) - al;
        kernel->set(0, (1.0+PhysCon::im*rr/2.0 + sig/2.0) + al*std::exp(-PhysCon::im*phi) * mu);
        kernel->set(1, al * std::exp(-PhysCon::im*2.0*phi) * 0.5 * (mu*mu - 1.0));
        for (int i = 2; i < order; i++)
            kernel->set(i, (2.0*i-1.0)/(i+1.0) * mu / lam * kernel->get(i-1) - (i-2.0)/(i+1.0) / (lam*lam) * kernel->get(i-2));
        
        for (int i = 0; i < order; i++)
            kernel->set(i, kernel->get(i) * std::exp(-PhysCon::im/PhysCon::hbar*(vb*i*dt)));
    }

    void HDTransparentBC::getRHS(std::complex<double>* psibd, std::complex<double>* psiad, double vb, std::complex<double>* res, int nElec){
        if (this->nElec != nElec)
            throw std::invalid_argument("Number of electrons in HDTransparentBC does not match the number of electrons in the system.");

        for (int i = 0; i < nElec; i++)
            res[i] = kernel->inner(psis[i]) - psiad[i];
    }
}

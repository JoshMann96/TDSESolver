#include "CORECommonHeader.h"
#include "PhysCon.h"
#include "FDBCs.h"

namespace FD_BCs{
    HDTransparentBC::HDTransparentBC(int order, double dx, double dt) : order(order), dx(dx / PhysCon::a0), dt(dt * PhysCon::hbar / PhysCon::auE_ha) {
        if (order < 1)
            throw std::invalid_argument("Order of HDTransparentBC must be greater than 0.");

        psis = new CyclicArray<std::complex<double>>(order);
        kernel = new CyclicArray<std::complex<double>>(order);
    }

    HDTransparentBC::~HDTransparentBC() {
        sq_free(psis);
        sq_free(kernel);
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
    }

}

int main(int argc, char** argv){
    FD_BCs::HDTransparentBC *bc = new FD_BCs::HDTransparentBC(10, 0.1*PhysCon::a0, 0.1*PhysCon::auE_ha / PhysCon::hbar);
    bc->calcKernel(3.0);
    bc->printKernel();
    // (0.921767,0.313862) (0.229096,0.196664) (-0.0452012,0.0583993) (0.0434235,0.0302205) (-0.0252054,0.0404439) (0.00784178,0.00435176) (-0.0135082,0.0275446) (-0.00605452,-0.00259745) (-0.00616852,0.0166615) (-0.0114493,-0.00359059) 
    delete bc;
}
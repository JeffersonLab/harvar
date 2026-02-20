//--------- code
// compile --- g++ -O2 -std=c++17 twobody.cc -I. `root-config --cflags --glibs` -o twobody
//
#include "TwoBodyDecay.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <iostream>

int main() {
    // 1) Define a parent four-vector (e.g., ρ meson in lab frame)
    //    Here we choose a ρ at rest with mass = 0.770 GeV
    TVector3 v(0.6,0.4,1.0);
  TLorentzVector rho;//(0.0, 1.0, 1.0, std::sqrt(0.770*0.770));
  rho.SetVectM(v,0.77);
    // 2) Daughter masses (e.g., pions)
    double m_pi = 0.13957;

    // 3) Reference vector to define x-axis in rest frame
    //    e.g., lab z-axis or some production plane normal
    TVector3 refVec(0.5, 0.5, 0.5);

    // 4) Instantiate and generate decay
    TwoBodyDecay decay(rho, m_pi, m_pi);
    decay.Generate(refVec);

    // 5) Output daughter four-vectors and angles
    std::cout << "Daughter 1 (π⁺): "; decay.daughter1.Print();
    std::cout << "Daughter 2 (π⁻): "; decay.daughter2.Print();
    std::cout << "Generated angles in rest frame:\n";
    std::cout << "  theta = " << decay.theta * 180.0/M_PI << " deg\n";
    std::cout << "  phi   = " << decay.phi   * 180.0/M_PI << " deg\n";

    // 6) Recompute angles to verify using CalculateAngles()
    decay.CalculateAngles(refVec);
    std::cout << "Recomputed angles via CalculateAngles():\n";
    std::cout << "  theta = " << decay.theta * 180.0/M_PI << " deg\n";
    std::cout << "  phi   = " << decay.phi   * 180.0/M_PI << " deg\n";
    
    std::cout << "  mass = " << rho.M() << std::endl;
    return 0;
}

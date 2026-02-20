#include "w_kernels.hpp"
#include <iostream>
using namespace Wkernels;

int main()
{
    Mat4 u{}, l{}, s{};          // Fill these with real data!
    double eps = 0.7;
    double phi = 0.5;            //   ϕ  in radians
    double kap = 0.3;            //   φ  in radians

    u[h('0')][h('0')][h('+')][h('+')] = {0.0123, 0.0000};   // real example
    u[h('0')][h('0')][h('0')][h('+')] = {0.0150, 0.0240};   // real example
    l[h('+')][h('-')][h('0')][h('+')] = {0.0,    -0.0071};  // imaginary
    
    for(double phi2 = -3.14; phi2<3.14; phi2+=0.2){
        phi = phi2;
        auto wUU = UU(u, eps, phi, kap);
        std::cout <<  "PHI = " <<  phi << "  W_UU^LL = " << wUU.LL << '\n';
    }

}

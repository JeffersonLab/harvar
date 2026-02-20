/**********************************************************************
 *  w_kernels.hpp   –   Angular kernels W^{LL,LT,TT}_{XY} calculator
 *
 *  Implements all real/imaginary combinations from Eqs. (15)–(19)
 *  in the “ρ-electroproduction” PDF.  Only the u, l (and optional s)
 *  3×3×3×3 complex matrices are required.
 *
 *  Usage:
 *      #include "w_kernels.hpp"
 *
 *      Mat4  u, l, s;          // fill these with your numbers
 *      double eps   = ...;     // ε(y,γ)   – Eq.(3)
 *      double phi   = ...;     // ϕ   (lepton–hadron plane)
 *      double kap   = ...;     // φ   (π⁺ decay plane)
 *
 *      auto wUU = Wkernels::UU(u, eps, phi, kap);
 *      std::cout << "W_UU^LL  = " << wUU.LL << '\n';
 *********************************************************************/
 #pragma once
 #include <complex>
 #include <cmath>
 #include <array>
 
 namespace Wkernels
 {
     using  C   = std::complex<double>;
     using Mat4 = std::array<std::array<std::array<std::array<C,3>,3>,3>,3>;
 
     // index helper: '+'=>0, '0'=>1, '-'=>2
     constexpr int h(char c)      { return c=='+'?0 : c=='0'?1 : 2; }
 
     // shorthand accessor -----------------------------------------------------
     inline const C& U(const Mat4& M, char nu, char nup, char mu, char mup)
     {   return M[h(nu)][h(nup)][h(mu)][h(mup)]; }
 
     struct Triple { double LL{}, LT{}, TT{}; };   // returned kernels
 
     // =======================================================================
     //  Unpolarised target, unpolarised beam  (Eq. 15)
     // =======================================================================
     inline Triple UU(const Mat4& u,
                      double eps, double phi, double kap)
     {
         double ce      = std::sqrt(eps*(1.0+eps));
         double cphi    = std::cos(phi);
         double c2phi   = std::cos(2.0*phi);
         double sphi    = std::sin(phi);            // needed later
         double cpk     = std::cos(kap);
         double cphi_k  = std::cos(phi+kap);
         double cphi_mk = std::cos(phi-kap);
         double c2phi_k = std::cos(2.0*phi+kap);
         double c2phi_mk= std::cos(2.0*phi-kap);
         double c2phi2k = std::cos(2.0*phi+2.0*kap);
         double cm2k    = std::cos(2.0*kap);
 
         Triple out;
 
         // ---------- LL ------------------------------------------------------
         out.LL  =  std::real(U(u,'0','0','+','+'))
                 + eps * std::real(U(u,'0','0','0','0'))
                 - 2.0 * cphi * ce * std::real(U(u,'0','0','0','+'))
                 - eps * c2phi * std::real(U(u,'0','0','-','+'));
 
         // ---------- LT ------------------------------------------------------
         out.LT  =  cphi_k * ce *
                    std::real( U(u,'0','+','0','+') - U(u,'-','0','0','+') )
                 - cpk *
                    std::real( U(u,'0','+','+','+')
                             - U(u,'-','0','+','+')
                             + 2.0*eps*U(u,'0','+','0','0') )
                 + c2phi_k * eps * std::real( U(u,'0','+','-','+') )
                 - cphi_mk * ce *
                    std::real( U(u,'0','-','0','+') - U(u,'+','0','0','+') )
                 + c2phi_mk * eps * std::real( U(u,'+','0','-','+') );
 
         // ---------- TT ------------------------------------------------------
         double cphi2k_p = std::cos(phi+2.0*kap);
         double cphi2k_m = std::cos(phi-2.0*kap);
         out.TT  = 0.5 * ( std::real(U(u,'+','+','+','+'))
                         + std::real(U(u,'-','-','+','+'))
                         + 2.0*eps*std::real(U(u,'+','+','0','0')) )
                 + 0.5 * c2phi2k * eps *
                     std::real( U(u,'-','+','-','+') )
                 - cphi * ce *
                     std::real( U(u,'+','+','0','+')
                              + U(u,'-','-','0','+') )
                 + cphi2k_p * ce *
                     std::real( U(u,'-','+','0','+') )
                 - cm2k *
                     std::real( U(u,'-','+','+','+')
                              + eps*U(u,'-','+','0','0') )
                 - c2phi * eps *
                     std::real( U(u,'+','+','-','+') )
                 + cphi2k_m * ce *
                     std::real( U(u,'+','-','0','+') )
                 + 0.5 * c2phi_mk * eps *
                     std::real( U(u,'+','-','-','+') );
 
         return out;
     }
 
     // =======================================================================
     //  Longitudinally polarised beam, unpolarised target  (Eq. 17)
     // =======================================================================
     inline Triple LU(const Mat4& u,
                      double eps, double phi, double kap)
     {
         double ce_m  = std::sqrt(eps*(1.0-eps));
         double sphi  = std::sin(phi);
         double spkp  = std::sin(phi+kap);
         double spkm  = std::sin(phi-kap);
         double sp2kp = std::sin(2.0*phi+kap);
 
         Triple out;
 
         out.LL = -2.0 * sphi * ce_m * std::imag( U(u,'0','0','0','+') );
 
         out.LT =  spkp * ce_m *
                     std::imag( U(u,'0','+','0','+') - U(u,'-','0','0','+') )
                 - std::sin(kap) * std::sqrt(1.0-eps*eps) *
                     std::imag( U(u,'0','+','+','+')
                              - U(u,'-','0','+','+') )
                 - spkm * ce_m *
                     std::imag( U(u,'0','-','0','+') - U(u,'+','0','0','+') );
 
         out.TT = - sphi * ce_m *
                     std::imag( U(u,'+','+','0','+') + U(u,'-','-','0','+') )
                 + sp2kp * ce_m *
                     std::imag( U(u,'-','+','0','+') )
                 - std::sin(2.0*kap) * std::sqrt(1.0-eps*eps) *
                     std::imag( U(u,'-','+','+','+') )
                 + std::sin(phi-2.0*kap) * ce_m *
                     std::imag( U(u,'+','-','0','+') );
 
         return out;
     }
 
     // =======================================================================
     //  Longitudinally polarised target, unpolarised beam  (Eq. 18)
     // =======================================================================
     inline Triple UL(const Mat4& l,
                      double eps, double phi, double kap)
     {
         double ce_p  = std::sqrt(eps*(1.0+eps));
         double sphi  = std::sin(phi);
         double spkp  = std::sin(phi+kap);
         double spkm  = std::sin(phi-kap);
         double sp2k  = std::sin(2.0*kap);
         double sp2kp = std::sin(2.0*phi+kap);
         double sp2km = std::sin(2.0*phi-kap);
 
         Triple out;
 
         out.LL = -2.0 * sphi * ce_p * std::imag( U(l,'0','0','0','+') )
                - eps * std::sin(2.0*phi) * std::imag( U(l,'0','0','-','+') );
 
         out.LT =  spkp * ce_p *
                     std::imag( U(l,'0','+','0','+') - U(l,'-','0','0','+') )
                 - std::sin(kap) *
                     std::imag( U(l,'0','+','+','+')
                              - U(l,'-','0','+','+')
                              + 2.0*eps*U(l,'0','+','0','0') )
                 + sp2kp * eps *
                     std::imag( U(l,'0','+','-','+') )
                 - spkm * ce_p *
                     std::imag( U(l,'0','-','0','+') - U(l,'+','0','0','+') )
                 + sp2km * eps *
                     std::imag( U(l,'+','0','-','+') );
 
         out.TT = 0.5 * sp2kp * eps *
                     std::imag( U(l,'-','+','-','+') )
                 - sphi * ce_p *
                     std::imag( U(l,'+','+','0','+') + U(l,'-','-','0','+') )
                 + sp2kp * ce_p *
                     std::imag( U(l,'-','+','0','+') )
                 - sp2k *
                     std::imag( U(l,'-','+','+','+')
                              + eps*U(l,'-','+','0','0') )
                 - eps * std::sin(2.0*phi) *
                     std::imag( U(l,'+','+','-','+') )
                 + sp2km * ce_p *
                     std::imag( U(l,'+','-','0','+') )
                 + 0.5 * sp2km * eps *
                     std::imag( U(l,'+','-','-','+') );
 
         return out;
     }
 
     // =======================================================================
     //  Longitudinal beam × Longitudinal target (Eq. 19)
     // =======================================================================
     inline Triple LL(const Mat4& l,
                      double eps, double phi, double kap)
     {
         double ce_m  = std::sqrt(eps*(1.0-eps));
         double cphi  = std::cos(phi);
         double cpk   = std::cos(kap);
         double cphi_k= std::cos(phi+kap);
         double cphi_mk=std::cos(phi-kap);
         double c2k   = std::cos(2.0*kap);
         double c2kp  = std::cos(phi+2.0*kap);
         double c2km  = std::cos(phi-2.0*kap);
 
         Triple out;
 
         out.LL = -2.0 * cphi * ce_m * std::real( U(l,'0','0','0','+') )
                + std::sqrt(1.0-eps*eps) *
                    std::real( U(l,'0','0','+','+') );
 
         out.LT =  cphi_k * ce_m *
                     std::real( U(l,'0','+','0','+') - U(l,'-','0','0','+') )
                 - cpk *
                     std::real( U(l,'0','+','+','+')
                              - U(l,'-','0','+','+')
                              + 2.0*eps*U(l,'0','+','0','0') )
                 + c2kp * eps *
                     std::real( U(l,'0','+','-','+') )
                 - cphi_mk * ce_m *
                     std::real( U(l,'0','-','0','+') - U(l,'+','0','0','+') )
                 + c2km * eps *
                     std::real( U(l,'+','0','-','+') );
 
         out.TT =  std::sqrt(1.0 - eps*eps) * 0.5 *
                     ( std::real(U(l,'+','+','+','+'))
                     + std::real(U(l,'-','-','+','+')) )
                 - cphi * ce_m *
                     std::real( U(l,'+','+','0','+') + U(l,'-','-','0','+') )
                 + c2kp * ce_m *
                     std::real( U(l,'-','+','0','+') )
                 - c2k *
                     std::real( U(l,'-','+','+','+')
                              + eps*U(l,'-','+','0','0') )
                 + c2km * ce_m *
                     std::real( U(l,'+','-','0','+') );
 
         return out;
     }

    // =======================================================================
    //  Unpolarized beam × Transverse target (Eq. 18)
    // =======================================================================
    inline Triple UT(const Mat4& n, const Mat4& s,
                     double eps, double phiS, double phi, double kap)
    {
        double ce_p = std::sqrt(eps*(1.0+eps));
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        double c2phi = std::cos(2.0*phi);
        double s2phi = std::sin(2.0*phi);
        double cpk = std::cos(kap);
        double spk = std::sin(kap);
        double c2k = std::cos(2.0*kap);
        double s2k = std::sin(2.0*kap);
        double cphi_k = std::cos(phi+kap);
        double sphi_k = std::sin(phi+kap);
        double cphi_mk = std::cos(phi-kap);
        double sphi_mk = std::sin(phi-kap);
        double c2phi_k = std::cos(2.0*phi+kap);
        double s2phi_k = std::sin(2.0*phi+kap);
        double c2phi_mk = std::cos(2.0*phi-kap);
        double s2phi_mk = std::sin(2.0*phi-kap);
        double cphi_2k = std::cos(phi+2.0*kap);
        double sphi_2k = std::sin(phi+2.0*kap);
        double cphi_m2k = std::cos(phi-2.0*kap);
        double sphi_m2k = std::sin(phi-2.0*kap);
        double c2phi_2k = std::cos(2.0*phi+2.0*kap);
        double s2phi_2k = std::sin(2.0*phi+2.0*kap);
        double c2phi_m2k = std::cos(2.0*phi-2.0*kap);
        double s2phi_m2k = std::sin(2.0*phi-2.0*kap);

        double s_phi_phiS = std::sin(phi - phiS);
        double c_phi_phiS = std::cos(phi - phiS);

        Triple out;

        // ---------- LL ------------------------------------------------------
        out.LL = s_phi_phiS * (
                   std::imag(U(n,'0','0','+','+'))
                 + eps * std::imag(U(n,'0','0','0','0'))
                 - 2.0*cphi*ce_p * std::imag(U(n,'0','0','0','+'))
                 - eps*c2phi * std::imag(U(n,'0','0','-','+'))
                 )
               + c_phi_phiS * (
                 -2.0*sphi*ce_p * std::imag(U(s,'0','0','0','+'))
                 - eps*s2phi * std::imag(U(s,'0','0','-','+'))
                 );

        // ---------- LT ------------------------------------------------------
        out.LT = s_phi_phiS * (
                   cphi_k*ce_p * std::imag(U(n,'0','+','0','+') - U(n,'-','0','0','+'))
                 - cpk * std::imag(U(n,'0','+','+','+') - U(n,'-','0','+','+')
                                 + 2.0*eps*U(n,'0','+','0','0'))
                 + c2phi_k*eps * std::imag(U(n,'0','+','-','+'))
                 - cphi_mk*ce_p * std::imag(U(n,'0','-','0','+') - U(n,'+','0','0','+'))
                 + c2phi_mk*eps * std::imag(U(n,'+','0','-','+'))
                 )
               + c_phi_phiS * (
                   sphi_k*ce_p * std::imag(U(s,'0','+','0','+') - U(s,'-','0','0','+'))
                 - spk * std::imag(U(s,'0','+','+','+') - U(s,'-','0','+','+')
                                 + 2.0*eps*U(s,'0','+','0','0'))
                 + s2phi_k*eps * std::imag(U(s,'0','+','-','+'))
                 - sphi_mk*ce_p * std::imag(U(s,'0','-','0','+') - U(s,'+','0','0','+'))
                 + s2phi_mk*eps * std::imag(U(s,'+','0','-','+'))
                 );

        // ---------- TT ------------------------------------------------------
        out.TT = s_phi_phiS * (
                   0.5 * std::imag(U(n,'+','+','+','+') + U(n,'-','-','+','+')
                                 + 2.0*eps*U(n,'+','+','0','0'))
                 + 0.5*c2phi_2k*eps * std::imag(U(n,'-','+','-','+'))
                 - cphi*ce_p * std::imag(U(n,'+','+','0','+') + U(n,'-','-','0','+'))
                 + cphi_2k*ce_p * std::imag(U(n,'-','+','0','+'))
                 - c2k * std::imag(U(n,'-','+','+','+') + eps*U(n,'-','+','0','0'))
                 - c2phi*eps * std::imag(U(n,'+','+','-','+'))
                 + cphi_m2k*ce_p * std::imag(U(n,'+','-','0','+'))
                 + 0.5*c2phi_m2k*eps * std::imag(U(n,'+','-','-','+'))
                 )
               + c_phi_phiS * (
                   0.5*s2phi_2k*eps * std::imag(U(s,'-','+','-','+'))
                 - sphi*ce_p * std::imag(U(s,'+','+','0','+') + U(s,'-','-','0','+'))
                 + sphi_2k*ce_p * std::imag(U(s,'-','+','0','+'))
                 - s2k * std::imag(U(s,'-','+','+','+') + eps*U(s,'-','+','0','0'))
                 - s2phi*eps * std::imag(U(s,'+','+','-','+'))
                 + sphi_m2k*ce_p * std::imag(U(s,'+','-','0','+'))
                 + 0.5*s2phi_m2k*eps * std::imag(U(s,'+','-','-','+'))
                 );

        return out;
    }

    // =======================================================================
    //  Longitudinally polarized beam × Transverse target (Eq. 19)
    // =======================================================================
    inline Triple LT(const Mat4& n, const Mat4& s,
                     double eps, double phiS, double phi, double kap)
    {
        double ce_m = std::sqrt(eps*(1.0-eps));
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);
        double cpk = std::cos(kap);
        double c2k = std::cos(2.0*kap);
        double cphi_k = std::cos(phi+kap);
        double sphi_k = std::sin(phi+kap);
        double cphi_mk = std::cos(phi-kap);
        double sphi_mk = std::sin(phi-kap);
        double cphi_2k = std::cos(phi+2.0*kap);
        double sphi_2k = std::sin(phi+2.0*kap);
        double cphi_m2k = std::cos(phi-2.0*kap);
        double sphi_m2k = std::sin(phi-2.0*kap);

        double s_phi_phiS = std::sin(phi - phiS);
        double c_phi_phiS = std::cos(phi - phiS);

        Triple out;

        // ---------- LL ------------------------------------------------------
        out.LL = s_phi_phiS * (
                   2.0*sphi*ce_m * std::real(U(n,'0','0','0','+'))
                 )
               + c_phi_phiS * (
                 -2.0*cphi*ce_m * std::real(U(s,'0','0','0','+'))
                 + std::sqrt(1.0-eps*eps) * std::real(U(s,'0','0','+','+'))
                 );

        // ---------- LT ------------------------------------------------------
        out.LT = s_phi_phiS * (
                 -sphi_k*ce_m * std::real(U(n,'0','+','0','+') - U(n,'-','0','0','+'))
                 + std::sin(kap)*std::sqrt(1.0-eps*eps) *
                   std::real(U(n,'0','+','+','+') - U(n,'-','0','+','+'))
                 + sphi_mk*ce_m * std::real(U(n,'0','-','0','+') - U(n,'+','0','0','+'))
                 )
               + c_phi_phiS * (
                   cphi_k*ce_m * std::real(U(s,'0','+','0','+') - U(s,'-','0','0','+'))
                 - cpk*std::sqrt(1.0-eps*eps) *
                   std::real(U(s,'0','+','+','+') - U(s,'-','0','+','+'))
                 - cphi_mk*ce_m * std::real(U(s,'0','-','0','+') - U(s,'+','0','0','+'))
                 );

        // ---------- TT ------------------------------------------------------
        out.TT = s_phi_phiS * (
                   sphi*ce_m * std::real(U(n,'+','+','0','+') + U(n,'-','-','0','+'))
                 - sphi_2k*ce_m * std::real(U(n,'-','+','0','+'))
                 + std::sin(2.0*kap)*std::sqrt(1.0-eps*eps) * std::real(U(n,'-','+','+','+'))
                 - sphi_m2k*ce_m * std::real(U(n,'+','-','0','+'))
                 )
               + c_phi_phiS * (
                   std::sqrt(1.0-eps*eps) * 0.5 *
                   (std::real(U(s,'+','+','+','+')) + std::real(U(s,'-','-','+','+')))
                 - cphi*ce_m * std::real(U(s,'+','+','0','+') + U(s,'-','-','0','+'))
                 + cphi_2k*ce_m * std::real(U(s,'-','+','0','+'))
                 - c2k*std::sqrt(1.0-eps*eps) * std::real(U(s,'-','+','+','+'))
                 + cphi_m2k*ce_m * std::real(U(s,'+','-','0','+'))
                 );

        return out;
    }
 }
 
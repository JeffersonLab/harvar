/********************************************************************
 *  xsection_v1.cpp  –  Formula-(4) cross-section evaluator
 *  ---------------------------------------------------------------
 *  Compile with:   g++ -std=c++17 -O2 xsection_v1.cpp -o xsec
 *******************************************************************/
 #include <iostream>
 #include <iomanip>
 #include <cmath>
 
 // ------------------------------------------------------------------
 // Physical constants
 // ------------------------------------------------------------------
 constexpr double alpha_em = 1.0 / 137.035999084;     // fine structure
 constexpr double M_N      = 0.9382720813;            // nucleon mass [GeV]
 
 // ------------------------------------------------------------------
 // Helper functions
 // ------------------------------------------------------------------
 double gamma_bj(double xB, double Q2)                     // Eq. (4)
 {
     return 2.0 * xB * M_N / std::sqrt(Q2);
 }
 
 double epsilon(double y, double gamma)                    // Eq. (3)
 {
     double y2g2 = y * y * gamma * gamma;
     return (1.0 - y - 0.25 * y2g2) /
            (1.0 - y + 0.5 * y * y + 0.25 * y2g2);
 }
 
 double dsigma_3fold(double xB, double Q2, double y,
                     double eps,
                     double dSigmaT_dt, double dSigmaL_dt) // Eq. (5)
 {
     double pref = alpha_em / (2.0 * M_PI);
     double kin  = (y * y) / (1.0 - eps) * (1.0 - xB) / xB / Q2;
     return pref * kin * (dSigmaT_dt + eps * dSigmaL_dt);
 }
 
 double dsigma_7fold(double WUU, double WLU, double WUL,
                     double WLL, double WUT, double WLT,
                     double Pl,  double SL,  double ST,
                     double dsigma3)                        // Eq. (4)
 {
     double S = WUU + Pl * WLU + SL * WUL + Pl * SL * WLL
                + ST * WUT + Pl * ST * WLT;
 
     return (1.0 / (4.0 * M_PI * M_PI)) * dsigma3 * S;      // (2π)² = 4π²
 }
 
 // ------------------------------------------------------------------
 // Main program
 // ------------------------------------------------------------------
 int main()
 {
     // ----------- user inputs --------------------------------------
     double xB, Q2,  y,  t;
     double dSigmaT_dt, dSigmaL_dt;
     double Pl, SL, ST;
     double WUU, WLU, WUL, WLL, WUT, WLT;
 
     std::cout << "Kinematics  (xB  Q2[GeV²]  y  t[GeV²]): ";
     if(!(std::cin >> xB >> Q2 >> y >> t)) return 0;
 
     std::cout << "Partial cross sections  (dσ_T/dt  dσ_L/dt) [nb/GeV²]: ";
     std::cin >> dSigmaT_dt >> dSigmaL_dt;
 
     std::cout << "Polarisations  (P_l  S_L  S_T): ";
     std::cin >> Pl >> SL >> ST;
 
     std::cout << "Six W-terms  (WUU WLU WUL WLL WUT WLT): ";
     std::cin >> WUU >> WLU >> WUL >> WLL >> WUT >> WLT;
 
     // ----------- calculations -------------------------------------
     double g      = gamma_bj(xB, Q2);
     double eps    = epsilon(y, g);
     double ds3    = dsigma_3fold(xB, Q2, y, eps, dSigmaT_dt, dSigmaL_dt);
     double ds7    = dsigma_7fold(WUU, WLU, WUL, WLL, WUT, WLT,
                                  Pl, SL, ST, ds3);
 
     // ----------- output -------------------------------------------
     std::cout << "\nComputed quantities\n"
               << "---------------------------------------------\n"
               << std::fixed << std::setprecision(6)
               << "gamma        = " << g   << '\n'
               << "epsilon       " << std::setprecision(8)
               << "= " << eps << '\n'
               << std::setprecision(6)
               << "dσ/dxBdQ²dt  = " << ds3 << "  [nb / GeV⁴]\n"
               << "d⁷σ (full)   = " << ds7 << "  [nb / (GeV⁴ sr²)]\n";
 }
 
// RhoLundIO.h ─────────────────────────────────────────────────────────
// Helper for printing a RhoEvent in LUND format (CLAS convention)

#ifndef RHO_CROSS_H
#define RHO_CROSS_H

#include "RhoEvent.h"
#include <iomanip>
#include <ostream>
#include <vector>

#include <iomanip>
#include <vector>
#include "w_kernels.hpp"
#include "AmplitudeLoader.h"

struct Candidate {
     RhoEvent  ev;
     double weight;
 };
 
 struct PhysicsInput {
     Wkernels::Mat4 u, l, n, s;
     double dsigmaT_dt, dsigmaL_dt;   // [nb / GeV²]
     double Pl = 0.0, SL = 0.0, ST = 0.0;   // polarisations
 };

class Cross {

  TRandom3 rng;
  double beamEnergy = 10.6;
  static constexpr double alem  = 1.0 / 137.035999084;

public:

  Cross(double beam = 10.6){ beamEnergy = beam; rng.SetSeed(0); }

  void generate(Candidate &cand, double Q2_min = 2.0, double Q2_max=2.05,
                double xB_min = 0.3, double xB_max=0.31) {
     double q2 = rng.Uniform(Q2_min,Q2_max);
     double xb = rng.Uniform(xB_min,xB_max);

     cand.ev.Q2  = q2;
     cand.ev.xB  = xb;

     // Set polarizations
     // For now: randomize beam, zero target (change as needed)
     double pl = rng.Uniform(-1.0,1.0);
     cand.ev.beam_pol = (pl > 0) ? 1.0 : -1.0;  // discretized ±1

     cand.ev.target_pol_long = 0.0;   // unpolarized target (set to ±1 for polarized)
     cand.ev.target_pol_trans = 0.0;  // no transverse polarization
     cand.ev.target_pol_azimuth = 0.0; // not used if ST=0

     // To enable target polarization, uncomment and modify:
     // cand.ev.target_pol_long = 1.0;
     // cand.ev.target_pol_trans = 0.5;
     // cand.ev.target_pol_azimuth = rng.Uniform(0, 2*M_PI);

     cand.ev.beamE = beamEnergy;
     cand.ev.BuildScatteredElectron();
     cand.ev.BuildVirtualPhoton();
     cand.ev.BuildRecoilProton();
     cand.ev.BuildRecoilPions();

     // CRITICAL: Compute PDF-defined angles (equations 20-24)
     cand.ev.ComputePDFAngles();
  }

   double dsigma_3fold(double xB, double Q2, double y,
                            double eps, double dsigmaT_dt, double dsigmaL_dt)
 {
     double pref = alem/(2.0*M_PI);
     double kin  = (y*y)/(1.0-eps)*(1.0-xB)/xB/Q2;
     return pref * kin * (dsigmaT_dt + eps*dsigmaL_dt);
 }

  double dsigma_7fold(double WUU,double WLU,double WUL,
                            double WLL,double WUT,double WLT,
                            double Pl,double SL,double ST,
                            double sigma3)
 {
     double S = WUU + Pl*WLU + SL*WUL + Pl*SL*WLL
                      + ST*WUT + Pl*ST*WLT;
     return sigma3 * S / (4.0*M_PI*M_PI);
 }


  void updateWeight(Candidate &cand) {

    double t = (cand.ev.p_out-cand.ev.p_in).M2();
    PhysicsInput ph = loadPhysicsInputs(cand.ev.xB, cand.ev.Q2, t);

    // Get polarizations from event
    double Pl = cand.ev.beam_pol;
    double SL = cand.ev.target_pol_long;
    double ST = cand.ev.target_pol_trans;
    double phiS = cand.ev.target_pol_azimuth;

    // Angles from PDF calculation
    double phi_prod = cand.ev.phi;
    double phi_decay = cand.ev.phiPi;
    double theta = cand.ev.thetaPi;

    double cosTheta = std::cos(theta);
    double sinTheta = std::sin(theta);
    double cos2 = cosTheta * cosTheta;
    double sin2 = sinTheta * sinTheta;
    double sqrt2_cos_sin = std::sqrt(2.0) * cosTheta * sinTheta;

    // Compute all 6 W kernel sets (equations 13, 15-19)
    auto WUU = Wkernels::UU(ph.u, cand.ev.eps, phi_prod, phi_decay);
    auto WLU = Wkernels::LU(ph.u, cand.ev.eps, phi_prod, phi_decay);
    auto WUL = Wkernels::UL(ph.l, cand.ev.eps, phi_prod, phi_decay);
    auto WLL = Wkernels::LL(ph.l, cand.ev.eps, phi_prod, phi_decay);
    auto WUT = Wkernels::UT(ph.n, ph.s, cand.ev.eps, phiS, phi_prod, phi_decay);
    auto WLT = Wkernels::LT(ph.n, ph.s, cand.ev.eps, phiS, phi_prod, phi_decay);

    // Apply theta structure (PDF equations 7-12)
    double W_UU_full = cos2*WUU.LL + sqrt2_cos_sin*WUU.LT + sin2*WUU.TT;
    double W_LU_full = cos2*WLU.LL + sqrt2_cos_sin*WLU.LT + sin2*WLU.TT;
    double W_UL_full = cos2*WUL.LL + sqrt2_cos_sin*WUL.LT + sin2*WUL.TT;
    double W_LL_full = cos2*WLL.LL + sqrt2_cos_sin*WLL.LT + sin2*WLL.TT;
    double W_UT_full = cos2*WUT.LL + sqrt2_cos_sin*WUT.LT + sin2*WUT.TT;
    double W_LT_full = cos2*WLT.LL + sqrt2_cos_sin*WLT.LT + sin2*WLT.TT;

    // 3-fold reduced cross section (PDF equation 5)
    double sigma3 = dsigma_3fold(cand.ev.xB, cand.ev.Q2, cand.ev.y,
                                  cand.ev.eps, ph.dsigmaT_dt, ph.dsigmaL_dt);

    // Full 7-fold cross section with ALL 6 polarization terms (PDF equation 4)
    double w = dsigma_7fold(W_UU_full, W_LU_full, W_UL_full, W_LL_full,
                            W_UT_full, W_LT_full,
                            Pl, SL, ST, sigma3);

    cand.weight = w;

    // Safeguard against negative weights
    if (w < 0) {
      std::cerr << "WARNING: negative weight " << w << " at Q2=" << cand.ev.Q2
                << ", xB=" << cand.ev.xB << ", t=" << t << "\n";
      std::cerr << "  Polarizations: Pl=" << Pl << ", SL=" << SL << ", ST=" << ST << "\n";
      std::cerr << "  W terms: UU=" << W_UU_full << ", LU=" << W_LU_full
                << ", UL=" << W_UL_full << ", LL=" << W_LL_full
                << ", UT=" << W_UT_full << ", LT=" << W_LT_full << "\n";
      cand.weight = 0.0;
    }

    // Check for NaN/Inf
    if (!std::isfinite(w)) {
      std::cerr << "ERROR: non-finite weight at Q2=" << cand.ev.Q2 << ", xB=" << cand.ev.xB << "\n";
      cand.weight = 0.0;
    }
  }
  
  PhysicsInput loadPhysicsInputs(double xB, double Q2,
                                double t /*GeV²*/)
  {
    PhysicsInput ph{};

    // Load amplitudes using AmplitudeLoader
    // Method: "SCHC" for simple test model, "FILE" for external data
    std::string method = "SCHC";  // Change to "FILE" when you have amplitude files
    std::string filename = "";     // Set to amplitude file path if method="FILE"

    try {
        AmplitudeSet amp = AmplitudeLoader::Load(method, xB, Q2, t, filename);

        // Copy amplitude matrices
        ph.u = amp.u;
        ph.l = amp.l;
        ph.n = amp.n;  // Now populated!
        ph.s = amp.s;  // Now populated!

        // Copy cross sections
        ph.dsigmaT_dt = amp.dsigmaT_dt;
        ph.dsigmaL_dt = amp.dsigmaL_dt;

        // Validate (with eps ~ 0.5 as typical value)
        // AmplitudeLoader::Validate(amp, 0.5, false);  // set true for debug output

    } catch (const std::exception& e) {
        std::cerr << "AmplitudeLoader error: " << e.what() << "\n";
        std::cerr << "Using fallback placeholder values\n";

        // Fallback: minimal non-zero values
        ph.dsigmaT_dt = 1.0;
        ph.dsigmaL_dt = 0.3;
        ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')][Wkernels::h('+')] = {0.5, 0.0};
        ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')] = {0.3, 0.0};
    }

    return ph;
 }
};


#endif /* RHO_LUND_IO_H */

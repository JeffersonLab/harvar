// RhoEvent.h
#ifndef RHO_EVENT_H
#define RHO_EVENT_H

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TwoBodyDecay.h"
#include "AngleCalculator.h"
#include <cmath>
#include <iostream>

class RhoEvent {

public:
    //--- inputs
  double Q2, xB, beamE;

  // Polarizations (PDF equation 4)
  double beam_pol;           // P_ℓ ∈ [-1, +1] (longitudinal)
  double target_pol_long;    // S_L ∈ [-1, +1] (longitudinal)
  double target_pol_trans;   // S_T ∈ [0, 1]   (transverse magnitude)
  double target_pol_azimuth; // ϕ_S ∈ [0, 2π)  (transverse direction)

  static constexpr double Mp = 0.93827;
  //--- stored four‑vectors
    TLorentzVector e_in, e_out;
    TLorentzVector p_in, p_out;
    TLorentzVector q,   v;        // virtual γ and ρ
    TLorentzVector piPlus, piMinus;

    double eps , y, enu;

    //--- PDF-defined angles (equations 20-24 from rho-v1.pdf)
    double phi;       // ϕ: production plane angle (eq 20-21)
    double phiPi;     // φ (kappa): decay plane angle (eq 22-23)
    double thetaPi;   // θ: decay polar angle in ρ rest frame (eq 24)

    //--- DEPRECATED: old angles from TwoBodyDecay (DO NOT USE FOR PHYSICS!)
    double theta_production_OLD;   // old production angle (WRONG frame)
    double phi_production_OLD;     // old production azimuth (WRONG frame)

    // ctor: set kinematics, beamE defaults to 10.6 GeV
    RhoEvent(double beamE_ = 10.6)
      : Q2(2.0), xB(0.3), beamE(beamE_),
        beam_pol(0.0), target_pol_long(0.0),
        target_pol_trans(0.0), target_pol_azimuth(0.0)
    {
        rng.SetSeed(0);
        // target proton at rest
        p_in.SetPxPyPzE(0,0,0, Mp);
        // incoming e⁻ along z
        e_in.SetPxPyPzE(0,0, beamE, beamE);
    }
    RhoEvent(double Q2_, double xB_, double beamE_ = 10.6)
      : Q2(Q2_), xB(xB_), beamE(beamE_),
        beam_pol(0.0), target_pol_long(0.0),
        target_pol_trans(0.0), target_pol_azimuth(0.0),
        phi(0), phiPi(0), thetaPi(0),
        theta_production_OLD(0), phi_production_OLD(0)
    {
        rng.SetSeed(0);
        // target proton at rest
        p_in.SetPxPyPzE(0,0,0, Mp);
        // incoming e⁻ along z
        e_in.SetPxPyPzE(0,0, beamE, beamE);
    }

    // build full event + decay + angles
    void Generate() {
        BuildVirtualPhoton();
        BuildScatteredElectron();
        BuildRecoilProton();
        BuildRhoMeson();
        DecayRho();
        ComputePDFAngles();  // CRITICAL: compute angles per PDF spec
    }

    // Compute ϕ, φ, θ angles per PDF equations (20-24)
    void ComputePDFAngles() {
        auto angles = AngleCalculator::ComputeAngles(
            e_in, e_out, p_in, p_out, v, piPlus);

        phi = angles.phi;       // ϕ (production plane)
        phiPi = angles.kappa;   // φ (decay plane)
        thetaPi = angles.theta; // θ (decay polar)
    }

    // print everything
    void Print() const {
        TLorentzVector t = p_out-p_in;
        std::cout << "--- RhoEvent ---\n";
        std::cout << "e_in:    "; e_in.Print();
        std::cout << "e_out:   "; e_out.Print();
        std::cout << "p_in:    "; p_in.Print();
        std::cout << "p_out:   "; p_out.Print();
        std::cout << "q:       "; q.Print();
        std::cout << "rho(v):  "; v.Print();
        std::cout << "pi+:     "; piPlus.Print();
        std::cout << "pi-:     "; piMinus.Print();
        std::cout << " Q2  = " << Q2 << std::endl;
        std::cout << " xB  = " << xB << std::endl;
        std::cout << " y   = " << y << std::endl;
        std::cout << " eps = " << eps << std::endl;
        std::cout << " -t  = " << t.M2() << std::endl;
        std::cout << "\n--- PDF Angles (eqs 20-24) ---\n";
        std::cout << "ϕ (phi):         " << phi*180/M_PI << "° (production plane)\n";
        std::cout << "φ (phiPi):       " << phiPi*180/M_PI << "° (decay plane)\n";
        std::cout << "θ (thetaPi):     " << thetaPi*180/M_PI << "° (decay polar)\n";
    }

public:
    TRandom3 rng;

    // 1) virtual photon: q = e_in - e_out
    void BuildVirtualPhoton() {
        q = e_in-e_out;
    }

    TLorentzVector ScatteredElectron(double Q2, double xB, double beamE) {
    // 1) energy transfer ν = Q2 / (2 Mp xB)
        double nu = Q2 / (2.0 * Mp * xB);
        enu = nu;
    // 2) scattered‐electron energy E' = E_beam − ν
     double Eprime = beamE - nu;
       y     = nu/beamE;
     double γ     = gamma_bj(xB,Q2);
      eps   = epsilon(y,γ);
    // 3) scattering angle θ_e from Q2 = 2 E E' (1 - cos θ)
        double cosTh = 1.0 - Q2 / (2.0 * beamE * Eprime);
    // guard against small numerical slip:
        if (cosTh > +1) cosTh = +1;
        if (cosTh < -1) cosTh = -1;
        double sinTh = std::sqrt(1.0 - cosTh*cosTh);

        // 4) build four‐vector (px, py, pz, E)
        //    assume scattering in the x–z plane: py = 0
        double px = Eprime * sinTh;
        double pz = Eprime * cosTh;

        return TLorentzVector(px, 0.0, pz, Eprime);
    }

    // 2) scattered electron (simple elastic kinematics)
    void BuildScatteredElectron() {
        e_out = ScatteredElectron(Q2,xB,beamE);
    }

    // 3) recoil proton from t
    void BuildRecoilProton() {
        TLorentzVector  cm = q + p_in;
        TVector3       ref = e_out.Vect();
        TwoBodyDecay decay (cm , Mp, 0.77526);
        decay.Generate(ref);
        p_out = decay.daughter1;
        v = decay.daughter2;

        // Store old (WRONG) angles for debugging/comparison
        theta_production_OLD = decay.theta;
        phi_production_OLD = decay.phi;
    }
    void BuildRecoilPions() {
        // Generate isotropic decay (for now - will be replaced with
        // amplitude-weighted decay in Phase 5)
        TLorentzVector  vrho = v;
        TVector3         ref = q.Vect();
        TwoBodyDecay decay (vrho , 0.139570, 0.139570);
        decay.Generate(ref);
        piPlus = decay.daughter1;
        piMinus = decay.daughter2;

        // PDF angles will be computed by ComputePDFAngles()
        // DO NOT use decay.theta, decay.phi - they are in wrong frame!
    }
    // 4) ρ meson from 4‑momentum conservation
    void BuildRhoMeson() {
        //v = q + p_in - e_out - p_out;
        v = q + p_in - p_out;
    }

    // 5) isotropic 2‑body decay in ρ rest frame
    //
    // NOTE: This generates ISOTROPIC decay for WEIGHTED event generation.
    // The PDF angles (ϕ,φ,θ) are computed from the generated kinematics
    // via ComputePDFAngles(), and events are weighted by W(ϕ,φ,θ) in Cross::updateWeight().
    //
    // This is physically correct for weighted events.
    //
    // For UNWEIGHTED events, you would need to:
    //   1. Sample (ϕ,φ,θ) from W distribution via rejection sampling
    //   2. Generate decay kinematics consistent with those specific angles
    //   3. Set all event weights to 1
    // This is more complex but eliminates event weights in analysis.
    //
    void DecayRho() {
        const double Mpi = 0.13957;
        double Mrho = v.M();
        double pMag = sqrt((Mrho*Mrho - 4*Mpi*Mpi))/2.;

        // Isotropic decay angles
        double cos_t = rng.Uniform(-1,1);
        double sin_t = sqrt(1 - cos_t*cos_t);
        double ph    = rng.Uniform(0, 2*M_PI);

        // momenta in ρ rest frame
        TLorentzVector p1(pMag*sin_t*cos(ph),
                          pMag*sin_t*sin(ph),
                          pMag*cos_t,
                          sqrt(pMag*pMag + Mpi*Mpi));
        TLorentzVector p2 = p1;
        p2.SetVect(-p1.Vect());  // back‑to‑back

        // boost back to lab
        TVector3 b = v.BoostVector();
        piPlus  = p1;  piPlus.Boost(b);
        piMinus = p2;  piMinus.Boost(b);
    }

    // 6) compute θ,φ in ρ rest frame with helicity axes

     double gamma_bj(double xB, double Q2)
    { return 2.0 * xB * Mp / std::sqrt(Q2); }

  double epsilon(double y, double γ)
 {
     double y2γ2 = y*y*γ*γ;
     return (1.0 - y - 0.25*y2γ2) /
            (1.0 - y + 0.5*y*y + 0.25*y2γ2);
 }
};

#endif // RHO_EVENT_H

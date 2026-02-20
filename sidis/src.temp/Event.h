// Event.h  --------------------------------------------------------------
#pragma once
#include <TLorentzVector.h>
#include <cmath>

class Event {
public:
    // mandatory four‑vectors (laboratory frame)
    TLorentzVector l;      // incoming lepton
    TLorentzVector lPrime; // scattered lepton
    TLorentzVector p;      // initial nucleon  (usually at rest)
    TLorentzVector pPrime; // recoil nucleon
    TLorentzVector v;      // ρ‑meson
    TLorentzVector k1;     // π+
    TLorentzVector k2;     // π−   (optional, not used here)

    // outputs -----------------------------------------------------------
    double cos_phiH = 0.0; //!<  cos ϕ   (production‑plane angle, Eq. 20)
    double sin_phiH = 0.0; //!<  sin ϕ   (Eq. 21)
    double   phiH   = 0.0; //!<  ϕ  in [0,2π)

    double cos_phi  = 0.0; //!<  cos φ   (decay angle, Eq. 22)
    double sin_phi  = 0.0; //!<  sin φ   (Eq. 23)
    double   phi    = 0.0; //!<  φ  in [0,2π)

    double cos_theta = 0.0; //!< cos θ    (Eq. 24)
    double   theta   = 0.0; //!< θ  in [0,π]

    //-------------------------------------------------------------------
    /** call after all four‑vectors are set */
    void computeAngles();

private:
    TVector3 boostToCM()  const;  ///< γ*p hadronic‑CM boost
    TVector3 boostToVM()  const;  ///< ρ rest‑frame boost
    static  double angleFromSinCos(double s, double c);
};

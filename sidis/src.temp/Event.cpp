// Event.cxx  ------------------------------------------------------------
#include "Event.h"

// Boost that brings γ*+p system to rest
TVector3 Event::boostToCM() const {
    TLorentzVector q = l - lPrime;
    TLorentzVector had = q + p;      // γ* + p
    return  –had.BoostVector();      // opposite to hadron motion
}

// Boost that brings ρ to rest
TVector3 Event::boostToVM() const { return –v.BoostVector(); }

// unified atan2 wrapper returning an angle in [0,2π)
double Event::angleFromSinCos(double s, double c) {
    double a = std::atan2(s, c);      // (–π,π]
    return (a < 0) ? a + 2*M_PI : a;  // → [0,2π)
}

void Event::computeAngles()
{
    // --- 1) go to hadronic CM -----------------------------------------
    TVector3 bCM = boostToCM();

    TLorentzVector l_CM  = l;      l_CM.Boost(bCM);
    TLorentzVector lp_CM = lPrime; lp_CM.Boost(bCM);
    TLorentzVector q_CM  = l_CM - lp_CM;
    TLorentzVector v_CM  = v;      v_CM.Boost(bCM);

    TVector3 q  = q_CM.Vect();
    TVector3 v3 = v_CM.Vect();
    TVector3 l3 = l_CM.Vect();
    TVector3 lp3= lp_CM.Vect();

    // Eq. 20–21  ── production‑plane angle ϕ (phiH)
    TVector3 n_prod = q.Cross(v3);
    TVector3 n_scat = l3.Cross(lp3);

    cos_phiH = n_prod.Dot(n_scat) /
               (n_prod.Mag() * n_scat.Mag());

    sin_phiH = n_prod.Cross(n_scat).Dot(q) /
               (n_prod.Mag() * n_scat.Mag() * q.Mag());

    phiH = angleFromSinCos(sin_phiH, cos_phiH);

    // Eq. 22–23  ── decay‑plane angle φ (phi)
    TVector3 n_decay = v3.Cross(k1.Vect().Boost(bCM)); // π+ in CM
    cos_phi = n_prod.Dot(n_decay) /
              (n_prod.Mag() * n_decay.Mag());

    //   [(q×v)×v] =  (n_prod)×v
    TVector3 cross_qv_v = n_prod.Cross(v3);
    sin_phi = cross_qv_v.Dot(k1.Vect().Boost(bCM)) /
              (cross_qv_v.Mag() * n_decay.Mag());

    phi = angleFromSinCos(sin_phi, cos_phi);

    // --- 2) go to ρ rest frame for θ ----------------------------------
    TLorentzVector pPrime_VM = pPrime; pPrime_VM.Boost(boostToVM());
    TLorentzVector k1_VM     = k1;     k1_VM.Boost(boostToVM());

    TVector3 pR = pPrime_VM.Vect().Unit();
    TVector3 kR = k1_VM.Vect().Unit();

    cos_theta = –pR.Dot(kR);     // Eq. 24  (minus sign!)
    theta     = std::acos(std::clamp(cos_theta,‑1.0,1.0));
}

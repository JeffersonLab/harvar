/**********************************************************************
 * AngleCalculator.h - Compute ϕ, φ, θ angles per rho-v1.pdf spec
 *
 * Implements equations (20-24) from the PDF specification:
 *   - ϕ: azimuthal angle between ρ production plane and lepton
 *        scattering plane in hadronic CM system
 *   - φ: azimuthal angle between ρ decay plane and production plane
 *   - θ: polar angle of π+ in vector-meson rest frame
 *        (z-axis opposite to p', y-axis along p'×q)
 *
 * CRITICAL: These are NOT the same as TwoBodyDecay angles!
 *********************************************************************/
#ifndef ANGLE_CALCULATOR_H
#define ANGLE_CALCULATOR_H

#include "TLorentzVector.h"
#include "TVector3.h"
#include <cmath>
#include <stdexcept>

class AngleCalculator {
public:
    struct Angles {
        double phi;      // ϕ: production plane angle [0, 2π)
        double kappa;    // φ: decay plane angle [0, 2π)
        double theta;    // θ: decay polar angle [0, π]
    };

    /**
     * Compute all PDF-defined angles from event kinematics.
     *
     * All four-vectors should be in the LAB frame.
     * The function handles necessary boosts internally.
     *
     * @param l       Incoming electron four-vector
     * @param l_prime Scattered electron four-vector
     * @param p       Initial proton four-vector (at rest)
     * @param p_prime Final proton four-vector
     * @param v       Rho meson four-vector
     * @param pi_plus Pi+ four-vector (lab frame)
     * @return Angles struct with ϕ, φ, θ
     */
    static Angles ComputeAngles(
        const TLorentzVector& l,
        const TLorentzVector& l_prime,
        const TLorentzVector& p,
        const TLorentzVector& p_prime,
        const TLorentzVector& v,
        const TLorentzVector& pi_plus);

private:
    static constexpr double EPSILON = 1e-10;

    /**
     * Compute ϕ using equations (20-21).
     * Angle between ρ production plane and lepton scattering plane
     * in hadronic CM system.
     */
    static double ComputePhi(const TVector3& l_cm,
                             const TVector3& l_prime_cm,
                             const TVector3& q_cm,
                             const TVector3& v_cm);

    /**
     * Compute φ using equations (22-23).
     * Angle between ρ decay plane and production plane
     * in hadronic CM system.
     */
    static double ComputeKappa(const TVector3& q_cm,
                               const TVector3& v_cm,
                               const TVector3& k1_cm);

    /**
     * Compute θ using equation (24).
     * Polar angle of π+ in ρ rest frame with specific axis definition:
     *   z-axis: opposite to outgoing nucleon p'
     *   y-axis: p' × q
     *
     * NOTE THE MINUS SIGN in the PDF equation (24)!
     */
    static double ComputeTheta(const TVector3& p_prime_rho,
                               const TVector3& q_rho,
                               const TVector3& pi_plus_rho);

    /**
     * Safe normalization: returns unit vector or throws if |v| < epsilon.
     */
    static TVector3 SafeUnit(const TVector3& v, const char* context);

    /**
     * Compute angle from sin and cos using atan2 (always in [0, 2π)).
     */
    static double AngleFromSinCos(double sin_val, double cos_val);

    /**
     * Clamp value to [-1, 1] for safe acos.
     */
    static double ClampForAcos(double x);
};

//=============================================================================
// Implementation
//=============================================================================

inline AngleCalculator::Angles AngleCalculator::ComputeAngles(
    const TLorentzVector& l,
    const TLorentzVector& l_prime,
    const TLorentzVector& p,
    const TLorentzVector& p_prime,
    const TLorentzVector& v,
    const TLorentzVector& pi_plus)
{
    Angles result;

    // Virtual photon
    TLorentzVector q = l - l_prime;

    // -------------------------------------------------------------------------
    // 1. Boost to hadronic CM frame: boost by -(q + p)
    // -------------------------------------------------------------------------
    TLorentzVector cm_system = q + p;
    TVector3 boost_to_cm = -cm_system.BoostVector();

    TLorentzVector l_cm = l;           l_cm.Boost(boost_to_cm);
    TLorentzVector l_prime_cm = l_prime; l_prime_cm.Boost(boost_to_cm);
    TLorentzVector q_cm = q;           q_cm.Boost(boost_to_cm);
    TLorentzVector v_cm = v;           v_cm.Boost(boost_to_cm);
    TLorentzVector pi_plus_cm = pi_plus; pi_plus_cm.Boost(boost_to_cm);

    // -------------------------------------------------------------------------
    // 2. Compute ϕ in hadronic CM (eqs 20-21)
    // -------------------------------------------------------------------------
    result.phi = ComputePhi(l_cm.Vect(), l_prime_cm.Vect(),
                             q_cm.Vect(), v_cm.Vect());

    // -------------------------------------------------------------------------
    // 3. Compute φ in hadronic CM (eqs 22-23)
    // -------------------------------------------------------------------------
    result.kappa = ComputeKappa(q_cm.Vect(), v_cm.Vect(),
                                 pi_plus_cm.Vect());

    // -------------------------------------------------------------------------
    // 4. Compute θ in ρ rest frame (eq 24)
    // -------------------------------------------------------------------------
    TLorentzVector p_prime_rho = p_prime;
    TLorentzVector q_rho = q;
    TLorentzVector pi_plus_rho = pi_plus;

    TVector3 boost_to_rho = -v.BoostVector();
    p_prime_rho.Boost(boost_to_rho);
    q_rho.Boost(boost_to_rho);
    pi_plus_rho.Boost(boost_to_rho);

    result.theta = ComputeTheta(p_prime_rho.Vect(), q_rho.Vect(),
                                 pi_plus_rho.Vect());

    return result;
}

inline double AngleCalculator::ComputePhi(
    const TVector3& l_cm,
    const TVector3& l_prime_cm,
    const TVector3& q_cm,
    const TVector3& v_cm)
{
    // Equation (20-21):
    // cos ϕ = [(q⃗ × v⃗) · (l⃗ × l⃗')] / [|q⃗ × v⃗| · |l⃗ × l⃗'|]
    // sin ϕ = {[(q⃗ × v⃗) × (l⃗ × l⃗')] · q⃗} / [|q⃗ × v⃗| · |l⃗ × l⃗'| · |q⃗|]

    TVector3 n_prod = q_cm.Cross(v_cm);        // production plane normal
    TVector3 n_lept = l_cm.Cross(l_prime_cm);  // lepton plane normal

    double n_prod_mag = n_prod.Mag();
    double n_lept_mag = n_lept.Mag();

    if (n_prod_mag < EPSILON || n_lept_mag < EPSILON) {
        // Degenerate geometry: planes are undefined
        // This can happen if q || v or l || l' (unphysical)
        return 0.0;  // default value
    }

    double cos_phi = n_prod.Dot(n_lept) / (n_prod_mag * n_lept_mag);

    // For sin ϕ: [(q⃗ × v⃗) × (l⃗ × l⃗')] · q⃗
    TVector3 temp = n_prod.Cross(n_lept);
    double q_mag = q_cm.Mag();

    if (q_mag < EPSILON) {
        return 0.0;  // degenerate
    }

    double sin_phi = temp.Dot(q_cm) / (n_prod_mag * n_lept_mag * q_mag);

    return AngleFromSinCos(sin_phi, cos_phi);
}

inline double AngleCalculator::ComputeKappa(
    const TVector3& q_cm,
    const TVector3& v_cm,
    const TVector3& k1_cm)
{
    // Equation (22-23):
    // cos φ = [(q⃗ × v⃗) · (v⃗ × k⃗₁)] / [|q⃗ × v⃗| · |v⃗ × k⃗₁|]
    // sin φ = {[(q⃗ × v⃗) × v⃗] · (k⃗₁ × v⃗)} / [|(q⃗ × v⃗) × v⃗| · |k⃗₁ × v⃗|]

    TVector3 n_prod = q_cm.Cross(v_cm);       // production plane normal
    TVector3 n_decay = v_cm.Cross(k1_cm);     // decay plane normal

    double n_prod_mag = n_prod.Mag();
    double n_decay_mag = n_decay.Mag();

    if (n_prod_mag < EPSILON || n_decay_mag < EPSILON) {
        return 0.0;  // degenerate
    }

    double cos_kappa = n_prod.Dot(n_decay) / (n_prod_mag * n_decay_mag);

    // For sin φ: [(q⃗ × v⃗) × v⃗] · (k⃗₁ × v⃗)
    TVector3 temp1 = n_prod.Cross(v_cm);
    TVector3 temp2 = k1_cm.Cross(v_cm);

    double temp1_mag = temp1.Mag();
    double temp2_mag = temp2.Mag();

    if (temp1_mag < EPSILON || temp2_mag < EPSILON) {
        return 0.0;
    }

    double sin_kappa = temp1.Dot(temp2) / (temp1_mag * temp2_mag);

    return AngleFromSinCos(sin_kappa, cos_kappa);
}

inline double AngleCalculator::ComputeTheta(
    const TVector3& p_prime_rho,
    const TVector3& q_rho,
    const TVector3& pi_plus_rho)
{
    // Equation (24):
    // cos θ = -[p⃗' · P⃗_π+] / [|p⃗'| · |P⃗_π+|]
    //
    // The z-axis in the ρ rest frame is defined as OPPOSITE to p⃗'
    // (hence the minus sign in the PDF equation)

    double p_prime_mag = p_prime_rho.Mag();
    double pi_plus_mag = pi_plus_rho.Mag();

    if (p_prime_mag < EPSILON || pi_plus_mag < EPSILON) {
        return 0.0;  // degenerate (shouldn't happen)
    }

    // NOTE: Minus sign as per PDF equation (24)
    double cos_theta = -(p_prime_rho.Dot(pi_plus_rho)) /
                       (p_prime_mag * pi_plus_mag);

    // Clamp to avoid numerical errors in acos
    cos_theta = ClampForAcos(cos_theta);

    return std::acos(cos_theta);
}

inline TVector3 AngleCalculator::SafeUnit(const TVector3& v,
                                           const char* context)
{
    double mag = v.Mag();
    if (mag < EPSILON) {
        throw std::runtime_error(
            std::string("AngleCalculator: degenerate vector in ") + context);
    }
    return v.Unit();
}

inline double AngleCalculator::AngleFromSinCos(double sin_val,
                                                double cos_val)
{
    // atan2 returns angle in [-π, π], convert to [0, 2π)
    double angle = std::atan2(sin_val, cos_val);
    if (angle < 0.0) angle += 2.0 * M_PI;
    return angle;
}

inline double AngleCalculator::ClampForAcos(double x)
{
    if (x > 1.0) return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}

#endif // ANGLE_CALCULATOR_H

# Implementation Plan: Rho-Meson Event Generator

## Overview

This plan implements the complete rho-meson electroproduction event generator according to `rho-v1.pdf`. The work is divided into phases with clear deliverables and verification steps.

---

## Phase 1: Critical Angle Calculation Fix (MUST DO FIRST)

### 1.1 Create AngleCalculator Class

**File**: `sidis/src/AngleCalculator.h`

**Purpose**: Implement PDF equations (20-24) for ϕ, φ, θ

**Interface**:
```cpp
class AngleCalculator {
public:
    struct Angles {
        double phi;      // ϕ: production plane angle
        double kappa;    // φ: decay plane angle
        double theta;    // θ: decay polar angle
    };

    // Compute all angles from event kinematics
    static Angles ComputeAngles(
        const TLorentzVector& l,      // incoming electron
        const TLorentzVector& l_prime, // scattered electron
        const TLorentzVector& p,      // initial proton
        const TLorentzVector& p_prime, // final proton
        const TLorentzVector& v,      // rho meson
        const TLorentzVector& pi_plus  // π+ (in rho rest frame)
    );

private:
    // Compute ϕ: between production and scattering planes (eqs 20-21)
    static double ComputePhi(const TVector3& l, const TVector3& l_prime,
                             const TVector3& q, const TVector3& v);

    // Compute φ: between production and decay planes (eqs 22-23)
    static double ComputeKappa(const TVector3& q, const TVector3& v,
                               const TVector3& k1_cm);  // k1 in hadronic CM

    // Compute θ: polar angle in rho rest frame (eq 24)
    static double ComputeTheta(const TVector3& p_prime_rho,
                               const TVector3& pi_plus_rho);
};
```

**Implementation Details**:

1. **ϕ Calculation (eqs 20-21)**:
   ```cpp
   TVector3 n_prod = q.Cross(v);           // production plane normal
   TVector3 n_lept = l.Cross(l_prime);     // lepton plane normal

   double cos_phi = n_prod.Dot(n_lept) / (n_prod.Mag() * n_lept.Mag());

   TVector3 temp = n_prod.Cross(n_lept);
   double sin_phi = temp.Dot(q) / (n_prod.Mag() * n_lept.Mag() * q.Mag());

   double phi = std::atan2(sin_phi, cos_phi);
   ```

2. **φ Calculation (eqs 22-23)**:
   ```cpp
   TVector3 n_prod = q.Cross(v);
   TVector3 n_decay = v.Cross(k1_cm);

   double cos_kappa = n_prod.Dot(n_decay) / (n_prod.Mag() * n_decay.Mag());

   TVector3 temp1 = n_prod.Cross(v);
   TVector3 temp2 = k1_cm.Cross(v);
   double sin_kappa = temp1.Dot(temp2) / (temp1.Mag() * temp2.Mag());

   double kappa = std::atan2(sin_kappa, cos_kappa);
   ```

3. **θ Calculation (eq 24)**:
   ```cpp
   // In rho rest frame:
   // z-axis: opposite to p_prime
   // y-axis: p_prime × q
   // Note the minus sign in the PDF equation!

   double cos_theta = -(p_prime_rho.Dot(pi_plus_rho)) /
                      (p_prime_rho.Mag() * pi_plus_rho.Mag());

   // Clamp to avoid acos domain errors
   if (cos_theta > 1.0) cos_theta = 1.0;
   if (cos_theta < -1.0) cos_theta = -1.0;

   double theta = std::acos(cos_theta);
   ```

**Numerical Safety**:
- Add checks: `if (mag < 1e-10) handle_degenerate_case()`
- Clamp acos arguments to [-1, 1]
- Document failure modes

**Unit Tests** (`tests/test_angles.cpp`):
```cpp
// Test 1: Known geometry (coplanar scattering, specific decay)
// Test 2: sin²+cos² = 1 for all angles
// Test 3: Angle ranges: ϕ,φ ∈ [0,2π), θ ∈ [0,π]
// Test 4: Boost invariance: angles same in any frame
```

---

### 1.2 Integrate into RhoEvent

**File**: `sidis/src/RhoEvent.h`

**Changes**:

1. Add hadronic CM frame support:
```cpp
TLorentzVector q_cm, v_cm, p_prime_cm;  // vectors in hadronic CM
TVector3 boost_to_cm;  // = (q + p_in).BoostVector()
```

2. Replace `DecayRho()` with two-stage approach:
```cpp
void DecayRho() {
    // Generate decay kinematics (temporarily isotropic)
    GenerateDecayKinematics();

    // Compute PDF-defined angles
    ComputePDFAngles();
}

void GenerateDecayKinematics() {
    // Same as current DecayRho but don't compute angles yet
}

void ComputePDFAngles() {
    // Use AngleCalculator
    auto angles = AngleCalculator::ComputeAngles(
        e_in, e_out, p_in, p_out, v, piPlus);
    phi = angles.phi;
    phiPi = angles.kappa;
    thetaPi = angles.theta;
}
```

3. Remove TwoBodyDecay angle calculations (keep only kinematics generation)

**Testing**:
- Generate 1000 events
- Verify all angles in correct ranges
- No NaNs/Infs
- Compare to old angles: should be DIFFERENT (old ones were wrong!)

---

## Phase 2: Complete Helicity Amplitude System

### 2.1 Create AmplitudeLoader

**File**: `sidis/src/AmplitudeLoader.h`

**Purpose**: Load/compute helicity amplitude matrices u, l, n, s

**Interface**:
```cpp
struct AmplitudeSet {
    Wkernels::Mat4 u;  // unpolarized
    Wkernels::Mat4 l;  // longitudinal beam
    Wkernels::Mat4 n;  // transverse target (sin term)
    Wkernels::Mat4 s;  // transverse target (cos term)
    double dsigmaT_dt;  // nb/GeV²
    double dsigmaL_dt;  // nb/GeV²
};

class AmplitudeLoader {
public:
    // Load from file for given kinematics
    static AmplitudeSet LoadFromFile(const std::string& filename,
                                      double xB, double Q2, double t);

    // Simple parameterization (for testing/defaults)
    static AmplitudeSet SimpleModel(double xB, double Q2, double t);

    // Validate amplitude set (hermiticity, normalization)
    static bool Validate(const AmplitudeSet& amp, double eps);
};
```

**File Format** (CSV example):
```
# xB, Q2, t, nu, nup, mu, mup, Re(u), Im(u), Re(l), Im(l), Re(n), Im(n), Re(s), Im(s)
0.3, 2.0, -0.5, +, +, +, +, 0.5, 0.0, 0.0, 0.1, ...
0.3, 2.0, -0.5, +, +, +, 0, 0.2, 0.05, ...
...
```

**Alternative**: Use JSON or HDF5 for easier structure

**SimpleModel** (for testing):
```cpp
// Implement s-channel helicity conservation (SCHC) as baseline
// or use Regge-inspired parameterization
// Must satisfy normalization constraint (PDF eq 14)
```

**Validation**:
```cpp
bool Validate(const AmplitudeSet& amp, double eps) {
    // Check equation (14):
    // u[++][++] + u[--][++] + 2ε u[++][00] = 1 - (u[00][++] + ε u[00][00])

    // Check positivity of SDME (if applicable)

    // Check ranges of dsigma/dt > 0

    return all_checks_passed;
}
```

---

### 2.2 Update Cross.h

**Changes to `loadPhysicsInputs`**:

```cpp
PhysicsInput loadPhysicsInputs(double xB, double Q2, double t) {
    PhysicsInput ph{};

    // Load from configuration
    std::string amp_file = config.amplitudeFile;

    if (config.useSimpleModel) {
        auto amp = AmplitudeLoader::SimpleModel(xB, Q2, t);
        ph.u = amp.u;
        ph.l = amp.l;
        ph.n = amp.n;
        ph.s = amp.s;
        ph.dsigmaT_dt = amp.dsigmaT_dt;
        ph.dsigmaL_dt = amp.dsigmaL_dt;
    } else {
        auto amp = AmplitudeLoader::LoadFromFile(amp_file, xB, Q2, t);
        if (!AmplitudeLoader::Validate(amp, eps)) {
            throw std::runtime_error("Invalid amplitude set");
        }
        ph.u = amp.u;
        // ... etc
    }

    // Set polarizations from configuration
    ph.Pl = config.beamPolarization;
    ph.SL = config.targetPolarizationLong;
    ph.ST = config.targetPolarizationTrans;

    return ph;
}
```

---

## Phase 3: Add Missing W Kernels (UT, LT)

**File**: `sidis/src/w_kernels.hpp`

### 3.1 Implement UT Kernel (PDF Equation 18)

```cpp
namespace Wkernels {
    // ================================================================
    // Unpolarized beam × Transverse target (Eq. 18)
    // ================================================================
    inline Triple UT(const Mat4& n, const Mat4& s,
                     double eps, double phiS, double phi, double kap)
    {
        double ce_p = std::sqrt(eps*(1.0+eps));
        double sphi = std::sin(phi);
        double cphi = std::cos(phi);
        double s_phi_phiS = std::sin(phi - phiS);
        double c_phi_phiS = std::cos(phi - phiS);
        // ... (implement all terms from eq 18)

        Triple out;

        // LL component
        out.LL = s_phi_phiS * ( /* n terms */ )
               + c_phi_phiS * ( /* s terms */ );

        // LT component
        out.LT = s_phi_phiS * ( /* n terms */ )
               + c_phi_phiS * ( /* s terms */ );

        // TT component
        out.TT = s_phi_phiS * ( /* n terms */ )
               + c_phi_phiS * ( /* s terms */ );

        return out;
    }
}
```

**Exact implementation**: Transcribe equation (18) term-by-term from PDF. Each term is labeled in the PDF, making this straightforward but tedious.

### 3.2 Implement LT Kernel (PDF Equation 19)

Similar structure to UT. Transcribe equation (19).

```cpp
inline Triple LT(const Mat4& n, const Mat4& s,
                 double eps, double phiS, double phi, double kap)
{
    // Implementation of eq (19)
    // ...
}
```

**Testing**:
- Create test amplitude sets where n,s are known simple values
- Hand-calculate expected W_UT, W_LT
- Compare against function output

---

## Phase 4: Extend Polarization Support

### 4.1 Add Polarization to RhoEvent

**File**: `sidis/src/RhoEvent.h`

```cpp
class RhoEvent {
public:
    // Polarizations
    double beam_pol;           // P_ℓ ∈ [-1, +1]
    double target_pol_long;    // S_L ∈ [-1, +1]
    double target_pol_trans;   // S_T ∈ [0, 1]
    double target_pol_azimuth; // ϕ_S ∈ [0, 2π)

    // ... rest of class
};
```

### 4.2 Update Cross::generate

```cpp
void generate(Candidate &cand, /* ... */) {
    // ... existing code ...

    // Sample polarizations from configuration
    if (config.polarization.randomizeBeam) {
        cand.ev.beam_pol = rng.Uniform(-1.0, 1.0);
    } else {
        cand.ev.beam_pol = config.polarization.beamPol;
    }

    cand.ev.target_pol_long = config.polarization.targetPolLong;
    cand.ev.target_pol_trans = config.polarization.targetPolTrans;
    cand.ev.target_pol_azimuth = rng.Uniform(0, 2*M_PI);  // or from config

    // ... rest
}
```

### 4.3 Update Cross::updateWeight

```cpp
void updateWeight(Candidate &cand) {
    // ... existing code to compute sigma3 ...

    // Compute all 6 W terms
    auto WUU = Wkernels::UU(ph.u, eps, phi, kappa);
    auto WLU = Wkernels::LU(ph.u, eps, phi, kappa);
    auto WUL = Wkernels::UL(ph.l, eps, phi, kappa);
    auto WLL = Wkernels::LL(ph.l, eps, phi, kappa);
    auto WUT = Wkernels::UT(ph.n, ph.s, eps, phiS, phi, kappa);
    auto WLT = Wkernels::LT(ph.n, ph.s, eps, phiS, phi, kappa);

    // Apply theta structure (eq 7-12)
    double c2 = cosTheta*cosTheta;
    double s2 = sinTheta*sinTheta;
    double cs = cosTheta*sinTheta*std::sqrt(2.0);

    double W_UU_full = c2*WUU.LL + cs*WUU.LT + s2*WUU.TT;
    double W_LU_full = c2*WLU.LL + cs*WLU.LT + s2*WLU.TT;
    double W_UL_full = c2*WUL.LL + cs*WUL.LT + s2*WUL.TT;
    double W_LL_full = c2*WLL.LL + cs*WLL.LT + s2*WLL.TT;
    double W_UT_full = c2*WUT.LL + cs*WUT.LT + s2*WUT.TT;
    double W_LT_full = c2*WLT.LL + cs*WLT.LT + s2*WLT.TT;

    // Full 6-term cross section (eq 4)
    double w = dsigma_7fold(
        W_UU_full, W_LU_full, W_UL_full, W_LL_full, W_UT_full, W_LT_full,
        ph.Pl, ph.SL, ph.ST, sigma3);

    cand.weight = w;
}
```

---

## Phase 5: Fix Decay Sampling

### Option A: Unweighted Events (Preferred)

Generate decay angles from the full angular distribution W(ϕ,φ,θ).

**Algorithm**:
1. For given Q², xB, t: compute max of W over all angles
2. Rejection sampling in (ϕ, φ, θ) space
3. Once accepted, generate decay with those specific angles
4. All events have weight = 1

**Implementation** (`RhoEvent.h`):
```cpp
void DecayRho(const Wkernels::Mat4& u, const Wkernels::Mat4& l,
              const Wkernels::Mat4& n, const Wkernels::Mat4& s,
              double eps, double Pl, double SL, double ST) {

    // Find maximum of W over angles (pre-computed or scanned)
    double W_max = FindMaxW(u, l, n, s, eps, Pl, SL, ST);

    // Rejection sampling
    while (true) {
        double phi_try = rng.Uniform(0, 2*M_PI);
        double kappa_try = rng.Uniform(0, 2*M_PI);
        double cos_theta_try = rng.Uniform(-1, 1);
        double theta_try = std::acos(cos_theta_try);

        double W_val = EvaluateW(u, l, n, s, eps, Pl, SL, ST,
                                 phi_try, kappa_try, theta_try);

        if (rng.Uniform(0, W_max) <= W_val) {
            // Accept: generate decay with these angles
            GenerateDecayWithAngles(theta_try, kappa_try);
            phi = phi_try;  // production angle
            phiPi = kappa_try;
            thetaPi = theta_try;
            break;
        }
    }
}

void GenerateDecayWithAngles(double theta, double kappa) {
    // In rho rest frame with PDF-defined axes:
    // z-axis opposite to p_prime, y-axis = p_prime × q

    // 1. Compute basis in rho rest frame
    TLorentzVector p_prime_rho = p_out;
    p_prime_rho.Boost(-v.BoostVector());

    TLorentzVector q_rho = q;
    q_rho.Boost(-v.BoostVector());

    TVector3 z_axis = -p_prime_rho.Vect().Unit();
    TVector3 y_axis = p_prime_rho.Vect().Cross(q_rho.Vect()).Unit();
    TVector3 x_axis = y_axis.Cross(z_axis);

    // 2. Pi+ momentum in this basis
    double Mpi = 0.13957;
    double Mrho = v.M();
    double p_mag = std::sqrt((Mrho*Mrho - 4*Mpi*Mpi)) / 2.0;

    TVector3 p_pi_rho = p_mag * (std::sin(theta)*std::cos(kappa)*x_axis
                               + std::sin(theta)*std::sin(kappa)*y_axis
                               + std::cos(theta)*z_axis);

    TLorentzVector pi_plus_rho(p_pi_rho, std::hypot(p_mag, Mpi));
    TLorentzVector pi_minus_rho(-p_pi_rho, std::hypot(p_mag, Mpi));

    // 3. Boost to lab
    piPlus = pi_plus_rho;  piPlus.Boost(v.BoostVector());
    piMinus = pi_minus_rho; piMinus.Boost(v.BoostVector());
}
```

### Option B: Weighted Events (Alternative)

- Keep isotropic decay
- Weight each event by W(ϕ,φ,θ) / W_max
- Simpler but requires carrying weights through analysis

**Recommendation**: Use Option A for cleaner physics analysis.

---

## Phase 6: Configuration System

### 6.1 Create Config Class

**File**: `sidis/src/Config.h`

```cpp
struct Config {
    // Beam
    double beamEnergy = 10.6;  // GeV

    // Kinematic ranges
    double Q2_min = 1.0, Q2_max = 10.0;
    double xB_min = 0.1, xB_max = 0.9;
    double t_min = -2.0, t_max = 0.0;

    // Polarization
    struct {
        bool randomizeBeam = false;
        double beamPol = 0.0;         // P_ℓ
        double targetPolLong = 0.0;   // S_L
        double targetPolTrans = 0.0;  // S_T
        bool randomizeTargetAzimuth = true;
    } polarization;

    // Amplitudes
    bool useSimpleModel = true;
    std::string amplitudeFile = "";

    // Output
    std::string outputFormat = "LUND";  // or "ROOT", "HepMC"
    bool writeWeights = false;

    // Load from YAML file
    static Config LoadFromFile(const std::string& filename);

    // Validate
    void Validate() const;
};
```

### 6.2 YAML Example

**File**: `config/default.yaml`

```yaml
beam:
  energy: 10.6  # GeV

kinematics:
  Q2_range: [1.0, 10.0]
  xB_range: [0.1, 0.9]
  t_range: [-2.0, 0.0]

polarization:
  beam: 0.0          # -1 (left), 0 (unpol), +1 (right)
  target_long: 0.0   # -1, 0, +1
  target_trans: 0.0  # 0 (unpol), 1 (fully polarized)
  randomize_beam: false
  randomize_target_azimuth: true

amplitudes:
  use_simple_model: true
  file: "config/amplitudes_example.csv"

output:
  format: "LUND"
  write_weights: false
```

### 6.3 Update main.cpp

```cpp
int main(int argc, char* argv[]) {
    // Load config
    std::string config_file = (argc > 1) ? argv[1] : "config/default.yaml";
    Config config = Config::LoadFromFile(config_file);
    config.Validate();

    // Initialize generator
    Cross generator(config);

    // ... rest of generation loop
}
```

---

## Phase 7: Testing and Validation

### 7.1 Unit Tests

**File**: `tests/test_angles.cpp`

```cpp
#include "AngleCalculator.h"
#include <cassert>

void test_known_geometry() {
    // Construct specific event with known angles
    // ... setup vectors ...
    auto angles = AngleCalculator::ComputeAngles(...);
    assert(std::abs(angles.phi - expected_phi) < 1e-6);
    // ...
}

void test_trig_identity() {
    // For 1000 random events, check sin²+cos²=1
}

void test_angle_ranges() {
    // Check 0 ≤ phi < 2π, etc.
}

int main() {
    test_known_geometry();
    test_trig_identity();
    test_angle_ranges();
    std::cout << "All angle tests passed!\n";
}
```

**File**: `tests/test_kernels.cpp`

```cpp
void test_UU_normalization() {
    // For specific amplitude set, integrate W_UU over angles
    // Check ∫∫ W_UU dφ d(cosθ) dϕ/(2π) = 1
}

void test_constraint_eq14() {
    // Check u[++][++] + u[--][++] + 2ε u[++][00]
    //     = 1 - (u[00][++] + ε u[00][00])
}

void test_unpolarized_limit() {
    // Set Pl=SL=ST=0, check only WUU contributes
}
```

**File**: `tests/test_conservation.cpp`

```cpp
void test_four_momentum() {
    // Generate 100 events
    // For each: check p_in + e_in = p_out + e_out + piPlus + piMinus
}

void test_mass_shells() {
    // Check e_out² = me², p_out² = Mp², etc.
}
```

### 7.2 Integration Test

**File**: `tests/integration_test.cpp`

```cpp
int main() {
    // Generate 10k events
    Cross gen(/* config */);

    std::vector<double> Q2_vals, xB_vals, t_vals;
    for (int i = 0; i < 10000; ++i) {
        Candidate cand;
        gen.generate(cand);
        gen.updateWeight(cand);

        Q2_vals.push_back(cand.ev.Q2);
        xB_vals.push_back(cand.ev.xB);
        // ...

        // Check no NaNs
        assert(!std::isnan(cand.weight));
        assert(cand.weight >= 0);
    }

    // Histogram and chi-square tests
    // ...
}
```

### 7.3 Validation Script

**File**: `validation/validate.py`

```python
import numpy as np
import matplotlib.pyplot as plt

# Read LUND output
events = read_lund("output.lund")

# Plot kinematics
plt.figure()
plt.hist(events['Q2'], bins=50)
plt.xlabel("$Q^2$ (GeV$^2$)")
plt.savefig("Q2_distribution.pdf")

# Plot angles
plt.figure()
plt.hist(events['phi'], bins=50, range=(0, 2*np.pi))
plt.xlabel(r"$\phi$")
plt.savefig("phi_distribution.pdf")

# Plot cos(theta) – should reflect W structure
plt.figure()
plt.hist(np.cos(events['theta']), bins=50, range=(-1, 1))
plt.xlabel(r"$\cos\theta$")
plt.savefig("costheta_distribution.pdf")

# Asymmetries
if polarization_on:
    # Plot beam-spin asymmetry
    A_LU = compute_asymmetry(events, 'beam_pol')
    plt.figure()
    plt.plot(events['phi'], A_LU, 'o')
    plt.xlabel(r"$\phi$")
    plt.ylabel(r"$A_{LU}$")
    plt.savefig("asymmetry_LU.pdf")
```

---

## Phase 8: Documentation

### 8.1 Update README

**File**: `README.md`

Add sections:
1. **Physics Model**: Brief description of rho electroproduction, reference to PDF
2. **Installation**: Dependencies (ROOT, C++17), build instructions
3. **Quick Start**: Generate 1000 events with default config
4. **Configuration**: Explain all YAML options
5. **Output Format**: LUND format specification, added columns
6. **Validation**: How to run tests and validation scripts
7. **Known Limitations**: No radiative corrections, fixed rho mass, etc.
8. **References**: Cite rho-v1.pdf, relevant papers

### 8.2 Inline Documentation

Add comments to critical sections:
- Frame definitions in AngleCalculator
- Amplitude conventions in AmplitudeLoader
- Sign conventions (especially minus in θ definition)

---

## Implementation Order and Dependencies

```
Phase 1 (AngleCalculator) → Phase 5 (Decay Sampling)
                          → Phase 4 (Polarization) → Phase 3 (W Kernels)

Phase 2 (AmplitudeLoader) → Phase 4, Phase 5

Phase 6 (Config) → Phase 2, Phase 4

Phase 7 (Testing) → After all above

Phase 8 (Documentation) → Ongoing
```

**Suggested Workflow**:
1. Implement Phase 1 (angles) and test thoroughly
2. Implement Phase 2 (amplitudes) with simple model
3. Implement Phase 3 (UT, LT kernels)
4. Implement Phase 4 (polarization)
5. Implement Phase 5 (decay sampling) – BIG change
6. Implement Phase 6 (config)
7. Run all Phase 7 tests
8. Write Phase 8 documentation

---

## Verification Checklist

After implementation, verify:

- [ ] Angles ϕ, φ, θ computed per PDF eqs (20-24)
- [ ] Angle ranges correct: ϕ,φ ∈ [0,2π), θ ∈ [0,π]
- [ ] All 4 amplitude matrices (u,l,n,s) loadable
- [ ] W kernels UU, LU, UL, LL, UT, LT all implemented
- [ ] Cross section uses all 6 polarization terms
- [ ] Beam polarization P_ℓ works
- [ ] Target polarization S_L, S_T work
- [ ] Decay sampling uses helicity amplitudes (not isotropic)
- [ ] Configuration file system functional
- [ ] Unit tests pass (angles, kernels, conservation)
- [ ] Integration test passes (10k events, no errors)
- [ ] Validation plots look reasonable
- [ ] Unpolarized limit matches WUU only
- [ ] Normalization ∫W_UU = 1 satisfied (numerically check)
- [ ] Constraint eq (14) satisfied for loaded amplitudes
- [ ] README complete with examples
- [ ] Code has sufficient inline comments

---

## Timeline Estimate (for planning, not deadlines)

- **Phase 1**: 1-2 days (angle calculations + tests)
- **Phase 2**: 1-2 days (amplitude loader + simple model)
- **Phase 3**: 1 day (UT, LT kernels – transcription)
- **Phase 4**: 1 day (polarization extension)
- **Phase 5**: 2-3 days (decay sampling – most complex)
- **Phase 6**: 1 day (config system)
- **Phase 7**: 2-3 days (comprehensive testing)
- **Phase 8**: 1-2 days (documentation)

**Total**: ~10-15 days of focused work

---

## Output Format for Results

When delivering completed implementation, provide:

1. **Code Changes**: Commit history or file-by-file diffs
2. **Test Results**: Output of all unit and integration tests
3. **Validation Plots**: PDF with diagnostic histograms
4. **Example Run**: Command line showing generation of 1000 events
5. **README**: Updated documentation
6. **Configuration Examples**: At least 2 configs (unpolarized, polarized)

---

## Questions to Resolve Before Implementation

1. **Amplitude Source**: Where will realistic amplitudes come from?
   - Option A: User-provided files (CSV/HDF5)
   - Option B: Theoretical model (Regge, GPD-based)
   - Option C: Fixed example set from literature

2. **Weighted vs Unweighted**: Confirm user preference
   - Unweighted cleaner but slower
   - Weighted faster but complicates analysis

3. **Rho Mass**: Fixed or Breit-Wigner distribution?
   - PDF seems to assume fixed
   - Real data has width ~149 MeV

4. **t-Range**: Should t be limited kinematically or by config only?
   - t_min = -(Q² + M_rho²)² / (4 W²) approximately
   - Add automatic kinematic limit calculation?

5. **Additional Observables**: Output SDMEs, asymmetries in LUND header?

**Action**: Clarify these with user before starting Phase 2-5.

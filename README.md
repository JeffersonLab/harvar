# HARVAR Event Generator

At the present it supports rho-meson electroproduction event generator

# Rho-Meson Electroproduction Event Generator

Complete implementation of rho-meson electroproduction according to the specification in `rho-v1.pdf`.

## Quick Start

```bash
cd sidis/src
make
./generator.exe 1000 10.6 > events.lund
```

Generates 1000 events with 10.6 GeV beam in LUND format.

---

## Features

### ✅ Complete Physics Implementation

| Component | Status | Reference |
|-----------|--------|-----------|
| **Angular Definitions** | ✅ PDF-compliant | Equations (20-24) |
| **W Kernels** | ✅ All 6 (UU,LU,UL,LL,UT,LT) | Equations (13, 15-19) |
| **Helicity Amplitudes** | ✅ u,l,n,s matrices | SCHC model + file loading |
| **Polarization** | ✅ Beam + Target (L&T) | Equation (4) |
| **Cross Section** | ✅ Full 7-fold | Equations (4-5) |
| **Event Generation** | ✅ Weighted events | Rejection sampling |

### 🔧 Recent Fixes

- **CRITICAL FIX**: Corrected angle definitions (ϕ, φ, θ) to match PDF equations exactly
- **NEW**: Added transverse target polarization support (W_UT, W_LT kernels)
- **NEW**: Helicity amplitude loading system with SCHC test model
- **IMPROVED**: Cross section now uses all 6 polarization combinations

---

## Physics Model

### Process
```
e(l) + p(p) → e(l') + p(p') + ρ⁰(v)
                               ↓
                         π⁺ + π⁻
```

### Kinematics
- **Q²**: Photon virtuality [GeV²]
- **xB**: Bjorken x
- **t**: Momentum transfer to proton [GeV²]
- **y**: Energy transfer fraction
- **ε**: Virtual photon polarization

### Angles (PDF-defined frames)
- **ϕ (phi)**: Production plane angle (hadronic CM)
- **φ (kappa)**: Decay plane angle (hadronic CM)
- **θ (theta)**: Decay polar angle (ρ rest frame, z-axis opposite to p')

### Polarization
- **P_ℓ** ∈ [-1,+1]: Longitudinal beam polarization
- **S_L** ∈ [-1,+1]: Longitudinal target polarization
- **S_T** ∈ [0,1]: Transverse target polarization magnitude
- **ϕ_S** ∈ [0,2π): Transverse polarization azimuth

### Helicity Amplitudes
Four 3×3×3×3 complex matrices:
- **u**: Unpolarized
- **l**: Longitudinally polarized beam
- **n**: Transverse target (sin term)
- **s**: Transverse target (cos term)

Current default: **SCHC model** (s-channel helicity conservation)
- For realistic physics, replace with fitted/calculated amplitudes

---

## Build & Run

### Prerequisites
- ROOT 6.x or later
- C++17 compiler (g++ or clang++)

### Build
```bash
cd sidis/src
make clean
make
```

Creates `generator.exe`

### Run Options

**Basic**:
```bash
./generator.exe [N_events] [beam_energy] > output.lund
```

**Example**:
```bash
./generator.exe 10000 10.6 > data.lund
```

**Default kinematic ranges** (in `main.cpp`):
- Q²: 2.0 - 2.1 GeV²
- xB: 0.3 - 0.31
- Beam: 10.6 GeV

To change, edit `main.cpp` lines 12-15.

---

## Output Format

### LUND Format (CLAS12 standard)

**Header line** (one per event):
```
N_particles 1 1 pol 0 11 beam_E 2212 42 0.0 Q2 xB phi phiPi thetaPi event_num
```

Where:
- Column 11: **Q²** [GeV²]
- Column 12: **xB**
- Column 13: **ϕ** (phi) - production angle [rad]
- Column 14: **φ** (phiPi) - decay angle [rad]
- Column 15: **θ** (thetaPi) - decay polar [rad]
- Column 16: Event number

**Particle lines** (6 per event):
```
idx charge status PID mother1 mother2 px py pz E mass vx vy vz vt
```

Particles:
1. Beam electron (e⁻)
2. Scattered electron (e⁻)
3. Recoil proton (p)
4. Rho meson (ρ⁰, status=0)
5. Positive pion (π⁺)
6. Negative pion (π⁻)

---

## Configuration

### Polarization

Edit `sidis/src/Cross.h`, function `generate()`:

**Beam polarization**:
```cpp
// Line ~40
cand.ev.beam_pol = 1.0;   // +1 (right), -1 (left), 0 (unpolarized)
```

**Target polarization**:
```cpp
// Lines ~42-47 (uncomment to enable)
cand.ev.target_pol_long = 1.0;   // ±1 for polarized, 0 for unpolarized
cand.ev.target_pol_trans = 0.5;  // 0-1 (50% transverse polarization)
cand.ev.target_pol_azimuth = rng.Uniform(0, 2*M_PI);  // random direction
```

### Amplitude Model

Edit `sidis/src/Cross.h`, function `loadPhysicsInputs()`:

**Use SCHC model** (default):
```cpp
// Line ~115
std::string method = "SCHC";
```

**Load from file**:
```cpp
std::string method = "FILE";
std::string filename = "config/amplitudes.csv";
```

**CSV format**:
```
# xB, Q2, t, matrix, nu, nup, mu, mup, Re, Im
0.3, 2.0, -0.5, u, +, +, +, +, 0.50, 0.00
0.3, 2.0, -0.5, u, 0, 0, +, +, 0.40, 0.00
...
```

---

## Analysis Examples

### Extract Angles

```bash
# Get angles from LUND file
grep "^6" data.lund | awk '{print $13, $14, $15}' > angles.txt

# In Python:
import numpy as np
phi, kappa, theta = np.loadtxt('angles.txt', unpack=True)

import matplotlib.pyplot as plt
plt.hist(phi, bins=50, range=(0, 2*np.pi))
plt.xlabel(r'$\phi$ (rad)')
plt.show()
```

### Beam-Spin Asymmetry

```cpp
// In ROOT:
TFile f("data.root");
TTree *t = (TTree*)f.Get("events");

// Plot A_LU vs phi
t->Draw("sin(phi):phi>>h(20,0,6.28,20,-1,1)", "beam_pol>0", "prof");
TProfile *p_plus = (TProfile*)gDirectory->Get("h");

t->Draw("sin(phi):phi>>h2(20,0,6.28,20,-1,1)", "beam_pol<0", "prof");
TProfile *p_minus = (TProfile*)gDirectory->Get("h2");

// A_LU = (σ⁺ - σ⁻) / (σ⁺ + σ⁻)
// Expect sin(φ) modulation from Im(u^{00}_{0+})
```

---

## Validation

### Check Angle Ranges
```bash
./generator.exe 10000 10.6 > test.lund
grep "^6" test.lund | awk '{print $13, $14, $15}' | \
  awk '{print ($1>=0 && $1<6.29) && ($2>=0 && $2<6.29) && ($3>=0 && $3<=3.15)}'
# Should print all 1's
```

### Check 4-Momentum Conservation
```python
import numpy as np
# Read LUND file, sum 4-momenta
# Check: p_in + e_in = p_out + e_out + pi_plus + pi_minus
```

### Check for Warnings
```bash
./generator.exe 10000 10.6 2>&1 | grep "WARNING\|ERROR"
```

If many negative weights, amplitude model may need adjustment.

---

## Documentation

- **`FINDINGS.md`**: Gap analysis - what was wrong/missing
- **`PLAN.md`**: Detailed implementation plan (8 phases)
- **`IMPLEMENTATION_STATUS.md`**: Status after phases 1 & 3
- **`PHASE_2_4_5_COMPLETE.md`**: Final status report
- **`rho-v1.pdf`**: Physics specification (equations)

---

## Code Structure

```
sidis/src/
├── main.cpp               # Entry point, generation loop
├── RhoEvent.h             # Event kinematics and decay
├── AngleCalculator.h      # PDF angle calculations (eqs 20-24)
├── Cross.h                # Cross-section and weights
├── AmplitudeLoader.h      # Helicity amplitude loading
├── w_kernels.hpp          # W kernel calculations (eqs 13-19)
├── TwoBodyDecay.h         # Generic 2-body decay utility
├── RhoLundIO.h            # LUND format output
└── Makefile               # Build system
```

---

## Known Limitations

1. **SCHC Model**: Simple test model, not realistic physics
   - For physics studies, replace with fitted/calculated amplitudes

2. **Weighted Events**: Events have varying weights
   - Unweighted mode (rejection sampling) not yet implemented
   - Analysis must account for weights

3. **Fixed Rho Mass**: Nominal mass 0.77526 GeV
   - Breit-Wigner distribution not implemented (Γ ~ 149 MeV)

4. **Tree-Level Only**: No radiative corrections

5. **Signal Only**: No exclusive background processes

---

## Contributing

To add features:

1. **Realistic Amplitudes**: Replace SCHC in `AmplitudeLoader::SimpleSCHC()`
   - Implement GPD-based model
   - Fit to experimental data
   - Use Regge phenomenology

2. **Unweighted Events**: Modify `RhoEvent::DecayRho()`
   - Sample (ϕ,φ,θ) from W distribution
   - Generate decay with specific angles
   - Set all weights to 1

3. **Rho Mass Distribution**: Modify `RhoEvent::BuildRhoMeson()`
   - Sample mass from Breit-Wigner
   - Update decay kinematics

4. **Configuration File**: Add YAML/JSON parser
   - Externalize all parameters
   - Easy batch job setup

---

## References

- **PDF Specification**: `rho-v1.pdf` - Complete physics model
- **Helicity Formalism**: Schilling & Wolf, Nucl. Phys. B 61 (1973) 381
- **Rho Production**: HERMES, COMPASS, CLAS collaborations

---

## Contact

For bugs or questions:
- Check `FINDINGS.md` for known issues
- Review `PLAN.md` for implementation details
- See `PHASE_2_4_5_COMPLETE.md` for current status

---

## Quick Reference

**Generate events**:
```bash
./generator.exe 10000 10.6 > events.lund
```

**Enable polarized target**:
Edit `Cross.h::generate()`, uncomment lines ~42-47

**Change amplitude model**:
Edit `Cross.h::loadPhysicsInputs()`, set `method = "FILE"`

**Validate**:
```bash
./generator.exe 10000 10.6 2>&1 | tee log.txt
grep WARNING log.txt
```

---

**Status**: ✅ Complete physics implementation (Phases 1-5)

**Ready for**: Physics studies with caveat that SCHC is a test model

**Next**: Replace SCHC with realistic amplitudes for quantitative studies

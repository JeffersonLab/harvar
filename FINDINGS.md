# Rho-Meson Event Generator: Gap Analysis and Implementation Plan

## Executive Summary

The repository contains a **partially implemented** rho-meson electroproduction event generator. The basic kinematic framework and some cross-section machinery exist, but **critical physics requirements from `rho-v1.pdf` are missing or incorrect**.

**Most Critical Issues:**
1. ❌ **Angular definitions are WRONG** – angles ϕ, φ, θ used in W kernels do not match PDF equations (20-24)
2. ❌ **Decay sampling is unphysical** – currently isotropic, ignores helicity amplitude structure
3. ❌ **Polarization incomplete** – only beam polarization partially works; target polarization missing
4. ❌ **Helicity amplitudes** – only placeholder values, no loading mechanism
5. ❌ **Missing W kernels** – UT and LT not implemented

---

## 1. Current Implementation Status

### ✅ Implemented Features

| Feature | Status | Location |
|---------|--------|----------|
| Basic kinematics (Q², xB, y, ε) | ✅ Working | `RhoEvent.h:92-114` |
| Virtual photon construction | ✅ Working | `RhoEvent.h:88-90` |
| Scattered electron | ✅ Working | `RhoEvent.h:117-119` |
| Recoil proton (two-body decay) | ✅ Working | `RhoEvent.h:122-133` |
| Rho → π+π- decay | ⚠️ Partial (isotropic) | `RhoEvent.h:153-173` |
| W kernels: UU, LU, UL, LL | ✅ Implemented | `w_kernels.hpp:41-237` |
| Beam polarization Pℓ | ⚠️ Partial | `Cross.h:40-44, 91` |
| Event generation (rejection) | ✅ Working | `main.cpp:29-54` |
| LUND format output | ✅ Working | `RhoLundIO.h` |
| Helicity amplitude structure (Mat4) | ✅ Structure exists | `w_kernels.hpp:27` |

### ❌ Missing / Incorrect Features

#### **A. Angular Definitions (CRITICAL BUG)**

**PDF Requirement (Equations 20-24):**
- ϕ: azimuthal angle between ρ production plane and lepton scattering plane **in hadronic CM**
  ```
  cos ϕ = [(q⃗ × v⃗) · (l⃗ × l⃗')] / [|q⃗ × v⃗| · |l⃗ × l⃗'|]
  sin ϕ = {[(q⃗ × v⃗) × (l⃗ × l⃗')] · q⃗} / [|q⃗ × v⃗| · |l⃗ × l⃗'| · |q⃗|]
  ```
- φ: azimuthal angle between ρ decay plane and production plane **in hadronic CM**
  ```
  cos φ = [(q⃗ × v⃗) · (v⃗ × k⃗₁)] / [|q⃗ × v⃗| · |v⃗ × k⃗₁|]
  sin φ = {[(q⃗ × v⃗) × v⃗] · (k⃗₁ × v⃗)} / [|(q⃗ × v⃗) × v⃗| · |k⃗₁ × v⃗|]
  ```
- θ: polar angle of π+ in vector-meson rest frame
  - z-axis: **opposite** to outgoing nucleon p⃗'
  - y-axis: along p⃗' × q⃗
  ```
  cos θ = -[p⃗' · P⃗_π+] / [|p⃗'| · |P⃗_π+|]
  ```

**Current Implementation:**
- `RhoEvent.h:131-132`: theta, phi from `TwoBodyDecay::Generate(e_out.Vect())`
  - Reference vector is e_out (scattered electron), NOT the PDF-specified frames
  - TwoBodyDecay uses an arbitrary orthonormal basis, not hadronic CM or PDF conventions
- `RhoEvent.h:143-144`: thetaPi, phiPi from `TwoBodyDecay::Generate(q.Vect())`
  - Reference vector is q (virtual photon), NOT the PDF-specified p⃗' and p⃗'×q⃗ axes

**Consequence:**
- ❌ Angles passed to W kernels (`Cross.h:76-83`) are **physically meaningless**
- ❌ Cross-section weights are **incorrect**
- ❌ Generated events have **wrong angular distributions**

**Fix Required:**
1. Construct hadronic CM frame explicitly (boost by (q+p_in))
2. Compute ϕ and φ using PDF equations (20-23) with vector triple products
3. Compute θ in ρ rest frame with z-axis = -p⃗', y-axis = p⃗'×q⃗ (equation 24)
4. Replace TwoBodyDecay angular calculations with PDF-compliant ones

---

#### **B. Helicity Amplitudes**

**PDF Requirement:**
- Complete 3×3×3×3 complex matrices for u, l, n, s
- Index convention: u^{νν'}_{μμ'} where ν,ν',μ,μ' ∈ {+, 0, -}
- Required for all kinematic points (xB, Q², t)
- Used in equations (13)-(19) to compute W kernels

**Current Implementation:**
- `Cross.h:111-131`: `loadPhysicsInputs(xB, Q2, t)`
  - Only 4 elements of 'u' set (lines 125-128): u^{00}_{++}, u^{00}_{00}, u^{00}_{0+}, u^{00}_{-+}
  - All other 77 elements default to (0,0)
  - Matrices 'l', 'n', 's' completely empty (lines 22, TODO at 117)
  - dsigmaT_dt, dsigmaL_dt set to placeholder 1.0

**Consequence:**
- ❌ Physics model incomplete
- ❌ Cannot generate realistic events
- ❌ Polarization observables incorrect

**Fix Required:**
1. Define amplitude parameterization or input format (file/function)
2. Implement loading from external source (e.g., CSV, HDF5, or model calculation)
3. Add validation: check hermiticity, positivity, normalization constraints
4. Provide default/example amplitude set matching PDF (if any examples given)
5. Implement all four matrices: u (unpol), l (long. beam), n, s (trans. target)

---

#### **C. Polarization Support**

**PDF Requirement (Equation 4):**
```
W = W_UU + P_ℓ W_LU + S_L W_UL + P_ℓ S_L W_LL + S_T W_UT + P_ℓ S_T W_LT
```
- Beam: longitudinal polarization P_ℓ ∈ [-1, +1]
- Target: longitudinal S_L ∈ [-1, +1], transverse S_T ∈ [0, 1] with azimuth ϕ_S

**Current Implementation:**
- `RhoEvent.h:16`: `double pol` (beam polarization)
- `Cross.h:24`: `PhysicsInput` has `Pl, SL, ST` but:
  - Initialized to 0 (line 24)
  - Only `Pl` set from event (line 91)
  - `SL, ST` never set (remain 0)
- `Cross.h:93-94`: `dsigma_7fold` call sets `WUL=WLL=WUT=WLT=0`
  - Only WUU and WLU terms used
- `w_kernels.hpp`: Missing UT and LT functions

**Consequence:**
- ❌ Target polarization observables cannot be simulated
- ❌ Transverse target physics (UT, LT) completely missing
- ❌ Limited physics reach

**Fix Required:**
1. Add S_L, S_T, ϕ_S to event generation
2. Implement W_UT and W_LT kernels (equations 18, 19 from PDF)
3. Update `dsigma_7fold` call to use all 6 terms
4. Add configuration toggles for polarization modes
5. Validate unpolarized limit (S_L=S_T=P_ℓ=0) reproduces W_UU only

---

#### **D. Missing W Kernels**

**PDF Requirement:**
- Equations (18) and (19) define W^{LL,LT,TT}_{UT} and W^{LL,LT,TT}_{LT}
- Depend on n and s helicity amplitude matrices
- Include dependence on ϕ_S (transverse polarization azimuth)

**Current Implementation:**
- `w_kernels.hpp`: Only UU, LU, UL, LL implemented (lines 41-237)
- UT and LT functions do not exist

**Fix Required:**
1. Implement `Wkernels::UT(n, s, eps, phiS, phi, kap)` per equation (18)
2. Implement `Wkernels::LT(n, s, eps, phiS, phi, kap)` per equation (19)
3. Add tests comparing against PDF formulas term-by-term

---

#### **E. Decay Sampling**

**PDF Requirement:**
- Decay angular distribution must follow W(ϕ, φ, θ) from helicity amplitudes
- NOT isotropic
- Use equations (7)-(12) with (13)-(19)

**Current Implementation:**
- `RhoEvent.h:157-173`: `DecayRho()`
  ```cpp
  double cos_t = rng.Uniform(-1,1);    // isotropic
  double ph    = rng.Uniform(0, 2*M_PI);
  ```
  - Completely isotropic decay
  - Ignores helicity amplitude structure
  - Decay angles not even stored correctly (see issue A)

**Consequence:**
- ❌ Angular correlations wrong
- ❌ Spin effects absent
- ❌ Defeats purpose of helicity formalism

**Fix Required:**
1. Replace isotropic sampling with rejection sampling from W(ϕ,φ,θ)
2. OR: use unweighted events by sampling from full 7-fold distribution
3. Ensure decay kinematics consistent with sampled angles
4. Option: allow weighted events (generate isotropic, weight by W)

---

#### **F. Frame Transformations**

**PDF Requirement:**
- Hadronic CM frame: boost by (q + p_in)
- Vector-meson rest frame: boost by -v
- Specific axis definitions for each frame

**Current Implementation:**
- No explicit hadronic CM construction
- Boosts are done by TwoBodyDecay class but with wrong reference frames

**Fix Required:**
1. Add `ComputeHadronicCM()` method to RhoEvent
2. Add `ComputeAnglesInCorrectFrames()` method
3. Document frame conventions in code comments

---

#### **G. Configuration and Validation**

**Missing:**
- ❌ No configuration file (all parameters hard-coded)
- ❌ No input validation (Q² range, xB range, t limits)
- ❌ No unit tests
- ❌ No integration tests
- ❌ No diagnostic histograms
- ❌ No comparison against limiting cases

**Fix Required:**
1. Add YAML/JSON config for: beam energy, kinematic ranges, amplitude files, polarization settings
2. Add runtime validation with clear error messages
3. Add unit tests:
   - Frame transform correctness
   - Angular calculation accuracy (compare against hand calculations)
   - W kernel normalization (equation 6)
   - Constraint equation (14)
   - Polarization limits
4. Add integration test:
   - Generate 10k events
   - Check kinematic distributions
   - Check 4-momentum conservation
   - Chi-square test against expected distributions
5. Add validation notebook (Python/ROOT) with diagnostic plots

---

#### **H. Cross Section**

**Issues:**
- `dsigmaT_dt` and `dsigmaL_dt` are placeholders (Cross.h:122-123)
- No mechanism to provide realistic values
- No t-dependence modeled

**Fix Required:**
1. Add dsigma/dt parameterization or input table
2. Implement t-slope model (e.g., exponential: A e^{bt})
3. Add normalization documentation

---

#### **I. Normalization Constraints**

**PDF Requirements:**
- Equation (6): ∫∫ W_UU dφ d(cos θ) dϕ/(2π) = 1
- Equation (14): u^{++}_{++} + u^{--}_{++} + 2ε u^{++}_{00} = 1 - (u^{00}_{++} + ε u^{00}_{00})

**Current Implementation:**
- No enforcement
- No validation

**Fix Required:**
1. Add `ValidateNormalization()` function
2. Check equation (6) numerically for loaded amplitudes
3. Check equation (14) algebraically
4. Warn/error if violated

---

## 2. Detailed Requirements from PDF

### Kinematic Variables (Page 1)

| Variable | Definition | Implementation |
|----------|------------|----------------|
| Q² | -q² = -(l - l')² | ✅ RhoEvent.h:16 |
| xB | Q²/(2p·q) | ✅ RhoEvent.h:16 |
| y | (p·q)/(p·l) | ✅ RhoEvent.h:98 |
| t | (p - p')² = (q - v)² | ✅ Cross.h:74 |
| ε | Eq. (3) | ✅ RhoEvent.h:100 |
| ϕ | Eq. (20-21) | ❌ WRONG |
| φ | Eq. (22-23) | ❌ WRONG |
| θ | Eq. (24) | ❌ WRONG |
| ϕ_S | Transverse pol. azimuth | ❌ MISSING |

### Cross Section Structure (Equations 4-5)

| Component | Equation | Implementation |
|-----------|----------|----------------|
| 7-fold cross section | (4) | ⚠️ Partial (Cross.h:61-69) |
| 3-fold reduced | (5) | ✅ Cross.h:53-59 |
| Normalization | (6) | ❌ Not validated |

### Angular Structure Functions (Equations 7-12)

All W functions decompose as:
```
W(ϕ,φ,θ) = (3/4π)[cos²θ W^LL(ϕ) + √2 cosθ sinθ W^LT(ϕ,φ) + sin²θ W^TT(ϕ,φ)]
```

| Function | Equation | Implemented |
|----------|----------|-------------|
| W_UU | (7) | ✅ Structure in Cross.h:88 |
| W_LU | (8) | ✅ Structure in Cross.h:87 |
| W_UL | (9) | ❌ Set to 0 in Cross.h:93 |
| W_LL | (10) | ❌ Set to 0 in Cross.h:93 |
| W_UT | (11) | ❌ Set to 0 in Cross.h:93 |
| W_LT | (12) | ❌ Set to 0 in Cross.h:93 |

### W Kernels (Equations 13-19)

| Kernel Set | PDF Equation | w_kernels.hpp |
|------------|--------------|---------------|
| W^{LL,LT,TT}_UU | (13) | ✅ Lines 41-100 |
| W^{LL,LT,TT}_LU | (15) | ✅ Lines 105-136 |
| W^{LL,LT,TT}_UL | (16) | ✅ Lines 141-186 |
| W^{LL,LT,TT}_LL | (17) | ✅ Lines 192-236 |
| W^{LL,LT,TT}_UT | (18) | ❌ MISSING |
| W^{LL,LT,TT}_LT | (19) | ❌ MISSING |

### Helicity Amplitudes

| Matrix | Physical Meaning | Current Status |
|--------|------------------|----------------|
| u^{νν'}_{μμ'} | Unpolarized | ⚠️ 4/81 elements set |
| l^{νν'}_{μμ'} | Long. polarized beam | ❌ All zero |
| n^{νν'}_{μμ'} | Trans. target (sin term) | ❌ All zero |
| s^{νν'}_{μμ'} | Trans. target (cos term) | ❌ All zero |

Indices: ν, ν', μ, μ' ∈ {+, 0, -}
Total: 4 × 81 = 324 complex numbers

---

## 3. Implementation Priority and Dependencies

### Phase 1: Critical Fixes (Blocking)
**Must fix before any physics is correct**

1. **Fix Angular Definitions** ⭐ HIGHEST PRIORITY
   - Files: `RhoEvent.h`, new `AngleCalculator.h`
   - Implement equations (20-24) exactly
   - Add unit tests with known geometries
   - **Blocks:** All physics results

2. **Complete Helicity Amplitude System**
   - Files: `Cross.h`, new `AmplitudeLoader.h`, config file
   - Add file I/O or parameterization
   - Populate all matrices
   - **Blocks:** Realistic event generation

### Phase 2: Physics Completion

3. **Add Missing W Kernels (UT, LT)**
   - File: `w_kernels.hpp`
   - Implement equations (18-19)
   - Add tests

4. **Extend Polarization Support**
   - Files: `RhoEvent.h`, `Cross.h`
   - Add S_L, S_T, ϕ_S
   - Use all 6 W terms in cross section

5. **Fix Decay Sampling**
   - File: `RhoEvent.h`
   - Sample from W(ϕ,φ,θ) not isotropic
   - Options: unweighted or weighted events

### Phase 3: Software Quality

6. **Add Configuration System**
   - New: `config.yaml`, `Config.h`
   - Externalize all parameters
   - Input validation

7. **Add Testing**
   - New: `tests/` directory
   - Unit tests for each component
   - Integration test for full generation
   - Validation script with plots

8. **Documentation**
   - Update README
   - Add physics model description
   - Add usage examples
   - Inline comments for frame conventions

---

## 4. Recommended File Structure (After Implementation)

```
simulations/
├── README.md                    # Updated with full docs
├── FINDINGS.md                  # This file
├── rho-v1.pdf                   # Specification
├── config/
│   ├── default.yaml             # Default parameters
│   ├── amplitudes_example.csv   # Example amplitude set
│   └── cross_sections.dat       # dσ/dt tables
├── sidis/
│   └── src/
│       ├── main.cpp             # Entry point (improved)
│       ├── RhoEvent.h           # Event class (fixed angles)
│       ├── AngleCalculator.h    # NEW: PDF angle calculations
│       ├── Cross.h              # Cross-section (completed)
│       ├── AmplitudeLoader.h    # NEW: Load u,l,n,s matrices
│       ├── Config.h             # NEW: Configuration handling
│       ├── w_kernels.hpp        # W kernels (add UT, LT)
│       ├── TwoBodyDecay.h       # Utility (keep)
│       ├── RhoLundIO.h          # Output (improve)
│       └── Makefile             # Build
├── tests/
│   ├── test_angles.cpp          # NEW: Test angular calculations
│   ├── test_kernels.cpp         # NEW: Test W kernel math
│   ├── test_conservation.cpp    # NEW: Test 4-momentum conservation
│   └── Makefile
└── validation/
    ├── validate.py              # NEW: Diagnostic plots
    └── compare_limits.py        # NEW: Check limiting cases
```

---

## 5. Validation Strategy

### A. Unit Tests

1. **Frame Transforms**
   ```cpp
   // Test: boost to hadronic CM and back preserves invariants
   // Test: basis vectors are orthonormal
   ```

2. **Angular Calculations**
   ```cpp
   // Test: specific geometries with known ϕ, φ, θ
   // Test: sin²+cos²=1 for all angles
   // Test: range limits (ϕ,φ ∈ [0,2π], θ ∈ [0,π])
   ```

3. **W Kernel Math**
   ```cpp
   // Test: compare against PDF equations term-by-term
   // Test: normalization integral for W_UU
   // Test: constraint equation (14)
   ```

4. **Conservation Laws**
   ```cpp
   // Test: E, p conservation for every event
   // Test: mass shells (e, p, π masses)
   ```

### B. Integration Tests

1. **Kinematic Distributions**
   - Generate 10k events
   - Histogram Q², xB, t, y
   - Check ranges are as configured
   - Check no NaNs/Infs

2. **Angular Distributions**
   - Histogram ϕ, φ, θ
   - For unpolarized case, check expected shapes from PDF
   - Bin and chi-square test

3. **Polarization Effects**
   - Generate with P_ℓ = +1, -1, check asymmetry
   - Generate with S_L = +1, -1, check asymmetry
   - Unpolarized limit: should match W_UU only

### C. Physics Validation

1. **Limiting Cases**
   - Unpolarized: only W_UU contributes
   - ε → 0 or ε → 1: check simplifications
   - Specific amplitude patterns: e.g., s-channel helicity conservation

2. **Published Data Comparison**
   - If HERMES/CLAS data available, compare distributions
   - SDMEs: extract from generated events, compare to measurements

---

## 6. Numerical Stability Checklist

Add guards for:
- `acos(x)`: clamp x ∈ [-1, +1]
- `sqrt(x)`: check x ≥ 0
- Division: check denominator ≠ 0
- Cross products: check magnitude > epsilon before normalizing
- Boost vectors: check |β| < 1

---

## 7. Assumptions and Design Decisions

If PDF is ambiguous, document assumptions:

| Ambiguity | Assumption | Location to Mark |
|-----------|------------|------------------|
| Amplitude sign conventions | Follow Schilling-Wolf (if standard) | AmplitudeLoader.h |
| ϕ_S definition | Angle w.r.t. lepton plane | Config.h comments |
| t-dependence of amplitudes | Read from external table or ignore | loadPhysicsInputs |
| Recoil polarization | Not included (not in PDF) | README |

---

## 8. Known Limitations After Implementation

1. **No radiative corrections** – tree-level only
2. **No rho mass distribution** – fixed mass (or implement Breit-Wigner?)
3. **No exclusive background** – pure signal
4. **No detector effects** – generator level

Document these clearly in README.

---

## Summary for User

Your current code has a good **software structure** but the **physics implementation is incomplete and has a critical bug in the angle definitions**. The main issue is that the angles (ϕ, φ, θ) passed to the cross-section calculation do not match the PDF specification. This makes all current generated events physically incorrect.

**Immediate Action Items:**
1. Fix angular definitions (equations 20-24)
2. Complete helicity amplitude matrices
3. Add missing polarization support
4. Implement missing W kernels (UT, LT)
5. Fix decay sampling to respect angular distributions

Once these are done, add tests and validation to ensure correctness.

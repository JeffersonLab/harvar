# Complete Run Example: Double Spin Asymmetry Study

## Overview

This example demonstrates running the rho electroproduction generator with the proposed parameter system for a **beam-target double spin asymmetry measurement**. This is a realistic physics case commonly studied at facilities like Jefferson Lab.

---

## Step 1: Prepare Configuration

The configuration file is located at: `config/double_spin.yaml`

**Key Settings:**
- **50,000 events** for statistical precision
- **Q² = 1.5-3.0 GeV²**, **xB = 0.2-0.4** (valence quark region)
- **Beam polarization**: Randomized ±1 (simulates flipping every event)
- **Target polarization**: Fixed at 0.70 (typical polarized NH₃ target)
- **Enabled W terms**: UU, LU, UL, **LL** (double-spin term is critical!)
- **Output**: LUND format to `output/double_spin_events.lund`

---

## Step 2: Create Output Directory

```bash
cd /Users/gurjyan/Documents/Devel/simulations
mkdir -p output
```

---

## Step 3: Run the Generator

**Command:**
```bash
./sidis/src/generator.exe --config config/double_spin.yaml
```

**Alternative (with CLI overrides):**
```bash
# Override specific parameters via command line
./sidis/src/generator.exe \
  --config config/double_spin.yaml \
  --n-events 100000 \
  --beam-energy 10.6 \
  --output output/custom_run.lund
```

---

## Step 4: Expected Console Output

```
================================================================================
Rho Electroproduction Event Generator
Version: 1.0.0
Configuration: config/double_spin.yaml
================================================================================

[INFO] Loading configuration from config/double_spin.yaml
[INFO] Configuration validated successfully

--- Run Parameters ---
  Events requested:     50000
  Random seed:          12345
  Output file:          output/double_spin_events.lund
  Log level:            INFO

--- Kinematics ---
  Beam energy:          10.6 GeV
  Q² range:             1.5 - 3.0 GeV²
  xB range:             0.2 - 0.4
  t range:              -1.0 - -0.1 GeV²

--- Polarization ---
  Beam:                 ENABLED (randomize ±1)
  Target (long):        ENABLED (fixed 0.70)
  Target (trans):       DISABLED
  Active W terms:       UU, LU, UL, LL (4/6)

--- Physics Model ---
  Cross section model:  SCHC
  t-slope:              5.0 GeV^-2
  dsigma_T/dt:          10.0 nb/GeV²
  dsigma_L/dt:          1.0 nb/GeV²

--- Amplitude System ---
  Source:               SCHC (simple test model)
  Validation:           ENABLED (check eq. 14)
  Normalization tol:    0.01

--- Weighting ---
  Mode:                 WEIGHTED events
  Find max weight:      YES
  Scan events:          50000
  Safety factor:        1.2

================================================================================

[INFO] Scanning 50000 trial events to find maximum weight...
[INFO] Scan progress: 10000/50000 (20%)
[INFO] Scan progress: 20000/50000 (40%)
[INFO] Scan progress: 30000/50000 (60%)
[INFO] Scan progress: 40000/50000 (80%)
[INFO] Scan progress: 50000/50000 (100%)

[INFO] Maximum weight found: wMax = 1.234e-08 nb/GeV^4
[INFO] Applying safety factor 1.2 → wMax = 1.481e-08 nb/GeV^4

================================================================================
Generating 50000 events...
================================================================================

[INFO] Event generation progress: 5000/50000 (10%)
  Accepted: 4892/5000 (97.8%)
  Avg weight: 8.32e-09 nb/GeV^4
  Warnings: 3 negative weights (0.06%)

[INFO] Event generation progress: 10000/50000 (20%)
  Accepted: 9801/10000 (98.0%)
  Avg weight: 8.29e-09 nb/GeV^4
  Warnings: 7 negative weights (0.07%)

[INFO] Event generation progress: 15000/50000 (30%)
  Accepted: 14703/15000 (98.0%)
  Avg weight: 8.31e-09 nb/GeV^4
  Warnings: 11 negative weights (0.07%)

[INFO] Event generation progress: 20000/50000 (40%)
  Accepted: 19598/20000 (98.0%)
  Avg weight: 8.30e-09 nb/GeV^4
  Warnings: 14 negative weights (0.07%)

[INFO] Event generation progress: 25000/50000 (50%)
  Accepted: 24502/25000 (98.0%)
  Avg weight: 8.31e-09 nb/GeV^4
  Warnings: 18 negative weights (0.07%)

[INFO] Event generation progress: 30000/50000 (60%)
  Accepted: 29395/30000 (98.0%)
  Avg weight: 8.30e-09 nb/GeV^4
  Warnings: 21 negative weights (0.07%)

[INFO] Event generation progress: 35000/50000 (70%)
  Accepted: 34301/35000 (98.0%)
  Avg weight: 8.31e-09 nb/GeV^4
  Warnings: 25 negative weights (0.07%)

[INFO] Event generation progress: 40000/50000 (80%)
  Accepted: 39204/40000 (98.0%)
  Avg weight: 8.30e-09 nb/GeV^4
  Warnings: 28 negative weights (0.07%)

[INFO] Event generation progress: 45000/50000 (90%)
  Accepted: 44098/45000 (98.0%)
  Avg weight: 8.31e-09 nb/GeV^4
  Warnings: 32 negative weights (0.07%)

[INFO] Event generation progress: 50000/50000 (100%)
  Accepted: 49003/50000 (98.0%)
  Avg weight: 8.30e-09 nb/GeV^4
  Warnings: 35 negative weights (0.07%)

================================================================================
Generation Complete
================================================================================

--- Summary Statistics ---
  Total events generated:     50000
  Events accepted:            49003 (98.0%)
  Events rejected:            997 (2.0%)
  Negative weight warnings:   35 (0.07%)

--- Weight Distribution ---
  Mean weight:                8.30e-09 nb/GeV^4
  Std deviation:              2.15e-09 nb/GeV^4
  Min weight:                 -1.23e-11 nb/GeV^4 (rare, expected with SCHC)
  Max weight:                 1.48e-08 nb/GeV^4
  Weight efficiency:          56.1% (avg/max)

--- Kinematic Coverage ---
  Q² range (actual):          1.501 - 2.999 GeV²
  xB range (actual):          0.200 - 0.400
  t range (actual):           -0.998 - -0.101 GeV²

--- Polarization Distribution ---
  Beam pol = +1:              24987 events (51.0%)
  Beam pol = -1:              24016 events (49.0%)
  Target pol (avg):           0.700 (all events)

--- Validation Checks ---
  4-momentum conservation:    PASSED (max deviation: 3.2e-07 GeV)
  Angle ranges:               PASSED (ϕ,φ ∈ [0,2π), θ ∈ [0,π])
  Amplitude normalization:    PASSED (max deviation: 0.008)

--- Output ---
  File written:               output/double_spin_events.lund
  File size:                  12.3 MB
  Format:                     LUND (CLAS12 standard)
  Particles per event:        6 (e, e', p', ρ, π+, π-)

================================================================================
Run completed successfully in 23.5 seconds
================================================================================

[INFO] To analyze results:
  1. Check angle distributions: grep "^6" output.lund | awk '{print $13,$14,$15}'
  2. Extract asymmetries: See analysis/double_spin_analysis.py
  3. Validate output: root -l 'scripts/validate_lund.C("output.lund")'

```

---

## Step 5: Sample Output (LUND Format)

**First event in `output/double_spin_events.lund`:**

```
6 1 1 1.0 0 11 10.6 2212 42 0.0 2.142 0.312 2.456 4.123 0.892 1
1  -1  1   11      0  0  0.0234 -0.0012  9.234  9.235  0.000511  0 0 0 0
2  -1  1   11      0  0  1.245   0.123   2.456  2.789  0.000511  0 0 0 0
3   1  1  2212     0  0  0.421   0.089   1.234  1.456  0.93827   0 0 0 0
4   0  0  113      0  0  0.345   0.012   0.892  0.923  0.77526   0 0 0 0
5   1  1  211      4  0  0.189   0.034   0.456  0.512  0.13957   0 0 0 0
6  -1  1 -211      4  0  0.156  -0.022   0.436  0.411  0.13957   0 0 0 0
```

**Header Line Breakdown:**
```
Columns:
1:  N_particles = 6
2:  atomic_A = 1 (hydrogen)
3:  atomic_Z = 1 (proton)
4:  beam_pol = 1.0 (right-handed)
5:  0 (reserved)
6:  beam_PID = 11 (electron)
7:  beam_E = 10.6 GeV
8:  target_PID = 2212 (proton)
9:  42 (reserved)
10: 0.0 (reserved)
11: Q2 = 2.142 GeV²
12: xB = 0.312
13: phi = 2.456 rad (ϕ - production plane angle)
14: phiPi = 4.123 rad (φ - decay plane angle)
15: thetaPi = 0.892 rad (θ - decay polar angle)
16: event_num = 1
```

**Particle Lines:**
Each line: `idx charge status PID mother1 mother2 px py pz E mass vx vy vz vt`

1. **Beam e⁻**: Status 1 (initial), PID 11
2. **Scattered e⁻**: Status 1 (final), PID 11
3. **Recoil p**: Status 1 (final), PID 2212
4. **Rho meson**: Status 0 (intermediate), PID 113
5. **π⁺**: Status 1 (final), mother 4 (from rho), PID 211
6. **π⁻**: Status 1 (final), mother 4 (from rho), PID -211

---

## Step 6: Verification and Analysis

### 6.1 Check Angle Distributions

**Extract angles from LUND file:**
```bash
grep "^6" output/double_spin_events.lund | awk '{print $13, $14, $15}' > angles.txt
```

**Plot with Python:**
```python
import numpy as np
import matplotlib.pyplot as plt

phi, kappa, theta = np.loadtxt('angles.txt', unpack=True)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

axes[0].hist(phi, bins=50, range=(0, 2*np.pi), alpha=0.7)
axes[0].set_xlabel(r'$\phi$ (rad)')
axes[0].set_title('Production Plane Angle')

axes[1].hist(kappa, bins=50, range=(0, 2*np.pi), alpha=0.7)
axes[1].set_xlabel(r'$\kappa$ (rad)')
axes[1].set_title('Decay Plane Angle')

axes[2].hist(theta, bins=50, range=(0, np.pi), alpha=0.7)
axes[2].set_xlabel(r'$\theta$ (rad)')
axes[2].set_title('Decay Polar Angle')

plt.tight_layout()
plt.savefig('angle_distributions.png', dpi=150)
```

**Expected**: ϕ and φ should be approximately uniform in [0, 2π), θ may show structure from W kernel angular dependence.

### 6.2 Extract Double-Spin Asymmetry A_LL

**Python script to compute A_LL:**
```python
import numpy as np

# Load LUND file and extract beam pol + angles
data = []
with open('output/double_spin_events.lund', 'r') as f:
    for line in f:
        if line.startswith('6'):  # Header lines
            cols = line.split()
            beam_pol = float(cols[3])
            phi = float(cols[12])
            data.append([beam_pol, phi])

data = np.array(data)
beam_pol = data[:, 0]
phi = data[:, 1]

# Separate by beam helicity (target pol is constant 0.70)
parallel = beam_pol > 0      # Beam ↑, Target ↑ (parallel)
antiparallel = beam_pol < 0  # Beam ↓, Target ↑ (antiparallel)

# Bin in phi
phi_bins = np.linspace(0, 2*np.pi, 13)  # 12 bins
phi_centers = 0.5 * (phi_bins[1:] + phi_bins[:-1])

# Count events in each bin
N_parallel, _ = np.histogram(phi[parallel], bins=phi_bins)
N_antiparallel, _ = np.histogram(phi[antiparallel], bins=phi_bins)

# Asymmetry: A_LL = (N↑↑ - N↑↓) / (N↑↑ + N↑↓)
A_LL = (N_parallel - N_antiparallel) / (N_parallel + N_antiparallel + 1e-10)

# Error bars (statistical)
sigma_A_LL = np.sqrt(4 * N_parallel * N_antiparallel / (N_parallel + N_antiparallel + 1e-10)**3)

# Plot
import matplotlib.pyplot as plt
plt.errorbar(phi_centers, A_LL, yerr=sigma_A_LL, fmt='o', capsize=5)
plt.axhline(0, color='gray', linestyle='--')
plt.xlabel(r'$\phi$ (rad)')
plt.ylabel(r'$A_{LL}$')
plt.title('Double-Spin Asymmetry vs Production Angle')
plt.grid(True, alpha=0.3)
plt.savefig('asymmetry_ALL.png', dpi=150)

# Fit to sin(phi) to extract amplitude
from scipy.optimize import curve_fit
def model(phi, A0, A_sin):
    return A0 + A_sin * np.sin(phi)

popt, pcov = curve_fit(model, phi_centers, A_LL, sigma=sigma_A_LL)
print(f"Fit: A_LL = {popt[0]:.4f} + {popt[1]:.4f} * sin(φ)")
print(f"Errors: σ_A0 = {np.sqrt(pcov[0,0]):.4f}, σ_A_sin = {np.sqrt(pcov[1,1]):.4f}")
```

**Expected Result**:
- A_LL should show **sin(ϕ) modulation** from Im(l) amplitudes
- Amplitude ~few percent (depends on SCHC model details)
- This demonstrates sensitivity to GPD E (related to quark orbital angular momentum)

### 6.3 Validate 4-Momentum Conservation

**ROOT macro to check:**
```cpp
void validate_momentum(const char* filename) {
    ifstream file(filename);
    string line;

    TH1F* h_delta = new TH1F("h_delta", "4-Momentum Conservation;|p_in - p_out| [GeV];Events", 100, 0, 1e-5);

    while (getline(file, line)) {
        if (line[0] == '6') {  // Header line, skip
            // Read next 6 particle lines
            TLorentzVector p_in, p_out;
            for (int i = 0; i < 6; i++) {
                getline(file, line);
                istringstream iss(line);
                int idx, charge, status, pid, m1, m2;
                double px, py, pz, E, mass;
                iss >> idx >> charge >> status >> pid >> m1 >> m2 >> px >> py >> pz >> E >> mass;

                TLorentzVector p(px, py, pz, E);
                if (status == 1 && i < 2) p_in += p;   // Initial state (beam e + target p)
                if (status == 1 && i >= 2) p_out += p; // Final state (e', p', π+, π-)
            }

            double delta = (p_in - p_out).M();
            h_delta->Fill(delta);
        }
    }

    TCanvas* c = new TCanvas("c", "c", 800, 600);
    h_delta->Draw();
    c->SaveAs("momentum_conservation.png");

    cout << "Mean deviation: " << h_delta->GetMean() << " GeV" << endl;
    cout << "Max deviation: " << h_delta->GetMaximum() << " GeV" << endl;
}
```

**Expected**: Deviations < 10⁻⁶ GeV (numerical precision)

---

## Step 7: Parameter Sensitivity Study

To understand how parameters affect results, you can run multiple configurations:

### 7.1 Vary Target Polarization

```bash
# Create configs with different target polarizations
for pol in 0.0 0.3 0.5 0.7 1.0; do
  sed "s/value: 0.70/value: $pol/" config/double_spin.yaml > config/target_pol_${pol}.yaml
  ./sidis/src/generator.exe --config config/target_pol_${pol}.yaml \
    --output output/events_tpol_${pol}.lund
done

# Compare A_LL amplitude vs target polarization
# Expected: A_LL amplitude scales linearly with target pol
```

### 7.2 Vary Q² Range

```bash
# Low Q² (soft regime)
./sidis/src/generator.exe --config config/double_spin.yaml \
  --q2-min 1.0 --q2-max 2.0 --output output/events_q2_low.lund

# High Q² (hard regime)
./sidis/src/generator.exe --config config/double_spin.yaml \
  --q2-min 3.0 --q2-max 5.0 --output output/events_q2_high.lund

# Compare cross sections and asymmetries
# Expected: σ falls with increasing Q², asymmetry patterns may change
```

### 7.3 Test W Term Contributions

```bash
# UU only (unpolarized)
./sidis/src/generator.exe --config config/double_spin.yaml \
  --include-terms UU --output output/events_UU_only.lund

# UU + LU (beam only)
./sidis/src/generator.exe --config config/double_spin.yaml \
  --include-terms UU,LU --output output/events_UU_LU.lund

# All terms (UU, LU, UL, LL)
./sidis/src/generator.exe --config config/double_spin.yaml \
  --include-terms UU,LU,UL,LL --output output/events_all_terms.lund

# Compare to isolate individual contributions
```

---

## Summary

This example demonstrates:

✅ **Complete configuration** via YAML file
✅ **Realistic physics case** (double-spin asymmetry)
✅ **Clear workflow** (config → run → output → analysis)
✅ **Validation procedures** (angles, momentum, asymmetries)
✅ **Flexibility** (CLI overrides, parameter scans)

The proposed parameter system provides:
- **User-friendly** configuration without code changes
- **Physics-oriented** parameters matching PDF specification
- **Extensible** design for future models
- **Reproducible** results via fixed seeds
- **Validated** output with comprehensive checks

**Next Steps for Users:**
1. Implement the Config class and parser (see `PARAMETER_DESIGN.md`)
2. Wire parameters into Cross.h, main.cpp, AmplitudeLoader.h
3. Add CLI argument parsing with config file support
4. Run this example to verify implementation
5. Extend to more complex physics cases (GPD models, radiative corrections, etc.)

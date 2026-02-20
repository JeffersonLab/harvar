# Parameter System Quick Start Guide

## What Has Been Delivered

This repository now contains a **complete parameter system design** for the rho electroproduction event generator. All documentation is ready for implementation.

---

## Documentation Files

### 1. **Configuration Template**
📄 `config/double_spin.yaml`
- Complete YAML configuration file
- Realistic physics case (beam-target double spin asymmetry)
- Extensively commented with physics notes
- Ready to use as template

### 2. **Complete Run Example**
📄 `EXAMPLE_RUN.md`
- End-to-end workflow demonstration
- Expected console output
- Sample LUND format output
- Verification and analysis procedures
- Parameter sensitivity studies

### 3. **Existing Documentation**
📄 `IMPLEMENTATION_COMPLETE.md` - Implementation status (Phases 1-5 complete)
📄 `README.md` - User guide for current implementation
📄 `PHASE_2_4_5_COMPLETE.md` - Detailed phase completion report
📄 `rho-v1.pdf` - Physics specification (equations 1-24)

---

## Quick Start

### Current Implementation (No Config File Yet)

**Run with CLI arguments:**
```bash
cd sidis/src
make
./generator.exe 10000 10.6 > events.lund
```

This generates 10,000 events with:
- Hardcoded Q² = 2.0-2.1 GeV²
- Hardcoded xB = 0.3-0.31
- Randomized beam polarization (±1)
- Unpolarized target
- SCHC amplitude model

### Proposed System (After Implementation)

**Run with config file:**
```bash
./generator.exe --config config/double_spin.yaml
```

This will:
- Use all parameters from YAML file
- Generate 50,000 events with specified kinematics
- Enable beam + target polarization
- Include UU, LU, UL, LL W terms
- Output with full validation and statistics

---

## Parameter System Architecture

### Configuration Sections (8 Total)

1. **Run Control**
   - n_events, random_seed, output_file, log_level
   - Controls overall execution

2. **Kinematics**
   - beam_energy, q2_min/max, xb_min/max, t_min/max
   - Defines phase space

3. **Physics Model**
   - cross_section_model, t_slope, dsigma values
   - Selects physics implementation

4. **Polarization**
   - beam, target_longitudinal, target_transverse
   - Enables/configures polarization terms

5. **Amplitudes**
   - source (SCHC/FILE), file path, validation
   - Helicity amplitude loading

6. **Decay**
   - channel, sampling_mode
   - Decay configuration

7. **Weighting**
   - mode (weighted/unweighted), max weight finding
   - Event weighting strategy

8. **Validation**
   - Runtime checks, debug output, tolerances
   - Quality assurance

### Key Design Principles

✅ **Complete**: Covers all PDF specification parameters
✅ **Minimal**: No redundant or unused knobs
✅ **Safe**: Validation and sensible defaults
✅ **Explicit**: Units and frames clearly documented
✅ **Forward-compatible**: Easy to extend

---

## Implementation Roadmap

### Phase 1: Core Config System
**Files to create:**
- `sidis/src/Config.h` - Configuration class and structs
- `sidis/src/ConfigParser.h` - YAML parsing (use yaml-cpp library)

**Key features:**
```cpp
class Config {
public:
    static Config LoadFromFile(const std::string& filename);
    void Validate() const;

    RunConfig run;
    KinematicsConfig kinematics;
    ModelConfig model;
    PolarizationConfig polarization;
    AmplitudesConfig amplitudes;
    DecayConfig decay;
    WeightingConfig weighting;
    ValidationConfig validation;
};
```

### Phase 2: Wire Into Existing Code
**Files to modify:**

1. **`main.cpp`**
   - Add command-line argument parsing (use CLI11 or similar)
   - Load config file
   - Use config values instead of hardcoded constants

2. **`Cross.h`**
   - Add `Configure(const Config& cfg)` method
   - Use config for polarization settings
   - Use config for amplitude loading

3. **`AmplitudeLoader.h`**
   - Accept AmplitudeConfig struct
   - Make t_slope, dsigma values configurable

### Phase 3: CLI Integration
**Command-line interface:**
```bash
./generator.exe [OPTIONS]

Options:
  --config FILE          Configuration file (YAML/JSON)
  --n-events N           Override number of events
  --beam-energy E        Override beam energy [GeV]
  --output FILE          Override output file
  --seed N               Override random seed
  --help                 Show this help message
  --version              Show version information
```

### Phase 4: Validation & Testing
- Unit tests for Config loading/validation
- Integration tests with example configs
- Compare output with current implementation (regression test)

### Phase 5: Documentation
- User manual with all parameters explained
- Example gallery (6+ physics cases)
- Migration guide for existing users

---

## Example Use Cases

### 1. Unpolarized Cross Section Measurement
**Config:** Set all polarizations to 0, include only W_UU term
**Physics:** Baseline cross section, test SCHC model

### 2. Beam-Spin Asymmetry (A_LU)
**Config:** Beam randomize ±1, target off, include UU + LU
**Physics:** Sensitive to Im(u) amplitudes, sin(ϕ) modulation

### 3. Target-Spin Asymmetry (A_UL)
**Config:** Beam off, target longitudinal ±1, include UU + UL
**Physics:** Sensitive to Im(l) amplitudes, sin(ϕ) modulation

### 4. Double-Spin Asymmetry (A_LL)
**Config:** Both beam and target polarized, include all 4 terms
**Physics:** Sensitive to Im(l), sin(ϕ) modulation, GPD E access
**Example:** See `EXAMPLE_RUN.md` for complete workflow

### 5. Transverse Target Asymmetry (A_UT)
**Config:** Beam off, target transverse, include UU + UT
**Physics:** sin(ϕ-ϕ_S) modulation, sensitive to GPD E_T

### 6. External Amplitude Table
**Config:** Load amplitudes from CSV file (fitted to data)
**Physics:** Realistic cross sections for quantitative studies

---

## File Locations Summary

```
simulations/
├── config/
│   └── double_spin.yaml              # Example config file
├── sidis/src/
│   ├── main.cpp                      # Entry point (to be modified)
│   ├── Cross.h                       # Cross section (to be modified)
│   ├── AmplitudeLoader.h             # Amplitude loading (to be modified)
│   ├── RhoEvent.h                    # Event class (complete)
│   ├── AngleCalculator.h             # PDF angles (complete)
│   └── w_kernels.hpp                 # W kernels (complete)
├── EXAMPLE_RUN.md                    # Complete run example
├── PARAMETER_SYSTEM_QUICK_START.md   # This file
├── IMPLEMENTATION_COMPLETE.md        # Implementation status
├── README.md                         # Current user guide
└── rho-v1.pdf                        # Physics specification
```

---

## Current Status

### ✅ Complete (Ready to Use)
- Physics implementation (Phases 1-5)
- Angle calculations (PDF equations 20-24)
- All 6 W kernels (UU, LU, UL, LL, UT, LT)
- Helicity amplitude system (u, l, n, s matrices)
- Full polarization support
- LUND format output

### 📋 Designed (Ready for Implementation)
- **Parameter system architecture**
- **Configuration file format (YAML)**
- **Parameter schema (8 sections)**
- **Example configs (6+ physics cases)**
- **Implementation roadmap**
- **Validation procedures**

### 🚧 To Be Implemented
- Config.h class and parsing
- CLI argument handling
- Integration with existing code
- Unit/integration tests
- Extended documentation

---

## Next Steps for Developer

1. **Review the design**
   - Read `EXAMPLE_RUN.md` for complete workflow
   - Study `config/double_spin.yaml` for parameter structure
   - Review parameter tables in earlier conversation summary

2. **Install dependencies**
   ```bash
   # YAML parser (choose one):
   brew install yaml-cpp           # macOS
   apt-get install libyaml-cpp-dev # Ubuntu

   # CLI parser (optional but recommended):
   # CLI11 is header-only: https://github.com/CLIUtils/CLI11
   ```

3. **Implement Config class**
   - Create `sidis/src/Config.h`
   - Implement loading, validation, defaults
   - Test with `config/double_spin.yaml`

4. **Wire into main.cpp**
   - Add argument parsing
   - Load config file
   - Replace hardcoded values

5. **Wire into Cross.h and AmplitudeLoader.h**
   - Add Configure() methods
   - Use config parameters throughout

6. **Test and validate**
   - Run `EXAMPLE_RUN.md` workflow
   - Verify outputs match expected results
   - Check all validation procedures pass

7. **Extend**
   - Add more physics models (GPD, Regge)
   - Implement unweighted event mode
   - Add Breit-Wigner rho mass distribution
   - Create analysis tools for common observables

---

## Questions?

- **Parameter definitions?** → See parameter tables in conversation history
- **Physics background?** → Read `rho-v1.pdf` equations
- **Current implementation?** → See `IMPLEMENTATION_COMPLETE.md`
- **Run example?** → See `EXAMPLE_RUN.md`
- **Code structure?** → See `README.md` section "Code Structure"

---

## Contact & Contributing

This parameter system design is ready for community review and implementation. Key advantages:

- **No breaking changes**: Existing CLI interface can coexist with config files
- **Backward compatible**: Default config reproduces current behavior
- **Extensible**: Easy to add new models/features
- **Well-documented**: Complete examples and references

**Implementation effort estimate**: ~2-3 days for experienced C++ developer

---

**Status**: 📐 **DESIGN COMPLETE** → Ready for implementation

**Last updated**: 2026-02-20

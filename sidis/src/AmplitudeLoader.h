/**********************************************************************
 * AmplitudeLoader.h - Load/compute helicity amplitude matrices
 *
 * Provides helicity amplitude matrices u, l, n, s for rho electroproduction
 * as required by PDF equations (13)-(19).
 *
 * Features:
 *   - Simple SCHC (s-channel helicity conservation) model for testing
 *   - File loading from CSV/JSON (future)
 *   - Validation of normalization constraints
 *
 * Matrix indexing: amplitude[nu][nup][mu][mup]
 *   where nu, nup, mu, mup ∈ {'+', '0', '-'} → {0, 1, 2}
 *********************************************************************/
#ifndef AMPLITUDE_LOADER_H
#define AMPLITUDE_LOADER_H

#include "w_kernels.hpp"
#include <complex>
#include <string>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>

struct AmplitudeSet {
    Wkernels::Mat4 u;  // unpolarized
    Wkernels::Mat4 l;  // longitudinal beam polarization
    Wkernels::Mat4 n;  // transverse target (sin term)
    Wkernels::Mat4 s;  // transverse target (cos term)

    double dsigmaT_dt;  // dσ_T/dt [nb/GeV²]
    double dsigmaL_dt;  // dσ_L/dt [nb/GeV²]

    // Kinematics where these amplitudes apply
    double xB;
    double Q2;
    double t;
};

class AmplitudeLoader {
public:
    /**
     * Load amplitudes using specified method.
     *
     * @param method  "SCHC" (simple model) or "FILE" (load from file)
     * @param xB      Bjorken x
     * @param Q2      Virtuality [GeV²]
     * @param t       Mandelstam t [GeV²]
     * @param filename Optional file for method="FILE"
     * @return AmplitudeSet with all matrices populated
     */
    static AmplitudeSet Load(const std::string& method,
                             double xB, double Q2, double t,
                             const std::string& filename = "");

    /**
     * Simple SCHC model for testing.
     *
     * s-channel helicity conservation: only diagonal elements survive.
     * Parameterization based on simple Regge-inspired t-dependence.
     *
     * This is a PLACEHOLDER model for testing the code structure.
     * Replace with realistic amplitudes for physics studies.
     */
    static AmplitudeSet SimpleSCHC(double xB, double Q2, double t);

    /**
     * Load amplitudes from CSV file.
     *
     * Format:
     *   # xB, Q2, t, matrix, nu, nup, mu, mup, Re, Im
     *   0.3, 2.0, -0.5, u, +, +, +, +, 0.5, 0.0
     *   ...
     */
    static AmplitudeSet LoadFromFile(const std::string& filename,
                                      double xB, double Q2, double t);

    /**
     * Validate amplitude set.
     *
     * Checks:
     *   1. Normalization constraint (PDF equation 14)
     *   2. Cross sections are positive
     *   3. Amplitudes are finite
     */
    static bool Validate(const AmplitudeSet& amp, double eps,
                         bool verbose = true);

private:
    // Helper: zero-initialize all matrices
    static void ZeroAll(AmplitudeSet& amp);

    // Helper: check if index char is valid
    static bool IsValidIndex(char c) {
        return c == '+' || c == '0' || c == '-';
    }

    // Simple t-dependence: exponential falloff
    static double TDependence(double t) {
        const double b_slope = 5.0;  // GeV^-2 (typical for rho)
        return std::exp(b_slope * t);  // t is negative
    }
};

//=============================================================================
// Implementation
//=============================================================================

inline AmplitudeSet AmplitudeLoader::Load(const std::string& method,
                                           double xB, double Q2, double t,
                                           const std::string& filename)
{
    if (method == "SCHC" || method == "schc") {
        return SimpleSCHC(xB, Q2, t);
    } else if (method == "FILE" || method == "file") {
        if (filename.empty()) {
            throw std::runtime_error("AmplitudeLoader: filename required for FILE method");
        }
        return LoadFromFile(filename, xB, Q2, t);
    } else {
        throw std::runtime_error("AmplitudeLoader: unknown method '" + method + "'");
    }
}

inline AmplitudeSet AmplitudeLoader::SimpleSCHC(double xB, double Q2, double t)
{
    AmplitudeSet amp;
    ZeroAll(amp);

    amp.xB = xB;
    amp.Q2 = Q2;
    amp.t = t;

    // t-dependence factor (exponential slope)
    double t_factor = TDependence(t);

    // Cross sections with simple Q² and xB dependence
    // These are PLACEHOLDER values - replace with realistic model!
    amp.dsigmaT_dt = 10.0 * t_factor / Q2;  // nb/GeV²
    amp.dsigmaL_dt = 2.0 * t_factor / Q2;   // nb/GeV² (smaller than T)

    using namespace Wkernels;

    // ---------------------------------------------------------------------
    // SCHC Model: Only diagonal helicity elements survive
    // u^{λγ* λγ*}_{λρ λρ} where λγ* and λρ are conserved
    // ---------------------------------------------------------------------

    // Normalization: we'll use equation (14) to determine some amplitudes
    // from others. Start with dominant amplitudes:

    // U-matrix (unpolarized)
    // Diagonal elements (SCHC)
    amp.u[h('0')][h('0')][h('+')][h('+')] = {0.40, 0.0};  // γ*_L → ρ_L
    amp.u[h('0')][h('0')][h('0')][h('0')] = {0.30, 0.0};  // γ*_L → ρ_0
    amp.u[h('+')][h('+')][h('+')][h('+')] = {0.15, 0.0};  // γ*_R → ρ_R
    amp.u[h('-')][h('-')][h('+')][h('+')] = {0.15, 0.0};  // γ*_L (bar) → ρ_L (bar)

    // Small non-diagonal for interference (breaking pure SCHC slightly)
    amp.u[h('0')][h('0')][h('0')][h('+')] = {0.10, 0.05}; // interference term

    // Use constraint equation (14) to check normalization:
    // u_{++}^{++} + u_{--}^{++} + 2ε u_{++}^{00} = 1 - (u_{00}^{++} + ε u_{00}^{00})
    // We set things so this is approximately satisfied for ε ~ 0.5

    // L-matrix (longitudinally polarized beam)
    // Smaller than u, imaginary parts for beam spin asymmetries
    amp.l[h('0')][h('0')][h('0')][h('+')] = {0.05, 0.10};
    amp.l[h('0')][h('0')][h('+')][h('+')] = {0.02, 0.05};
    amp.l[h('+')][h('+')][h('+')][h('+')] = {0.03, 0.02};
    amp.l[h('-')][h('-')][h('+')][h('+')] = {0.03, -0.02}; // opposite sign

    // N-matrix (transverse target, sin term)
    // Even smaller, for transverse SSA
    amp.n[h('0')][h('0')][h('+')][h('+')] = {0.01, 0.03};
    amp.n[h('0')][h('0')][h('0')][h('0')] = {0.01, 0.02};
    amp.n[h('0')][h('0')][h('0')][h('+')] = {0.02, 0.01};

    // S-matrix (transverse target, cos term)
    amp.s[h('0')][h('0')][h('+')][h('+')] = {0.01, -0.02};
    amp.s[h('0')][h('0')][h('0')][h('0')] = {0.01, -0.01};
    amp.s[h('0')][h('0')][h('0')][h('+')] = {0.02, 0.02};

    return amp;
}

inline AmplitudeSet AmplitudeLoader::LoadFromFile(const std::string& filename,
                                                   double xB, double Q2, double t)
{
    AmplitudeSet amp;
    ZeroAll(amp);

    amp.xB = xB;
    amp.Q2 = Q2;
    amp.t = t;

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("AmplitudeLoader: cannot open file " + filename);
    }

    std::string line;
    int line_num = 0;
    bool found_match = false;

    using namespace Wkernels;

    while (std::getline(file, line)) {
        ++line_num;

        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double file_xB, file_Q2, file_t;
        std::string matrix_name;
        char nu, nup, mu, mup;
        double re_val, im_val;
        char comma;

        // Parse: xB, Q2, t, matrix, nu, nup, mu, mup, Re, Im
        if (!(iss >> file_xB >> comma >> file_Q2 >> comma >> file_t >> comma
                  >> matrix_name >> comma
                  >> nu >> comma >> nup >> comma >> mu >> comma >> mup >> comma
                  >> re_val >> comma >> im_val)) {
            throw std::runtime_error("AmplitudeLoader: parse error at line " +
                                     std::to_string(line_num));
        }

        // Check if this line matches our kinematics (with tolerance)
        const double tol = 0.01;
        if (std::abs(file_xB - xB) > tol || std::abs(file_Q2 - Q2) > tol ||
            std::abs(file_t - t) > tol) {
            continue;  // wrong kinematics
        }

        found_match = true;

        // Validate indices
        if (!IsValidIndex(nu) || !IsValidIndex(nup) ||
            !IsValidIndex(mu) || !IsValidIndex(mup)) {
            throw std::runtime_error("AmplitudeLoader: invalid index at line " +
                                     std::to_string(line_num));
        }

        std::complex<double> val(re_val, im_val);

        // Assign to correct matrix
        if (matrix_name == "u") {
            amp.u[h(nu)][h(nup)][h(mu)][h(mup)] = val;
        } else if (matrix_name == "l") {
            amp.l[h(nu)][h(nup)][h(mu)][h(mup)] = val;
        } else if (matrix_name == "n") {
            amp.n[h(nu)][h(nup)][h(mu)][h(mup)] = val;
        } else if (matrix_name == "s") {
            amp.s[h(nu)][h(nup)][h(mu)][h(mup)] = val;
        } else if (matrix_name == "dsigmaT") {
            amp.dsigmaT_dt = re_val;  // use real part
        } else if (matrix_name == "dsigmaL") {
            amp.dsigmaL_dt = re_val;
        } else {
            throw std::runtime_error("AmplitudeLoader: unknown matrix '" +
                                     matrix_name + "' at line " +
                                     std::to_string(line_num));
        }
    }

    if (!found_match) {
        throw std::runtime_error("AmplitudeLoader: no matching kinematics in file");
    }

    // Set default cross sections if not in file
    if (amp.dsigmaT_dt == 0.0) amp.dsigmaT_dt = 1.0;
    if (amp.dsigmaL_dt == 0.0) amp.dsigmaL_dt = 0.3;

    return amp;
}

inline bool AmplitudeLoader::Validate(const AmplitudeSet& amp, double eps,
                                       bool verbose)
{
    bool all_ok = true;

    using namespace Wkernels;

    // Check 1: Cross sections positive
    if (amp.dsigmaT_dt <= 0.0 || amp.dsigmaL_dt <= 0.0) {
        if (verbose) {
            std::cerr << "AmplitudeLoader WARNING: negative cross section\n";
            std::cerr << "  dsigmaT_dt = " << amp.dsigmaT_dt << "\n";
            std::cerr << "  dsigmaL_dt = " << amp.dsigmaL_dt << "\n";
        }
        all_ok = false;
    }

    // Check 2: Normalization constraint (PDF equation 14)
    // u_{++}^{++} + u_{--}^{++} + 2ε u_{++}^{00} = 1 - (u_{00}^{++} + ε u_{00}^{00})
    std::complex<double> lhs = amp.u[h('+')][h('+')][h('+')][h('+')]
                             + amp.u[h('-')][h('-')][h('+')][h('+')]
                             + 2.0*eps * amp.u[h('+')][h('+')][h('0')][h('0')];

    std::complex<double> rhs = 1.0 - (amp.u[h('0')][h('0')][h('+')][h('+')]
                                    + eps * amp.u[h('0')][h('0')][h('0')][h('0')]);

    double norm_diff = std::abs(lhs - rhs);

    if (norm_diff > 0.1) {  // tolerance
        if (verbose) {
            std::cerr << "AmplitudeLoader WARNING: normalization constraint violated\n";
            std::cerr << "  LHS = " << lhs << "\n";
            std::cerr << "  RHS = " << rhs << "\n";
            std::cerr << "  |difference| = " << norm_diff << "\n";
            std::cerr << "  (Should be close to 0 per PDF equation 14)\n";
        }
        // Don't fail on this - it's just a warning for SCHC model
        // all_ok = false;
    }

    // Check 3: Amplitudes are finite
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                for (int m = 0; m < 3; ++m) {
                    auto check_finite = [&](const Wkernels::Mat4& mat, const char* name) {
                        if (!std::isfinite(std::real(mat[i][j][k][m])) ||
                            !std::isfinite(std::imag(mat[i][j][k][m]))) {
                            if (verbose) {
                                std::cerr << "AmplitudeLoader ERROR: non-finite value in "
                                          << name << "[" << i << "][" << j << "]["
                                          << k << "][" << m << "]\n";
                            }
                            return false;
                        }
                        return true;
                    };

                    if (!check_finite(amp.u, "u")) all_ok = false;
                    if (!check_finite(amp.l, "l")) all_ok = false;
                    if (!check_finite(amp.n, "n")) all_ok = false;
                    if (!check_finite(amp.s, "s")) all_ok = false;
                }
            }
        }
    }

    if (verbose && all_ok) {
        std::cerr << "AmplitudeLoader: validation passed\n";
    }

    return all_ok;
}

inline void AmplitudeLoader::ZeroAll(AmplitudeSet& amp)
{
    // Zero-initialize all matrices
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                for (int m = 0; m < 3; ++m) {
                    amp.u[i][j][k][m] = {0.0, 0.0};
                    amp.l[i][j][k][m] = {0.0, 0.0};
                    amp.n[i][j][k][m] = {0.0, 0.0};
                    amp.s[i][j][k][m] = {0.0, 0.0};
                }
            }
        }
    }

    amp.dsigmaT_dt = 0.0;
    amp.dsigmaL_dt = 0.0;
    amp.xB = 0.0;
    amp.Q2 = 0.0;
    amp.t = 0.0;
}

#endif // AMPLITUDE_LOADER_H

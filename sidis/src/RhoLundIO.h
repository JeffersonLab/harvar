// RhoLundIO.h ─────────────────────────────────────────────────────────
// Helper for printing a RhoEvent in LUND format (CLAS convention)

#ifndef RHO_LUND_IO_H
#define RHO_LUND_IO_H

#include "RhoEvent.h"
#include <iomanip>
#include <ostream>
#include <vector>

#include <iomanip>
#include <vector>

/// Write one event in CLAS12 LUND format.
/// @param ev        the filled RhoEvent
/// @param eventNo   running event number (starts at 1)
/// @param weight    event weight (1 = unweighted)
/// @param vtx       production vertex [cm] – defaults to (0,0,0)
/// @param out       output stream   (std::cout by default)
inline void printLundCLAS12(const RhoEvent&   ev,
                            int               eventNo,
                            double            weight = 1.0,
                            const TVector3&   vtx    = TVector3(0,0,0),
                            std::ostream&     out    = std::cout)
{
    using P4 = const TLorentzVector&;

    struct Rec { int q; int status; int pid; P4 p4; double m; };
    std::vector<Rec> particles = {
        { -1, 0, 11, ev.e_in,        0.000511 },   // beam e⁻  (optional, tag status=0 if desired)
        { -1, 1, 11, ev.e_out,   0.000511 },   // scattered e⁻
        { +1, 1, 2212, ev.p_out,  0.938272 },   // recoil p
        {  0, 0, 113, ev.v,       ev.v.M()  },  // ρ
        { +1, 1, 211, ev.piPlus,      0.13957  },   // π⁺
        { -1, 1, -211, ev.piMinus,      0.13957  }    // π⁻
    };

    /* ---------------- header (10 numbers) ---------------- */
    out << std::fixed << std::setprecision(6);
    out << particles.size()             << ' '
        << 1                       << ' '
        << 1                        << ' '
        << ev.pol             << ' '
        << 0             << ' '
        << 11             << ' '
        << ev.beamE             << ' '
        << 2212             << ' '
        << 42             << ' '
        << 0.0             << ' '
        << ev.Q2                         << ' '
        << ev.xB                         << ' '
        << ev.theta                         << ' '
        << ev.phi                        << ' '
        << ev.thetaPi                         << ' '
        << ev.phiPi             << ' '   // ν = y·Ebeam

        << eventNo                      << '\n';

    /* -------------- particle records (15 cols) ----------- */
    int idx = 1;
    for (const auto& r : particles)
    {
        int status = 1;            // 1 = final state
        int m1 = 0, m2 = 0;        // mother indices not used here

        out << std::setw(3)  << idx           << ' '
            << std::setw(3)  << r.q           << ' '
            << std::setw(3)  << r.status        << ' '
            << std::setw(6)  << r.pid         << ' '
            << std::setw(3)  << m1            << ' '
            << std::setw(3)  << m2            << ' '
            << std::setw(12) << r.p4.Px()     << ' '
            << std::setw(12) << r.p4.Py()     << ' '
            << std::setw(12) << r.p4.Pz()     << ' '
            << std::setw(12) << r.p4.E()      << ' '
            << std::setw(10) << r.m           << ' '
            << std::setw(9)  << vtx.X()       << ' '
            << std::setw(9)  << vtx.Y()       << ' '
            << std::setw(9)  << vtx.Z()       << ' '
            << std::setw(9)  << 0.0           << '\n';   // production time
        ++idx;
    }
}



#endif /* RHO_LUND_IO_H */

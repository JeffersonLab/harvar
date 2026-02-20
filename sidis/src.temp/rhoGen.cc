//  rhoGen.cc  --  simple ρ0 electro‑production Monte‑Carlo
// 
//  needs:  ROOT >= 6  (TLorentzVector, TRandom3)
// compile with - 
// g++ -O2 rhoGen.cc `root-config --cflags --glibs` -o rhoGen
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <iostream>
#include <iomanip>
#include <cmath>

constexpr double Mp   = 0.938272;     // proton   mass [GeV]
constexpr double Mrho = 0.77526;      // rho0     mass [GeV]
constexpr double Mpi  = 0.13957;      // charged‑pion mass [GeV]
constexpr double me   = 0.000511;     // electron mass [GeV]

struct Event {
    TLorentzVector l, lp;             // e beam & scattered e′
    TLorentzVector p, pp;             // p beam & recoil p′
    TLorentzVector q;                 // virtual photon
    TLorentzVector v;                 // ρ0
    TLorentzVector k1, k2;            // π+  π–
};

/// Helper: isotropic direction
TVector3 isotropicDir(TRandom3& R)  {
    const double cosTh = R.Uniform(-1., 1.);
    const double sinTh = std::sqrt(1. - cosTh*cosTh);
    const double phi   = R.Uniform(0., 2.*M_PI);
    return TVector3( sinTh*std::cos(phi),
                     sinTh*std::sin(phi),
                     cosTh );
}

/// One Monte‑Carlo event -----------------------------------------------------
Event generate(TRandom3& R,
               double Ebeam = 11.0,        // beam energy [GeV]
               double  xmin = 0.05, double  xmax = 0.7,
               double Q2min = 1.0,  double Q2max = 5.0,
               double  tmin =-1.0,  double  tmax =-0.05 )
{
    //-----------------------------------------------------------------------
    // 1)  pick inclusive kinematics  (flat → replace by PDF, structure fn…)
    //-----------------------------------------------------------------------
    double xB = R.Uniform(xmin, xmax);
    double Q2 = R.Uniform(Q2min, Q2max);                // [GeV²]
    double nu = Q2 / (2.0*Mp*xB);                       // γ* energy
    double Ep = Ebeam;                                  // shorthand

    // scattered‑electron energy  E′ = E − ν   (Eq.(10))
    double EpPrime = Ep - nu;
    if (EpPrime <= me) return generate(R,Ebeam,xmin,xmax,Q2min,Q2max,tmin,tmax);

    // θ_e from  Q² = 4 E E′ sin²(θ/2)
    double sinHalf = std::sqrt(Q2 /(4.0*Ep*EpPrime));
    double theta_e = 2.0*std::asin( std::min(1.0, sinHalf) );
    double phi_e   = R.Uniform(0., 2.*M_PI);

    //-----------------------------------------------------------------------
    // 2)  build lepton 4‑vectors  (lab z‑axis = beam)
    //-----------------------------------------------------------------------
    double plep = std::sqrt(Ep*Ep - me*me);
    TLorentzVector l(0, 0,  plep,  Ep);

    double plp  = std::sqrt(EpPrime*EpPrime - me*me);
    TLorentzVector lp( plp*std::sin(theta_e)*std::cos(phi_e),
                       plp*std::sin(theta_e)*std::sin(phi_e),
                       plp*std::cos(theta_e), EpPrime );

    TLorentzVector q = l - lp;               // γ*  (Eq.(4))

    //-----------------------------------------------------------------------
    // 3)  γ* p  →  ρ p′   in γ*p c.m.                                       
    //-----------------------------------------------------------------------
    TLorentzVector p(0,0,0,Mp);              // target at rest

    double W2 = Mp*Mp + 2.*Mp*nu - Q2;       // Eq.(9)
    double W  = std::sqrt(W2);
    if (W < Mp + Mrho) return generate(R,Ebeam,xmin,xmax,Q2min,Q2max,tmin,tmax);

    // two‑body c.m. momentum
    double pcm = std::sqrt( (W2 - std::pow(Mp+Mrho,2.)) *
                            (W2 - std::pow(Mp-Mrho,2.)) ) / (2.*W);

    // choose production angles of ρ in c.m.
    TVector3 dir_cm = isotropicDir(R);
    TLorentzVector v_cm( pcm*dir_cm, std::sqrt(pcm*pcm + Mrho*Mrho) );
    TLorentzVector pp_cm( -pcm*dir_cm, std::sqrt(pcm*pcm + Mp*Mp) );

    // boost γ*p c.m. → lab   (β⃗ = q⃗ /(ν+Mp) )
    //TVector3 beta_cm = q.Vect() / ( q.E() + Mp );
    // --- replacement: scale with 1/(Eγ*+Mp) -----------
    TVector3 beta_cm = q.Vect();                // copy the 3‑vector
    beta_cm *= 1.0 / (q.E() + Mp);              // scale in place
    v_cm.Boost(beta_cm);
    pp_cm.Boost(beta_cm);

    //-----------------------------------------------------------------------
    // 4)  enforce t  = (p−p′)²  (Eq.(6))  by rejecting until match
    //-----------------------------------------------------------------------
    double t   = (p - pp_cm).Mag2();
    if ( t > tmax || t < tmin )
        return generate(R,Ebeam,xmin,xmax,Q2min,Q2max,tmin,tmax);

    //-----------------------------------------------------------------------
    // 5)  decay ρ0 → π+π−  isotropically in its rest frame  (Eq.(17))
    //-----------------------------------------------------------------------
    TVector3 beta_rho = v_cm.BoostVector();

    double ppi = std::sqrt( std::pow(Mrho*Mrho/4. - Mpi*Mpi, 1.0) );
    TVector3 k1_dir = isotropicDir(R);
    TLorentzVector k1_rf(  ppi*k1_dir,  std::sqrt(ppi*ppi + Mpi*Mpi) );
    TLorentzVector k2_rf( -ppi*k1_dir,  std::sqrt(ppi*ppi + Mpi*Mpi) );

    TLorentzVector k1 = k1_rf;  k1.Boost(beta_rho);
    TLorentzVector k2 = k2_rf;  k2.Boost(beta_rho);

    // ---------------------------------------------------------------------
    Event ev;
    ev.l   = l;      ev.lp  = lp;
    ev.p   = p;      ev.pp  = pp_cm;
    ev.q   = q;      ev.v   = v_cm;
    ev.k1  = k1;     ev.k2  = k2;
    return ev;
}
//------------------------------------------------------------------
// ρ‑event → stdout in LUND format  (see GEMC spec) :contentReference[oaicite:0]{index=0}
//
//   Header (10 numbers) : Npart  A  Z  Tpol  Bpol  BeamType  Ebeam
//                         NucleonID  ProcID  Weight
//   One line per particle: idx  τ  type  PDG  parent  firstD
//                          px  py  pz  E  m  vx  vy  vz
//------------------------------------------------------------------
#include <fstream>
#include <iomanip>

void write_lund(const Event& ev, std::ostream& os = std::cout,
                int A=1, int Z=1, double beamPol=0.0, double targetPol=0.0,
                int procID=1, double weight=1.0, double Ebeam=11.0)
{
    constexpr int BeamType    = 11;   // electron
    constexpr int NucleonID   = 2212; // struck proton
    const int     Npart       = 4;    // e'  p'  π+  π‑

    //--- HEADER ----------------------------------------------------
    os << std::fixed << std::setprecision(6);
    os << Npart << " " << A << " " << Z << " "
       << targetPol << " " << beamPol << " "
       << BeamType << " " << Ebeam << " "
       << NucleonID << " " << procID << " " << weight << "\n";

    // helper lambda: print one particle line
    auto pline = [&os](int idx, const TLorentzVector& p4,
                       int pdg, double mass, double tau = -1.,
                       int type=1, int parent=0, int firstD=0)
    {
        os << idx << " " << tau << " " << type << " "
           << pdg << " " << parent << " " << firstD << " "
           << p4.Px() << " " << p4.Py() << " " << p4.Pz() << " "
           << p4.E()  << " " << mass       << " "
           << 0.0 << " " << 0.0 << " " << 0.0 << "\n";   // vertex (cm)
    };

    //--- PARTICLE LINES -------------------------------------------
    pline(1, ev.lp ,  11 , me   );        // scattered e⁻
    pline(2, ev.pp , 2212, Mp   );        // recoil p
    pline(3, ev.k1 ,  211, Mpi  );        // π⁺
    pline(4, ev.k2 , -211, Mpi  );        // π⁻
}

/// ------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    TRandom3 R(0);                      // seed with /dev/urandom
    const int Nevt = 200000;                // quick demo

    std::cout << std::fixed << std::setprecision(6);

    for (int i = 0; i < Nevt; ++i) {
        Event ev = generate(R);
        write_lund(ev);          // → stdout
    }
    /*
    for (int i=0;i<Nevt;++i) {
        Event ev = generate(R);

        std::cout << "EVENT " << i+1 << "\n";
        auto pr = [](const char* n, const TLorentzVector& v){
            std::cout << std::setw(4) << n << " "
                      << v.Px() << "  " << v.Py() << "  "
                      << v.Pz() << "  " << v.E()  << "\n";
        };
        pr("e  ", ev.l);   pr("e' ", ev.lp);
        pr("p  ", ev.p);   pr("p' ", ev.pp);
        pr("rho", ev.v);   pr("pi+", ev.k1);  pr("pi-", ev.k2);
        std::cout << std::string(60,'-') << "\n";
    }*/
}

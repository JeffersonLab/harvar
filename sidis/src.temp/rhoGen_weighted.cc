/**************************************************************
 *  rhoGen_weighted.cc  –  ρ0 electro‑production generator
 *                       distributed according to d⁷σ
 *  -----------------------------------------------------------
 *  Compile:  g++ -O2 -std=c++17 rhoGen_weighted.cc `root-config --cflags --glibs`  -o rhoGen_weighted
 *
 *  Output:   LUND records to stdout (one event per line)
 *************************************************************/
 #include <TVector3.h>
 #include <TLorentzVector.h>
 #include <TRandom3.h>
 #include <TMath.h>
 
 #include <complex>
 #include <array>
 #include <iostream>
 #include <iomanip>
 #include <cmath>
 #include <limits>
 
 #include "w_kernels.hpp"           // gives Wkernels::UU…()
 
 //--------------------------------------------------------------------
 // 1.  Constants
 //--------------------------------------------------------------------
 constexpr double Mp   = 0.9382720813;   // proton mass   [GeV]
 constexpr double Mrho = 0.77526;        // ρ⁰ mass       [GeV]
 constexpr double Mpi  = 0.13957;        // π^{±} mass    [GeV]
 constexpr double me   = 0.000511;       // e⁻  mass      [GeV]
 constexpr double αem  = 1.0 / 137.035999084;
 
 //--------------------------------------------------------------------
 // 2.  Utility structures
 //--------------------------------------------------------------------
 struct Event {
     TLorentzVector  l, lp;          // e beam & scattered e′
     TLorentzVector  p, pp;          // p beam & recoil p′
     TLorentzVector  q;              // virtual photon
     TLorentzVector  v;              // ρ⁰
     TLorentzVector  k1, k2;         // π⁺ π⁻
 };
 
 using Mat4 = std::array<std::array<std::array<std::array<
              std::complex<double>,3>,3>,3>,3>;
 
 //--------------------------------------------------------------------
 // 3.  Physics helpers copied from cross_section.cc
 //--------------------------------------------------------------------
 inline double gamma_bj(double xB, double Q2)
 { return 2.0 * xB * Mp / std::sqrt(Q2); }
 
 inline double epsilon(double y, double γ)
 {
     double y2γ2 = y*y*γ*γ;
     return (1.0 - y - 0.25*y2γ2) /
            (1.0 - y + 0.5*y*y + 0.25*y2γ2);
 }
 
 inline double dsigma_3fold(double xB, double Q2, double y,
                            double eps, double dσT_dt, double dσL_dt)
 {
     double pref = αem/(2.0*M_PI);
     double kin  = (y*y)/(1.0-eps)*(1.0-xB)/xB/Q2;
     return pref * kin * (dσT_dt + eps*dσL_dt);
 }
 
 inline double dsigma_7fold(double WUU,double WLU,double WUL,
                            double WLL,double WUT,double WLT,
                            double Pl,double SL,double ST,
                            double σ3)
 {
     double S = WUU + Pl*WLU + SL*WUL + Pl*SL*WLL
                      + ST*WUT + Pl*ST*WLT;
     return σ3 * S / (4.0*M_PI*M_PI);
 }
 
 //--------------------------------------------------------------------
 // 4.  Kinematic utilities
 //--------------------------------------------------------------------
 TVector3 isotropicDir(TRandom3& R)
 {
     const double cosθ = R.Uniform(-1.,1.);
     const double sinθ = std::sqrt(1.-cosθ*cosθ);
     const double φ    = R.Uniform(0.,2.*M_PI);
     return { sinθ*std::cos(φ), sinθ*std::sin(φ), cosθ };
 }
 
 /// Trento‑convention azimuth between two planes (n₁ × n₂ orientation)
 double trentoPhi(const TVector3& kIn , const TVector3& kOut,
                  const TVector3& q   , const TVector3& pPr)
 {
     TVector3 nL = kIn.Cross(kOut).Unit();    // leptonic
     TVector3 nH = q.Cross(pPr ).Unit();      // hadronic
     double   φ  = std::acos( nL.Dot(nH) );
     // orientation (sign of sinφ):
     TVector3 z  = q.Unit();
     if ( z.Dot( nL.Cross(nH) ) < 0.0 ) φ = 2.*M_PI - φ;
     return φ;
 }
 
 /// Decay angles of π⁺ in the ρ rest frame.
 /// Returns (cosθ_h, φ_h) following Schilling‑Wolf conventions.
 std::pair<double,double> decayAngles(const TLorentzVector& k1,
                                      const TLorentzVector& v,
                                      const TLorentzVector& q)
 {
     // Boost to ρ rest frame
     TLorentzVector k1_rf = k1;
     k1_rf.Boost( -v.BoostVector() );
 
     // Define z′ along –q* in ρ rest frame
     TLorentzVector q_rf = q;
     q_rf.Boost( -v.BoostVector() );
     TVector3 zPrime = (-q_rf.Vect()).Unit();
 
     // y′ perpendicular to production plane
     TVector3 yPrime = v.Vect().Cross(q.Vect()).Unit();
 
     // x′ = y′ × z′
     TVector3 xPrime = yPrime.Cross(zPrime);
 
     TVector3 kVec   = k1_rf.Vect().Unit();
     double cosθh    =  kVec.Dot(zPrime);
     double φh       =  std::atan2( kVec.Dot(yPrime),
                                    kVec.Dot(xPrime) );
     if(φh<0) φh += 2.*M_PI;
     return {cosθh, φh};
 }
 
 //--------------------------------------------------------------------
 // 5.  Physics input (fill this for your model)
 //--------------------------------------------------------------------
 struct PhysicsInput {
     Mat4 u, l, s;
     double dσT_dt, dσL_dt;   // [nb / GeV²]
     double Pl = 0.0, SL = 0.0, ST = 0.0;   // polarisations
 };
 
 PhysicsInput loadPhysicsInputs(double /*xB*/, double /*Q2*/,
                                double t /*GeV²*/)
 {
     PhysicsInput ph{};
 
     /* ---------------------------------------------------------------
      *  TODO:
      *  1.  Set the u,l,s helicity matrices for the chosen (xB,Q²,t)
      *      or read them from a file.
      *  2.  Provide dσT/dt and dσL/dt (nb/GeV²).
      * --------------------------------------------------------------*/
     ph.dσT_dt = 1.0;   // placeholder
     ph.dσL_dt = 1.0;   // placeholder
     ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')][Wkernels::h('+')] = {0.0123, 0.0000};   // real example
     //ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')] = {0.0150, 0.0240};   // real example
     ph.l[Wkernels::h('+')][Wkernels::h('-')][Wkernels::h('0')][Wkernels::h('+')] = {0.0,    -0.0071};  // imaginary
     return ph;
 }
 
 //--------------------------------------------------------------------
 // 6.  One candidate event (+ its weight)
 //--------------------------------------------------------------------
 struct Candidate {
     Event  ev;
     double weight;
 };
 
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
 Candidate generateCandidate(TRandom3& R,
                             double Ebeam = 11.0,
                             double xmin  = 0.05, double xmax  = 0.7,
                             double Q2min = 1.0 , double Q2max = 5.0,
                             double tmin  =-1.0 , double tmax  =-0.05 )
 {
     /* --------  step 1: inclusive kinematics  ------------------- */
     double xB = R.Uniform(xmin,xmax);
     double Q2 = R.Uniform(Q2min,Q2max);
     double ν  = Q2 / (2.*Mp*xB);
     double Ee = Ebeam;
     double EeP= Ee - ν;
     if(EeP<=me) return generateCandidate(R,Ebeam,xmin,xmax,Q2min,Q2max,tmin,tmax);
 
     double sinHalf = std::sqrt( Q2 /(4.*Ee*EeP) );
     double θe      = 2.*std::asin( std::min(1.0, sinHalf) );
     double φe      = R.Uniform(0.,2.*M_PI);
 
     /* --------  step 2: lepton 4‑vectors ------------------------ */
     double plep = std::sqrt(Ee*Ee - me*me);
     TLorentzVector l ( 0, 0,  plep, Ee );
     double plp  = std::sqrt(EeP*EeP - me*me);
     TLorentzVector lp( plp*std::sin(θe)*std::cos(φe),
                        plp*std::sin(θe)*std::sin(φe),
                        plp*std::cos(θe), EeP );
     TLorentzVector q = l - lp;
 
     /* --------  step 3: γ* p → ρ p′  ---------------------------- */
     double W2 = Mp*Mp + 2.*Mp*ν - Q2;
     double W  = std::sqrt(W2);
     if (W < Mp + Mrho)
         return generateCandidate(R,Ebeam,xmin,xmax,Q2min,Q2max,tmin,tmax);
 
     double pcm = std::sqrt( (W2 - std::pow(Mp+Mrho,2.)) *
                             (W2 - std::pow(Mp-Mrho,2.)) ) / (2.*W);
     TVector3 dir_cm = isotropicDir(R);
     TLorentzVector v_cm ( pcm*dir_cm, std::sqrt(pcm*pcm+Mrho*Mrho) );
     TLorentzVector pp_cm( -pcm*dir_cm, std::sqrt(pcm*pcm+Mp*Mp) );
 
     // boost c.m. → lab
     TVector3 β_cm = q.Vect() * (1.0/(q.E()+Mp));
     v_cm .Boost(β_cm);
     pp_cm.Boost(β_cm);
 
     double t = (TLorentzVector(0,0,0,Mp) - pp_cm).Mag2();
     if(t>tmax || t<tmin)
         return generateCandidate(R,Ebeam,xmin,xmax,Q2min,Q2max,tmin,tmax);
 
     /* --------  step 4: ρ → π⁺π⁻ ------------------------------- */
     TVector3 βρ = v_cm.BoostVector();
     double pπ   = std::sqrt( Mrho*Mrho/4. - Mpi*Mpi );
     TVector3 kDir = isotropicDir(R);
     TLorentzVector k1_rf(  pπ*kDir,  std::sqrt(pπ*pπ + Mpi*Mpi) );
     TLorentzVector k2_rf( -pπ*kDir,  std::sqrt(pπ*pπ + Mpi*Mpi) );
     TLorentzVector k1 = k1_rf;  k1.Boost(βρ);
     TLorentzVector k2 = k2_rf;  k2.Boost(βρ);
 
     /* --------  step 5: physics weight -------------------------- */
     double y     = ν/Ee;
     double γ     = gamma_bj(xB,Q2);
     double eps   = epsilon(y,γ);
 
     PhysicsInput ph = loadPhysicsInputs(xB,Q2,t);
     
     auto WUU = Wkernels::UU(ph.u,eps, trentoPhi(l.Vect(),lp.Vect(),
                                                 q.Vect(),pp_cm.Vect()),
                                       decayAngles(k1,v_cm,q).second);
    double WLU = 0.0, WUL = 0.0, WLL = 0.0, WUT = 0.0, WLT = 0.0;
     //auto WLU = Wkernels::LU(ph.u,eps, ... );   // fill as needed
     /* If you only need the unpolarised cross section you can skip
        the other  five kernels and set them to zero.               */
 
     double σ3 = dsigma_3fold(xB,Q2,y,eps,ph.dσT_dt,ph.dσL_dt);
     //std::cout << " WUU  " << WUU.LL << std::endl;
     double cos2theta = std::cos(θe)*std::cos(θe);
     double W_UU = WUU.LL*cos2theta;
     std::cout << " WUU  " << W_UU << "   " << W_UU*3./(4*3.14) << " " << WUU.LL << std::endl;
     double w  = dsigma_7fold(W_UU,0,0,0,0,0,
                              ph.Pl,ph.SL,ph.ST, σ3);
     //std::cout << xB << "  " << Q2 << "  " << t << "  " <<  w << std::endl;
     /* --------  pack / return ---------------------------------- */
     Candidate cand;
     cand.ev = { l, lp, TLorentzVector(0,0,0,Mp), pp_cm,
                 q, v_cm, k1, k2 };
     cand.weight = w;
     return cand;
 }
 
 //--------------------------------------------------------------------
 // 7.  Main program – un‑weighting by accept/reject
 //--------------------------------------------------------------------
 int main(int argc,char* argv[])
 {
     const long     Nrequest = (argc>1) ? std::stol(argv[1]) : 1000000;
     const double   Ebeam    = (argc>2) ? std::stod(argv[2]) : 10.6;
 
     TRandom3 R(0);
 
     /* -- quick pre‑scan to get a safe upper bound for w ---------- */
     const int Nscan = 50000;
     double wMax = 0.0;
     for(int i=0;i<Nscan;++i) {
         //double w = generateCandidate(R,Ebeam).weight;
         double w = generateCandidate(R,Ebeam,0.299,0.301,2.49,2.51,-0.51,-0.49).weight;
         if(w>wMax) wMax = w;
     }
     wMax *= 1.2;                                // safety margin
 
     std::cout << std::fixed << std::setprecision(6);
     long Nkept=0, Nattempt=0;
     while(Nkept < Nrequest) {
        auto cand = generateCandidate(R,Ebeam,0.299,0.301,2.49,2.51,-0.51,-0.49);
         //auto cand = generateCandidate(R,Ebeam);
         ++Nattempt;
         if( R.Uniform(0.,wMax) <= cand.weight ) {
             ++Nkept;
             const auto& e = cand.ev;
             write_lund(e);
             /* -------  write LUND  --------------------------------- *
              * You can of course switch to ROOT/hipo here.           */
             /*std::cout << 6 << ' ' << 1 << '\n'   // #part, event# (dummy)
                       << " 1  -1  "  << e.l.Px()  << ' ' << e.l.Py()  << ' ' << e.l.Pz()  << ' ' << e.l.E()  << '\n'
                       << " 1  -1  "  << e.lp.Px() << ' ' << e.lp.Py() << ' ' << e.lp.Pz() << ' ' << e.lp.E() << '\n'
                       << " 1   1  "  << e.k1.Px() << ' ' << e.k1.Py() << ' ' << e.k1.Pz() << ' ' << e.k1.E() << '\n'
                       << " 1  -1  "  << e.k2.Px() << ' ' << e.k2.Py() << ' ' << e.k2.Pz() << ' ' << e.k2.E() << '\n'
                       << " 1 2212  "<< e.pp.Px() << ' ' << e.pp.Py() << ' ' << e.pp.Pz() << ' ' << e.pp.E() << '\n'
                       << " 1 2212  0 0 0 " << Mp << '\n';  // target
                       */
         }
     }
     std::cerr << "Generated " << Nkept << " events in "
               << Nattempt << " attempts  (ε = "
               << 100.*Nkept/Nattempt << " %)\n";
     return 0;
 }
 
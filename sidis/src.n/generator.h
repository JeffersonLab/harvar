#include "reaction.h"
#include "w_kernels.hpp"
#include "cross.h"

#ifndef __GENERATOR__
#define __GENERATOR__

namespace sim {
  
  struct physicsInput {
    Wkernels::Mat4 u, l, s;
    double dsigmaT_dt, dsigmaL_dt;   // [nb / GeV²]
    double Pl = 0.0, SL = 0.0, ST = 0.0;   // polarisations
  };
  
class generator {
 private:

  static constexpr double Mp = 0.93827;
  static constexpr double alem  = 1.0 / 137.035999084;
  candidate *cand;
  TRandom3   rand;
  
  double q2_min, q2_max;
  double xb_min, xb_max;
  double E;

  int generated, accepted;
  
 public:
  generator(candidate *c){
    cand = c; E = cand->react.E();
    generated = 0; accepted = 0;
  }
  
  void setRange(double __q2min, double __q2max, double __xbmin, double __xbmax){
    q2_min = __q2min; q2_max = __q2max;
    xb_min = __xbmin; xb_max = __xbmax;
  }

  int  is_valid(double q2, double xb){
    double upper = 2.0*Mp*E*xb/(1+Mp*xb/E);
    double lower = 2.73*xb/(1.0-xb);
    if(q2>1.0&&q2>lower&&q2<upper) return 1;
      return 0;
  }
  
  void generate(){
    double q2 = rand.Uniform(q2_min,q2_max);
    double xb = rand.Uniform(xb_min,xb_max);
    int counter = 0;
    while(is_valid(q2,xb)==0){
      q2 = rand.Uniform(q2_min,q2_max);
      xb = rand.Uniform(xb_min,xb_max);
      counter++;
    }
    cand->react.generate(q2,xb);
    
    //printf("random (Q2,xb) : %4d %8.5f %8.5f\n",counter,q2,xb);
  }
  
  double dsigma_3fold(double xB, double Q2, double y,
		      double eps, double dsigmaT_dt, double dsigmaL_dt)
  {
    double pref = alem/(2.0*M_PI);
    double kin  = (y*y)/(1.0-eps)*(1.0-xB)/xB/Q2;
    return pref * kin * (dsigmaT_dt + eps*dsigmaL_dt);
  }

  double dsigma_7fold(double WUU,double WLU,double WUL,
		      double WLL,double WUT,double WLT,
		      double Pl,double SL,double ST,
		      double sigma3)
  {
    double S = WUU + Pl*WLU + SL*WUL + Pl*SL*WLL
      + ST*WUT + Pl*ST*WLT;
    return sigma3 * S / (4.0*M_PI*M_PI);
  }

  double scan(int count){
    double maxWeight = 0.0;
    for(int i = 0; i < count; i++){
      generate();
      updateWeight();
      //printf("weight = %f\n",cand->weight);
      if(cand->weight>maxWeight) maxWeight = cand->weight;
    }
    //printf("max weight = %f\n",maxWeight);
    return maxWeight*1.2;
  }

  void generate(double maxWeight){
    int misses = 0;
    generate();
    updateWeight();
    generated++;
    double accept = rand.Uniform(0,maxWeight);
    //printf("try 0 : %f %f \n",accept, cand->weight);
    while(accept>cand->weight){
      misses++;
      //printf("missed = %d : %f %f \n",misses,accept, cand->weight);
      generate();
      updateWeight();
      accept = rand.Uniform(0,maxWeight);
      generated++;
    }
    accepted++;
    //printf("misses = %d\n",misses);
  }
  void stats(){
    printf("\ngenerator %d %d ~ %f\n\n",generated, accepted,
	   ((double) accepted)/generated  );
  }

  double scan2(int count){
    double maxWeight = 0.0;
    for(int i = 0; i < count; i++){
      generate();
      updateWeight2();
      if(cand->weight>maxWeight) maxWeight = cand->weight;
    }
    return maxWeight*1.2;
  }
  
  void updateWeight2(){
    
    cross c;
    double t = (cand->react.p_out-cand->react.p_in).m2();
    double w = c.weight(cand->react.Q2(),cand->react.xB(),t,
			cand->react.prodTheta, cand->react.prodPhi,
			cand->react.decayTheta,
			cand->react.decayPhi,1.0);
    cand->weight = w;
  }
  
  void updateWeight(){
     double t = (cand->react.p_out-cand->react.p_in).m2();
     physicsInput ph = loadPhysicsInputs(cand->react.xB(),cand->react.Q2(),t);

     double  Q2 = cand->react.Q2();
     double  xB = cand->react.xB();
     double   y = cand->react.Y();
     double eps = cand->react.Eps();
     double sigma3 = dsigma_3fold(xB,Q2,y,eps,ph.dsigmaT_dt,ph.dsigmaL_dt);
     
     auto WUU = Wkernels::UU(ph.u,eps, cand->react.prodPhi ,
			     cand->react.decayPhi);
     
     auto WLU = Wkernels::LU(ph.u,eps, cand->react.prodPhi,
			     cand->react.decayPhi);

     double cosTheta = std::cos(cand->react.decayTheta);
     double sinTheta = std::sin(cand->react.decayTheta);

     double W_LU = cosTheta*cosTheta*WLU.LL+std::sqrt(2)*cosTheta*sinTheta*WLU.LT+sinTheta*sinTheta*WLU.TT;
     double W_UU = cosTheta*cosTheta*WUU.LL+std::sqrt(2)*cosTheta*sinTheta*WUU.LT+sinTheta*sinTheta*WUU.TT;

     double pol = cand->react.rpol<0?-1:1;
     
     //printf("-t = %f %f %f %f\n",t,y,eps,sigma3);
     double w  = dsigma_7fold(W_UU,W_LU,0,0,0,0,
			     pol,ph.SL,ph.ST, sigma3);
     //printf("sigmas = %f %f %f %f\n",W_LU,W_UU,sigma3, w);
     cand->weight = w;
  }

  
  physicsInput loadPhysicsInputs(double xB, double Q2,
  				 double t /*GeV²*/)
  {
    physicsInput ph{};
    
    /* ---------------------------------------------------------------
     *  TODO:
     *  1.  Set the u,l,s helicity matrices for the chosen (xB,Q²,t)
     *      or read them from a file.
     *  2.  Provide dσT/dt and dσL/dt (nb/GeV²).
     * --------------------------------------------------------------*/
    ph.dsigmaT_dt = 1.0;   // placeholder
    ph.dsigmaL_dt = 1.0;   // placeholder
    
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')][Wkernels::h('+')] = {1.0*exp(2*t), 0.0000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')] = {0.00, 0.0000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')] = {0.00, 0.0000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('-')][Wkernels::h('+')] = {0.00, 0.0000};
    
    return ph;
  }
};
}
#endif

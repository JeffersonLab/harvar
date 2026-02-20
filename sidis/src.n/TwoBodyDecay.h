#ifndef TWO_BODY_DECAY_H
#define TWO_BODY_DECAY_H

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include <cmath>

class TwoBodyDecay {
public:
    // Parent four-vector and daughter masses
    TLorentzVector parent;
    double m1, m2;

    // Outputs: daughters and decay angles in parent rest frame
    TLorentzVector daughter1, daughter2;
    double theta;  // polar angle (radians)
    double phi;    // azimuthal angle (radians)

    // Constructor: initializes parent & masses; RNG seed optional
    TwoBodyDecay(const TLorentzVector& parent_,
                 double mass1,
                 double mass2,
                 UInt_t seed = 0)
      : parent(parent_), m1(mass1), m2(mass2), rng(seed)
    {}

    // Generate decay given a reference vector 'refVec' that defines
    // the x-axis in the parent rest frame (must not be colinear with parent).
    // 'refVec' is provided in the lab frame.
    void Generate(const TVector3& refVec) {
        // --- 1) Compute momentum magnitude in rest frame ---
        double M = parent.M();
        double M2 = M*M;
        double term = (M2 - (m1 + m2)*(m1 + m2))
                    * (M2 - (m1 - m2)*(m1 - m2));
        double p = (term > 0 ? std::sqrt(term)/(2.0*M) : 0.0);

        //std::cout << "********  momentum = " << p << std::endl;
        // --- 2) Build orthonormal basis in lab: uz along parent momentum ---
        TVector3 uz = parent.Vect().Unit();
        TVector3 tmp = refVec.Unit();

        /*if (tmp.Cross(uz).Mag2() < 1e-6) {
            tmp = TVector3(1,0,0);
            if (std::fabs(uz.Dot(tmp)) > 0.9)
                tmp = TVector3(0,1,0);
        }*/
        TVector3 ydir = tmp.Cross(uz);
        TVector3 uy = ydir.Unit();
        TVector3 ux = uy.Cross(uz);
        //ux.Print();
        //uy.Print();
        //uz.Print();

        // --- 3) Randomly sample angles ---
        double cosTh = rng.Uniform(-1.0, 1.0);
        theta = std::acos(cosTh);
        phi   = rng.Uniform(0.0, 2*M_PI);
        double sinTh = std::sqrt(1 - cosTh*cosTh);

        // --- 4) Momentum vector in rest-frame basis ---
        //ux*(p * sinTh * std::cos(phi))
         //              + uy*(p * sinTh * std::sin(phi))
         //              + uz*(p * cosTh);
        TVector3 pVec =ux*(p * sinTh * std::cos(phi))
                       + uy*(p * sinTh * std::sin(phi))
                       + uz*(p * cosTh);
/*        pVec.SetXYZ(ux*(p * sinTh * std::cos(phi))
                       , uy*(p * sinTh * std::sin(phi))
                       , uz*(p * cosTh));
                       TVector3 pVec (ux*(p * sinTh * std::cos(phi))
                       , uy*(p * sinTh * std::sin(phi))
                       , uz*(p * cosTh));*/
        //std::cout << "  p = " << p  <<  "  " ;pVec.Print() ;
        // --- 5) Build daughters in parent rest frame ---
        TLorentzVector p1_rest(pVec,  std::hypot(p, m1));
        TLorentzVector p2_rest(-pVec, std::hypot(p, m2));

        //std::cout << "####  theta " << theta*180./3.14 << " phi " << phi*180./3.14 << " masses " << p1_rest.M() << "  " << p2_rest.M() << std::endl;
        // --- 6) Boost back to lab frame ---
        TVector3 boost = parent.BoostVector();
        daughter1 = p1_rest;  daughter1.Boost(boost);
        daughter2 = p2_rest;  daughter2.Boost(boost);
    }

    // Recompute theta, phi for existing daughters using refVec defining x-axis
    void CalculateAngles(const TVector3& refVec) {
        // --- 1) Build orthonormal basis as in Generate ---
        
        TVector3 uz = parent.Vect().Unit();
        TVector3 tmp = refVec.Unit();

        /*if (tmp.Cross(uz).Mag2() < 1e-6) {
            tmp = TVector3(1,0,0);
            if (std::fabs(uz.Dot(tmp)) > 0.9)
                tmp = TVector3(0,1,0);
        }*/
        TVector3 uy = tmp.Cross(uz).Unit();
        TVector3 ux = uy.Cross(uz);
        /*        TVector3 tmp = refVec;
        if (tmp.Cross(uz).Mag2() < 1e-6) {
            tmp = TVector3(1,0,0);
            if (std::fabs(uz.Dot(tmp)) > 0.9)
                tmp = TVector3(0,1,0);
        }
        TVector3 ux = (tmp - uz*(tmp.Dot(uz))).Unit();
        TVector3 uy = uz.Cross(ux);
        */
        // --- 2) Boost daughter1 to parent rest frame ---
        TLorentzVector d1 = daughter1;
        d1.Boost(-parent.BoostVector());
        TVector3 p = d1.Vect();

        // --- 3) Compute angles in this basis ---
        double pMag = p.Mag();
        if (pMag > 0) {
            theta = std::acos(p.Dot(uz) / pMag);
            phi   = std::atan2(p.Dot(uy), p.Dot(ux));
        } else {
            theta = 0.0;
            phi   = 0.0;
        }
    }

private:
    TRandom3 rng;
};

#endif // TWO_BODY_DECAY_H

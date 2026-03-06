// ROOTCompat.h - Minimal replacement for ROOT classes
// This header provides basic implementations of TLorentzVector, TVector3, and TRandom3
// to allow compilation without ROOT framework

#ifndef ROOT_COMPAT_H
#define ROOT_COMPAT_H

#include <cmath>
#include <iostream>
#include <random>

// Basic type compatibility
typedef unsigned int UInt_t;

//=============================================================================
// TVector3 - 3D vector class
//=============================================================================
class TVector3 {
private:
    double fX, fY, fZ;

public:
    TVector3() : fX(0), fY(0), fZ(0) {}
    TVector3(double x, double y, double z) : fX(x), fY(y), fZ(z) {}

    void SetXYZ(double x, double y, double z) { fX = x; fY = y; fZ = z; }
    double X() const { return fX; }
    double Y() const { return fY; }
    double Z() const { return fZ; }

    double Mag2() const { return fX*fX + fY*fY + fZ*fZ; }
    double Mag() const { return std::sqrt(Mag2()); }

    TVector3 Unit() const {
        double mag = Mag();
        return (mag > 0) ? TVector3(fX/mag, fY/mag, fZ/mag) : TVector3();
    }

    double Dot(const TVector3& v) const {
        return fX*v.fX + fY*v.fY + fZ*v.fZ;
    }

    TVector3 Cross(const TVector3& v) const {
        return TVector3(fY*v.fZ - fZ*v.fY,
                       fZ*v.fX - fX*v.fZ,
                       fX*v.fY - fY*v.fX);
    }

    TVector3 operator+(const TVector3& v) const {
        return TVector3(fX + v.fX, fY + v.fY, fZ + v.fZ);
    }

    TVector3 operator-(const TVector3& v) const {
        return TVector3(fX - v.fX, fY - v.fY, fZ - v.fZ);
    }

    TVector3 operator-() const {
        return TVector3(-fX, -fY, -fZ);
    }

    TVector3 operator*(double s) const {
        return TVector3(fX*s, fY*s, fZ*s);
    }

    void Print() const {
        std::cout << "(" << fX << ", " << fY << ", " << fZ << ")" << std::endl;
    }
};

//=============================================================================
// TLorentzVector - 4-vector class for relativistic kinematics
//=============================================================================
class TLorentzVector {
private:
    double fP[4];  // (px, py, pz, E)

public:
    TLorentzVector() { fP[0] = fP[1] = fP[2] = fP[3] = 0; }
    TLorentzVector(double px, double py, double pz, double e) {
        fP[0] = px; fP[1] = py; fP[2] = pz; fP[3] = e;
    }
    TLorentzVector(const TVector3& v, double e) {
        fP[0] = v.X(); fP[1] = v.Y(); fP[2] = v.Z(); fP[3] = e;
    }

    void SetPxPyPzE(double px, double py, double pz, double e) {
        fP[0] = px; fP[1] = py; fP[2] = pz; fP[3] = e;
    }

    double Px() const { return fP[0]; }
    double Py() const { return fP[1]; }
    double Pz() const { return fP[2]; }
    double E() const { return fP[3]; }

    double P2() const { return fP[0]*fP[0] + fP[1]*fP[1] + fP[2]*fP[2]; }
    double P() const { return std::sqrt(P2()); }

    double M2() const { return fP[3]*fP[3] - P2(); }
    double M() const {
        double m2 = M2();
        return (m2 >= 0) ? std::sqrt(m2) : -std::sqrt(-m2);
    }

    TVector3 Vect() const { return TVector3(fP[0], fP[1], fP[2]); }

    TVector3 BoostVector() const {
        double e = E();
        return (e > 0) ? TVector3(fP[0]/e, fP[1]/e, fP[2]/e) : TVector3();
    }

    void Boost(const TVector3& beta) {
        double b2 = beta.Mag2();
        if (b2 == 0) return;

        double gamma = 1.0 / std::sqrt(1.0 - b2);
        double bp = beta.Dot(Vect());
        double gamma2 = (b2 > 0) ? (gamma - 1.0) / b2 : 0.0;

        TVector3 p = Vect() + beta * (gamma2 * bp + gamma * fP[3]);
        fP[0] = p.X();
        fP[1] = p.Y();
        fP[2] = p.Z();
        fP[3] = gamma * (fP[3] + bp);
    }

    void SetVect(const TVector3& v) {
        fP[0] = v.X(); fP[1] = v.Y(); fP[2] = v.Z();
    }

    TLorentzVector operator+(const TLorentzVector& v) const {
        return TLorentzVector(fP[0] + v.fP[0], fP[1] + v.fP[1],
                             fP[2] + v.fP[2], fP[3] + v.fP[3]);
    }

    TLorentzVector operator-(const TLorentzVector& v) const {
        return TLorentzVector(fP[0] - v.fP[0], fP[1] - v.fP[1],
                             fP[2] - v.fP[2], fP[3] - v.fP[3]);
    }

    void Print() const {
        std::cout << "(" << fP[0] << ", " << fP[1] << ", " << fP[2]
                  << "; " << fP[3] << ")" << std::endl;
    }
};

//=============================================================================
// TRandom3 - Random number generator
//=============================================================================
class TRandom3 {
private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dist;

public:
    TRandom3(UInt_t seed = 0) : gen(seed), dist(0.0, 1.0) {}

    void SetSeed(UInt_t seed) {
        gen.seed(seed);
    }

    double Uniform(double min, double max) {
        return min + (max - min) * dist(gen);
    }
};

#endif // ROOT_COMPAT_H

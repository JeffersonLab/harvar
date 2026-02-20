// ----
// compile g++ -O2 -std=c++17 main.cpp -I. `root-config --cflags --glibs` -o main
#include "RhoEvent.h"
#include "RhoLundIO.h"
#include "Cross.h"

int main(int argc,char* argv[]) {

     const long     Nrequest = (argc>1) ? std::stol(argv[1]) : 200;
     const double   Ebeam    = (argc>2) ? std::stod(argv[2]) : 10.6;
    // choose your kinematics:
     double qmin = 2.0;
     double qmax = 2.1;
     double xmin = 0.3;
     double xmax = 0.31;
     
     double Q2 = 2.0;      // GeV^2
    double xB = 0.3;      
    //double t  = -0.35;    // GeV^2
    Cross cross;
    Candidate cand;
TRandom3 R(0);
    cross.generate(cand);
    //cand.ev.Print();

    cross.updateWeight(cand);


    const int Nscan = 50000;
     double wMax = 0.0;
     for(int i=0;i<Nscan;++i) {
         //double w = generateCandidate(R,Ebeam).weight;
       cross.generate(cand,qmin,qmax,xmin,xmax);
         cross.updateWeight(cand);
         double w = cand.weight;
         if(w>wMax) wMax = w;
     }
     wMax *= 1.2;                                // safety margin

     std::cerr << "Cross section max value = " << wMax << std::endl;
     //std::cout << " Max Cross = " << wMax << std::endl;
     //return 1;
     std::cout << std::fixed << std::setprecision(6);
     long Nkept=0, Nattempt=0;
     while(Nkept < Nrequest) {
       //auto cand = generateCandidate(R,Ebeam,0.299,0.301,2.49,2.51,-0.51,-0.49);
       //auto cand = generateCandidate(R,Ebeam);
       cross.generate(cand,qmin,qmax,xmin,xmax);
       cross.updateWeight(cand);
       ++Nattempt;
       if( R.Uniform(0.,wMax) <= cand.weight ) {
	 ++Nkept;
	 printLundCLAS12(cand.ev,12);
       }
     }
     std::cerr << "Generated " << Nkept << " events in "
               << Nattempt << " attempts  (ε = "
               << 100.*Nkept/Nattempt << " %)\n";

    return 0;
}

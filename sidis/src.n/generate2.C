#include "generator.h"

void generate2(){
  // define a reaction to decay to proton and rho
  // where rho decays to two pions
  sim::candidate cr(10.6,0.77,0.13957,0.13957);
  sim::generator gen(&cr);
  sim::event     event;
  
  gen.setRange(1.5,2.5,0.05,1.0);

  //gen.generate();
  //gen.updateWeight();
  
  double weight = gen.scan(150000);
  //printf("maximum wiegh scan = %f\n",weight);  
  // set range for generation, Q2 min, Q2 max, Xb min, Xb max
  
  
  for(int j = 0; j < 12000; j++) {
    gen.generate(weight);
    cr.react.getEvent(event);
    event.show();
    //cr.react.getEvent(event);
    //event.show();
  }
}

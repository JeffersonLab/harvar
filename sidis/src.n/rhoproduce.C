#include "generator.h"

void rhoproduce(){
  // define a reaction to decay to proton and rho
  // where rho decays to two pions
  sim::candidate cr(10.6,0.77,0.13957,0.13957);
  sim::generator gen(&cr);
  sim::event     event;
  
  // set range for generation, Q2 min, Q2 max, Xb min, Xb max
  gen.setRange(1.45,4.5,0.1,0.7);
  double weight = gen.scan(150000);
  
  for(int j = 0; j < 24000; j++) {
    gen.generate(weight);
    //cr.react.show();
    cr.react.getEvent(event);
    // printf(" %d\n", event.beamPol);
    event.show();
  }

  //gen.stats();
  
}

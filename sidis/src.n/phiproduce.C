#include "generator.h"

void phiproduce(){
  // define a reaction to decay to proton and rho
  // where rho decays to two pions

  sim::candidate cr(10.6,1.02,0.49368,0.49368);
  sim::generator gen(&cr);
  sim::event     event;
  
  cr.react.setDecayIds(333,321,-321);
  // set range for generation, Q2 min, Q2 max, Xb min, Xb max
  gen.setRange(2.45,2.5,0.35,0.36);
  double weight = gen.scan(150000);
  
  for(int j = 0; j < 24000; j++) {
    gen.generate(weight);
    //cr.react.show();
    cr.react.getEvent(event);
    // printf(" %d\n", event.beamPol);
    if(event.hasNaN()==false)
      event.show();
  }

  //gen.stats();
  
}

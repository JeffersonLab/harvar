#include "reaction.h"

void generate(){
  // define a reaction to decay to proton and rho
  // where rho decays to two pions
  sim::candidate cr(10.6,0.77,0.13957,0.13957);
  for(int i = 0; i < 25; i++){
    cr.react.generate(2.3, 0.3);
    cr.react.show();
  }

  printf("*****************************\n");
  printf("generate phi meson decay\n");
  printf("*****************************\n");

  sim::candidate cf(10.6,1.02,0.49368,0.49368);
  for(int i = 0; i < 25; i++){
    cf.react.generate(2.3, 0.3);
    cf.react.show();
  }
}

#include "Simulator.h"

int main(int argc, char *argv[]){

  assert_ge(argc,2);
  Simulator simulator;
  simulator.init(argv[1]);
  simulator.print();
  simulator.run();
  return 0;
}

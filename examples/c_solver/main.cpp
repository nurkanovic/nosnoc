#include"HomotopySolver.hpp"

int main(){
  HomotopySolver solver;

  std::map<std::string, casadi::DM> nlp_arg = {{"lbx", -casadi::inf},
                                               {"ubx",  casadi::inf},
                                               {"lbg",  0},
                                               {"ubg",  casadi::inf},
                                               {"x0",   0}};
  solver.solve(nlp_arg);
}

#include<HomotopySolver.hpp>
#include<iostream>

HomotopySolver::HomotopySolver()
{
  nlp_solver = casadi::nlpsol("solver", "ipopt", "nosnoc_solver_nlp.casadi");
  complemetarity_function = casadi::external("comp_res", "nosnoc_solver_comp.casadi");
}

uint32_t HomotopySolver::solve(std::map<std::string, casadi::DM> arg)
{
  p0[0] = 1;
  
  std::map<std::string, casadi::DM> nlp_arg = {{"lbx", lbw},
                                               {"ubx", ubw},
                                               {"lbg", lbg},
                                               {"ubg", ubg},
                                               {"x0", x0},
                                               {"p", p0}};
  

  auto nlpsolin = casadi::nlpsol_in();
  for(auto str : nlpsolin)
  {
    std::cout << str << std::endl;
  }
  for(auto p : p0)
  {
    std::cout << p << ", ";
  }

  auto res = nlp_solver(arg);
  return 0;
}
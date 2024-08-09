#include<HomotopySolver.hpp>
#include<iostream>

HomotopySolver::HomotopySolver()
{
  nlp_solver = casadi::nlpsol("solver", "ipopt", "{{opts.solver_name}}_nlp.casadi");
  complemetarity_function = casadi::external("comp_res", "{{opts.solver_name}}_comp.casadi");
}

uint32_t HomotopySolver::solve(std::map<std::string, casadi::DM> arg)
{
  p0[0] = 1;
  {% raw %}
  std::map<std::string, casadi::DM> nlp_arg = {{"lbx", lbw},
                                               {"ubx", ubw},
                                               {"lbg", lbg},
                                               {"ubg", ubg},
                                               {"x0", x0},
                                               {"p", p0}};
  {% endraw %}

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

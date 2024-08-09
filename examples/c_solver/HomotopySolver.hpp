#include <casadi/casadi.hpp>
#include <stdint.h>
#include <vector>
using namespace casadi;

class HomotopySolver{
 public:
  HomotopySolver();
  uint32_t solve(std::map<std::string, DM> arg);
 private:
  std::vector<double> lbw = {-inf, -inf};
  std::vector<double> ubw = {inf, inf};
  std::vector<double> lbg = {-inf};
  std::vector<double> ubg = {0};
  std::vector<double> p0 = {0, 1};
  std::vector<double> x0 = {0, 0};
  const std::vector<int> ind_mpcc = {0, 1};
  Function nlp_solver;
  Function complemetarity_function;
};
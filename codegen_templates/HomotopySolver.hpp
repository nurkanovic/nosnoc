#include <casadi/casadi.hpp>
#include <stdint.h>
#include <vector>
using namespace casadi;

class HomotopySolver{
 public:
  HomotopySolver();
  uint32_t solve(std::map<std::string, DM> arg);
 private:
  std::vector<double> lbw = {{'{' ~ nlp_lbw|join(', ') ~ '}' }};
  std::vector<double> ubw = {{'{' ~ nlp_ubw|join(', ') ~ '}' }};
  std::vector<double> lbg = {{'{' ~ nlp_lbg|join(', ') ~ '}' }};
  std::vector<double> ubg = {{'{' ~ nlp_ubg|join(', ') ~ '}' }};
  std::vector<double> p0 = {{'{' ~ nlp_p0|join(', ') ~ '}' }};
  std::vector<double> x0 = {{'{' ~ nlp_x0|join(', ') ~ '}' }};
  const std::vector<int> ind_mpcc = {{'{' ~ ind_mpcc|join(', ') ~ '}' }};
  Function nlp_solver;
  Function complemetarity_function;
};

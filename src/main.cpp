#include "lib.hpp"
using namespace std;

int main(int argc, char **argv) {
  // Semente para números aleatórios.
  srand(time(NULL));

  // set cout precision globally
  std::cout << std::fixed << std::setprecision(6);

  solve_it(FA, Matrix(vector<double>{3.0,1.0}), 4, 1e-7, 1e-2, GRADIENT);

  return 0;
}



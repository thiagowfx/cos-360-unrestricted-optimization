#include "lib.hpp"
using namespace std;

int main(int argc, char **argv) {
  // Semente para números aleatórios.
  srand(time(NULL));

  // cout precision globally
  // std::cout << std::fixed << std::setprecision(6);

  Matrix x0sub(2,1);
  x0sub.set(1, 3.0);
  x0sub.set(2, 1.0);

  solve_it(FA, x0sub, 10, 1e-6);
  // gradient_method(fa, gradfa, x, 1e-7);

  return 0;
}

#include "lib.hpp"
using namespace std;

int main(int argc, char **argv) {
  // Semente para números aleatórios.
  srand(time(NULL));

  double x1 = rand_int(10);
  double x2 = rand_int(10);

  Matrix x(2,1);
  x.set(1, x1);
  x.set(2, x2);
  x.debug();

  gradient_method(fa, gradfa, x, 1e-7);

  return 0;
}

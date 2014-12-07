#include "lib.hpp"
using namespace std;

int main(int argc, char **argv) {
  // Semente para números aleatórios.
  // Iniciar com uma semente diferente a cada execução do programa.
  // Essa semente é obtida em função do número de segundos desde 1o de janeiro de 1970.
  srand(time(NULL));

  // set cout precision globally
  DEFAULT_PRECISION = 7;
  std::cout << std::fixed << std::setprecision(DEFAULT_PRECISION);


  // Exemplo com o método de Newton com busca de armijo.
  /*
  solve_it(
      FA,
      Matrix(vector<double>{1.0,0.0}),
      4,
      1e-7,
      1e-2,
      NEWTON
      );
    */
   

  // Exemplo com o método de Newton Puro.
  /*
  solve_it(
      FA,
      Matrix(vector<double>{1.0,0.0}),
      4,
      1e-7,
      1e-2,
      NEWTONPURE
      );
   */


  // Exemplo com o método do gradiente.
  /*
  solve_it(
      FA,
      Matrix(vector<double>{3.0,1.0}),
      4,
      1e-7,
      1e-2,
      GRADIENT
      );
   */

  return 0;
}



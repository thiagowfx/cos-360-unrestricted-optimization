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

  // Resposta: o ponto encontrado pelo método de otimização.
  Matrix ans;

  /* 
     Como declarar uma Matrix 2 x 1?
    
     Método #1: use-o quando tiver que usar a matriz várias vezes
        Matrix m(2,1);
        m.set(1, 3.0);
        m.set(2, 1.0);

     Método #2: use-o quando só precisar usar a matriz uma vez, para passar como parâmetro para alguma função
        Matrix(vector<double>{3.0, 1.0});
  */


  // Exemplo com o método de Newton com busca de armijo.
  
  /*
  ans = solve_it(
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
  ans = solve_it(
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
  ans = solve_it(
      FA,
      Matrix(vector<double>{3.0,1.0}),
      4,
      1e-7,
      1e-2,
      GRADIENT
      );
   */

  // Exemplo com o método de Quasi Newton (BFGS)
  
  ans = solve_it(
      FA,
      Matrix(vector<double>{2.0,1.0}),
      4,
      1e-7,
      1e-1,
      QUASINEWTON
      );

  std::cout << "Error #1: " << (ans - Matrix(vector<double>{0.0,1.0})).mod() << std::endl;
  std::cout << "Error #2: " << fa(ans) - fa(Matrix(vector<double>{0.0,1.0})) << std::endl;

  return 0;
}



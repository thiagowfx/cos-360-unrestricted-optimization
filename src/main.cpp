#include "lib.hpp"
using namespace std;

int sum(int a, int b) {
  return a + b;
}

int diff(int a, int b) {
  return a - b;
}

int apply(std::function<int(int,int)> f, int x, int y) {
  return f(x,y);
}

int main(int argc, char **argv) {
  //cout << sum(1,2) << endl << apply(sum,1,2) << endl;
  //for (int i = 0; i < 5; ++i) {
  //cout << apply( [i](int x, int y) -> int {return x + y + i;} , 1, 2) << endl; }

  // Semente para números aleatórios.
  srand(time(NULL));

  // set cout precision globally
  // std::cout << std::fixed << std::setprecision(6);

  solve_it(FA, Matrix(vector<double>{3.0,1.0}), 4, 1e-5);
  // gradient_method(fa, gradfa, Matrix(vector<double>{3.0,1.0}), 1e-7);

  return 0;
}



#include "lib.hpp"
using namespace std;

int main() {
  Matrix m(2,2,1.0);
  m.set(2,1,2.0);
  m.debug();

  return 0;
}

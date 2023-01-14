#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>


using namespace std;

int main() {
  for (int j = 2; j < 6; j++) {
    for (int i = 1; i < 9; i++) {
      ofstream mon_flux2;
      mon_flux2.open("Resultats/plot_" + to_string(j) + "_test" + to_string(i), ios::out);

      mon_flux2 << "plot 'tf_Rusanov_test_" + to_string(i) +".txt' u 1:" + to_string(j) + " w lp" << endl;
      mon_flux2 << "replot 'tf_HLL_test_" + to_string(i) +".txt' u 1:" + to_string(j) + " w lp" << endl;
      mon_flux2 << "replot 'tf_HLLC_test_" + to_string(i) +".txt' u 1:" + to_string(j) + " w lp" << endl;

      mon_flux2.close();
    }

  }

  vector<vector<vector<double>>> U(2,vector<double> (2));

  U[0][0][0] = 1.;
  cout << U[0][0][0] << endl;


  return 0;
}

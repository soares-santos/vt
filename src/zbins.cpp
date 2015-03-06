#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

void usage(void) {
  cerr << "Usage: zbins <zmin> <zmax> <dzmin> <zzerr_file> \n" ;
  exit(8);
} 


// w(theta) = A * theta^(1-gamma)
// The parameters A, gamma vary with redshift.
// This relationship was determined with simulations.
// The relevant file is wtheta_fit_DES_mock_v1.07.txt
double gamma(double z){ return 1.71-(0.32*z); }
double A(double z){ return 0.02; }



int main(int argc, char *argv[]) {

  if(argc != 5) usage();

  double zi=atof(argv[1]);
  double zf=atof(argv[2]);
  double dz=atof(argv[3]);

  double z0=zi, z1=zi+dz;

  int n=0;

  double z=0,zerr=0,sigmaz=0;

  ifstream zzerr (argv[4]);

  if (zzerr.is_open()){
    while(!zzerr.eof()){
      zzerr >> z >> zerr;
      n++;
      sigmaz+=zerr/(1+z);
    }
    zzerr.close();
    sigmaz/=n;
  }

  sigmaz=(sigmaz>dz)?sigmaz:dz;

  cout << "# zmin \t\t zmax \t\t A \t gamma\n";

  while (z0<zf) {
    z1=z0+2*sigmaz*(1+z0);
    cout << "  " << z0 << " \t " << z1 << " \t " << A(z0) << " \t " << gamma(z0) << "\n";
    z0=z1;
  }

  return 0;
}


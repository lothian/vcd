#ifndef PERTURBATION_H
#define PERTURBATION_H

#include "psi4/libmints/mintshelper.h"

// Allowed operators:
//
// Mu = electric dipole = -r
// P  = velocity-gauge electric dipole = -i Nabla
// P* = complex conjugate of P
// L  = magnetic dipole = -(1/2) r x P
// L* = complex conjugate of L
// Q  = traceless quadrupole
// RR = quadrupole
// 

using namespace std;

namespace psi {

class Perturbation {
public:
  std::string operator_; // perturbation name
  Perturbation(std::string op, shared_ptr<Wavefunction> ref, 
               shared_ptr<MintsHelper> mints, bool full_virtual_space);
  ~Perturbation();
  double **prop_p(int i) { return prop_[i]; }
  double **prop_p(int i, int j) // for quarupolar peturbations
  { 
    int ij = ((i) > (j) ? (i)*((i)+1)/2 + (j) : (j)*((j)+1)/2 + (i));
    return prop_[ij];
  }
  void print(std::string out);
  void print(int i, std::string out);
  void print(int i, int j, std::string out);
  void print();
  void print(int i);
  void print(int i, int j);

protected:
  int nact_;
  double ***prop_;

private:
  bool allowed(std::string op);
  bool dipole(std::string op);
  bool quadrupole(std::string op);
};

} // psi

#endif // PERTURBATION_H

#ifndef MP2_H
#define MP2_H

#include "hamiltonian.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"

using namespace std;

namespace psi {

class MP2Wfn: public Wavefunction {
public:
  MP2Wfn(shared_ptr<Wavefunction> reference, shared_ptr<Hamiltonian> H, Options &options);
  virtual ~MP2Wfn();
  double energy();

protected:
  std::string wfn_;     // wfn type (CCSD, CCSD_T, etc.)
  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  shared_ptr<Hamiltonian> H_; // integrals and Fock matrix

  // Energy denominators
  double ****D2_;
  // First-order wave function
  double ****t2_;
}; // 

} // psi

#endif // MP2_H

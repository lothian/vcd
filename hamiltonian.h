#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "psi4/libmints/mintshelper.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.h"

using namespace std;

namespace psi {

class Hamiltonian {
public:
  Hamiltonian(shared_ptr<PSIO>, shared_ptr<Wavefunction>, std::vector<shared_ptr<MOSpace> >);
  virtual ~Hamiltonian();

protected:
  int nmo_;
  int nso_;
  int nact_;
  int nfzc_;
  int nfzv_;
  double efzc_;

  double **fock_;
  double ****ints_;
  double ****L_;

  friend class MP2Wfn;
}; // Hamiltonian

} // psi

#endif // HAMILTONIAN_H

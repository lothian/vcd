#include "mp2wfn.h"
#include "array.h"
#include "hamiltonian.h"
#include <psi4/libqt/qt.h>
#include <psi4/libciomr/libciomr.h>
#include <cmath>

namespace psi {

MP2Wfn::MP2Wfn(shared_ptr<Wavefunction> reference, shared_ptr<Hamiltonian> H,
             Options &options) : Wavefunction(options)
{
  wfn_ = options.get_str("WFN");

  // What does this do?
  set_reference_wavefunction(reference);
  shallow_copy(reference);

  H_ = H;

  int nfrzv = 0;
  no_ = nv_ = 0;
  for(int i=0; i < nirrep_; i++) {
    no_ += doccpi()[i] - frzcpi_[i];
    nv_ += nmopi_[i] - doccpi()[i] - frzvpi_[i];
    nfrzv += frzvpi_[i];
  }

  std::vector<std::string> labels = molecule_->irrep_labels();

  outfile->Printf("\n\tReference Wfn Parameters:\n");
  outfile->Printf("\t---------------------------\n");
  outfile->Printf("\tNumber of irreps        = %d\n", nirrep_);
  outfile->Printf("\tNumber of MOs           = %d\n", nmo_);
  outfile->Printf("\tNumber of active MOs    = %d\n", no_+nv_);
  outfile->Printf("\tNumber of active occ    = %d\n", no_);
  outfile->Printf("\tNumber of active vir    = %d\n", nv_);
  outfile->Printf("\tNumber of frozen occ    = %d\n", nfrzc_);
  outfile->Printf("\tNumber of frozen vir    = %d\n\n", nfrzv);
  outfile->Printf("\tLabel\t# MOs\t# FZDC\t# DOCC\t# VIRT\t# FZVR\n");
  outfile->Printf("\t-----\t-----\t------\t------\t------\t------\n");
  for(int i=0; i < nirrep_; i++) {
      outfile->Printf("\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\n",
              labels[i].c_str(),nmopi_[i],frzcpi_[i],doccpi()[i],nmopi_[i]-doccpi()[i],frzvpi_[i]);
    }
  outfile->Printf("\n\tNuclear Repulsion Energy    = %20.15f\n", molecule_->nuclear_repulsion_energy(reference_wavefunction_->get_dipole_field_strength()));
  outfile->Printf( "\tFrozen Core Energy          = %20.15f\n", H_->efzc_);
  outfile->Printf( "\tTotal SCF Energy (ref)      = %20.15f\n", reference_wavefunction_->energy());

  // Prepare energy denominators
  int no = no_;
  int nv = nv_;
  double **fock = H_->fock_;
  double ****ints = H_->ints_;

  D2_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          D2_[i][j][a][b] = fock[i][i] + fock[j][j] - fock[a+no][a+no] - fock[b+no][b+no];

  t2_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          t2_[i][j][a][b] = ints[i][j][a+no][b+no]/D2_[i][j][a][b];
}

MP2Wfn::~MP2Wfn()
{
  int no = no_;
  int nv = nv_;
  free_4d_array(D2_, no, no, nv);
  free_4d_array(t2_, no, no, nv);
}

double MP2Wfn::energy() {
  double emp2=0.0;

  int no = no_;
  int nv = nv_;
  double ****L = H_->L_;
  double ****t2 = t2_; 

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          emp2 += t2[i][j][a][b]*L[i][j][a+no][b+no];

  return emp2;
}

} // psi


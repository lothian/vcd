/*
 * @BEGIN LICENSE
 *
 * ugacc by T. Daniel Crawford, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"

#include <psi4/libciomr/libciomr.h>
#include <psi4/libqt/qt.h>
#include "psi4/libmints/mintshelper.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/dipole.h"
#include <map>

#include "hamiltonian.h"
#include "perturbation.h"
#include "mp2wfn.h"

#include "array.h"


using namespace std;

namespace psi {

double levi(int a, int b, int c);

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
  if(name == "VCD" || options.read_globals()) {
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "MP2", "RHF MP2 CCSD CCSD(T)");
    options.add_str("DERTYPE", "NONE");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("R_CONVERGENCE", 1e-7);
    options.add_str("AAT_TYPE", "ANALYTIC", "NUMERICAL");
  }

  return true;
}

extern "C" PSI_API
SharedWavefunction vcd(SharedWavefunction ref, Options& options)
{
  outfile->Printf("\t*************************\n");
  outfile->Printf("\t*                       *\n");
  outfile->Printf("\t*          VCD          *\n");
  outfile->Printf("\t*                       *\n");
  outfile->Printf("\t*************************\n");
  outfile->Printf("\n");

  outfile->Printf("\tWave function  = %s\n", options.get_str("WFN").c_str());
  outfile->Printf("\tDertype        = %s\n", options.get_str("DERTYPE").c_str());

  // Error trapping – need closed-shell SCF in place
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");
  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references.");
  for(int h=0; h < ref->nirrep(); h++)
    if(ref->soccpi()[h]) throw PSIEXCEPTION("VCD is for closed-shell systems only.");

  // Error trapping – no frozen core allowed for now
  if(ref->nfrzc() != 0) throw PSIEXCEPTION("VCD must be run without frozen core for now.");

  // Error trapping – need no symmetry for now
  if(ref->nirrep() != 1) throw PSIEXCEPTION("VCD is for C1 symmetry only for now.");
  
  // Set up I/O object
  shared_ptr<PSIO> psio(_default_psio_lib_);

  // Prepare MO space vector that runs over all orbitals
  std::vector<shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);


  // Prepare Hamiltonian (transform the integrals and sort them into member arrays)
  shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces));

  // shared_ptr<MP2Wfn> mp2(new MP2Wfn(ref, H, options));
  // double emp2 = mp2->energy(); // MP2 correlation energy
  // outfile->Printf("\tMP2 Correlation Energy      = %20.14f\n", emp2);
  // outfile->Printf("\tMP2 Total Energy            = %20.14f\n", emp2+ref->energy());

  outfile->Printf("\n\tReference Wfn Parameters:\n");
  outfile->Printf("\t---------------------------\n");
  outfile->Printf("\tNumber of irreps        = %d\n", ref->nirrep());
  int nmo = ref->nmo();
  outfile->Printf("\tNumber of MOs           = %d\n", nmo);
  int no = ref->doccpi()[0];
  int nv = ref->nmopi()[0] - ref->doccpi()[0];
  outfile->Printf("\tNumber of active MOs    = %d\n", no+nv);
  outfile->Printf("\tNumber of active occ    = %d\n", no);
  outfile->Printf("\tNumber of active vir    = %d\n", nv);
  outfile->Printf("\tNumber of frozen occ    = %d\n", ref->nfrzc());
  outfile->Printf("\tNumber of frozen vir    = %d\n\n", ref->frzvpi()[0]);
  int natom = ref->molecule()->natom();
  outfile->Printf("\tNumber of atoms         = %d\n", natom);

  shared_ptr<MintsHelper> mints(new MintsHelper(ref->basisset(), options, 0));

  // Prepare spin-adapted TEIs
  SharedMatrix TEI = mints->mo_eri(ref->Ca(), ref->Ca(), ref->Ca(), ref->Ca());
  SharedMatrix L(new Matrix("Spin-Adapted TEIs", nmo*nmo, nmo*nmo));
  for(int p=0; p < nmo; p++)
    for(int q=0; q < nmo; q++) {
      int pq = p * nmo + q;
      for(int r=0; r < nmo; r++) {
        int pr = p * nmo + r;
        for(int s=0; s < nmo; s++) {
          int rs = r * nmo + s;
          int qs = q * nmo + s;
          int ps = p * nmo + s;
          int qr = q * nmo + r;
          L->set(pq,rs, 2.0 * TEI->get(pr,qs) - TEI->get(ps,qr));
        }
      }
    }

  // =============
  // RHF Energy 
  // =============

  SharedMatrix h = mints->ao_kinetic();
  h->add(mints->ao_potential());
  h->set_name("Core Hamiltonian");
  h->transform(ref->Ca());
  double e1 = 0.0;
  double e2 = 0.0;
  for(int i=0; i < no; i++) {
    e1 += 2.0 * h->get(i,i);
    for(int j=0; j < no; j++) {
      int ij = i * nmo + j;
      e2 += L->get(ij, ij);
    }
  }
  double enuc = ref->molecule()->nuclear_repulsion_energy(ref->get_dipole_field_strength());
  outfile->Printf("\tRHF one-electron energy     = %20.12f\n", e1);
  outfile->Printf("\tRHF two-electron energy     = %20.12f\n", e2);
  outfile->Printf("\tRHF total electronic energy = %20.12f\n", e1 + e2);
  outfile->Printf("\tRHF total energy            = %20.12f\n", e1 + e2 + enuc);

  // ==================
  // RHF Dipole Moment
  // ==================
  std::string cart = "XYZ";
  std::vector<SharedMatrix> dipole = mints->ao_dipole();
  SharedMatrix mu(new Matrix("MO Basis Electric Dipole Integrals", no, no));
  outfile->Printf("\n\tElectric Dipole Moment:\n");
  for(int coord=0; coord < 3; coord++) {
    mu->transform(ref->Ca_subset("AO", "OCC"), dipole[coord], ref->Ca_subset("AO", "OCC"));   
    double dipmom_e=0.0;
    for(int i=0; i < no; i++) dipmom_e += mu->get(i,i);
    dipmom_e *= 2.0;

    double dipmom_n=0.0;
    for(int atom=0; atom < natom; atom++) {
      double geom = ref->molecule()->geometry().get(atom, coord);
      double z = ref->molecule()->Z(atom);
      dipmom_n += geom * z;
    }

    outfile->Printf("\tmu_e(%c) = %20.14f \t mu_n(%c) = %20.14f \t mu(%c) = %20.14f\n", 
                    cart[coord], dipmom_e, cart[coord], dipmom_n, cart[coord], dipmom_e + dipmom_n);
  }
  outfile->Printf("\n");

  // ====================
  // Derivative integrals
  // ====================
  std::vector<SharedMatrix> S_deriv1;
  std::vector<SharedMatrix> h_deriv1;
  std::vector<SharedMatrix> F_deriv1;
  std::vector<SharedMatrix> L_deriv1;
  {
    std::vector<SharedMatrix> dS;
    std::vector<SharedMatrix> dh;
    std::vector<SharedMatrix> dV;
    std::vector<SharedMatrix> dTEI;
    SharedMatrix dL(new Matrix("Spin-Adapted TEI Derivatives", nmo*nmo, nmo*nmo));
    SharedMatrix dF(new Matrix("Skeleton Fock Derivative", nmo, nmo));
    for(int atom=0; atom < natom; atom++) {
      dS = mints->mo_oei_deriv1("OVERLAP", atom, ref->Ca(), ref->Ca());
      dh = mints->mo_oei_deriv1("KINETIC", atom, ref->Ca(), ref->Ca());
      dV = mints->mo_oei_deriv1("POTENTIAL", atom, ref->Ca(), ref->Ca());
      dTEI = mints->mo_tei_deriv1(atom, ref->Ca(), ref->Ca(), ref->Ca(), ref->Ca());
      for(int coord=0; coord < 3; coord++) {
	int R_coord = atom * 3 + coord;

        std::string s = "Overlap Integral Derivative (" + to_string(atom) + ", " + to_string(coord) + ")";
        dS[coord]->set_name(s);
        S_deriv1.push_back(dS[coord]->clone());
	
        s = "Core Hamiltonian Integral Derivative (" + to_string(atom) + ", " + to_string(coord) + ")";
        dh[coord]->set_name(s);
        dh[coord]->add(dV[coord]);
	h_deriv1.push_back(dh[coord]->clone());

        // Build spin adapted TEI derivs for current coordinate
        for(int p=0; p < nmo; p++) {
          for(int q=0; q < nmo; q++) {
            int pq = p * nmo + q;
            for(int r=0; r < nmo; r++) {
              int pr = p * nmo + r;
              for(int s=0; s < nmo; s++) {
                int rs = r * nmo + s;
                int qs = q * nmo + s;
                int ps = p * nmo + s;
                int qr = q * nmo + r;
                dL->set(pq,rs, 2.0 * dTEI[coord]->get(pr,qs) - dTEI[coord]->get(ps,qr));
              }
            }
          }
	}
        s = "Spin-Adapted TEI Derivative (" + to_string(atom) + ", " + to_string(coord) + ")";
	dL->set_name(s);
	L_deriv1.push_back(dL->clone());

        // Build skeleton Fock derivs
        for(int p=0; p < nmo; p++) {
          for(int q=0; q < nmo; q++) {
            double val = h_deriv1[R_coord]->get(p,q);
            for(int i=0; i < no; i++) {
              int pi = p * nmo + i;
              int qi = q * nmo + i;
              val += L_deriv1[R_coord]->get(pi,qi);
            }
            dF->set(p, q, val);
          } // q
        } // p
        s = "Skeleton Fock Derivative (" + to_string(atom) + ", " + to_string(coord) + ")";
        dF->set_name(s);
        F_deriv1.push_back(dF->clone());

      } // coord
    } // atom
  } // derivative integral preparation 


  // =============
  // RHF Gradient
  // =============

  SharedMatrix grad_total(new Matrix("RHF Gradient", natom, 3));
  SharedMatrix grad_one(new Matrix("Core Hamiltonian Gradient", natom, 3));
  SharedMatrix grad_S(new Matrix("Overlap Gradient", natom, 3));
  SharedMatrix grad_two(new Matrix("Coulomb Gradient", natom, 3));
  SharedMatrix f = ref->Fa(); // Fock matrix (already transformed by Hamiltonian object)
  for(int atom=0; atom < natom; atom++) {
    for(int coord=0; coord < 3; coord++) { // Cartesian coordinate
      int R_coord = atom * 3 + coord;
      double grad1 = 0.0;
      double gradS = 0.0;
      double grad2 = 0.0;
      for(int i=0; i < no; i++) {
        grad1 += 2.0 * h_deriv1[R_coord]->get(i,i);
        gradS -= 2.0 * S_deriv1[R_coord]->get(i,i) * f->get(i,i);
        int ii = i * nmo + i;
        for(int j=0; j < no; j++) {
          int jj = j * nmo + j;
          int ij = i * nmo + j;
          int ji = j * nmo + i;
          grad2 += L_deriv1[R_coord]->get(ij,ij);
        } // j
      } // i
      grad_one->set(atom, coord, grad1);
      grad_S->set(atom, coord, gradS);
      grad_two->set(atom, coord, grad2);
      grad_total->set(atom, coord, grad1 + gradS + grad2);
    } // p
  } // atom

  Matrix nucgrad = ref->molecule()->nuclear_repulsion_energy_deriv1(ref->get_dipole_field_strength());
  nucgrad.print("outfile");
  grad_one->print("outfile");
  grad_S->print("outfile");
  grad_two->print("outfile");
  grad_total->add(nucgrad);
  grad_total->print("outfile");

  std::vector<SharedMatrix> U_R; // CPHF coefficients for nuclear coordinate perturbations
  std::vector<SharedMatrix> U_F; // CPHF coefficients for electric field perturbations
  {
    // ================================================
    // CPHF Equations (real perturbations)
    // ================================================

    // MO Hessian for real perturbations
    SharedMatrix G(new Matrix("MO Hessian (Real)", no*nv, no*nv));
    for(int a=0; a < nv; a++) {
      for(int i=0; i < no; i++) {
        int ai = a * no + i;
        for(int e=0; e < nv; e++) {
          int ae = (a + no) * nmo + (e + no);
          for(int m=0; m < no; m++) {
            int em = e * no + m;
            int im = i * nmo + m;
            int am = (a + no) * nmo + m;
            int ie = i * nmo + (e + no);
            double val = (f->get(e+no,e+no) - f->get(m, m))*(e==a)*(m==i);
            val += L->get(ae,im) + L->get(am,ie);
            G->set(ai, em, val);
          }
        }
      }
    }
    double **Gp = G->pointer();
  
    SharedMatrix Gclone = G->clone();
  
    // CPHF for nuclear-coordinate and electric-field perturbations
    SharedMatrix B(new Matrix("CPHF RHS Nuclear Coordinates", nv, no));
    double **Bp = B->pointer();
    int *ipiv = init_int_array(no*nv);
    for(int atom=0; atom < natom; atom++) {
      for(int coord=0; coord < 3; coord++) {
        int R_coord = atom * 3 + coord;
        // CPHF B vector
        for(int a=0; a < nv; a++) {
          for(int i=0; i < no; i++) {
            double val = F_deriv1[R_coord]->get(a+no, i);
            val -= S_deriv1[R_coord]->get(a+no, i) * f->get(i, i);
            for(int m=0; m < no; m++)
              for(int n=0; n < no; n++) {
                int am = (a + no) * nmo + m;
                int an = (a + no) * nmo + n;
                int im = i * nmo + m;
                int in = i * nmo + n;
                val -= 0.5 * S_deriv1[R_coord]->get(m,n) * (L->get(am,in) + L->get(an,im));
              }
            B->set(a, i, -val);
          }
        }
  
        // Solve CPHF Equations
        for(int ai=0; ai < no*nv; ai++) ipiv[ai] = 0.0;
        int errcod = C_DGESV(nv*no, 1, Gp[0], nv*no, ipiv, Bp[0], nv*no);
        G->copy(Gclone); // restore the MO Hessian
        SharedMatrix U = B->clone();
        std::string s = "CPHF Coefficients for R(" + to_string(atom) + ", " + to_string(coord) + ")";
        U->set_name(s);
        U->print();
        U_R.push_back(U);

      } // coord
    } // atom

    // Electric-Field Response
    SharedMatrix BF(new Matrix("CPHF RHS Electric Field", nv, no));
    double **BFp = BF->pointer();
    for(int coord=0; coord < 3; coord++) {
      BF->transform(ref->Ca_subset("AO","VIR"), dipole[coord], ref->Ca_subset("AO","OCC"));
      // BF->scale(-1.0);
      for(int ai=0; ai < no*nv; ai++) ipiv[ai] = 0.0;
      int errcod = C_DGESV(nv*no, 1, Gp[0], nv*no, ipiv, BFp[0], nv*no);
      G->copy(Gclone); // restore the MO Hessian
      SharedMatrix U = BF->clone();
      std::string s = "CPHF Coefficients for F(";
      s.push_back(cart[coord]);
      s += ")";
      U->set_name(s);
      U->print();
      U_F.push_back(U);
    }

  } // CPHF(R, F)

  std::vector<SharedMatrix> U_B; // CPHF coefficients for magnetic field perturbations
  {
    // ================================================
    // CPHF Equations (magnetic field perturbation)
    // ================================================

    // MO Hessian for imaginary perturbations
    SharedMatrix G(new Matrix("MO Hessian (Imaginary)", no*nv, no*nv));
    for(int a=0; a < nv; a++) {
      for(int i=0; i < no; i++) {
        int ai = a * no + i;
        for(int e=0; e < nv; e++) {
          int ae = (a + no) * nmo + (e + no);
          for(int m=0; m < no; m++) {
            int em = e * no + m;
            int im = i * nmo + m;
            int mi = m * nmo + i;
            int am = (a + no) * nmo + m;
            int ie = i * nmo + (e + no);
            int me = m * nmo + (e + no);
            int ei = (e + no) * nmo + i;
            double val = (f->get(e+no,e+no) - f->get(m, m))*(e==a)*(m==i);
            val += L->get(am,ie) - L->get(ae,im);
            G->set(ai, em, val);
          }
        }
      }
    }
    double **Gp = G->pointer();

    SharedMatrix Gclone = G->clone();

    // CPHF for magnetic field perturbations
    SharedMatrix B(new Matrix("CPHF RHS Magnetic Field Components", nv, no));
    double **Bp = B->pointer();
    int *ipiv = init_int_array(no*nv);

    std::vector<SharedMatrix> angmom = mints->so_angular_momentum();
    SharedMatrix dm(new Matrix("Magnetic Moment Integrals", nv, no));
    for(int coord=0; coord < 3; coord++) {
      dm->transform(ref->Ca_subset("AO","VIR"), angmom[coord], ref->Ca_subset("AO","OCC"));
      dm->scale(-0.5);

      for(int a=0; a < nv; a++)
        for(int i=0; i < no; i++)
          B->set(a,i, dm->get(a,i));

      // Solve CPHF Equations
      for(int ai=0; ai < no*nv; ai++) ipiv[ai] = 0.0;
      int errcod = C_DGESV(nv*no, 1, Gp[0], nv*no, ipiv, Bp[0], nv*no);
      G->copy(Gclone); // restore the MO Hessian
      SharedMatrix U = B->clone();
      std::string s = "CPHF Coefficients for B(" + to_string(coord) + ")";
      U->set_name(s);
      U->print();
      U_B.push_back(U);
    } // coord
  } // CPHF(B) 

    // ================================================
    // HF APTs
    // ================================================
  {
    SharedMatrix APT(new Matrix("APT Total", natom*3, 3));

    // Contribution of skeleton derivative integrals
    SharedMatrix APT_elec = mints->dipole_grad(ref->Da_subset("AO"));
    APT_elec->scale(2.0);

    std::vector<SharedMatrix> Mu = mints->ao_dipole();
    for(int coord=0; coord < 3; coord++) Mu[coord]->transform(ref->Ca());

    for(int atom=0; atom < natom; atom++) {
      for(int coord=0; coord < 3; coord++) {
        int R_coord = atom * 3 + coord;
        for(int dip_coord=0; dip_coord < 3; dip_coord++) {

          // Contribution of overlap derivatives
          double val=0;
          for(int i=0; i < no; i++) {
            for(int j=0; j < no; j++) {
              val += -2.0 * S_deriv1[R_coord]->get(i,j) * Mu[dip_coord]->get(i,j);
            } 
          }

          // Contribution of CPHF(R) coefficients
          for(int a=0; a < nv; a++) {
            for(int i=0; i < no; i++) {
              val += +4.0 * U_R[R_coord]->get(a,i) * Mu[dip_coord]->get(a+no,i);
            }
          }
          APT_elec->add(R_coord, dip_coord, val);

        } // dip_coord
      } // atom coord
    } // atom

    // Gradient of nuclear dipole moment
    SharedMatrix APT_n(new Matrix("APT Nuclear Component", natom*3, 3));
    for(int atom=0; atom < natom; atom++) {
      double z = ref->molecule()->Z(atom);
      APT_n->set(atom * 3 + 0, 0, z);
      APT_n->set(atom * 3 + 1, 1, z);
      APT_n->set(atom * 3 + 2, 2, z);
    }
    APT_n->print();

    APT->add(APT_elec);
    APT->add(APT_n);
    APT->print();
  }

  {
    // ================================================
    // HF AATs
    // ================================================

    double psi_dipmom_au2si = 8.478353552E-30;
    double psi_hbar = 1.054571800E-34;
    double psi_m2angstroms = 1.0E10;
    double aatconv = psi_dipmom_au2si * (1 / psi_hbar) * (1 / psi_m2angstroms) * 1.0E6;

    std::vector<SharedMatrix> halfS_deriv;
    SharedMatrix halfdS(new Matrix("Half-Derivative Overlap", no, nv));
    SharedMatrix AAT_elec(new Matrix("AAT Electronic Component", natom*3, 3));
    SharedMatrix AAT_nuc(new Matrix("AAT Nuclear Component", natom*3, 3));
    ref->molecule()->geometry().print();
    for(int atom = 0; atom < natom; atom++) {
      halfS_deriv = mints->ao_overlap_half_deriv1("LEFT", atom);
      for(int coord=0; coord < 3; coord++) {
        halfdS->transform(ref->Ca_subset("AO","OCC"), halfS_deriv[coord], ref->Ca_subset("AO","VIR"));
//	halfdS->scale(-1.0);
	int R_coord = atom * 3 + coord;
        for(int B_coord = 0; B_coord < 3; B_coord++) {

          double val = 0.0;
	  for(int a=0; a < nv; a++)
	    for(int i=0; i < no; i++) {
              val += U_R[R_coord]->get(a,i) * U_B[B_coord]->get(a,i);
              val += halfdS->get(i,a) * U_B[B_coord]->get(a,i);
	    }
          AAT_elec->set(R_coord, B_coord, 2.0 * val);

	  val = 0.0;
	  for(int gamma=0; gamma < 3; gamma++) {
	    double geom = ref->molecule()->geometry().get(atom, gamma);
	    double z = ref->molecule()->Z(atom);
	    val += levi(coord, B_coord, gamma) * geom * z / 4.0;
	  }
          AAT_nuc->set(R_coord, B_coord, val);

        } // B_coord
      } // coord (R)
    } // atom

    // AAT_elec * psi_dipmom_au2si * (1 / psi_hbar) * (1 / psi_m2angstroms) * 10**6

    outfile->Printf("\n\t*******************\n");
    outfile->Printf(  "\t** Atomic Units: **\n");
    outfile->Printf(  "\t*******************\n\n");
    AAT_elec->print();
    AAT_nuc->print();
    SharedMatrix AAT(new Matrix("AAT Total", natom*3, 3));
    AAT->add(AAT_elec);
    AAT->add(AAT_nuc);
    AAT->print();

    outfile->Printf("\n\t**********************\n");
    outfile->Printf(  "\t** T^-1 A^_1 Units: **\n");
    outfile->Printf(  "\t**********************\n\n");
    AAT_elec->scale(aatconv);
    AAT_elec->print();
    AAT_nuc->scale(aatconv);
    AAT_nuc->print();
    AAT->scale(aatconv);
    AAT->print();

  }
  
  return ref;
}

/* A stupid Levi-Civita evaluator */
double levi(int a, int b, int c) {
  double val=0;
  int x=0, y=1, z=2;

  if(a==x && b==y && c==z) val=1;
  else if(a==y && b==z && c==x) val=1;
  else if(a==z && b==x && c==y) val=1;
  else if(a==x && b==z && c==y) val=-1;
  else if(a==y && b==x && c==z) val=-1;
  else if(a==z && b==y && c==x) val=-1;
  else val=0;

  return val;
}


} // End namespaces

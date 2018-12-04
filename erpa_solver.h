/*
 *@BEGIN LICENSE
 *
 * ERPA, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
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
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 *
 *@END LICENSE
 *
 */

#ifndef ERPA_SOLVER_H
#define ERPA_SOLVER_H


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <psi4/libplugin/plugin.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libpsio/psio.hpp>
//#include <../bin/fnocc/blas.h>
//#include "blas.h"
//#include <psi4/libqt/qt.h>

#include<psi4/libiwl/iwl.h>
#include <psi4/psifiles.h>

#include<psi4/libmints/wavefunction.h>
#include<psi4/libmints/matrix.h>
#include<psi4/libmints/vector.h>

#include "fortran.h"

// TODO: move to psifiles.h
#define PSIF_DCC_QMO          268
#define PSIF_V2RDM_CHECKPOINT 269
#define PSIF_V2RDM_D2AA       270
#define PSIF_V2RDM_D2AB       271
#define PSIF_V2RDM_D2BB       272
#define PSIF_V2RDM_D3AAA      273
#define PSIF_V2RDM_D3AAB      274
#define PSIF_V2RDM_D3BBA      275
#define PSIF_V2RDM_D3BBB      276
#define PSIF_V2RDM_D1A        277
#define PSIF_V2RDM_D1B        278

namespace psi{ namespace erpa{

class ERPASolver: public Wavefunction{
  public:
    ERPASolver(SharedWavefunction reference_wavefunction,Options & options);
    ~ERPASolver();
    void common_init();
    //double compute_energy();
    virtual bool same_a_b_orbs() const { return same_a_b_orbs_; }
    virtual bool same_a_b_dens() const { return same_a_b_dens_; }

    // public methods
    void EvaluateSigma(long int N,long int maxdim, long int L, double ** sigma, double **b);

    /// extended RPA
    void IPExtendedRPA();
    void OldExtendedRPA(std::string type);
    void ExtendedRPA(std::string type);

    double DOCIAlphaERPA(double alpha,bool is_singlet);
    double AlphaERPA(double alpha,bool is_singlet,bool is_excited_state_computation);

    /// efficient ERPA implementation
    void NewExtendedRPA(std::string type);
    void SpinAdaptedExtendedRPA();

    /// ERPA 2
    void ExtendedRPA2();

  protected:

    void ReconstructD3(double * D3aaa, double * D3aab, double * D3bba, double * D3bbb, double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b);

    /// build excited-state D1 using ground-state 1-,2-,3-RDM and ERPA wave function
    void BuildExcitedStateD1(int h,double*c,double*D1a,double*D1b,double*D2aa,double*D2ab,double*D2bb,double*D3aaa,double*D3aab,double*D3abb,double*D3bbb,int print);
    void BuildExcitedStateD2(int h,double*c,double*D1a,double*D1b,double*D2aa,double*D2ab,double*D2bb,double*D3aaa,double*D3aab,double*D3bba,double*D3bbb,double*D4aaaa,double*D4aaab,double*D4aabb,double*D4bbba,double*D4bbbb,int print);

    /// read FCI density
    void FCIDensity(double*D1a,double*D1b,double*D2aa,double*D2bb,double*D2ab);

    /// read density from detci, called by the function FCIDensity
    void ReadDetciDensity(double * d2aa, double * d2bb, double * d2ab, double * d1a, double * d1b, int nmo); 

    /// symmetry product table:
    int * table;

    /// returns symmetry product for two orbitals
    int SymmetryPair(int i, int j);

    /// returns symmetry product for four orbitals
    int TotalSym(int i, int j,int k, int l);

    int * symmetry;
    int * symmetry_full;
    int * symmetry_really_full;
    int * energy_to_pitzer_order;
    int * energy_to_pitzer_order_really_full;
    int * symmetry_energy_order;
    int * pitzer_offset;
    int * pitzer_offset_full;
    int * pitzer_to_energy_order;

    /// orbitals that were active in the reference computation 
    bool * active_reference_orbitals_;

    /// orbitals that were inactive in the reference computation 
    bool * inactive_reference_orbitals_;

    /// mapping from full set of orbitals to set excluding external orbitals
    int * full_to_no_ext_;

    /// the number of non-external orbitals
    int nmo_no_ext_;


    /// active-space geminals for each symmetry:
    std::vector < std::vector < std::pair<int,int> > > gems;

    /// total number of active molecular orbitals
    int amo_;

    /// total number of frozen core orbitals
    int nfrzc_;

    /// total number of frozen virtual orbitals
    int nfrzv_;

    /// total number of restricted doubly occupied orbitals
    int nrstc_;

    /// total number of restricted unoccupied orbitals
    int nrstv_;

    /// active molecular orbitals per irrep
    int * amopi_;

    /// restricted core orbitals per irrep.  these will be optimized optimized
    int * rstcpi_;

    /// restricted virtual orbitals per irrep.  these will be optimized optimized
    int * rstvpi_;

    /// number of auxilliary basis functions
    int nQ_;

    /// read three-index integrals and transform them to MO basis
    void ThreeIndexIntegrals();

    /// three-index integral buffer
    double * Qmo_;

    // mapping arrays with abelian symmetry
    void BuildBasis();
    int * full_basis;

    /// mapping arrays with symmetry:
    int * gems_erpa;
    int * gems_ab;
    int * gems_full;
    int * gems_plus_core;
    int *** bas_erpa_sym;
    int *** bas_ab_sym;
    int *** bas_full_sym;
    int *** bas_really_full_sym;
    int *** ibas_erpa_sym;
    int *** ibas_ab_sym;
    int *** ibas_full_sym;

    /// grab one- and two-electron integrals
    void GetIntegrals();

    /// read two-electron integrals from disk
    void GetTEIFromDisk();

    /// grab a specific two-electron integral
    double TEI(int i, int j, int k, int l, int h);

    /// SCF energy
    double escf_;

    /// nuclear repulsion energy
    double enuc_;

    // read teis from disk:
    void ReadIntegrals(double * tei,long int nmo);

    // multiplicity
    int multiplicity_;

    /// full space of integrals blocked by symmetry
    double * tei_full_sym_;
    double * oei_full_sym_;

    long int tei_full_dim_;
    int oei_full_dim_;

    /// compute frozen core energy and adjust oeis
    void FrozenCoreEnergy();

    /// are we using 3-index integrals?
    bool is_df_;

    /// reference rdm type (e.g. CI, CCSD, ...)
    std::string ref_rdm_;

    /// read v2RDM 1-RDM from disk
    void ReadOPDM(double * D1a, double * D1b);

    /// read v2RDM-CASSCF 2RDM from disk
    void ReadTPDM(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b);


    /// read v2RDM-CASSCF 2RDM from disk (low memory)
    void ReadTPDMLowMemory(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b);

    /// determine number of non external orbitals
    void CountNonExternalOrbitals();

    /// read full CI 2-RDM from disk
    void ReadTPDM_CI(double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b);

    /// read 2-RDM from disk (IWL format)
    void ReadTPDM_IWL(iwlbuf *Buf,double*d2,int nmo);


    /// read 3-RDM from disk
    void Read3PDM(double * D3aaa, double * D3aab, double * D3bba, double * D3bbb);

    /// buffers for dipole integrals
 
    SharedMatrix mux_moa_;
    SharedMatrix muy_moa_;
    SharedMatrix muz_moa_;
    SharedMatrix mux_mob_;
    SharedMatrix muy_mob_;
    SharedMatrix muz_mob_;
   
    /// symmetries of dipole integrals
    int hx_;
    int hy_;
    int hz_;

    double TransitionDipoleMoment(double * C, double * D1a, double * D1b, int nh, int h, std::string component);
    double TransitionDipoleMoment_SpinAdapted(double * C, double * D1a, double * D1b, int nh, int h, std::string component);

    void TransitionDensity(double * C, double * D1a, double * D1b, double * TDa, double * TDb, int nh, int h);
    void TransitionDensity_SpinAdapted(double * C, double * D1a, double * D1b, double * TDa, double * TDb, int nh, int h);
    void TransitionDensity_SpinAdaptedl(double * Cl, double * D1a, double * D1b, double * TDa, double * TDb, int nh, int h);

    double pruning_threshold;
    int number_coefficients; 
    double coefficient_threshold;
    
};

}}
#endif


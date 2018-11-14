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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<psi4/libmints/writer.h>
#include<psi4/libmints/writer_file_prefix.h>

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>

#include<psi4/libplugin/plugin.h>
#include<psi4/psi4-dec.h>
#include<psi4/liboptions/liboptions.h>

#include<psi4/libpsio/psio.hpp>
#include<psi4/libmints/wavefunction.h>
#include<psi4/psifiles.h>
#include<psi4/libmints/mintshelper.h>

#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>

#include<psi4/libtrans/integraltransform.h>
#include<psi4/libtrans/mospace.h>
#include<psi4/libiwl/iwl.h>

#include<time.h>

#include <psi4/libmints/molecule.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/basisset.h>

#include<psi4/libqt/qt.h>

#include "erpa_solver.h"
#include "fortran.h"

#include <psi4/libqt/lapack_intfc_mangle.h>
#include "blas.h"

#ifdef _OPENMP
    #include<omp.h>
#else
    #define omp_get_wtime() ( (double)clock() / CLOCKS_PER_SEC )
    #define omp_get_max_threads() 1
#endif

using namespace psi;
using namespace fnocc;
using namespace std;

long int INDEX999(long int n, long int i, long int j) {
  return i * n + j;
}

long int UIndex(long int i,long int j,long int k,long int l,long int ijtype,long int kltype,long int n) {
    return (i + n*j + ijtype * n*n) * 4L * n*n + (k + n*l + kltype * n*n);
}

namespace psi{ namespace erpa{

ERPASolver::ERPASolver(SharedWavefunction reference_wavefunction,Options & options):
    Wavefunction(options){
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

ERPASolver::~ERPASolver()
{
}

void  ERPASolver::common_init(){

    is_df_ = false;
    if ( options_.get_str("SCF_TYPE") == "DF" ) {
        is_df_ = true;
    }else if ( options_.get_str("SCF_TYPE") == "CD" ) {
        throw PsiException("plugin erpa does not yet work with cholesky-decomposed integrals",__FILE__,__LINE__);
    }

    ref_rdm_ = options_.get_str("REF_RDM");

    shallow_copy(reference_wavefunction_);

    escf_     = reference_wavefunction_->reference_energy();
    nalpha_   = reference_wavefunction_->nalpha();
    nbeta_    = reference_wavefunction_->nbeta();
    nalphapi_ = reference_wavefunction_->nalphapi();
    nbetapi_  = reference_wavefunction_->nbetapi();
    doccpi_   = reference_wavefunction_->doccpi();
    soccpi_   = reference_wavefunction_->soccpi();
    frzcpi_   = reference_wavefunction_->frzcpi();
    frzvpi_   = reference_wavefunction_->frzvpi();
    nmopi_    = reference_wavefunction_->nmopi();
    nirrep_   = reference_wavefunction_->nirrep();
    nso_      = reference_wavefunction_->nso();
    nmo_      = reference_wavefunction_->nmo();
    nsopi_    = reference_wavefunction_->nsopi();
    molecule_ = reference_wavefunction_->molecule();
    enuc_     = molecule_->nuclear_repulsion_energy({0.0,0.0,0.0});

    // need somewhere to store gradient, if required
    gradient_ =  reference_wavefunction_->matrix_factory()->create_shared_matrix("Total gradient", molecule_->natom(), 3);

    // restricted doubly occupied orbitals per irrep (optimized)
    rstcpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstcpi_,'\0',nirrep_*sizeof(int));

    // restricted unoccupied occupied orbitals per irrep (optimized)
    rstvpi_   = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)rstvpi_,'\0',nirrep_*sizeof(int));

    // active orbitals per irrep:
    amopi_    = (int*)malloc(nirrep_*sizeof(int));
    memset((void*)amopi_,'\0',nirrep_*sizeof(int));

    // multiplicity:
    multiplicity_ = Process::environment.molecule()->multiplicity();

    if (options_["FROZEN_DOCC"].has_changed()) {
        throw PsiException("FROZEN_DOCC is currently disabled.",__FILE__,__LINE__);

        if (options_["FROZEN_DOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_double();
        }
    }
    // when restarting from casscf, the number of frozen core orbitals is wrong
    for (int h = 0; h < nirrep_; h++) {
        if ( frzcpi_[h] > 0 ) {
            outfile->Printf("    <<< Warning >>> Reducing number of frozen core orbitals in irrep %5i to zero.\n",h);
            outfile->Printf("                    Check that your active space is correct.\n");
        }
        frzcpi_[h] = 0;
        if ( frzvpi_[h] > 0 ) {
            outfile->Printf("    <<< Warning >>> Reducing number of frozen virtual obitals in irrep %5i to zero.\n",h);
            outfile->Printf("                    Check that your active space is correct.\n");
        }
        frzvpi_[h] = 0;
    }
    if (options_["RESTRICTED_DOCC"].has_changed()) {
        if (options_["RESTRICTED_DOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_DOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstcpi_[h] = options_["RESTRICTED_DOCC"][h].to_double();
        }
    }
    if (options_["RESTRICTED_UOCC"].has_changed()) {
        if (options_["RESTRICTED_UOCC"].size() != nirrep_) {
            throw PsiException("The RESTRICTED_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            rstvpi_[h] = options_["RESTRICTED_UOCC"][h].to_double();
        }
    }
    if (options_["FROZEN_UOCC"].has_changed()) {

        //if ( !is_df_ ) {
        //    throw PsiException("FROZEN_UOCC is currently enabled only for SCF_TYPE CD and DF.",__FILE__,__LINE__);
        //}
        if (options_["FROZEN_UOCC"].size() != nirrep_) {
            throw PsiException("The FROZEN_UOCC array has the wrong dimensions_",__FILE__,__LINE__);
        }
        for (int h = 0; h < nirrep_; h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_double();
        }
    }

    // user could specify active space with ACTIVE array
    if ( options_["ACTIVE"].has_changed() ) {
        //throw PsiException("The ACTIVE array is not yet enabled.",__FILE__,__LINE__);
        if (options_["ACTIVE"].size() != nirrep_) {
            throw PsiException("The ACTIVE array has the wrong dimensions_",__FILE__,__LINE__);
        }

        // warn user that active array takes precedence over restricted_uocc array
        if (options_["RESTRICTED_UOCC"].has_changed()) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING!! >>>\n");
            outfile->Printf("\n");
            outfile->Printf("    The ACTIVE array takes precedence over the RESTRICTED_UOCC array.\n");
            outfile->Printf("    Check below whether your active space was correctly specified.\n");
            outfile->Printf("\n");
        }

        // overwrite rstvpi_ array using the information in the frozen_docc,
        // restricted_docc, active, and frozen_virtual arrays.  start with nso total
        // orbitals and let the linear dependency check below adjust the spaces as needed
        for (int h = 0; h < nirrep_; h++) {
            amopi_[h]  = options_["ACTIVE"][h].to_double();
            rstvpi_[h] = nsopi_[h] - frzcpi_[h] - rstcpi_[h] - frzvpi_[h] - amopi_[h];
        }
    }

    // were there linear dependencies in the primary basis set?
    if ( nmo_ != nso_ ) {

        // which irreps lost orbitals?
        int * lost = (int*)malloc(nirrep_*sizeof(int));
        memset((void*)lost,'\0',nirrep_*sizeof(int));
        bool active_space_changed = false;
        for (int h = 0; h < factory_->nirrep(); h++){
            lost[h] = nsopi_[h] - nmopi_[h];
            if ( lost[h] > 0 ) {
                active_space_changed = true;
            }

            // eliminate frozen virtual orbitals first
            if ( frzvpi_[h] > 0 && lost[h] > 0 ) {
                frzvpi_[h] -= ( frzvpi_[h] < lost[h] ? frzvpi_[h] : lost[h] );
                lost[h]    -= ( frzvpi_[h] < lost[h] ? frzvpi_[h] : lost[h] );
            }
            // if necessary, eliminate restricted virtual orbitals next
            if ( rstvpi_[h] > 0 && lost[h] > 0 ) {
                rstvpi_[h] -= ( rstvpi_[h] < lost[h] ? rstvpi_[h] : lost[h] );
            }
        }
        if ( active_space_changed ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING!! >>>\n");
            outfile->Printf("\n");
            outfile->Printf("    Your basis set may have linear dependencies.\n");
            outfile->Printf("    The number of restricted or frozen virtual orbitals per irrep may have changed.\n");
            outfile->Printf("\n");
            outfile->Printf("    No. orbitals removed per irrep: [");
            for (int h = 0; h < nirrep_; h++)
                outfile->Printf("%4i",nsopi_[h] - nmopi_[h]);
            outfile->Printf(" ]\n");
            //outfile->Printf("    No. frozen virtuals per irrep:  [");
            //for (int h = 0; h < nirrep_; h++)
            //    outfile->Printf("%4i",frzvpi_[h]);
            //outfile->Printf(" ]\n");
            //outfile->Printf("\n");
            outfile->Printf("    Check that your active space is still correct.\n");
            outfile->Printf("\n");
        }
    }


    Ca_ = SharedMatrix(reference_wavefunction_->Ca());
    Cb_ = SharedMatrix(reference_wavefunction_->Cb());

    S_  = SharedMatrix(reference_wavefunction_->S());

    Fa_ = SharedMatrix(reference_wavefunction_->Fa());
    Fb_ = SharedMatrix(reference_wavefunction_->Fb());

    Da_ = SharedMatrix(reference_wavefunction_->Da());
    Db_ = SharedMatrix(reference_wavefunction_->Db());

    // Lagrangian matrix
    Lagrangian_ = SharedMatrix(reference_wavefunction_->Lagrangian());

    epsilon_a_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_a_->copy(reference_wavefunction_->epsilon_a().get());
    epsilon_b_= std::shared_ptr<Vector>(new Vector(nirrep_, nmopi_));
    epsilon_b_->copy(reference_wavefunction_->epsilon_b().get());

    amo_      = 0;
    nfrzc_    = 0;
    nfrzv_    = 0;
    nrstc_    = 0;
    nrstv_    = 0;

    int ndocc = 0;
    int nvirt = 0;
    for (int h = 0; h < nirrep_; h++){
        nfrzc_   += frzcpi_[h];
        nrstc_   += rstcpi_[h];
        nrstv_   += rstvpi_[h];
        nfrzv_   += frzvpi_[h];
        amo_   += nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
        ndocc    += doccpi_[h];
        amopi_[h] = nmopi_[h]-frzcpi_[h]-rstcpi_[h]-rstvpi_[h]-frzvpi_[h];
    }

    int ndoccact = ndocc - nfrzc_ - nrstc_;
    nvirt    = amo_ - ndoccact;

    // sanity check for orbital occupancies:
    for (int h = 0; h < nirrep_; h++) {
        int tot = doccpi_[h] + soccpi_[h] + rstvpi_[h] + frzvpi_[h];
        if (doccpi_[h] + soccpi_[h] + rstvpi_[h] + frzvpi_[h] > nmopi_[h] ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING >>> irrep %5i has too many orbitals:\n",h);
            outfile->Printf("\n");
            outfile->Printf("                    docc = %5i\n",doccpi_[h]);
            outfile->Printf("                    socc = %5i\n",soccpi_[h]);
            outfile->Printf("                    rstu = %5i\n",rstvpi_[h]);
            outfile->Printf("                    frzv = %5i\n",frzvpi_[h]);
            outfile->Printf("                    tot  = %5i\n",doccpi_[h] + soccpi_[h] + rstvpi_[h] + frzvpi_[h]);
            outfile->Printf("\n");
            outfile->Printf("                    total no. orbitals should be %5i\n",nmopi_[h]);
            outfile->Printf("\n");
            throw PsiException("at least one irrep has too many orbitals",__FILE__,__LINE__);
        }
        if (frzcpi_[h] + rstcpi_[h] > doccpi_[h] ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< WARNING >>> irrep %5i has too many frozen and restricted core orbitals:\n",h);
            outfile->Printf("                    frzc = %5i\n",frzcpi_[h]);
            outfile->Printf("                    rstd = %5i\n",rstcpi_[h]);
            outfile->Printf("                    docc = %5i\n",doccpi_[h]);
            outfile->Printf("\n");
            throw PsiException("at least one irrep has too many frozen core orbitals",__FILE__,__LINE__);
        }
    }

    // memory is from process::environment
    memory_ = Process::environment.get_memory();

    // set the wavefunction name
    name_ = "ERPA";

    // build mapping arrays and determine the number of geminals per block
    BuildBasis();

    int ms = (multiplicity_ - 1)/2;

    // print orbitals per irrep in each space
    outfile->Printf("  ==> Active space details <==\n");
    outfile->Printf("\n");
    //outfile->Printf("        Freeze core orbitals?                   %5s\n",nfrzc_ > 0 ? "yes" : "no");
    outfile->Printf("        Number of frozen core orbitals:         %5i\n",nfrzc_);
    outfile->Printf("        Number of restricted occupied orbitals: %5i\n",nrstc_);
    outfile->Printf("        Number of active occupied orbitals:     %5i\n",ndoccact);
    outfile->Printf("        Number of active virtual orbitals:      %5i\n",nvirt);
    outfile->Printf("        Number of restricted virtual orbitals:  %5i\n",nrstv_);
    outfile->Printf("        Number of frozen virtual orbitals:      %5i\n",nfrzv_);
    outfile->Printf("\n");

    std::vector<std::string> labels = reference_wavefunction_->molecule()->irrep_labels();
    outfile->Printf("        Irrep:           ");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4s",labels[h].c_str());
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" \n");
    outfile->Printf(" \n");

    outfile->Printf("        frozen_docc     [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",frzcpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        restricted_docc [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",rstcpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        active          [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",amopi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        restricted_uocc [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",rstvpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("        frozen_uocc     [");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("%4i",frzvpi_[h]);
        if ( h < nirrep_ - 1 ) {
            outfile->Printf(",");
        }
    }
    outfile->Printf(" ]\n");
    outfile->Printf("\n");

    
    // eigevalue threshold below which A matrix is pruned 
    pruning_threshold = options_.get_double("PRUNING_THRESHOLD");
    coefficient_threshold = options_.get_double("COEFFICIENT_THRESHOLD");
    number_coefficients = options_.get_double("NUMBER_COEFFICIENTS");
    
    outfile->Printf("  ==> Eigenvalue threshold below which A matrix is pruned: %3.2le\n", pruning_threshold);
    outfile->Printf("  ==> Number of expansion coeffients to be printed if all are above threshold: %5i\n", number_coefficients);
    outfile->Printf("  ==> Threshold above which excitation coefficients are printed: %3.2le\n", coefficient_threshold);
    outfile->Printf("\n");
  
    // TODO
    //outfile->Printf("\n");
    //outfile->Printf("  ==> Memory requirements <==\n");
    //outfile->Printf("\n");

    long int nn1o2 = nmo_*(nmo_+1)/2;

    if ( is_df_ ) {
        // storage requirements for df integrals
        nQ_ = Process::environment.globals["NAUX (SCF)"];
        if ( options_.get_str("SCF_TYPE") == "DF" ) {
            std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");
            nQ_ = auxiliary->nbf();
            Process::environment.globals["NAUX (SCF)"] = nQ_;
        }
    }else {
        // four-index integrals will be factorized using SVD
        nQ_ = nn1o2;
    }

    long int maxgem = 0;
    for (int h = 0; h < nirrep_; h++) {
        if ( gems_ab[h] > maxgem ) maxgem = gems_ab[h];
    }
    double required_memory   = 4.0 * nmo_*nmo_;
    required_memory         += 3.0 * nmo_*nmo_*nmo_*nmo_;
    required_memory         += ((1.0 * nn1o2*nn1o2) > (2.0 * nmo_*nmo_*nQ_) ? (1.0 * nn1o2*nn1o2) : (2.0 * nmo_*nmo_*nQ_));
    required_memory         += 7.0 * 4 * maxgem * maxgem;
    required_memory         += 2.0 * nQ_ * maxgem;

    outfile->Printf("\n");
    outfile->Printf("        Total memory requirements:     %7.2lf mb\n",required_memory * 8.0 / 1024.0 / 1024.0);
    outfile->Printf("        (note this estimate might be wrong when using four-index integrals ...)\n");
    outfile->Printf("\n");

    if ( required_memory * 8.0 > (double)memory_ ) {
        outfile->Printf("\n");
        outfile->Printf("        Not enough memory!\n");
        outfile->Printf("\n");
        //throw PsiException("Not enough memory",__FILE__,__LINE__);
    }

    // if using 3-index integrals, transform them before allocating any memory integrals, transform
    if ( is_df_ ) {
        outfile->Printf("    ==> Transform three-electron integrals <==\n");
        outfile->Printf("\n");

        double start = omp_get_wtime();
        ThreeIndexIntegrals();
        double end = omp_get_wtime();

        outfile->Printf("\n");
        outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        outfile->Printf("\n");
    } else {
        // transform integrals
        outfile->Printf("    ==> Transform two-electron integrals <==\n");
        outfile->Printf("\n");

        double start = omp_get_wtime();
        std::vector<shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::all);
        std::shared_ptr<IntegralTransform> ints(new IntegralTransform(reference_wavefunction_, spaces, IntegralTransform::TransformationType::Restricted,
            				      IntegralTransform::OutputType::IWLOnly, IntegralTransform::MOOrdering::PitzerOrder, IntegralTransform::FrozenOrbitals::None, false));
        ints->set_dpd_id(0);
        ints->set_keep_iwl_so_ints(true);
        ints->set_keep_dpd_so_ints(true);
        ints->initialize();
        ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
        double end = omp_get_wtime();
        outfile->Printf("\n");
        outfile->Printf("        Time for integral transformation:  %7.2lf s\n",end-start);
        outfile->Printf("\n");

    }

    GetIntegrals();

    // even if we use rhf/rohf reference, we need same_a_b_orbs_=false
    // to trigger the correct integral transformations in deriv.cc
    same_a_b_orbs_ = false;
    same_a_b_dens_ = false;

    // dipole integrals in symmetry orbital basis

    std::shared_ptr<MintsHelper> mints(new MintsHelper(reference_wavefunction_));
    std::vector<std::shared_ptr<Matrix> > dipole = mints->so_dipole();

    // Transforming dipole integrals from the symmetry orbital to the spin orbital basis
    
    mux_moa_ = std::shared_ptr<Matrix>(new Matrix(dipole[0]));
    muy_moa_ = std::shared_ptr<Matrix>(new Matrix(dipole[1]));
    muz_moa_ = std::shared_ptr<Matrix>(new Matrix(dipole[2]));
    mux_mob_ = std::shared_ptr<Matrix>(new Matrix(dipole[0]));
    muy_mob_ = std::shared_ptr<Matrix>(new Matrix(dipole[1]));
    muz_mob_ = std::shared_ptr<Matrix>(new Matrix(dipole[2]));

    // Getting the symmetry info. of the dipole integrals
    hx_ = dipole[0]->symmetry();
    hy_ = dipole[1]->symmetry();
    hz_ = dipole[2]->symmetry();

    for (int h = 0; h < nirrep_; h++) {
         
        int h2 = SymmetryPair(h,hx_);
        for (int i = 0; i < nmopi_[h]; i++) {
            for (int j = 0; j < nmopi_[h2]; j++) {
                double duma = 0.0;
                double dumb = 0.0;
                for (int mu = 0; mu < nsopi_[h]; mu++) {
                    for (int nu = 0; nu < nsopi_[h2]; nu++) {
                        duma += dipole[0]->pointer(h)[mu][nu] * Ca_->pointer(h)[mu][i] * Ca_->pointer(h2)[nu][j];
                        dumb += dipole[0]->pointer(h)[mu][nu] * Cb_->pointer(h)[mu][i] * Cb_->pointer(h2)[nu][j];
                    }
                }
                mux_moa_->pointer(h)[i][j] = duma;
                mux_mob_->pointer(h)[i][j] = dumb;
            }
        }

        h2 = SymmetryPair(h,hy_);
        for (int i = 0; i < nmopi_[h]; i++) {
            for (int j = 0; j < nmopi_[h2]; j++) {
                double duma = 0.0;
                double dumb = 0.0;
                for (int mu = 0; mu < nmopi_[h]; mu++) {
                    for (int nu = 0; nu < nmopi_[h2]; nu++) {
                        duma += dipole[1]->pointer(h)[mu][nu] * Ca_->pointer(h)[mu][i] * Ca_->pointer(h2)[nu][j];
                        dumb += dipole[1]->pointer(h)[mu][nu] * Cb_->pointer(h)[mu][i] * Cb_->pointer(h2)[nu][j];
                    }
                }
                muy_moa_->pointer(h)[i][j] = duma;
                muy_mob_->pointer(h)[i][j] = dumb;
            }
        }

        h2 = SymmetryPair(h,hz_);
        for (int i = 0; i < nmopi_[h]; i++) {
            for (int j = 0; j < nmopi_[h2]; j++) {
                double duma = 0.0;
                double dumb = 0.0;
                for (int mu = 0; mu < nmopi_[h]; mu++) {
                    for (int nu = 0; nu < nmopi_[h2]; nu++) {
                        duma += dipole[2]->pointer(h)[mu][nu] * Ca_->pointer(h)[mu][i] * Ca_->pointer(h2)[nu][j];
                        dumb += dipole[2]->pointer(h)[mu][nu] * Cb_->pointer(h)[mu][i] * Cb_->pointer(h2)[nu][j];
                    }
                }
                muz_moa_->pointer(h)[i][j] = duma;
                muz_mob_->pointer(h)[i][j] = dumb;
            }
        }
    }

}

int ERPASolver::SymmetryPair(int i,int j) {
    return table[i*8+j];
}

int ERPASolver::TotalSym(int i,int j,int k, int l) {
    return SymmetryPair(SymmetryPair(symmetry[i],symmetry[j]),SymmetryPair(symmetry[k],symmetry[l]));
}

void ERPASolver::ExtendedRPA(std::string type) {
    throw PsiException("this is dead code. only a matter of time before it is removed.",__FILE__,__LINE__);
}


// variant 2 of the extended RPA as described in
// Chatterjee and Pernal, JCP 137, 204109 (2012)

void ERPASolver::ExtendedRPA2() {

    throw PsiException("ERPA2 is not implemented.",__FILE__,__LINE__);
/*

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals
    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // build two-electron integrals

    // note from AED (9/11/17): since the ERPA code below was designed for use 
    // with 3-index integrals, I had to play some tricks to make it work with
    // 4-index integrals.  If you want to use 4-index integrals (which is a bad
    // idea), the code will perform an SVD of the ERI tensor to yield factorized 
    // integrals that are appropriate for the 3-index code.  The factorization 
    // isn't symmetric, so there are now two buffers for 3-index integrals: Qmo_ 
    // and Qmo2.  In the case of DF, Qmo2 = Qmo_.  I say using the 4-index code
    // is a bad idea because it requires more storage.

    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));

    int info;
    double * Qmo2;
    if ( is_df_ ) {
        F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
        Qmo2 = Qmo_;
    }else { 
        // unpack two-electron integrals (which were initially blocked by symmetry)
        long int offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        // now, compute SVD of tei_full_sym_ and store result in Qmo_, Qmo2
        double * U  = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * VT = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * S  = (double*)malloc(nn1o2*sizeof(double));
        memset((void*)U,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)VT,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)S,'\0',nn1o2*sizeof(double));

        Qmo_ = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        memset((void*)Qmo_,'\0',nn1o2*nn1o2*sizeof(double));

        int lwork      = 5 * nn1o2;
        double * work  = (double*)malloc(lwork*sizeof(double));

        //info = C_DGESVD('A','A',nn1o2,nn1o2,tei,nn1o2,S,VT,nn1o2,U,nn1o2,work,lwork);
        char jobu  = 'A';
        char jobvt = 'A';
        int info;
        int dim = (int)nn1o2;
        info = F_DGESVD(&jobvt,&jobu,&dim,&dim,tei,&dim,S,U,&dim,VT,&dim,work,&lwork,&info);
        if ( info != 0 ) {
            throw PsiException("something went wrong in the SVD of the 4-index ERI tensor.",__FILE__,__LINE__);
        }

        // repack two-electron integrals (which were destroyed by SVD)
        memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));
        offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        double err = 0.0;
        for (int Q = 0; Q < nn1o2; Q++) {
            for (int kl = 0; kl < nn1o2; kl++) {
                Qmo_[kl*nn1o2 + Q] = U[Q*nn1o2 + kl] * S[Q];
            }
        }
        free(U);
        free(S);
        Qmo2 = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        C_DCOPY(nn1o2*nn1o2,VT,1,Qmo2,1);
        free(VT);

        //for (int ij = 0; ij < nn1o2; ij++) {
        //    for (int kl = 0; kl < nn1o2; kl++) {
        //        double dum = 0.0;
        //        for (int Q = 0; Q < nn1o2; Q++) {
        //            //dum += U[ij*nn1o2 + Q] * VT[kl*nn1o2 + Q] * S[Q];
        //            //dum += U[Q*nn1o2 + ij] * VT[Q*nn1o2 + kl] * S[Q];
        //            //dum += VT[ij*nn1o2 + Q] * U[kl*nn1o2 + Q] * S[Q];
        //            //dum += VT[Q*nn1o2 + ij] * U[Q*nn1o2 + kl] * S[Q];

        //            //dum += U[ij*nn1o2 + Q] * VT[Q*nn1o2 + kl] * S[Q];
        //            // this one:
        //            //dum += VT[ij*nn1o2 + Q] * U[Q*nn1o2 + kl] * S[Q];
        //            //dum += U[Q*nn1o2 + ij] * VT[Q*nn1o2 + kl] * S[Q];

        //            dum += VT[ij*nn1o2 + Q] * Qmo_[kl*nn1o2 + Q];
        //        }
        //        double diff = fabs(tei[ij*nn1o2+kl] - dum);
        //        if ( diff > 1e-6 ) {
        //            printf("%5i %5i %20.12lf %20.12lf\n",ij,kl,tei[ij*nn1o2+kl],dum);
        //        }
        //        err += diff*diff;
        //    }
        //}
        //err = sqrt(err);
        //printf("%20.12lf\n",err);
        nQ_ = nn1o2;
        F_DGEMM('t','n',nn1o2,nn1o2,nQ_,1.0,Qmo_,nQ_,Qmo2,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
    }

    // read TDPM from disk
    ReadTPDM(D2aa,D2bb,D2ab,D1a,D1b,tei);

    // build intermediates for funky 3-index sums:
    double * tempa  = (double*)malloc(nmo_*nmo_*nQ_*sizeof(double));
    double * tempb  = (double*)malloc(nmo_*nmo_*nQ_*sizeof(double));
    double * Ia = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * Ib = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)tempa,'\0',nmo_*nmo_*nQ_*sizeof(double));
    memset((void*)tempb,'\0',nmo_*nmo_*nQ_*sizeof(double));
    memset((void*)Ia,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Ib,'\0',nmo_*nmo_*sizeof(double));

    for (int l = 0; l < amo_; l++) {

        int hl = symmetry[l];

        for (int p = 0; p < amo_; p++) {

            int hp = symmetry[p];

            int hpl = hp ^ hl;

            for (int Q = 0; Q < nQ_; Q++) {

                double duma = 0.0;
                double dumb = 0.0;

                for (int qs = 0; qs < gems_ab[hpl]; qs++) {

                    long int q = bas_ab_sym[hpl][qs][0];
                    long int s = bas_ab_sym[hpl][qs][1];

                    long int pq = p*nmo_+q;
                    long int ls = l*nmo_+s;

                    long int qp = q*nmo_+p;
                    long int sl = s*nmo_+l;

                    duma += Qmo_[INDEX(q,s)*nQ_+Q] * ( D2aa[pq*nmo_*nmo_+ls] + D2ab[pq*nmo_*nmo_+ls] );
                    dumb += Qmo_[INDEX(q,s)*nQ_+Q] * ( D2bb[pq*nmo_*nmo_+ls] + D2ab[qp*nmo_*nmo_+sl] );

                }

                tempa[l*nmo_*nQ_ + p*nQ_ + Q] = duma;
                tempb[l*nmo_*nQ_ + p*nQ_ + Q] = dumb;
            }
        }
    }
    for (int j = 0; j < amo_; j++) {
        for (int l = 0; l < amo_; l++) {
            double duma = 0.0;
            double dumb = 0.0;
            for (int p = 0; p < amo_; p++) {
                duma += C_DDOT(nQ_,tempa + l*nmo_*nQ_ + p*nQ_,1,Qmo2 + INDEX(p,j)*nQ_,1);
                dumb += C_DDOT(nQ_,tempb + l*nmo_*nQ_ + p*nQ_,1,Qmo2 + INDEX(p,j)*nQ_,1);
            }
            Ia[j*nmo_+l] = duma;
            Ib[j*nmo_+l] = dumb;
        }
    }
    free(tempa);
    free(tempb);

    long int n = nmo_*nmo_;

    // check energy
    double e2 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2 +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];
                }
            }
        }
    }
    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }

    printf("%20.12lf %20.12lf %20.12lf %20.12lf\n",e1,e2,e1+e2,e1+e2+enuc_);
    //printf("%5i\n",info);
    exit(0);


    // build and diagonalize ERPA2 matrix ... note that we are
    // only concerned with the totally symmetric irrep in ERPA2 
    for (int h = 0; h < 1; h++) {

        long int nh  = gems_erpa[h];
        //long int nh2 = gems_erpa[h];
        //long int nh2 = gems_ab[h];
        //long int nh2 = gems_aa[h];
        long int nh2 = gems_00[h];
        //if ( nh == 0 && nh2 == 0 ) continue;

        // dimension: (ca, cb, d) => 3 nh
        long int totdim = 2 * nh + nh2;

        double * newB = (double*)malloc(totdim*totdim*sizeof(double));
        double * newA = (double*)malloc(totdim*totdim*sizeof(double));
        double * cc   = (double*)malloc(totdim*totdim*sizeof(double));
        double * eig  = (double*)malloc(totdim*sizeof(double));

        memset((void*)newA,'\0',totdim*totdim*sizeof(double));
        memset((void*)newB,'\0',totdim*totdim*sizeof(double));
        memset((void*)cc,'\0',totdim*totdim*sizeof(double));
        memset((void*)eig,'\0',totdim*sizeof(double));

        outfile->Printf("\n");
        outfile->Printf("    Symmetry: %5i\n",h);
        outfile->Printf("\n");

        // build B parts specific to erpa 2
        for (long int ij = 0; ij < nh2; ij++) {
//continue;

            long int i = bas_00_sym[h][ij][0];
            long int j = bas_00_sym[h][ij][1];
            //long int i = bas_erpa_sym[h][ij][0];
            //long int j = bas_erpa_sym[h][ij][1];

            long int ij1 = ibas_erpa_sym[h][i][j];

            for (long int kl = 0; kl < nh2; kl++) {

                long int k = bas_00_sym[h][kl][0];
                long int l = bas_00_sym[h][kl][1];
                //long int k = bas_erpa_sym[h][kl][0];
                //long int l = bas_erpa_sym[h][kl][1];

                long int kl1 = ibas_erpa_sym[h][k][l];

                double dum_a_d = 0.0;
                double dum_b_d = 0.0;
                double dum_d_a = 0.0;
                double dum_d_b = 0.0;
                double dum_d_d = 0.0;

                if ( j == l ) {

                    dum_a_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+l];
                    dum_b_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+l*nmo_+i];
                    dum_d_a += D2ab[k*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+i];
                    dum_d_b += D2ab[j*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+i];
                    dum_d_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+i];

                }
                if ( i == k ) {

                    dum_a_d -= D2ab[j*nmo_*nmo_*nmo_+k*nmo_*nmo_+l*nmo_+l];
                    dum_b_d -= D2ab[k*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+l];
                    dum_d_a -= D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+i];
                    dum_d_b -= D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+l];
                    dum_d_d -= D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+l];

                }

                long int a_d_sym = (ij1 + 0 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                long int b_d_sym = (ij1 + 1 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                long int d_a_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl1 + 0 * nh);
                long int d_b_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl1 + 1 * nh);
                long int d_d_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                    
                if ( i != j ) newB[a_d_sym] = dum_a_d;
                if ( i != j ) newB[b_d_sym] = dum_b_d;
                if ( k != l ) newB[d_a_sym] = dum_d_a;
                if ( k != l ) newB[d_b_sym] = dum_d_b;
                newB[d_d_sym] = dum_d_d;

//if ( i!=j ) printf("aa-d %5i\n",a_d_sym);
//if ( i!=j ) printf("bb-d %5i\n",b_d_sym);
//if ( k!=l ) printf("d-aa %5i\n",d_a_sym);
//if ( k!=l ) printf("d-bb %5i\n",d_b_sym);
//printf("d-d  %5i\n",d_d_sym);fflush(stdout);

            }
        }

        // build A1 parts specific to erpa 2
        for (long int ij = 0; ij < nh2; ij++) {
//continue;

            long int i = bas_00_sym[h][ij][0];
            long int j = bas_00_sym[h][ij][1];
            //long int i = bas_erpa_sym[h][ij][0];
            //long int j = bas_erpa_sym[h][ij][1];

            long int ij1 = ibas_erpa_sym[h][i][j];

            for (long int kl = 0; kl < nh2; kl++) {

                long int k = bas_00_sym[h][kl][0];
                long int l = bas_00_sym[h][kl][1];
                //long int k = bas_erpa_sym[h][kl][0];
                //long int l = bas_erpa_sym[h][kl][1];

                long int kl1 = ibas_erpa_sym[h][k][l];

                double dum_a_d = 0.0;
                double dum_b_d = 0.0;
                double dum_d_a = 0.0;
                double dum_d_b = 0.0;
                double dum_d_d = 0.0;

                dum_a_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+l] * h1[l*nmo_+j] + D2ab[j*nmo_*nmo_*nmo_+k*nmo_*nmo_+l*nmo_+l] * h1[i*nmo_+k];
                dum_b_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+l*nmo_+i] * h1[l*nmo_+j] + D2ab[k*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+l] * h1[i*nmo_+k];

                dum_d_a += D2ab[k*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+i] * h1[l*nmo_+j] + D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+i] * h1[i*nmo_+k];
                dum_d_b += D2ab[j*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+i] * h1[l*nmo_+j] + D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+l] * h1[i*nmo_+k];

                if ( j == l ) {

                    dum_d_d += 2.0 * D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+i] * h1[l*nmo_+j];

                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];

                            dum_a_d -= D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+q*nmo_+l] * h1[i*nmo_+q];

                            dum_b_d -= D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+l*nmo_+q] * h1[i*nmo_+q];

                            dum_d_a -= D2ab[k*nmo_*nmo_*nmo_+j*nmo_*nmo_+q*nmo_+i] * h1[i*nmo_+q];
                            dum_d_a -= D2ab[k*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+q] * h1[i*nmo_+q];
                            dum_d_a += D2ab[k*nmo_*nmo_*nmo_+q*nmo_*nmo_+i*nmo_+i] * h1[q*nmo_+j]; // "q" -> "p" in aed notes

                            dum_d_b -= D2ab[j*nmo_*nmo_*nmo_+k*nmo_*nmo_+q*nmo_+i] * h1[i*nmo_+q];
                            dum_d_b -= D2ab[j*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+q] * h1[i*nmo_+q];
                            dum_d_b += D2ab[q*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+i] * h1[q*nmo_+j]; // "q" -> "p" in aed notes

                            dum_d_d -= D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+q] * h1[i*nmo_+q];
                            dum_d_d -= D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+q*nmo_+i] * h1[i*nmo_+q];
                        }
                    }

                }
                if ( i == k ) {

                    dum_d_d += 2.0 * D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+l] * h1[i*nmo_+k];

                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];

                            dum_a_d -= D2ab[p*nmo_*nmo_*nmo_+k*nmo_*nmo_+l*nmo_+l] * h1[p*nmo_+j];

                            dum_b_d -= D2ab[k*nmo_*nmo_*nmo_+p*nmo_*nmo_+l*nmo_+l] * h1[p*nmo_+j];

                            dum_d_a -= D2ab[p*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+i] * h1[p*nmo_+j];
                            dum_d_a -= D2ab[j*nmo_*nmo_*nmo_+p*nmo_*nmo_+l*nmo_+i] * h1[p*nmo_+j];
                            dum_d_a += D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+p] * h1[i*nmo_+p]; // "p" -> "q" in aed notes

                            dum_d_b -= D2ab[p*nmo_*nmo_*nmo_+j*nmo_*nmo_+i*nmo_+l] * h1[p*nmo_+j];
                            dum_d_b -= D2ab[j*nmo_*nmo_*nmo_+p*nmo_*nmo_+i*nmo_+l] * h1[p*nmo_+j];
                            dum_d_b += D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+p*nmo_+l] * h1[i*nmo_+p]; // "p" -> "q" in aed notes

                            dum_d_d -= D2ab[p*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+l] * h1[p*nmo_+j];
                            dum_d_d -= D2ab[j*nmo_*nmo_*nmo_+p*nmo_*nmo_+l*nmo_+l] * h1[p*nmo_+j];
                        }
                    }
                }

                long int a_d_sym = (ij1 + 0 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                long int b_d_sym = (ij1 + 1 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                long int d_a_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl1 + 0 * nh);
                long int d_b_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl1 + 1 * nh);
                long int d_d_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                    
                if ( i != j ) newA[a_d_sym] = dum_a_d;
                if ( i != j ) newA[b_d_sym] = dum_b_d;
                if ( k != l ) newA[d_a_sym] = dum_d_a;
                if ( k != l ) newA[d_b_sym] = dum_d_b;
                newA[d_d_sym] = dum_d_d;

            }
        }

        // build A2 parts specific to erpa 2
        for (long int ij = 0; ij < nh2; ij++) {
//continue;

            long int i = bas_00_sym[h][ij][0];
            long int j = bas_00_sym[h][ij][1];
            //long int i = bas_erpa_sym[h][ij][0];
            //long int j = bas_erpa_sym[h][ij][1];

            long int ij1 = ibas_erpa_sym[h][i][j];

            for (long int kl = 0; kl < nh2; kl++) {

                long int k = bas_00_sym[h][kl][0];
                long int l = bas_00_sym[h][kl][1];
                //long int k = bas_erpa_sym[h][kl][0];
                //long int l = bas_erpa_sym[h][kl][1];

                long int kl1 = ibas_erpa_sym[h][k][l];

                double dum_a_d = 0.0;
                double dum_b_d = 0.0;
                double dum_d_a = 0.0;
                double dum_d_b = 0.0;
                double dum_d_d = 0.0;

                // a,d

                for (int hs = 0; hs < nirrep_; hs++) {
                    for (int ss = 0; ss < amopi_[hs]; ss++) {
                        long int s = ss + pitzer_offset_full[hs];
                        double e_ljls = tei[INDEX(l,j)*nn1o2+INDEX(l,s)];
                        dum_a_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+s] * e_ljls;
                    }
                }
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        double e_ikqk = tei[INDEX(i,k)*nn1o2+INDEX(q,k)];
                        dum_a_d += D2ab[j*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+l] * e_ikqk;
                    }
                }
                if ( j == l ) {
                    for (int hr = 0; hr < nirrep_; hr++) {
                        for (int rr = 0; rr < amopi_[hr]; rr++) {
                            long int r = rr + pitzer_offset_full[hr];
                            for (int hs = 0; hs < nirrep_; hs++) {
                                for (int ss = 0; ss < amopi_[hs]; ss++) {
                                    long int s = ss + pitzer_offset_full[hs];
                                    double e_irls = tei[INDEX(i,r)*nn1o2+INDEX(l,s)];
                                    dum_a_d -= D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+r*nmo_+s] * e_irls;
                                }
                            }
                        }
                    }
                }
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    double e_pjqk = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                                    dum_a_d -= D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+l] * e_pjqk;
                                }
                            }
                        }
                    }
                }

                // b,d

                for (int hr = 0; hr < nirrep_; hr++) {
                    for (int rr = 0; rr < amopi_[hr]; rr++) {
                        long int r = rr + pitzer_offset_full[hr];
                        double e_lrlj = tei[INDEX(l,r)*nn1o2+INDEX(l,j)];
                        dum_b_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+r*nmo_+i] * e_lrlj;
                    }
                }
                for (int hp = 0; hp < nirrep_; hp++) {
                    for (int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        double e_pkik = tei[INDEX(p,k)*nn1o2+INDEX(i,k)];
                        dum_b_d += D2ab[p*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+l] * e_pkik;
                    }
                }
                if ( j == l ) {
                    for (int hr = 0; hr < nirrep_; hr++) {
                        for (int rr = 0; rr < amopi_[hr]; rr++) {
                            long int r = rr + pitzer_offset_full[hr];
                            for (int hs = 0; hs < nirrep_; hs++) {
                                for (int ss = 0; ss < amopi_[hs]; ss++) {
                                    long int s = ss + pitzer_offset_full[hs];
                                    double e_lris = tei[INDEX(l,r)*nn1o2+INDEX(i,s)];
                                    dum_b_d -= D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+r*nmo_+s] * e_lris;
                                }
                            }
                        }
                    }
                }
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    double e_pkqj = tei[INDEX(p,k)*nn1o2+INDEX(q,j)];
                                    dum_b_d -= D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+l] * e_pkqj;
                                }
                            }
                        }
                    }
                }

                // d,a

                for (int hq = 0; hq < nirrep_; hq++) {
                    for (int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        double e_ljqj = tei[INDEX(l,j)*nn1o2+INDEX(q,j)];
                        dum_d_a += D2ab[k*nmo_*nmo_*nmo_+q*nmo_*nmo_+i*nmo_+i] * e_ljqj;
                    }
                }
                for (int hs = 0; hs < nirrep_; hs++) {
                    for (int ss = 0; ss < amopi_[hs]; ss++) {
                        long int s = ss + pitzer_offset_full[hs];
                        double e_ikis = tei[INDEX(i,k)*nn1o2+INDEX(i,s)];
                        dum_d_a += D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+s] * e_ikis;
                    }
                }
                if ( j == l ) {
                    for (int hr = 0; hr < nirrep_; hr++) {
                        for (int rr = 0; rr < amopi_[hr]; rr++) {
                            long int r = rr + pitzer_offset_full[hr];
                            for (int hs = 0; hs < nirrep_; hs++) {
                                for (int ss = 0; ss < amopi_[hs]; ss++) {
                                    long int s = ss + pitzer_offset_full[hs];
                                    double e_iris = tei[INDEX(i,r)*nn1o2+INDEX(i,s)];
                                    dum_d_a -= D2ab[k*nmo_*nmo_*nmo_+j*nmo_*nmo_+r*nmo_+s] * e_iris;
                                }
                            }
                        }
                    }
                }
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    double e_pjqj = tei[INDEX(p,j)*nn1o2+INDEX(q,j)];
                                    dum_d_a -= D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i] * e_pjqj;
                                }
                            }
                        }
                    }
                }

                // d,b

                for (int hp = 0; hp < nirrep_; hp++) {
                    for (int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        double e_pjlj = tei[INDEX(p,j)*nn1o2+INDEX(l,j)];
                        dum_d_b += D2ab[p*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+i] * e_pjlj;
                    }
                }
                for (int hr = 0; hr < nirrep_; hr++) {
                    for (int rr = 0; rr < amopi_[hr]; rr++) {
                        long int r = rr + pitzer_offset_full[hr];
                        double e_irik = tei[INDEX(i,r)*nn1o2+INDEX(i,k)];
                        dum_d_b += D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+r*nmo_+l] * e_irik;
                    }
                }
                if ( j == l ) {
                    for (int hr = 0; hr < nirrep_; hr++) {
                        for (int rr = 0; rr < amopi_[hr]; rr++) {
                            long int r = rr + pitzer_offset_full[hr];
                            for (int hs = 0; hs < nirrep_; hs++) {
                                for (int ss = 0; ss < amopi_[hs]; ss++) {
                                    long int s = ss + pitzer_offset_full[hs];
                                    double e_iris = tei[INDEX(i,r)*nn1o2+INDEX(i,s)];
                                    dum_d_b -= D2ab[j*nmo_*nmo_*nmo_+k*nmo_*nmo_+r*nmo_+s] * e_iris;
                                }
                            }
                        }
                    }
                }
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    double e_pjqj = tei[INDEX(p,j)*nn1o2+INDEX(q,j)];
                                    dum_d_b -= D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+i*nmo_+l] * e_pjqj;
                                }
                            }
                        }
                    }
                }


                // d,d

                dum_d_d += D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+i*nmo_+i] * tei[INDEX(l,j)*nn1o2+INDEX(l,j)];
                dum_d_d += D2ab[j*nmo_*nmo_*nmo_+j*nmo_*nmo_+l*nmo_+l] * tei[INDEX(i,k)*nn1o2+INDEX(i,k)];

                if ( j == l ) {
                    for (int hr = 0; hr < nirrep_; hr++) {
                        for (int rr = 0; rr < amopi_[hr]; rr++) {
                            long int r = rr + pitzer_offset_full[hr];
                            for (int hs = 0; hs < nirrep_; hs++) {
                                for (int ss = 0; ss < amopi_[hs]; ss++) {
                                    long int s = ss + pitzer_offset_full[hs];
                                    double e_iris = tei[INDEX(i,r)*nn1o2+INDEX(i,s)];
                                    dum_d_d -= D2ab[k*nmo_*nmo_*nmo_+k*nmo_*nmo_+r*nmo_+s] * e_iris;
                                }
                            }
                        }
                    }
                }
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    double e_pjqj = tei[INDEX(p,j)*nn1o2+INDEX(q,j)];
                                    dum_d_d -= D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+l] * e_pjqj;
                                }
                            }
                        }
                    }
                }

                long int a_d_sym = (ij1 + 0 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                long int b_d_sym = (ij1 + 1 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                long int d_a_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl1 + 0 * nh);
                long int d_b_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl1 + 1 * nh);
                long int d_d_sym = (ij  + 2 * nh) * (2 * nh + nh2) + (kl  + 2 * nh);
                    
                if ( i != j ) newA[a_d_sym] += dum_a_d;
                if ( i != j ) newA[b_d_sym] += dum_b_d;
                if ( k != l ) newA[d_a_sym] += dum_d_a;
                if ( k != l ) newA[d_b_sym] += dum_d_b;
                newA[d_d_sym] += dum_d_d;

            }
        }

        // build B(ij,kl) = <| [ k*l, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {
            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];
            for (long int kl = 0; kl < nh; kl++) {
                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;

                if ( j == l ) {

                    duma    += D1a[k*nmo_+i];
                    dumb    += D1b[k*nmo_+i];

                }
                if ( i == k ) {

                    duma    -= D1a[j*nmo_+l];
                    dumb    -= D1b[j*nmo_+l];

                }

                long int aaaa_sym = (ij + 0 * nh) * (2 * nh + nh2) + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * (2 * nh + nh2) + (kl + 1 * nh);
//printf("aaaa %5i\n",aaaa_sym);
//printf("bbbb %5i\n",bbbb_sym);

                newB[aaaa_sym] = duma;
                newB[bbbb_sym] = dumb;

            }
        }

        // build A1(ij,kl) = <| [ k*l, [H1, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;

                duma    += D1a[k*nmo_+i] * h1[l*nmo_+j];
                duma    += D1a[j*nmo_+l] * h1[i*nmo_+k];

                dumb    += D1b[k*nmo_+i] * h1[l*nmo_+j];
                dumb    += D1b[j*nmo_+l] * h1[i*nmo_+k];

                if ( j == l ) {


                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                            dumb    -= D1b[k*nmo_+q] * h1[i*nmo_+q];
                        }
                    }

                }
                if ( i == k ) {

                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                            dumb    -= D1b[q*nmo_+l] * h1[q*nmo_+j];
                        }
                    }

                }

                long int aaaa_sym = (ij + 0 * nh) * (2 * nh + nh2) + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * (2 * nh + nh2) + (kl + 1 * nh);

                newA[aaaa_sym] = duma;
                newA[bbbb_sym] = dumb;

            }
        }
        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>

        // N^5 version of A2 coulomb-like terms (1) and (2)
        // 
        // 1: (ik|qs) D(jq;ls)
        // 2: (lj|qs) D(kq;is)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            if ( gems_ab[hjl] == 0 ) continue;

            double * tempa = (double*)malloc(nQ_*gems_ab[hjl]*sizeof(double));
            double * tempb = (double*)malloc(nQ_*gems_ab[hjl]*sizeof(double));

            memset((void*)tempa,'\0',nQ_*gems_ab[hjl]*sizeof(double));
            memset((void*)tempb,'\0',nQ_*gems_ab[hjl]*sizeof(double));

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                for (int Q = 0; Q < nQ_; Q++) {

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {

                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        long int jq = j*nmo_+q;
                        long int ls = l*nmo_+s;
                        
                        long int qj = q*nmo_+j;
                        long int sl = s*nmo_+l;
                                    
                        duma += Qmo_[INDEX(q,s)*nQ_+Q] * (D2aa[jq * nmo_*nmo_ + ls] + D2ab[jq * nmo_*nmo_ + ls]);
                        dumb += Qmo_[INDEX(q,s)*nQ_+Q] * (D2bb[jq * nmo_*nmo_ + ls] + D2ab[qj * nmo_*nmo_ + sl]);

                    }

                    tempa[Q*gems_ab[hjl] + jl] = duma;
                    tempb[Q*gems_ab[hjl] + jl] = dumb;
                }
            }

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {
                
                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;
                
                    int hi = symmetry[i];
                
                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int Q = 0; Q < nQ_; Q++) {
                        duma += tempa[Q*gems_ab[hjl] + jl] * Qmo2[INDEX(i,k)*nQ_+Q];
                        dumb += tempb[Q*gems_ab[hjl] + jl] * Qmo2[INDEX(i,k)*nQ_+Q];
                    }

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    long int aaaa_sym = (ij + 0 * nh) * (2 * nh + nh2) + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * (2 * nh + nh2) + (kl + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                    // AED: I think term (2) above is just a funny transpose of term (1)
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aaaa_sym = (lk + 0 * nh) * (2 * nh + nh2) + (ji + 0 * nh);
                    bbbb_sym = (lk + 1 * nh) * (2 * nh + nh2) + (ji + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;
                }
            }
            free(tempa);
            free(tempb);
        }

        // A2 exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int kq = k*nmo_+q;
                        long int is = i*nmo_+s;

                        duma -= eint * D2aa[kq*n+is];
                        dumb -= eint * D2bb[kq*n+is];

                        //eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                        //long int jq = j*nmo_+q;
                        //long int ls = l*nmo_+s;

                        //duma -= eint * D2aa[jq*n+ls];
                        //dumb -= eint * D2bb[jq*n+ls];
                    }

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    long int aaaa_sym = (ij + 0 * nh) * (2 * nh + nh2) + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * (2 * nh + nh2) + (kl + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aaaa_sym = (lk + 0 * nh) * (2 * nh + nh2) + (ji + 0 * nh);
                    bbbb_sym = (lk + 1 * nh) * (2 * nh + nh2) + (ji + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                }
            }
        }

        // A2 exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int qk = q*nmo_+k;
                        long int is = i*nmo_+s;

                        long int kq = k*nmo_+q;
                        long int si = s*nmo_+i;

                        dumab += eint * D2ab[qk*n+is];
                        dumba += eint * D2ab[kq*n+si];

                        //dumabab -= eint * D2ab[kq*n+is];
                        //dumbaba -= eint * D2ab[qk*n+si];

                    }

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    long int aabb_sym = (ij + 0 * nh) * (2 * nh + nh2) + (kl + 1 * nh);
                    long int bbaa_sym = (ij + 1 * nh) * (2 * nh + nh2) + (kl + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aabb_sym = (ji + 0 * nh) * (2 * nh + nh2) + (lk + 1 * nh);
                    bbaa_sym = (ji + 1 * nh) * (2 * nh + nh2) + (lk + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                }
            }
        }

        // N^5 version of A2 funky sums over 3 indices:
        // 
        // 1: -dik (qs|pj) D(pq;ls)
        // 2: -djl (qs|ip) D(kq;ps)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 
        // Also, the intermediates for this guy can be built
        // at N^4 cost and can be done outside of these loops (see above)

        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            int hj = symmetry[j];
            for (int ll = 0; ll < amopi_[hj]; ll++) {

                int l = ll + pitzer_offset[hj];

                if ( i == l ) continue;

                long int il = ibas_erpa_sym[h][i][l];

                long int aaaa_sym = (ij + 0 * nh) * (2 * nh + nh2) + (il + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * (2 * nh + nh2) + (il + 1 * nh);

                newA[aaaa_sym] -= Ia[j*nmo_+l];
                newA[bbbb_sym] -= Ib[j*nmo_+l];

                // AED: I think second funky term is a transpose of the first
                // A(ij,il)(2) = A(ji,li)(1)

                long int ji = ibas_erpa_sym[h][j][i];
                long int li = ibas_erpa_sym[h][l][i];

                aaaa_sym = (ji + 0 * nh) * (2 * nh + nh2) + (li + 0 * nh);
                bbbb_sym = (ji + 1 * nh) * (2 * nh + nh2) + (li + 1 * nh);

                newA[aaaa_sym] -= Ia[j*nmo_+l];
                newA[bbbb_sym] -= Ib[j*nmo_+l];

            }
        }

        // A2: last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        duma += eint * D2aa[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];
                        dumb += eint * D2bb[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];

                    }

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    long int aaaa_sym = (ij + 0 * nh) * (2 * nh + nh2) + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * (2 * nh + nh2) + (kl + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aaaa_sym = (ji + 0 * nh) * (2 * nh + nh2) + (lk + 0 * nh);
                    bbbb_sym = (ji + 1 * nh) * (2 * nh + nh2) + (lk + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                }

            }

        }

        // A2: last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)

        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        dumab -= eint * D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+i*nmo_+l];
                        dumba -= eint * D2ab[q*nmo_*nmo_*nmo_+p*nmo_*nmo_+l*nmo_+i];

                    }

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    long int aabb_sym = (ij + 0 * nh) * (2 * nh + nh2) + (kl + 1 * nh);
                    long int bbaa_sym = (ij + 1 * nh) * (2 * nh + nh2) + (kl + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aabb_sym = (ji + 0 * nh) * (2 * nh + nh2) + (lk + 1 * nh);
                    bbaa_sym = (ji + 1 * nh) * (2 * nh + nh2) + (lk + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                }

            }

        }

        std::shared_ptr<Matrix> Amat ( new Matrix(totdim,totdim) );
        std::shared_ptr<Matrix> Bmat ( new Matrix(totdim,totdim) );
        for (int ij = 0; ij < totdim; ij++) {
            for (int kl = 0; kl < totdim; kl++) {
                Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                Bmat->pointer()[ij][kl] = newB[ij*totdim+kl];
            }
        }

        // project out null (and negative) space of A matrix
        std::shared_ptr<Matrix> eigvec ( new Matrix(totdim,totdim) );
        std::shared_ptr<Vector> eigval ( new Vector(totdim) );

        Amat->diagonalize(eigvec,eigval,descending);
        //eigval->print();

        bool prune = false;
        long int newdim = 0;

        // get rid of small/negative eigenvalues of A?
        if ( prune ) {

            newdim = 0;
            int * map = (int*)malloc(totdim*sizeof(int));
            memset((void*)map,'\0',totdim*sizeof(int));
            int * skip = (int*)malloc(totdim*sizeof(int));
            memset((void*)skip,'\0',totdim*sizeof(int));
            for (long int i = 0; i < totdim; i++) {
                double val = eigval->pointer()[i];
                skip[i] = 1;
                if ( val > pruning_threshold ) {
                    map[newdim] = i;
                    skip[i] = 0;
                    newdim++;
                }
            }
            //if ( newdim % 2 == 1 ) { // eigenvalues should come in pairs
            //    newdim--;
            //}
            for (int ij = 0; ij < totdim; ij++) {
                for (int kl = 0; kl < totdim; kl++) {
                    Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                }
            }
            Amat->transform(eigvec);
            Bmat->transform(eigvec);
        }else {
            newdim = totdim;
        }

        for (long int i = 0; i < newdim; i++) {
            for (long int j = 0; j < newdim; j++) {
                newA[i*newdim+j] = Amat->pointer()[i][j];
                newB[i*newdim+j] = Bmat->pointer()[i][j];
            }
        }

        outfile->Printf("\n");
        if ( newdim < totdim ) {
            outfile->Printf("    <<< warning >>> ");
            outfile->Printf("\n");
            outfile->Printf("    reducing dimension of [ k*l, [H, j*i ] ] from %5i to %5i\n",totdim,newdim);
            outfile->Printf("\n");
        }

        int info = 0;
        if ( prune ) {
            info = SymmetricGeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }else {
            GeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }

        // back transform wave functions to original untruncated basis
        std::shared_ptr<Matrix>  Cmat (new Matrix(totdim,totdim));
        for (int i = 0; i < newdim; i++) {
            for (int j = 0; j < totdim; j++) {
                double dum = 0.0;
                for (int k = 0; k < newdim; k++) {
                    dum += eigvec->pointer()[j][k] * newB[i * newdim + k];
                }
                Cmat->pointer()[i][j] = dum;
            }
        }

        // oscillator strengths: 
        // symmetry of state is h, consider only cartesian component with same symmetry
       
        if ( info != 0 ) {
            outfile->Printf("\n");
            outfile->Printf("    Error: diagonalization failed.  info = %5i\n",info);
            outfile->Printf("\n");
        }

        outfile->Printf("    state");
        outfile->Printf("          energy (Eh)");
        outfile->Printf("       ex energy (Eh)");
        outfile->Printf("       ex energy (eV)");
        outfile->Printf("       f, osc. strength\n");

        for (long int state = newdim-1; state >= 0; state--){
          
            if ( eig[state] > 0.0 ) {

                if ( prune ) eig[state] = 1.0 / eig[state];

                 double val = 0.0;//2./3. * eig[state] * (dumx*dumx+dumy*dumy+dumz*dumz);
              
                 outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138, val);
                 //outfile->Printf("                Transition Moment: X: %10.6lf    Y: %10.6lf   Z: %10.6lf\n", dumx, dumy, dumz);
               
            } 
        }
   
        free(newA);
        free(newB);
        free(cc);
        free(eig);
    }

    free(Ia);
    free(Ib);
    free(tei);
    free(h1);
    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);
    if ( !is_df_ ) {
        free(Qmo2);
    }
*/
}

double ERPASolver::DOCIAlphaERPA(double alpha,bool is_singlet) {

    //outfile->Printf("\n");
    //outfile->Printf("    ==> Alpha ERPA <==\n");
    //outfile->Printf("\n");
    //outfile->Printf("    alpha:    %20.12lf\n",alpha);
    //outfile->Printf("\n");
    outfile->Printf("    %20.12lf",alpha);

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals
    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // build two-electron integrals

    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));

    int info;
    if ( is_df_ ) {
        F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
    }else { 
        // unpack two-electron integrals (which were initially blocked by symmetry)
        long int offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }
    }

    // read TDPM from disk
    ReadTPDM(D2aa,D2bb,D2ab,D1a,D1b,tei);

    //ReadOPDM(D1a,D1b);

    long int n = nmo_*nmo_;

    // check energy
    double e2_aa = 0.0;
    double e2_ab = 0.0;
    double e2_bb = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2_ab +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2_aa += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2_bb += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];
                }
            }
        }
    }
    double e2 = e2_aa + e2_ab + e2_bb;

    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }
    Process::environment.globals["v2RDM TOTAL ENERGY"] = e1 + e2 + enuc_;

    //printf("%20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",e1,e2_ab,e2_aa,e2_bb,e1+e2,e1+e2+enuc_);
    //printf("%5i\n",info);
    //exit(0);

    // scale one-electron integrals ha = h(0) + ah'
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = 0; j < nmo_; j++) {

            //if ( i == j ) {
            //    h1[i*nmo_+j] = h1[i*nmo_+j];
            //}else {
            //    h1[i*nmo_+j] = alpha * h1[i*nmo_+j];
            //}
            if ( i != j ) h1[i*nmo_+j] *= alpha;
        }
    }

    // correlated (ERPA) ground-state 2-RDM
    double * D2aa_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
    double * D2bb_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
    double * D2ab_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));

    memset((void*)D2aa_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));
    memset((void*)D2bb_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));
    memset((void*)D2ab_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));

    // Gamma(pqrs) = gamma(pr)gamma(qs) + sum_v gamma_0v(pr)gamma_0v(qs) - gamma(qr)delta(ps)

    // D2aa(prqs) = gamma_a(pr)gamma_a(qs) + sum_v gamma_a_0v(pr)gamma_a_0v(qs) - gamma_a(qr)delta(ps)
    // D2bb(prqs) = gamma_b(pr)gamma_b(qs) + sum_v gamma_b_0v(pr)gamma_b_0v(qs) - gamma_b(qr)delta(ps)
    // D2ab(prqs) = gamma_a(pr)gamma_b(qs) + sum_v gamma_a_0v(pr)gamma_b_0v(qs)

    if ( is_singlet ) {
        for (int p = 0; p < nso_; p++) {
            for (int q = 0; q < nso_; q++) {
                for (int r = 0; r < nso_; r++) {
                    for (int s = 0; s < nso_; s++) {

                        double ga_pr = D1a[p*nso_+r];
                        double ga_qr = D1a[q*nso_+r];
                        double ga_qs = D1a[q*nso_+s];

                        double gb_pr = D1b[p*nso_+r];
                        double gb_qr = D1b[q*nso_+r];
                        double gb_qs = D1b[q*nso_+s];

                        D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_pr * gb_qs;
                        D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_pr * ga_qs - ga_qr * (p==s);// * 0.5;
                        D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gb_pr * gb_qs - gb_qr * (p==s);// * 0.5;

                    }
                }
            }
        }
    }

    double * alpha_tei_aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * alpha_tei_ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)alpha_tei_aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)alpha_tei_ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // scale two-electron integrals
    for (long int i = 0; i < nmo_; i++) {
        for (long int k = 0; k < nmo_; k++) {
            for (long int j = 0; j < nmo_; j++) {
                for (long int l = 0; l < nmo_; l++) {

                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];

                    // ab
                    if ( i == k && j == l ) {
                        // D2(ij,ij)
                        alpha_tei_ab[i*nmo_*nmo_*nmo_ + k*nmo_*nmo_ + j*nmo_ + l] = eint;
                    }else if ( i == j && k == l ) {
                        // D2(ii,jj)
                        alpha_tei_ab[i*nmo_*nmo_*nmo_ + k*nmo_*nmo_ + j*nmo_ + l] = eint;
                    }else {
                        // D2(ij,kl)
                        alpha_tei_ab[i*nmo_*nmo_*nmo_ + k*nmo_*nmo_ + j*nmo_ + l] = alpha * eint;
                    }

                    // aa
                    if ( i == k && j == l ) {
                        // D2(ij,ij)
                        alpha_tei_aa[i*nmo_*nmo_*nmo_ + k*nmo_*nmo_ + j*nmo_ + l] = eint;
                    }else if ( i == l && j == k ) {
                        // D2(ij,ji)
                        alpha_tei_aa[i*nmo_*nmo_*nmo_ + k*nmo_*nmo_ + j*nmo_ + l] = eint;
                    //}else if ( i == j && k == l ) {
                    //    // D2(ii,jj)
                    //    alpha_tei_aa[i*nmo_*nmo_*nmo_ + k*nmo_*nmo_ + j*nmo_ + l] = 0.0;
                    }else {
                        // D2(ij,kl)
                        alpha_tei_aa[i*nmo_*nmo_*nmo_ + k*nmo_*nmo_ + j*nmo_ + l] = alpha * eint;
                    }
                }
            }
        }
    }

    for (int h = 0; h < nirrep_; h++) {

        long int nh = gems_erpa[h];
        //long int nh = gems_ab[h];
        if ( nh == 0 ) continue;

        long int totdim = nh;
        //if ( type == "AA" ) totdim = 2*nh;

        double * newB = (double*)malloc(totdim*totdim*sizeof(double));
        double * newA = (double*)malloc(totdim*totdim*sizeof(double));
        double * cc   = (double*)malloc(totdim*totdim*sizeof(double));
        double * eig  = (double*)malloc(totdim*sizeof(double));

        memset((void*)newA,'\0',totdim*totdim*sizeof(double));
        memset((void*)newB,'\0',totdim*totdim*sizeof(double));
        memset((void*)cc,'\0',totdim*totdim*sizeof(double));
        memset((void*)eig,'\0',totdim*sizeof(double));

        //outfile->Printf("\n");
        //outfile->Printf("    Symmetry: %5i\n",h);
        //outfile->Printf("\n");

        // build B(ij,kl) = <| [ k*l, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {
            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];
            for (long int kl = 0; kl < nh; kl++) {
                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;
                double dumabab = 0.0;
                double dumbaba = 0.0;

                if ( j == l ) {
                    duma    += D1a[k*nmo_+i];
                    dumb    += D1b[k*nmo_+i];
                }
                if ( i == k ) {
                    duma    -= D1a[j*nmo_+l];
                    dumb    -= D1b[j*nmo_+l];
                }

                long int id = ij * nh + kl;
                newB[id] = 0.5 * (duma + dumb);
            }
        }

        // build A1(ij,kl) = <| [ k*l, [H1, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;

                duma    += D1a[k*nmo_+i] * h1[l*nmo_+j];
                duma    += D1a[j*nmo_+l] * h1[i*nmo_+k];

                dumb    += D1b[k*nmo_+i] * h1[l*nmo_+j];
                dumb    += D1b[j*nmo_+l] * h1[i*nmo_+k];

                if ( j == l ) {
                    //for (int hq = 0; hq < nirrep_; hq++) {
                    //    for (int qq = 0; qq < amopi_[hq]; qq++) {
                    //        long int q = qq + pitzer_offset_full[hq];
                    //        duma    -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                    //        dumb    -= D1b[k*nmo_+q] * h1[i*nmo_+q];
                    //    }
                    //}
                    duma    -= D1a[k*nmo_+k] * h1[i*nmo_+k]; // DOCI D1 is diagonal
                    dumb    -= D1b[k*nmo_+k] * h1[i*nmo_+k]; // DOCI D1 is diagonal
                }
                if ( i == k ) {
                    //for (int hq = 0; hq < nirrep_; hq++) {
                    //    for (int qq = 0; qq < amopi_[hq]; qq++) {
                    //        long int q = qq + pitzer_offset_full[hq];
                    //        duma    -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                    //        dumb    -= D1b[q*nmo_+l] * h1[q*nmo_+j];
                    //    }
                    //}
                    duma    -= D1a[l*nmo_+l] * h1[l*nmo_+j]; // DOCI D1 is diagonal
                    dumb    -= D1b[l*nmo_+l] * h1[l*nmo_+j]; // DOCI D1 is diagonal
                }

                long int id = ij * nh + kl;
                newA[id] += 0.5 * (duma + dumb);
            }
        }

//#if 0
        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>


        // 2 coulomb-like terms:
        // 
        // 1: (ik|qs) D(jq;ls)
        // 2: (lj|qs) D(kq;is)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma = 0.0;
                double dumb = 0.0;

                // 1: (ik|qs) D(jq;ls)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                double eint_aa = alpha_tei_aa[INDEX999(nmo_,i,k)*nmo_*nmo_+INDEX999(nmo_,q,s)];
                                double eint_ab = alpha_tei_ab[INDEX999(nmo_,i,k)*nmo_*nmo_+INDEX999(nmo_,q,s)];

                                long int jq = j*nmo_+q;
                                long int ls = l*nmo_+s;

                                long int qj = q*nmo_+j;
                                long int sl = s*nmo_+l;

                                duma    += eint_aa * D2aa[jq*n+ls] + eint_ab * D2ab[jq*n+ls];
                                dumb    += eint_aa * D2bb[jq*n+ls] + eint_ab * D2ab[qj*n+sl];
                            }
                        }
                    }
                }
                // testing DOCI-specific structure of D2
/*
                long int jj = j * nmo_ + j;
                long int ll = l * nmo_ + l;

                double eint = alpha_tei_ab[INDEX999(nmo_,i,k)*nmo_*nmo_+INDEX999(nmo_,j,l)];
                duma += eint * D2ab[jj*n+ll];
                dumb += eint * D2ab[jj*n+ll];
                if ( j == l ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (long int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            if ( j == q ) continue;
                            eint = alpha_tei_ab[INDEX999(nmo_,i,k)*nmo_*nmo_+INDEX999(nmo_,q,q)];
                            long int jq = j * nmo_ + q;
                            long int qj = q * nmo_ + j;
                            duma += eint * D2ab[jq*n+jq];
                            dumb += eint * D2ab[qj*n+qj];
                        }
                    }
                }
*/

                // 2: (lj|qs) D(kq;is)

                long int id = ij * nh + kl;
                newA[id] += 0.5 * ( duma + dumb );

                // AED: I think term (2) above is just a funny transpose of term (1)
                // A(ij,kl)(2) = A(lk,ji)(1)

                long int ji = ibas_erpa_sym[h][j][i];
                long int lk = ibas_erpa_sym[h][l][k];

                id = ji * nh + lk;
                newA[id] += 0.5 * (duma + dumb);

            }
        }

        // funky sums over 3 indices:
        // 
        // 1: -dik (qs|pj) D(pq;ls)
        // 2: -djl (qs|ip) D(kq;ps)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma = 0.0;
                double dumb = 0.0;

                // funky sums over 3 indices:
                // - dik (qs|pj) D(pq;ls) + djl (qs|ip) D(kq;ps)
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    for (int hs = 0; hs < nirrep_; hs++) {
                                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                            long int s = ss + pitzer_offset_full[hs];
                                            double eint_aa = alpha_tei_aa[INDEX999(nmo_,q,s)*nmo_*nmo_+INDEX999(nmo_,p,j)];
                                            double eint_ab = alpha_tei_ab[INDEX999(nmo_,q,s)*nmo_*nmo_+INDEX999(nmo_,p,j)];

                                            long int pq = p*nmo_+q;
                                            long int ls = l*nmo_+s;

                                            long int qp = q*nmo_+p;
                                            long int sl = s*nmo_+l;

                                            duma -= eint_aa * D2aa[pq*n+ls] + eint_ab * D2ab[pq*n+ls];
                                            dumb -= eint_aa * D2bb[pq*n+ls] + eint_ab * D2ab[qp*n+sl];
                                        }
                                    }
                                }
                            }
                        }
                    }

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * ( duma + dumb );

                    // AED: I think second funky term is a transpose of the first
                    // A(ij,il)(2) = A(ji,li)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int li = ibas_erpa_sym[h][l][i];

                    id = ji * nh + li;
                    newA[id] += 0.5 * ( duma + dumb );

                }

                //if ( j == l ) {
                //    for (int hp = 0; hp < nirrep_; hp++) {
                //        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                //            long int p = pp + pitzer_offset_full[hp];
                //            for (int hq = 0; hq < nirrep_; hq++) {
                //                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                //                    long int q = qq + pitzer_offset_full[hq];
                //                    for (int hs = 0; hs < nirrep_; hs++) {
                //                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                //                            long int s = ss + pitzer_offset_full[hs];
                //                            //double eint = TEI(q,s,i,p,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,s),1,Qmo_+nQ_*INDEX(i,p),1);
                //                            double eint = alpha_tei[INDEX(q,s)*nn1o2+INDEX(i,p)];

                //                            long int kq = k*nmo_+q;
                //                            long int ps = p*nmo_+s;

                //                            long int qk = q*nmo_+k;
                //                            long int sp = s*nmo_+p;

                //                            duma -= eint * ( D2aa[kq*n+ps] + D2ab[kq*n+ps] );
                //                            dumb -= eint * ( D2bb[kq*n+ps] + D2ab[qk*n+sp] );
                //                        }
                //                    }
                //                }
                //            }
                //        }
                //    }
                //}

                // AED: I think second funky term is a transpose of the first (accounted for above)
                // A(ij,il)(2) = A(ji,li)(1)

                //long int ji = ibas_erpa_sym[h][j][i];
                //long int li = ibas_erpa_sym[h][l][i];

                //id = ji * nh + li;
                //newA[id] += 0.5 * ( duma + dumb );

            }
        }
/*
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma = 0.0;
                double dumb = 0.0;

                double dumab = 0.0;
                double dumba = 0.0;

                // exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                double eint = alpha_tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                long int kq = k*nmo_+q;
                                long int is = i*nmo_+s;

                                duma -= eint * D2aa[kq*n+is];
                                dumb -= eint * D2bb[kq*n+is];

                                eint = alpha_tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                long int jq = j*nmo_+q;
                                long int ls = l*nmo_+s;

                                duma -= eint * D2aa[jq*n+ls];
                                dumb -= eint * D2bb[jq*n+ls];
                            }
                        }
                    }
                }

                // exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];

                                double eint = alpha_tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                long int qk = q*nmo_+k;
                                long int is = i*nmo_+s;

                                long int kq = k*nmo_+q;
                                long int si = s*nmo_+i;

                                dumab += eint * D2ab[qk*n+is];
                                dumba += eint * D2ab[kq*n+si];

                                eint = alpha_tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                long int jq = j*nmo_+q;
                                long int sl = s*nmo_+l;

                                long int qj = q*nmo_+j;
                                long int ls = l*nmo_+s;

                                dumab += eint * D2ab[jq*n+sl];
                                dumba += eint * D2ab[qj*n+ls];
                            }
                        }
                    }
                }

                // last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
                for (int hp = 0; hp < nirrep_; hp++) {
                    for (long int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        for (int hq = 0; hq < nirrep_; hq++) {
                            for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                long int q = qq + pitzer_offset_full[hq];

                                //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                double eint = alpha_tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                long int pq = p*nmo_+q;
                                long int li = l*nmo_+i;

                                duma += eint * D2aa[pq*n+li];
                                dumb += eint * D2bb[pq*n+li];

                                //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                eint = alpha_tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                long int kj = k*nmo_+j;

                                duma += eint * D2aa[kj*n+pq];
                                dumb += eint * D2bb[kj*n+pq];
                            }
                        }
                    }
                }

                // last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)
                for (int hp = 0; hp < nirrep_; hp++) {
                    for (long int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        for (int hq = 0; hq < nirrep_; hq++) {
                            for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                long int q = qq + pitzer_offset_full[hq];

                                //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                double eint = alpha_tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                long int pq = p*nmo_+q;
                                long int il = i*nmo_+l;

                                long int qp = q*nmo_+p;
                                long int li = l*nmo_+i;

                                dumab -= eint * D2ab[pq*n+il];
                                dumba -= eint * D2ab[qp*n+li];

                                //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                eint = alpha_tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                long int kj = k*nmo_+j;
                                long int jk = j*nmo_+k;

                                dumab -= eint * D2ab[jk*n+pq];
                                dumba -= eint * D2ab[kj*n+qp];
                            }
                        }
                    }
                }

                long int id = ij * nh + kl;

                newA[id] += 0.5 * ( duma + dumab + dumba + dumb );

            }
        }
*/

        // A2 exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)
        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint_aa = alpha_tei_aa[INDEX999(nmo_,q,j)*nmo_*nmo_+INDEX999(nmo_,l,s)];

                        long int kq = k*nmo_+q;
                        long int is = i*nmo_+s;

                        duma -= eint_aa * D2aa[kq*n+is];
                        dumb -= eint_aa * D2bb[kq*n+is];

                        //eint = alpha_tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                        //long int jq = j*nmo_+q;
                        //long int ls = l*nmo_+s;

                        //duma -= eint * D2aa[jq*n+ls];
                        //dumb -= eint * D2bb[jq*n+ls];
                    }

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * (duma + dumb);

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    newA[id] += 0.5 * (duma + dumb);
                }
            }
        }

        // A2 exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint_ab = alpha_tei_ab[INDEX999(nmo_,q,j)*nmo_*nmo_+INDEX999(nmo_,l,s)];

                        long int qk = q*nmo_+k;
                        long int is = i*nmo_+s;

                        long int kq = k*nmo_+q;
                        long int si = s*nmo_+i;

                        dumab += eint_ab * D2ab[qk*n+is];
                        dumba += eint_ab * D2ab[kq*n+si];

                    }

                    long int id = ij * nh + kl;
                    if ( is_singlet ) {
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    if ( is_singlet ) {
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                }
            }
        }

        // A2: last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint_aa = alpha_tei_aa[INDEX999(nmo_,p,j)*nmo_*nmo_+INDEX999(nmo_,q,k)];
                        
                        duma += eint_aa * D2aa[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];
                        dumb += eint_aa * D2bb[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];

                    }

                    long int id = ij * nh + kl;

                    newA[id] += 0.5 * (duma + dumb);

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitaitons, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    newA[id] += 0.5 * (duma + dumb);

                }

            }

        }

        // A2: last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint_ab = alpha_tei_ab[INDEX999(nmo_,p,j)*nmo_*nmo_+INDEX999(nmo_,q,k)];
                        
                        dumab -= eint_ab * D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+i*nmo_+l];
                        dumba -= eint_ab * D2ab[q*nmo_*nmo_*nmo_+p*nmo_*nmo_+l*nmo_+i];

                    }

                    long int id = ij * nh + kl;
                    if ( is_singlet ) { 
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    if ( is_singlet ) { 
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                }

            }

        }

        std::shared_ptr<Matrix> Amat ( new Matrix(totdim,totdim) );
        std::shared_ptr<Matrix> Bmat ( new Matrix(totdim,totdim) );
        for (int ij = 0; ij < totdim; ij++) {
            for (int kl = 0; kl < totdim; kl++) {
                Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                Bmat->pointer()[ij][kl] = newB[ij*totdim+kl];
            }
        }
   //      outfile->Printf("A is");
   //     Amat->print();
   //     outfile->Printf("B is");
   //     Bmat->print();

        // symmetrize
        // A(kl,ji) = A(ij,kl)
        // A(kl,ji) = A(lk,ij)

        // project out null (and negative) space of A matrix
        std::shared_ptr<Matrix> eigvec ( new Matrix(totdim,totdim) );
        std::shared_ptr<Vector> eigval ( new Vector(totdim) );
//        outfile->Printf("A is");
//        Amat->print();
//        outfile->Printf("B is");
//        Bmat->print();
       // eigval->print();

        bool prune = false;
        long int newdim = 0;

        // get rid of small/negative eigenvalues of A?
        if ( prune ) {

            Amat->diagonalize(eigvec,eigval,descending);

            newdim = 0;
            int * map = (int*)malloc(totdim*sizeof(int));
            memset((void*)map,'\0',totdim*sizeof(int));
            int * skip = (int*)malloc(totdim*sizeof(int));
            memset((void*)skip,'\0',totdim*sizeof(int));
            for (long int i = 0; i < totdim; i++) {
                double val = eigval->pointer()[i];
                skip[i] = 1;
                if ( val > pruning_threshold ) {
                    map[newdim] = i;
                    skip[i] = 0;
                    newdim++;
                }
            }
            //if ( newdim % 2 == 1 ) { // eigenvalues should come in pairs
            //    newdim--;
            //}
            for (int ij = 0; ij < totdim; ij++) {
                for (int kl = 0; kl < totdim; kl++) {
                    Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                }
            }
            Amat->transform(eigvec);
            Bmat->transform(eigvec);
            for (long int i = 0; i < newdim; i++) {
                for (long int j = 0; j < newdim; j++) {
                    newA[i*newdim+j] = Amat->pointer()[i][j];
                    newB[i*newdim+j] = Bmat->pointer()[i][j];
                }
            }

        }else {
            newdim = totdim;
        }

        if ( newdim < totdim ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< warning >>> ");
            outfile->Printf("\n");
            outfile->Printf("    reducing dimension of [ k*l, [H, j*i ] ] from %5i to %5i\n",totdim,newdim);
            outfile->Printf("\n");
        }

        std::shared_ptr<Matrix>  Cmat (new Matrix(totdim,totdim));
        std::shared_ptr<Matrix>  Cmatl (new Matrix(totdim,totdim));
        bool * energy_is_real = (bool*)malloc(totdim*sizeof(bool));

        int info = 0;
        if ( prune ) {
            info = SymmetricGeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
             outfile->Printf("Symmetric_solver\n");
            // back transform wave functions to original untruncated basis
            for (int i = 0; i < newdim; i++) {
                for (int j = 0; j < totdim; j++) {
                    double dum = 0.0;
                    for (int k = 0; k < newdim; k++) {
                        dum += eigvec->pointer()[j][k] * newB[i * newdim + k];
                    }
                    Cmat->pointer()[i][j] = dum;
                }
            }

        }else {
            // symmetrize A Nam, remove!
            for (int i = 0; i < newdim; i++) {
                for (int j = i; j < newdim; j++) {
                    double dum = newA[i*newdim +j] + newA[j*newdim+i];
                    newA[i*newdim+j] = 0.5 * dum;
                    newA[j*newdim+i] = 0.5 * dum;
	        }
	    }
            GeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig,energy_is_real);
              outfile->Printf("General_solver\n");
            // no need to back transform wave functions to original untruncated basis in this case
            // NOTE: this is either the left or right eigenfunction ... still need the other one
            for (int i = 0; i < newdim; i++) {
                for (int j = 0; j < totdim; j++) {
                    Cmat->pointer()[i][j] = newB[i * newdim + j];
                    Cmatl->pointer()[i][j] = newA[i * newdim + j];
                }
            }
                //  outfile->Printf("Cmat is");
                //  Cmat->print();
                //  outfile->Printf("Cmat left is");
                //  Cmatl->print();
        }

        // oscillator strengths: 
        // symmetry of state is h, consider only cartesian component with same symmetry
       
        if ( info != 0 ) {
            outfile->Printf("\n");
            outfile->Printf("    Error: diagonalization failed.  info = %5i\n",info);
            outfile->Printf("\n");
            throw PsiException("Error: diagonalization failed.",__FILE__,__LINE__);
        }
        //Cmatl->transpose();
        //outfile->Printf("    state");
        //outfile->Printf("          energy (Eh)");
        //outfile->Printf("       ex energy (Eh)");
        //outfile->Printf("       ex energy (eV)");
        //outfile->Printf("       f, osc. strength\n");

        for (long int state = newdim-1; state >= 0; state--){
          
            //printf("%5i %20.12lf\n",state,eig[state]);
            if ( eig[state] > 0.0 ) {

                if ( prune ) eig[state] = 1.0 / eig[state];
		else if ( !energy_is_real[state] ) {
		    continue;
		}

                // avoid considering nonsensical states
                if ( fabs(eig[state] ) > 10000.0) continue;
                outfile->Printf("\neigenvalue %20.12lf", eig[state]);
//printf("hey excited state: %20.12lf\n",eig[state]+e2+e1+enuc_);
               // normalize C://origin 11/01/
          //     double nrm = 0.0;
          //     double * C = Cmat->pointer()[state];
          //     for (int ij = 0; ij < nh; ij++) {
          //         int i  = bas_erpa_sym[h][ij][0];
          //         int j  = bas_erpa_sym[h][ij][1];
          //         int hi = symmetry[i];
          //         int hj = symmetry[j];
          //         for( int kk = 0; kk < amopi_[hi]; kk++){
          //             int k   = kk + pitzer_offset[hi];
          //             int kj  = ibas_erpa_sym[h][k][j];

          //             // check if kj is a proper ERPA excitations
          //             if ( kj == -999 ) continue;

          //             nrm    += 0.5 * C[ij] * C[kj] * D1a[k*amo_+i];
          //             nrm    += 0.5 * C[ij] * C[kj] * D1b[k*amo_+i];
          //         }
          //         for( int ll = 0; ll < amopi_[hj]; ll++){
          //             int l  = ll + pitzer_offset[hj];
          //             int il = ibas_erpa_sym[h][i][l];

          //             // check if il is a proper ERPA excitations
          //             if ( il == -999 ) continue;

          //             nrm   -= 0.5 * C[ij] * C[il] * D1a[j*amo_+l];
          //             nrm   -= 0.5 * C[ij] * C[il] * D1b[j*amo_+l];
          //         }
          //     }


          //     if ( nrm < 0.0 ) continue;

          //     nrm = sqrt(nrm);
          //     outfile->Printf("\n norm is %20lf",nrm);
          //      outfile->Printf("\n");
          //     for (int ij = 0; ij < nh; ij++) {
          //         C[ij] /= nrm;
          //     }
          //    outfile->Printf("normalized eigenvector\n");
          //    Cmat->print();
            // normalize C with both right and left eigenvectors
            // Ckl(l)<0|[k*l,j*i|0>Cij(r)
            // outfile->Printf("\nnormalize both left and right eigevectors\n");
             double nrm = 0.0;
             double * C = Cmat->pointer()[state];
             double * Cl = Cmatl->pointer()[state];
             for (int ij = 0; ij < nh; ij++) {
                 int i  = bas_erpa_sym[h][ij][0];
                 int j  = bas_erpa_sym[h][ij][1];
                 int hi = symmetry[i];
                 int hj = symmetry[j];
                 for( int kk = 0; kk < amopi_[hi]; kk++){
                     int k   = kk + pitzer_offset[hi];
                     int kj  = ibas_erpa_sym[h][k][j];

                     // check if kj is a proper ERPA excitations
                     if ( kj == -999 ) continue;
                     nrm    += 0.5 * C[ij] * Cl[kj] * D1a[k*amo_+i];
                     nrm    += 0.5 * C[ij] * Cl[kj] * D1b[k*amo_+i];
                 }
                 for( int ll = 0; ll < amopi_[hj]; ll++){
                     int l  = ll + pitzer_offset[hj];
                     int il = ibas_erpa_sym[h][i][l];

                     // check if il is a proper ERPA excitations
                     if ( il == -999 ) continue;

                     nrm   -= 0.5 * C[ij] * Cl[il] * D1a[j*amo_+l];
                     nrm   -= 0.5 * C[ij] * Cl[il] * D1b[j*amo_+l];
                 }
             }
            
             if ( nrm < 0.0 ) continue;

             nrm = sqrt(nrm);
             outfile->Printf("\n norm is %20lf\n\n\n", nrm);
             for (int ij = 0; ij < nh; ij++) {
                 C[ij] /= nrm;
                 Cl[ij]/= nrm;
             }
 
                double * TDa = (double*)malloc(nso_*nso_*sizeof(double));//right eigvec
                double * TDb = (double*)malloc(nso_*nso_*sizeof(double));
                double * TDla = (double*)malloc(nso_*nso_*sizeof(double));//for left eigenvectors
                double * TDlb = (double*)malloc(nso_*nso_*sizeof(double));

                memset((void*)TDa,'\0',nso_*nso_*sizeof(double));
                memset((void*)TDb,'\0',nso_*nso_*sizeof(double));
                memset((void*)TDla,'\0',nso_*nso_*sizeof(double));
                memset((void*)TDlb,'\0',nso_*nso_*sizeof(double));

                // build transition density matrix:
                TransitionDensity_SpinAdapted(Cmat->pointer()[state], D1a, D1b, TDa,TDb,nh,h);
                TransitionDensity_SpinAdaptedl(Cmatl->pointer()[state], D1a, D1b, TDla,TDlb,nh,h);
                if ( !is_singlet ) {
                   C_DSCAL(nso_*nso_,-1.0,TDb,1);
                   C_DSCAL(nso_*nso_,-1.0,TDlb,1);
                }
                // use transition density matrix to update ground-state 2-RDM
       //         for (int p = 0; p < nso_; p++) {
       //             for (int q = 0; q < nso_; q++) {
       //                 for (int r = 0; r < nso_; r++) {
       //                     for (int s = 0; s < nso_; s++) {
       //                         double ga_rp = TDa[p*nso_+r]; // careful with transpose ...
       //                         double gb_qs = TDb[s*nso_+q];
       //                         // Eq(10)
       //                         D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_rp * gb_qs;// * 0.5;
       //                         //if ( p != q && r != s) {

       //                             double ga_qs = TDa[s*nso_+q];
       //                             D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_rp * ga_qs;// * 0.5;

       //                             double gb_rp = TDb[p*nso_+r]; // careful with transpose ...
       //                             D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gb_rp * gb_qs;// * 0.5;
       //                         //}

       //                     }
       //                 }
       //             }
       //         }
       //         free(TDa);
       //         free(TDb);
               for (int p = 0; p < nso_; p++) {
                    for (int q = 0; q < nso_; q++) {
                        for (int r = 0; r < nso_; r++) {
                            for (int s = 0; s < nso_; s++) {
                                double gra_pr = TDa[r*nso_+p]; //right
                                double grb_pr = TDb[r*nso_+p]; 
                                double gla_qs = TDla[s*nso_+q];//left
                                double glb_qs = TDlb[s*nso_+q];
                                // Eq(10)       <0|r*p|v><v|s*q|0>
                                D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gra_pr * glb_qs;// * 0.5;
                                //if ( p != q && r != s) {

                                    D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gra_pr * gla_qs;// * 0.5;

                                    D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += grb_pr * glb_qs;// * 0.5;
                                //}

                            }
                        }
                    }
                }
                free(TDa);
                free(TDb);
                free(TDla);
                free(TDlb);

                //outfile->Printf(" NEWWAY   %5i %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138);

                //double val = 2./3. * eig[state] * (dumx*dumx+dumy*dumy+dumz*dumz);
                //
                //outfile->Printf(" NEWWAY   %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138, val);
                //outfile->Printf("                Transition Moment: X: %10.6lf    Y: %10.6lf   Z: %10.6lf\n", dumx, dumy, dumz);
                
                //bool * skip = (bool*)malloc(nh*sizeof(bool));
                //memset((void*)skip,'\0',nh*sizeof(bool));
                //
                //for (int k = 0; k < number_coefficients; k++) {

                //    double max = -999.0e99;
                //    int maxij = 0;
                //    for (int ij = 0; ij < nh; ij++) {
                //        if ( skip[ij] ) continue;
                //        if ( fabs(Cmat->pointer()[state][ij]) > max ) {
                //            max = fabs(Cmat->pointer()[state][ij]);
                //            maxij = ij;
                //        }
                //    }
                //    skip[maxij] = true;
                //    int id = maxij;
                //    int i = bas_erpa_sym[h][id][0];
                //    int j = bas_erpa_sym[h][id][1];
                //    int hi = symmetry[i];
                //    int hj = symmetry[j];
                //    //outfile->Printf("                %s %3i, symm %3i -> %3i, symm %3i:  %20.6lf\n",(id == maxij) ? "alpha" : "beta ",i, hi, j, hj, Cmat->pointer()[state][maxij]);
                //    //if ( fabs(Cmat->pointer()[state][maxij]) > coefficient_threshold ) {
                //    //    outfile->Printf("                %3i -> %3i  %20.6lf\n",i+1, j+1, Cmat->pointer()[state][maxij]);
                //    //}
                //}

                //free(skip);
                
                //outfile->Printf("\n");                   
            } 
        }
   
        free(newA);
        free(newB);
        free(cc);
        free(eig);
	free(energy_is_real);
    }

    // check energy

    double ac_e2 = 0.0;

    // recompute zeroth-order energy
    e2 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {

                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];

                    e2 +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];//D2ba=D2ab??
                    e2 += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];

                    // ab part ... exclude D2(ij,ij), D2(ii,jj)
                    if ( i == k && j == l ) {
                    }else if ( i == j && k == l ) {
                    }else {
                      ac_e2 += eint * D2ab_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                    }

                    // aa/bb part ... exclude D2(ij,ij), D2(ij,ji)
                    if ( i == k && j == l ) {
                    }else if ( i == l && j == k ) {
                    }else {
                      ac_e2 += 0.5 * eint * D2aa_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                      ac_e2 += 0.5 * eint * D2bb_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                    }

                }
            }
        }
    }

    e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }

    free(alpha_tei_aa);
    free(alpha_tei_ab);
    free(tei);
    free(h1);
    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);
    free(D2aa_prime);
    free(D2bb_prime);
    free(D2ab_prime);

    //outfile->Printf("    w(%lf): %20.12lf %20.12lf\n",alpha,ac_e2,e2);
    outfile->Printf("%20.12lf\n",ac_e2);

    //if ( !is_singlet ) 
    //    return ac_e2;
    //return ac_e2 + e1 + enuc_;
    return ac_e2;
}
double ERPASolver::AlphaERPA(double alpha,bool is_singlet,bool is_excited_state_computation) {

    //outfile->Printf("\n");
    //outfile->Printf("    ==> Alpha ERPA <==\n");
    //outfile->Printf("\n");
    //outfile->Printf("    alpha:    %20.12lf\n",alpha);
    //outfile->Printf("\n");

    if ( !is_excited_state_computation ) {
        outfile->Printf("    %20.12lf",alpha);
    }

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals
    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // build two-electron integrals

    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));

    int info;
    if ( is_df_ ) {
        F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
    }else { 
        // unpack two-electron integrals (which were initially blocked by symmetry)
        long int offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }
    }
    //double * tei_unscaled = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    //C_DCOPY(nn1o2*nn1o2,tei,1,tei_unscaled,1);

    CountNonExternalOrbitals();
    int nne = nmo_no_ext_;

    // read TDPM from disk
    double * D2aa = (double*)malloc(nne*nne*nne*nne*sizeof(double));
    double * D2bb = (double*)malloc(nne*nne*nne*nne*sizeof(double));
    double * D2ab = (double*)malloc(nne*nne*nne*nne*sizeof(double));

    ReadTPDMLowMemory(D2aa,D2bb,D2ab,D1a,D1b);

    // make list of geminals that does not include external orbitals
    int * gems_no_ext   = (int*)malloc(nirrep_*sizeof(int));
    int *** bas_no_ext  = (int***)malloc(nirrep_*sizeof(int**));
    int *** ibas_no_ext = (int***)malloc(nirrep_*sizeof(int**));
    for (int h = 0; h < nirrep_; h++) {
        ibas_no_ext[h]      = (int**)malloc(amo_*sizeof(int*));
        bas_no_ext[h]       = (int**)malloc(amo_*amo_*sizeof(int*));
        for (int i = 0; i < amo_; i++) {
            ibas_no_ext[h][i] = (int*)malloc(amo_*sizeof(int));
            for (int j = 0; j < amo_; j++) {
                ibas_no_ext[h][i][j] = -999;
            }
        }
        for (int i = 0; i < amo_*amo_; i++) {
            bas_no_ext[h][i] = (int*)malloc(2*sizeof(int));
            for (int j = 0; j < 2; j++) {
                bas_no_ext[h][i][j] = -999;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        int n_gems_no_ext = 0;
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
            if ( i_external ) continue;

            bool j_external = ( !active_reference_orbitals_[j] && !inactive_reference_orbitals_[j] );
            if ( j_external ) continue;

            bas_no_ext[h][n_gems_no_ext][0] = i;
            bas_no_ext[h][n_gems_no_ext][1] = j;
            ibas_no_ext[h][i][j] = n_gems_no_ext;

            n_gems_no_ext++;
        }
        gems_no_ext[h] = n_gems_no_ext;
    }

    long int n = nmo_*nmo_;

    // check energy
    double e2_aa = 0.0;
    double e2_ab = 0.0;
    double e2_bb = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        bool k_external = ( !active_reference_orbitals_[k] && !inactive_reference_orbitals_[k] );
        if ( k_external ) continue;
        int kk = full_to_no_ext_[k];
        for (long int l = 0; l < nmo_; l++) {
            bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
            if ( l_external ) continue;
            int ll = full_to_no_ext_[l];
            for (long int i = 0; i < nmo_; i++) {
                bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
                if ( i_external ) continue;
                int ii = full_to_no_ext_[i];
                for (long int j = 0; j < nmo_; j++) {
                    bool j_external = ( !active_reference_orbitals_[j] && !inactive_reference_orbitals_[j] );
                    if ( j_external ) continue;
                    int jj = full_to_no_ext_[j];
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2_ab +=       eint * D2ab[(ii*nne+jj)*nne*nne+(kk*nne+ll)];
                    e2_aa += 0.5 * eint * D2aa[(ii*nne+jj)*nne*nne+(kk*nne+ll)];
                    e2_bb += 0.5 * eint * D2bb[(ii*nne+jj)*nne*nne+(kk*nne+ll)];
                }
            }
        }
    }
    double e2 = e2_aa + e2_ab + e2_bb;

    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }
    Process::environment.globals["v2RDM TOTAL ENERGY"] = e1 + e2 + enuc_;

    //printf("%20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",e1,e2_ab,e2_aa,e2_bb,e1+e2,e1+e2+enuc_);
    //exit(0);

    // build effective hamiltonian
    double * heff = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)heff,'\0',nmo_*nmo_*sizeof(double));
    for (long int i = 0; i < nmo_; i++) {
        int i_space;
        if  ( inactive_reference_orbitals_[i] ){ 
            i_space = 0;
        }else if ( active_reference_orbitals_[i] ){  
            i_space = 1;
        }else{                                      
            i_space = 2;
        }
        for (long int j = 0; j < nmo_; j++) {
            int j_space;
            if  ( inactive_reference_orbitals_[j] ){ 
                j_space = 0;
            }else if ( active_reference_orbitals_[j] ){  
                j_space = 1;
            }else{                                      
                j_space = 2;
            }
            if ( i_space != j_space ) continue;

            heff[i*nmo_+j] = h1[i*nmo_+j];

            double dum = 0.0;
            for (long int k = 0; k < nmo_; k++) {

                int k_space;
                if  ( inactive_reference_orbitals_[k] ){ 
                    k_space = 0;
                }else if ( active_reference_orbitals_[k] ){  
                    k_space = 1;
                }else{                                      
                    k_space = 2;
                }
                if ( i_space == k_space ) continue;

                // TODO: this will only work for closed shells ... if it works at all
                double nk_a = D1a[k*nso_+k];
                double nk_b = D1b[k*nso_+k];

                double dum1 = tei[INDEX(i,j)*nn1o2+INDEX(k,k)];
                double dum2 = tei[INDEX(i,k)*nn1o2+INDEX(j,k)];
                //dum += (nk_a + nk_b)/2.0 * (2.0 * dum1 - dum2);
                dum += nk_a * (dum1 - dum2);
                dum += nk_b *  dum1;

            }
            heff[i*nmo_+j] += dum;
        }
    }

    // scale one-electron integrals
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = 0; j < nmo_; j++) {
            int nactive = 0;
            if ( active_reference_orbitals_[i] ) nactive++;
            if ( active_reference_orbitals_[j] ) nactive++;

            int ninactive = 0;
            if ( inactive_reference_orbitals_[i] ) ninactive++;
            if ( inactive_reference_orbitals_[j] ) ninactive++;

            double factor;
            // if all in same space:   factor = 1.0 - alpha
            // if in different spaces: factor = 0.0
            if ( nactive != 2 && ninactive != 2 && ( nactive + ninactive ) != 0  ){
               factor = 0.0;
	    }else {
               factor = 1.0 - alpha;
	    }

            h1[i*nmo_+j] = alpha * h1[i*nmo_+j] + factor * heff[i*nmo_+j];
        }
    }

    // Gamma(pqrs) = gamma(pr)gamma(qs) + sum_v gamma_0v(pr)gamma_0v(qs) - gamma(qr)delta(ps)

    // D2aa(prqs) = gamma_a(pr)gamma_a(qs) + sum_v gamma_a_0v(pr)gamma_a_0v(qs) - gamma_a(qr)delta(ps)
    // D2bb(prqs) = gamma_b(pr)gamma_b(qs) + sum_v gamma_b_0v(pr)gamma_b_0v(qs) - gamma_b(qr)delta(ps)
    // D2ab(prqs) = gamma_a(pr)gamma_b(qs) + sum_v gamma_a_0v(pr)gamma_b_0v(qs)

    // compute ground-state contribution to AC energy

    double ac_e2 = 0.0;

    if ( is_singlet ) {
        for (int p = 0; p < nso_; p++) {

            int pact   = 0;
            int pinact = 0;
            if ( active_reference_orbitals_[p] ) pact++;
            if ( inactive_reference_orbitals_[p] ) pinact++;

            for (int q = 0; q < nso_; q++) {

                int pqact   = pact;
                int pqinact = pinact;
                if ( active_reference_orbitals_[q] ) pqact++;
                if ( inactive_reference_orbitals_[q] ) pqinact++;

                for (int r = 0; r < nso_; r++) {

                    int pqract   = pqact;
                    int pqrinact = pqinact;
                    if ( active_reference_orbitals_[r] ) pqract++;
                    if ( inactive_reference_orbitals_[r] ) pqrinact++;

                    for (int s = 0; s < nso_; s++) {

                        int nactive = pqract;
                        if ( active_reference_orbitals_[s] ) nactive++;
                        if ( nactive == 4 ) continue;

                        int ninactive = pqrinact;
                        if ( inactive_reference_orbitals_[s] ) ninactive++;
                        if ( ninactive == 4 ) continue;

                        if ( nactive + ninactive == 0 ) continue;

                        double np_a = D1a[p*nso_+p];
                        double nq_a = D1a[q*nso_+q];
        
                        double np_b = D1b[p*nso_+p];
                        double nq_b = D1b[q*nso_+q];

                        double eint = tei[INDEX(p,r)*nn1o2+INDEX(q,s)];

                        ac_e2 += 0.5 * eint * (np_a - 1.0) * nq_a * (r==q)*(p==s);
                        ac_e2 += 0.5 * eint * (np_b - 1.0) * nq_b * (r==q)*(p==s);
                    }
                }
            }
        }
    }

    // transition density matrices

    double * TDa = (double*)malloc(nso_*nso_*sizeof(double));
    double * TDb = (double*)malloc(nso_*nso_*sizeof(double));

    for (int h = 0; h < nirrep_; h++) {

        // scale two-electron integrals
        for (long int i = 0; i < nmo_; i++) {
            for (long int k = i; k < nmo_; k++) {
                for (long int j = 0; j < nmo_; j++) {
                    for (long int l = j; l < nmo_; l++) {
    
                        long int ij = i*nmo_+j;
                        long int kl = k*nmo_+l;
    
                        double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                        int nactive = 0;
                        if ( active_reference_orbitals_[i] ) nactive++;
                        if ( active_reference_orbitals_[j] ) nactive++;
                        if ( active_reference_orbitals_[k] ) nactive++;
                        if ( active_reference_orbitals_[l] ) nactive++;
    
                        int ninactive = 0;
                        if ( inactive_reference_orbitals_[i] ) ninactive++;
                        if ( inactive_reference_orbitals_[j] ) ninactive++;
                        if ( inactive_reference_orbitals_[k] ) ninactive++;
                        if ( inactive_reference_orbitals_[l] ) ninactive++;
    
                        double factor;
                        // if all in same space:   factor = 1.0
                        // if in different spaces: factor = alpha
                        if ( nactive != 4 && ninactive != 4 && ( nactive + ninactive ) != 0  ){
                           factor = alpha;
    		            }else {
                           factor = 1.0;
    		            }
                        tei[INDEX(i,k)*nn1o2+INDEX(j,l)] = factor * eint;
                    }
                }
            }
        }

        long int nh = gems_erpa[h];
        //long int nh = gems_ab[h];
        if ( nh == 0 ) continue;

        long int totdim = nh;
        //if ( type == "AA" ) totdim = 2*nh;

        double * newB = (double*)malloc(totdim*totdim*sizeof(double));
        double * newA = (double*)malloc(totdim*totdim*sizeof(double));
        double * cc   = (double*)malloc(totdim*totdim*sizeof(double));
        double * eig  = (double*)malloc(totdim*sizeof(double));

        memset((void*)newA,'\0',totdim*totdim*sizeof(double));
        memset((void*)newB,'\0',totdim*totdim*sizeof(double));
        memset((void*)cc,'\0',totdim*totdim*sizeof(double));
        memset((void*)eig,'\0',totdim*sizeof(double));

        if ( is_excited_state_computation ) {
            outfile->Printf("\n");
            outfile->Printf("    Symmetry: %5i\n",h);
            outfile->Printf("\n");
        }

        // build B(ij,kl) = <| [ k*l, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {
            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];
            for (long int kl = 0; kl < nh; kl++) {
                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;
                double dumabab = 0.0;
                double dumbaba = 0.0;

                if ( j == l ) {
                    duma    += D1a[k*nmo_+i];
                    dumb    += D1b[k*nmo_+i];
                }
                if ( i == k ) {
                    duma    -= D1a[j*nmo_+l];
                    dumb    -= D1b[j*nmo_+l];
                }

                long int id = ij * nh + kl;
                newB[id] = 0.5 * (duma + dumb);
            }
        }

        // build A1(ij,kl) = <| [ k*l, [H1, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;

                duma    += D1a[k*nmo_+i] * h1[l*nmo_+j];
                duma    += D1a[j*nmo_+l] * h1[i*nmo_+k];

                dumb    += D1b[k*nmo_+i] * h1[l*nmo_+j];
                dumb    += D1b[j*nmo_+l] * h1[i*nmo_+k];

                if ( j == l ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                            dumb    -= D1b[k*nmo_+q] * h1[i*nmo_+q];
                        }
                    }
                }
                if ( i == k ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                            dumb    -= D1b[q*nmo_+l] * h1[q*nmo_+j];
                        }
                    }
                }

                long int id = ij * nh + kl;
                newA[id] += 0.5 * (duma + dumb);
            }
        }

//#if 0
        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>

        // 2 coulomb-like terms:
        // 
        // 1: (ik|qs) D(jq;ls)
        // 2: (lj|qs) D(kq;is)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            bool j_external = ( !active_reference_orbitals_[j] && !inactive_reference_orbitals_[j] );
            if ( j_external ) continue;

            int jne = full_to_no_ext_[j];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
                if ( l_external ) continue;

                int lne = full_to_no_ext_[l];

                double duma = 0.0;
                double dumb = 0.0;

                // 1: (ik|qs) D(jq;ls)
                
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];

                        bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                        if ( q_external ) continue;
                        int qne = full_to_no_ext_[q];

                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];

                                bool s_external = ( !active_reference_orbitals_[s] && !inactive_reference_orbitals_[s] );
                                if ( s_external ) continue;
                                int sne = full_to_no_ext_[s];

                                double eint = tei[INDEX(i,k)*nn1o2+INDEX(q,s)];

                                long int jq = jne*nne+qne;
                                long int ls = lne*nne+sne;

                                long int qj = qne*nne+jne;
                                long int sl = sne*nne+lne;

                                duma    += eint * ( D2aa[jq*nne*nne+ls] + D2ab[jq*nne*nne+ls] );
                                dumb    += eint * ( D2bb[jq*nne*nne+ls] + D2ab[qj*nne*nne+sl] );
                            }
                        }
                    }
                }
                // 2: (lj|qs) D(kq;is)

                long int id = ij * nh + kl;
                newA[id] += 0.5 * ( duma + dumb );

                // AED: I think term (2) above is just a funny transpose of term (1)
                // A(ij,kl)(2) = A(lk,ji)(1)

                long int ji = ibas_erpa_sym[h][j][i];
                long int lk = ibas_erpa_sym[h][l][k];

                id = ji * nh + lk;
                newA[id] += 0.5 * (duma + dumb);

            }
        }

        // funky sums over 3 indices:
        // 
        // 1: -dik (qs|pj) D(pq;ls)
        // 2: -djl (qs|ip) D(kq;ps)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
                if ( l_external ) continue;
                int lne = full_to_no_ext_[l];

                double duma = 0.0;
                double dumb = 0.0;

                // funky sums over 3 indices:
                // - dik (qs|pj) D(pq;ls) + djl (qs|ip) D(kq;ps)
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];

                            bool p_external = ( !active_reference_orbitals_[p] && !inactive_reference_orbitals_[p] );
                            if ( p_external ) continue;
                            int pne = full_to_no_ext_[p];

                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];

                                    bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                                    if ( q_external ) continue;
                                    int qne = full_to_no_ext_[q];

                                    for (int hs = 0; hs < nirrep_; hs++) {
                                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                            long int s = ss + pitzer_offset_full[hs];

                                            bool s_external = ( !active_reference_orbitals_[s] && !inactive_reference_orbitals_[s] );
                                            if ( s_external ) continue;
                                            int sne = full_to_no_ext_[s];

                                            double eint = tei[INDEX(q,s)*nn1o2+INDEX(p,j)];

                                            long int pq = pne*nne+qne;
                                            long int ls = lne*nne+sne;

                                            long int qp = qne*nne+pne;
                                            long int sl = sne*nne+lne;

                                            duma -= eint * ( D2aa[pq*nne*nne+ls] + D2ab[pq*nne*nne+ls] );
                                            dumb -= eint * ( D2bb[pq*nne*nne+ls] + D2ab[qp*nne*nne+sl] );
                                        }
                                    }
                                }
                            }
                        }
                    }

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * ( duma + dumb );
    
                    // AED: I think second funky term is a transpose of the first
                    // A(ij,il)(2) = A(ji,li)(1)
    
                    long int ji = ibas_erpa_sym[h][j][i];
                    long int li = ibas_erpa_sym[h][l][i];
    
                    id = ji * nh + li;
                    newA[id] += 0.5 * ( duma + dumb );
                }

            }
        }

        // A2 exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)
        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_no_ext[hjl]; ik++) {

                    long int i = bas_no_ext[hjl][ik][0];
                    long int k = bas_no_ext[hjl][ik][1];

                    long int ine = full_to_no_ext_[i];
                    long int kne = full_to_no_ext_[k];

                    //long int i = bas_ab_sym[hjl][ik][0];
                    //long int k = bas_ab_sym[hjl][ik][1];

                    //bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
                    //bool k_external = ( !active_reference_orbitals_[k] && !inactive_reference_orbitals_[k] );
                    //if ( k_external ) continue;
                    //if ( i_external ) continue;

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    for (int qs = 0; qs < gems_no_ext[hjl]; qs++ ) {
                        
                        long int q = bas_no_ext[hjl][qs][0];
                        long int s = bas_no_ext[hjl][qs][1];

                        long int qne = full_to_no_ext_[q];
                        long int sne = full_to_no_ext_[s];

                        //long int q = bas_ab_sym[hjl][qs][0];
                        //bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                        //if ( q_external ) continue;

                        //long int s = bas_ab_sym[hjl][qs][1];
                        //bool s_external = ( !active_reference_orbitals_[s] && !inactive_reference_orbitals_[s] );
                        //if ( s_external ) continue;

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int kq = kne*nne+qne;
                        long int is = ine*nne+sne;

                        duma -= eint * D2aa[kq*nne*nne+is];
                        dumb -= eint * D2bb[kq*nne*nne+is];

                        //eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                        //long int jq = j*nmo_+q;
                        //long int ls = l*nmo_+s;

                        //duma -= eint * D2aa[jq*n+ls];
                        //dumb -= eint * D2bb[jq*n+ls];
                    }

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * (duma + dumb);

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    newA[id] += 0.5 * (duma + dumb);
                }
            }
        }

        // A2 exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_no_ext[hjl]; ik++) {

                    long int i = bas_no_ext[hjl][ik][0];
                    long int k = bas_no_ext[hjl][ik][1];

                    long int ine = full_to_no_ext_[i];
                    long int kne = full_to_no_ext_[k];

                    //long int i = bas_ab_sym[hjl][ik][0];
                    //long int k = bas_ab_sym[hjl][ik][1];

                    //bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
                    //bool k_external = ( !active_reference_orbitals_[k] && !inactive_reference_orbitals_[k] );
                    //if ( k_external ) continue;
                    //if ( i_external ) continue;

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int qs = 0; qs < gems_no_ext[hjl]; qs++ ) {
                        
                        long int q = bas_no_ext[hjl][qs][0];
                        long int s = bas_no_ext[hjl][qs][1];

                        long int qne = full_to_no_ext_[q];
                        long int sne = full_to_no_ext_[s];
                        
                        //long int q = bas_ab_sym[hjl][qs][0];
                        //bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                        //if ( q_external ) continue;

                        //long int s = bas_ab_sym[hjl][qs][1];
                        //bool s_external = ( !active_reference_orbitals_[s] && !inactive_reference_orbitals_[s] );
                        //if ( s_external ) continue;

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int qk = qne*nne+kne;
                        long int is = ine*nne+sne;

                        long int kq = kne*nne+qne;
                        long int si = sne*nne+ine;

                        dumab += eint * D2ab[qk*nne*nne+is];
                        dumba += eint * D2ab[kq*nne*nne+si];

                    }

                    long int id = ij * nh + kl;
                    if ( is_singlet ) {
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    if ( is_singlet ) {
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                }
            }
        }

        // A2: last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_no_ext[hli]; li++) {

                long int l = bas_no_ext[hli][li][0];
                long int i = bas_no_ext[hli][li][1];

                long int lne = full_to_no_ext_[l];
                long int ine = full_to_no_ext_[i];

                //long int l = bas_ab_sym[hli][li][0];
                //long int i = bas_ab_sym[hli][li][1];

                //bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
                //bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
                //if ( l_external ) continue;
                //if ( i_external ) continue;

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int pq = 0; pq < gems_no_ext[hli]; pq++ ) {

                        long int p = bas_no_ext[hli][pq][0];
                        long int q = bas_no_ext[hli][pq][1];

                        long int pne = full_to_no_ext_[p];
                        long int qne = full_to_no_ext_[q];

                        //long int p = bas_ab_sym[hli][pq][0];
                        //bool p_external = ( !active_reference_orbitals_[p] && !inactive_reference_orbitals_[p] );
                        //if ( p_external ) continue;

                        //long int q = bas_ab_sym[hli][pq][1];
                        //bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                        //if ( q_external ) continue;

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        duma += eint * D2aa[pne*nne*nne*nne+qne*nne*nne+lne*nne+ine];
                        dumb += eint * D2bb[pne*nne*nne*nne+qne*nne*nne+lne*nne+ine];

                    }

                    long int id = ij * nh + kl;

                    newA[id] += 0.5 * (duma + dumb);

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitaitons, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    newA[id] += 0.5 * (duma + dumb);

                }

            }
        }

        // A2: last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_no_ext[hli]; li++) {

                long int l = bas_no_ext[hli][li][0];
                long int i = bas_no_ext[hli][li][1];

                long int lne = full_to_no_ext_[l];
                long int ine = full_to_no_ext_[i];

                //long int l = bas_ab_sym[hli][li][0];
                //long int i = bas_ab_sym[hli][li][1];

                //bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
                //bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
                //if ( l_external ) continue;
                //if ( i_external ) continue;

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];


                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int pq = 0; pq < gems_no_ext[hli]; pq++ ) {

                        long int p = bas_no_ext[hli][pq][0];
                        long int q = bas_no_ext[hli][pq][1];

                        long int pne = full_to_no_ext_[p];
                        long int qne = full_to_no_ext_[q];

                        //long int p = bas_ab_sym[hli][pq][0];
                        //bool p_external = ( !active_reference_orbitals_[p] && !inactive_reference_orbitals_[p] );
                        //if ( p_external ) continue;

                        //long int q = bas_ab_sym[hli][pq][1];
                        //bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                        //if ( q_external ) continue;

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        dumab -= eint * D2ab[pne*nne*nne*nne+qne*nne*nne+ine*nne+lne];
                        dumba -= eint * D2ab[qne*nne*nne*nne+pne*nne*nne+lne*nne+ine];

                    }

                    long int id = ij * nh + kl;
                    if ( is_singlet ) {
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    if ( is_singlet ) {
                        newA[id] += 0.5 * (dumab + dumba);
                    }else {
                        newA[id] -= 0.5 * (dumab + dumba);
                    }

                }

            }

        }

        std::shared_ptr<Matrix> Amat ( new Matrix(totdim,totdim) );
        std::shared_ptr<Matrix> Bmat ( new Matrix(totdim,totdim) );
        for (int ij = 0; ij < totdim; ij++) {
            for (int kl = 0; kl < totdim; kl++) {
                Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                Bmat->pointer()[ij][kl] = newB[ij*totdim+kl];
            }
        }

        // symmetrize
        // A(kl,ji) = A(ij,kl)
        // A(kl,ji) = A(lk,ij)

        // project out null (and negative) space of A matrix
        std::shared_ptr<Matrix> eigvec ( new Matrix(totdim,totdim) );
        std::shared_ptr<Vector> eigval ( new Vector(totdim) );

        Amat->diagonalize(eigvec,eigval,descending);
        //eigval->print();

        bool prune = true;
        long int newdim = 0;

        // get rid of small/negative eigenvalues of A?
        if ( prune ) {

            newdim = 0;
            int * map = (int*)malloc(totdim*sizeof(int));
            memset((void*)map,'\0',totdim*sizeof(int));
            int * skip = (int*)malloc(totdim*sizeof(int));
            memset((void*)skip,'\0',totdim*sizeof(int));
            for (long int i = 0; i < totdim; i++) {
                double val = eigval->pointer()[i];
                skip[i] = 1;
                if ( val > pruning_threshold ) {
                    map[newdim] = i;
                    skip[i] = 0;
                    newdim++;
                }
            }
            //if ( newdim % 2 == 1 ) { // eigenvalues should come in pairs
            //    newdim--;
            //}
            for (int ij = 0; ij < totdim; ij++) {
                for (int kl = 0; kl < totdim; kl++) {
                    Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                }
            }
            Amat->transform(eigvec);
            Bmat->transform(eigvec);
        }else {
            newdim = totdim;
        }

        for (long int i = 0; i < newdim; i++) {
            for (long int j = 0; j < newdim; j++) {
                newA[i*newdim+j] = Amat->pointer()[i][j];
                newB[i*newdim+j] = Bmat->pointer()[i][j];
            }
        }

        if ( is_excited_state_computation  ){
            if ( newdim < totdim ) {
                outfile->Printf("\n");
                outfile->Printf("    <<< warning >>> ");
                outfile->Printf("\n");
                outfile->Printf("    reducing dimension of [ k*l, [H, j*i ] ] from %5i to %5i\n",totdim,newdim);
                outfile->Printf("\n");
            }
        }

        int info = 0;
        if ( prune ) {
            info = SymmetricGeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }else {
            bool * energy_is_real = (bool*)malloc(totdim*sizeof(bool));
            GeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig,energy_is_real);
            free(energy_is_real);
        }

        // back transform wave functions to original untruncated basis
        std::shared_ptr<Matrix>  Cmat (new Matrix(totdim,totdim));
        for (int i = 0; i < newdim; i++) {
            for (int j = 0; j < totdim; j++) {
                double dum = 0.0;
                for (int k = 0; k < newdim; k++) {
                    dum += eigvec->pointer()[j][k] * newB[i * newdim + k];
                }
                Cmat->pointer()[i][j] = dum;
            }
        }

        // oscillator strengths: 
        // symmetry of state is h, consider only cartesian component with same symmetry
       
        if ( info != 0 ) {
            outfile->Printf("\n");
            outfile->Printf("    Error: diagonalization failed.  info = %5i\n",info);
            outfile->Printf("\n");
            throw PsiException("Error: diagonalization failed.",__FILE__,__LINE__);
        }

        if ( is_excited_state_computation ) {
            outfile->Printf("    state");
            outfile->Printf("          energy (Eh)");
            outfile->Printf("       ex energy (Eh)");
            outfile->Printf("       ex energy (eV)");
            outfile->Printf("       f, osc. strength\n");
        }

        // rebuild unscaled two-electron integrals
        if ( is_df_ ) {
            F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
        }else {
            // unpack two-electron integrals (which were initially blocked by symmetry)
            long int offset = 0;
            for (int h = 0; h < nirrep_; h++) {
                for (long int ij = 0; ij < gems_full[h]; ij++) {
                    long int i = bas_really_full_sym[h][ij][0];
                    long int j = bas_really_full_sym[h][ij][1];
                    for (long int kl = 0; kl < gems_full[h]; kl++) {
                        long int k = bas_really_full_sym[h][kl][0];
                        long int l = bas_really_full_sym[h][kl][1];
                        tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                    }
                }
                offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
            }
        }

        for (long int state = newdim-1; state >= 0; state--){
          
            if ( eig[state] > 0.0 ) {

                if ( prune ) eig[state] = 1.0 / eig[state];

                // avoid considering nonsensical states
                if ( fabs(eig[state] ) > 10000.0 ) continue;

                //printf("hey excited state: %20.12lf\n",eig[state]+e2+e1+enuc_);

                // normalize C:
                double nrm = 0.0;

                double * C = Cmat->pointer()[state];
                for (int ij = 0; ij < nh; ij++) {
                    int i  = bas_erpa_sym[h][ij][0];
                    int j  = bas_erpa_sym[h][ij][1];
                    int hi = symmetry[i];
                    int hj = symmetry[j];
                    for( int kk = 0; kk < amopi_[hi]; kk++){
                        int k   = kk + pitzer_offset[hi];
                        int kj  = ibas_erpa_sym[h][k][j];

                        // check if kj is a proper ERPA excitations
                        if ( kj == -999 ) continue;

                        nrm    += 0.5 * C[ij] * C[kj] * D1a[k*amo_+i];
                        nrm    += 0.5 * C[ij] * C[kj] * D1b[k*amo_+i];
                    }
                    for( int ll = 0; ll < amopi_[hj]; ll++){
                        int l  = ll + pitzer_offset[hj];
                        int il = ibas_erpa_sym[h][i][l];

                        // check if il is a proper ERPA excitations
                        if ( il == -999 ) continue;

                        nrm   -= 0.5 * C[ij] * C[il] * D1a[j*amo_+l];
                        nrm   -= 0.5 * C[ij] * C[il] * D1b[j*amo_+l];
                    }
                }

                if ( nrm < 0.0 ) continue;

                nrm = sqrt(nrm);
                for (int ij = 0; ij < nh; ij++) {
                    C[ij] /= nrm;
                }


                // compute transition dipole moments:

                double dumx = 0.0;
                double dumy = 0.0;
                double dumz = 0.0;

                if ( h == hx_ ) {

                    dumx = TransitionDipoleMoment_SpinAdapted(C, D1a, D1b, nh,h,"X");

                }
                if ( h == hy_ ) {

                    dumy = TransitionDipoleMoment_SpinAdapted(C, D1a, D1b, nh,h,"Y");

                }
                if ( h == hz_ ) {

                    dumz = TransitionDipoleMoment_SpinAdapted(C, D1a, D1b, nh,h,"Z");

                }
                
                // use transition density matrix to update ground-state 2-RDM
                if ( !is_excited_state_computation ) {

                    memset((void*)TDa,'\0',nso_*nso_*sizeof(double));
                    memset((void*)TDb,'\0',nso_*nso_*sizeof(double));

                    // build transition density matrix:
                    TransitionDensity_SpinAdapted(C, D1a, D1b, TDa,TDb,nh,h);
                    if ( !is_singlet ) {
                       C_DSCAL(nso_*nso_,-1.0,TDb,1);
                    }

                    for (int p = 0; p < nso_; p++) {
            
                        int pact   = 0;
                        int pinact = 0;
                        if ( active_reference_orbitals_[p] ) pact++;
                        if ( inactive_reference_orbitals_[p] ) pinact++;
            
                        for (int q = 0; q < nso_; q++) {
            
                            int pqact   = pact;
                            int pqinact = pinact;
                            if ( active_reference_orbitals_[q] ) pqact++;
                            if ( inactive_reference_orbitals_[q] ) pqinact++;
            
                            for (int r = 0; r < nso_; r++) {
            
                                int pqract   = pqact;
                                int pqrinact = pqinact;
                                if ( active_reference_orbitals_[r] ) pqract++;
                                if ( inactive_reference_orbitals_[r] ) pqrinact++;
            
                                for (int s = 0; s < nso_; s++) {
            
                                    int nactive = pqract;
                                    if ( active_reference_orbitals_[s] ) nactive++;
                                    if ( nactive == 4 ) continue;
            
                                    int ninactive = pqrinact;
                                    if ( inactive_reference_orbitals_[s] ) ninactive++;
                                    if ( ninactive == 4 ) continue;
            
                                    if ( nactive + ninactive == 0 ) continue;

                                    double ga_rp = TDa[p*nso_+r]; // careful with transpose ...
                                    double gb_qs = TDb[s*nso_+q];
                                    double ga_qs = TDa[s*nso_+q];
                                    double gb_rp = TDb[p*nso_+r]; // careful with transpose ...

                                    //double eint = tei_unscaled[INDEX(p,r)*nn1o2+INDEX(q,s)];
                                    double eint = tei[INDEX(p,r)*nn1o2+INDEX(q,s)];

                                    ac_e2 +=       eint * ga_rp * gb_qs;
                                    ac_e2 += 0.5 * eint * ga_rp * ga_qs;
                                    ac_e2 += 0.5 * eint * gb_rp * gb_qs;

                                }
                            }
                        }
                    }
                }else {

                    //outfile->Printf(" NEWWAY   %5i %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138);

                    // evaluate violations in the consistency condition [1-d(tu) ]n(t) c(ut) + sum_{v!=w} [ 2D(tv;uw) + 2D(t\bar{v};u\bar{w}) ] c(vw).
                    double violation = 0.0;
                    int nact = 0;
                    for (int t = 0; t < nso_; t++) {
                        if ( active_reference_orbitals_[t] ) nact++;
                    }
                    for (int t = 0; t < nso_; t++) {
                        if ( !active_reference_orbitals_[t] ) continue;
                        for (int u = 0; u < nso_; u++) {
                            if ( !active_reference_orbitals_[u] ) continue;

                            double sum = 0.0;
                            if ( t != u ) {
                                int ut = ibas_erpa_sym[h][u][t];
                                sum += D1a[t*nso_+t] * C[ut];
                            }
                            for (int v = 0; v < nso_; v++) {
                                if ( !active_reference_orbitals_[v] ) continue;
                                for (int w = 0; w < nso_; w++) {
                                    if ( !active_reference_orbitals_[w] ) continue;
                                    if ( w == v ) continue;

                                    long int tne = full_to_no_ext_[t];
                                    long int vne = full_to_no_ext_[v];
                                    long int une = full_to_no_ext_[u];
                                    long int wne = full_to_no_ext_[w];

                                    double d2aa = D2aa[tne*nne*nne*nne+vne*nne*nne+une*nne+wne];
                                    double d2ab = D2ab[tne*nne*nne*nne+vne*nne*nne+une*nne+wne];
                                    int vw = ibas_erpa_sym[h][v][w];
                                    sum += (d2aa + d2ab) * C[vw];
                                }
                            }
                            violation += sum * sum;
                        }
                    }

                    double val = 2./3. * eig[state] * (dumx*dumx+dumy*dumy+dumz*dumz);

                    outfile->Printf(" LOWMEM   %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138, val);
                    outfile->Printf("                Transition Moment: X: %10.6lf    Y: %10.6lf   Z: %10.6lf\n", dumx, dumy, dumz);
                    outfile->Printf("                Consistency violation: %20.12lf\n",sqrt(violation/(nact*nact)));

                    bool * skip = (bool*)malloc(nh*sizeof(bool));
                    memset((void*)skip,'\0',nh*sizeof(bool));

                    for (int k = 0; k < number_coefficients; k++) {

                        double max = -999.0e99;
                        int maxij = 0;
                        for (int ij = 0; ij < nh; ij++) {
                            if ( skip[ij] ) continue;
                            if ( fabs(C[ij]) > max ) {
                                max = fabs(C[ij]);
                                maxij = ij;
                            }
                        }
                        skip[maxij] = true;
                        int id = maxij;
                        int i = bas_erpa_sym[h][id][0];
                        int j = bas_erpa_sym[h][id][1];
                        int hi = symmetry[i];
                        int hj = symmetry[j];
                        //outfile->Printf("                %s %3i, symm %3i -> %3i, symm %3i:  %20.6lf\n",(id == maxij) ? "alpha" : "beta ",i, hi, j, hj, Cmat->pointer()[state][maxij]);
                        if ( fabs(C[maxij]) > coefficient_threshold ) {
                            outfile->Printf("                %3i -> %3i  %20.6lf\n",i+1, j+1, C[maxij]);
                        }
                    }

                    free(skip);

                    outfile->Printf("\n");
                }
            } 
        }
   
        free(newA);
        free(newB);
        free(cc);
        free(eig);
    }

    free(TDa);
    free(TDb);
    free(tei);
    //free(tei_unscaled);
    free(heff);
    free(h1);
    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);

    //outfile->Printf("    w(%lf): %20.12lf\n",alpha,ac_e2);
    if ( !is_excited_state_computation ) {
        outfile->Printf("%20.12lf\n",ac_e2);
    }

    return ac_e2;
}
/*void ERPASolver::SpinFlipExtendedRPA() {

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals
    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // build two-electron integrals

    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));

    int info;
    if ( is_df_ ) {
        F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
    }else { 
        // unpack two-electron integrals (which were initially blocked by symmetry)
        long int offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }
    }
    //double * tei_unscaled = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    //C_DCOPY(nn1o2*nn1o2,tei,1,tei_unscaled,1);

    CountNonExternalOrbitals();
    int nne = nmo_no_ext_;

    // read TDPM from disk
    double * D2aa = (double*)malloc(nne*nne*nne*nne*sizeof(double));
    double * D2bb = (double*)malloc(nne*nne*nne*nne*sizeof(double));
    double * D2ab = (double*)malloc(nne*nne*nne*nne*sizeof(double));

    ReadTPDMLowMemory(D2aa,D2bb,D2ab,D1a,D1b);

    // make list of geminals that does not include external orbitals
    int * gems_no_ext   = (int*)malloc(nirrep_*sizeof(int));
    int *** bas_no_ext  = (int***)malloc(nirrep_*sizeof(int**));
    int *** ibas_no_ext = (int***)malloc(nirrep_*sizeof(int**));
    for (int h = 0; h < nirrep_; h++) {
        ibas_no_ext[h]      = (int**)malloc(amo_*sizeof(int*));
        bas_no_ext[h]       = (int**)malloc(amo_*amo_*sizeof(int*));
        for (int i = 0; i < amo_; i++) {
            ibas_no_ext[h][i] = (int*)malloc(amo_*sizeof(int));
            for (int j = 0; j < amo_; j++) {
                ibas_no_ext[h][i][j] = -999;
            }
        }
        for (int i = 0; i < amo_*amo_; i++) {
            bas_no_ext[h][i] = (int*)malloc(2*sizeof(int));
            for (int j = 0; j < 2; j++) {
                bas_no_ext[h][i][j] = -999;
            }
        }
    }
    for (int h = 0; h < nirrep_; h++) {
        int n_gems_no_ext = 0;
        for (int ij = 0; ij < gems_ab[h]; ij++) {
            int i = bas_ab_sym[h][ij][0];
            int j = bas_ab_sym[h][ij][1];

            bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
            if ( i_external ) continue;

            bool j_external = ( !active_reference_orbitals_[j] && !inactive_reference_orbitals_[j] );
            if ( j_external ) continue;

            bas_no_ext[h][n_gems_no_ext][0] = i;
            bas_no_ext[h][n_gems_no_ext][1] = j;
            ibas_no_ext[h][i][j] = n_gems_no_ext;

            n_gems_no_ext++;
        }
        gems_no_ext[h] = n_gems_no_ext;
    }

    long int n = nmo_*nmo_;

    // check energy
    double e2_aa = 0.0;
    double e2_ab = 0.0;
    double e2_bb = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        bool k_external = ( !active_reference_orbitals_[k] && !inactive_reference_orbitals_[k] );
        if ( k_external ) continue;
        int kk = full_to_no_ext_[k];
        for (long int l = 0; l < nmo_; l++) {
            bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
            if ( l_external ) continue;
            int ll = full_to_no_ext_[l];
            for (long int i = 0; i < nmo_; i++) {
                bool i_external = ( !active_reference_orbitals_[i] && !inactive_reference_orbitals_[i] );
                if ( i_external ) continue;
                int ii = full_to_no_ext_[i];
                for (long int j = 0; j < nmo_; j++) {
                    bool j_external = ( !active_reference_orbitals_[j] && !inactive_reference_orbitals_[j] );
                    if ( j_external ) continue;
                    int jj = full_to_no_ext_[j];
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2_ab +=       eint * D2ab[(ii*nne+jj)*nne*nne+(kk*nne+ll)];
                    e2_aa += 0.5 * eint * D2aa[(ii*nne+jj)*nne*nne+(kk*nne+ll)];
                    e2_bb += 0.5 * eint * D2bb[(ii*nne+jj)*nne*nne+(kk*nne+ll)];
                }
            }
        }
    }
    double e2 = e2_aa + e2_ab + e2_bb;

    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }
    Process::environment.globals["v2RDM TOTAL ENERGY"] = e1 + e2 + enuc_;

    for (int h = 0; h < nirrep_; h++) {

        long int nh = gems_erpa[h];
        //long int nh = gems_ab[h];
        if ( nh == 0 ) continue;

        long int totdim = nh;
        //if ( type == "AA" ) totdim = 2*nh;

        double * newB = (double*)malloc(totdim*totdim*sizeof(double));
        double * newA = (double*)malloc(totdim*totdim*sizeof(double));
        double * cc   = (double*)malloc(totdim*totdim*sizeof(double));
        double * eig  = (double*)malloc(totdim*sizeof(double));

        memset((void*)newA,'\0',totdim*totdim*sizeof(double));
        memset((void*)newB,'\0',totdim*totdim*sizeof(double));
        memset((void*)cc,'\0',totdim*totdim*sizeof(double));
        memset((void*)eig,'\0',totdim*sizeof(double));

        outfile->Printf("\n");
        outfile->Printf("    Symmetry: %5i\n",h);
        outfile->Printf("\n");

        // build B(ij,kl) = <| [ k*l, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {
            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];
            for (long int kl = 0; kl < nh; kl++) {
                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double dumabab = 0.0;
                double dumbaba = 0.0;

                if ( j == l ) {
                    dumabab += D1a[k*nmo_+i];
                    dumbaba += D1b[k*nmo_+i];
                }
                if ( i == k ) {
                   dumabab -= D1b[j*nmo_+l];
                   dumbaba -= D1a[j*nmo_+l];
                }

                long int id = ij * nh + kl;

                newB[id] = dumabab;
                //newB[id] = dumbaba;
            }
        }

        // build A1(ij,kl) = <| [ k*l, [H1, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double dumbaba    = 0.0;
                double dumabab    = 0.0;

                dumbaba += D1b[k*nmo_+i] * h1[l*nmo_+j];
                dumbaba += D1a[j*nmo_+l] * h1[i*nmo_+k];

                dumabab += D1a[k*nmo_+i] * h1[l*nmo_+j];
                dumabab += D1b[j*nmo_+l] * h1[i*nmo_+k];

                if ( j == l ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            dumabab -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                            dumbaba -= D1b[k*nmo_+q] * h1[i*nmo_+q];
                        }
                    }
                }
                if ( i == k ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            dumabab -= D1b[q*nmo_+l] * h1[q*nmo_+j];
                            dumbaba -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                        }
                    }
                }

                long int id = ij * nh + kl;
                newA[id] += dumabab;
                //newA[id] += dumbaba;
            }
        }

        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>

        // 2 coulomb-like terms:
        // 
        // 1: (ik|qs) D(jq;ls)
        // 2: (lj|qs) D(kq;is)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            bool j_external = ( !active_reference_orbitals_[j] && !inactive_reference_orbitals_[j] );
            if ( j_external ) continue;

            int jne = full_to_no_ext_[j];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
                if ( l_external ) continue;

                int lne = full_to_no_ext_[l];

                double dumbaba = 0.0;
                double dumabab = 0.0;

                // 1: (ik|qs) D(jq;ls)
                
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];

                        bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                        if ( q_external ) continue;
                        int qne = full_to_no_ext_[q];

                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];

                                bool s_external = ( !active_reference_orbitals_[s] && !inactive_reference_orbitals_[s] );
                                if ( s_external ) continue;
                                int sne = full_to_no_ext_[s];

                                double eint = tei[INDEX(i,k)*nn1o2+INDEX(q,s)];

                                long int jq = jne*nne+qne;
                                long int ls = lne*nne+sne;

                                long int qj = qne*nne+jne;
                                long int sl = sne*nne+lne;

                                dumabab += eint * ( D2bb[jq*nne*nne+ls] + D2ab[qj*nne*nne+sl] );
                                dumbaba += eint * ( D2aa[jq*nne*nne+ls] + D2ab[jq*nne*nne+ls] );
                            }
                        }
                    }
                }
                // 2: (lj|qs) D(kq;is)

                long int id = ij * nh + kl;
                newA[id] += dumabab;
                //newA[id] += dumbaba;

                // AED: I think term (2) above is just a funny transpose of term (1)
                // A(ij,kl)(2) = A(lk,ji)(1)

                long int ji = ibas_erpa_sym[h][j][i];
                long int lk = ibas_erpa_sym[h][l][k];

                id = ji * nh + lk;
                newA[id] += dumabab;
                //newA[id] += dumbaba;

            }
        }

        // funky sums over 3 indices:
        // 
        // 1: -dik (qs|pj) D(pq;ls)
        // 2: -djl (qs|ip) D(kq;ps)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                bool l_external = ( !active_reference_orbitals_[l] && !inactive_reference_orbitals_[l] );
                if ( l_external ) continue;
                int lne = full_to_no_ext_[l];

                double dumbaba = 0.0;
                double dumabab = 0.0;

                // funky sums over 3 indices:
                // - dik (qs|pj) D(pq;ls) + djl (qs|ip) D(kq;ps)
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];

                            bool p_external = ( !active_reference_orbitals_[p] && !inactive_reference_orbitals_[p] );
                            if ( p_external ) continue;
                            int pne = full_to_no_ext_[p];

                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];

                                    bool q_external = ( !active_reference_orbitals_[q] && !inactive_reference_orbitals_[q] );
                                    if ( q_external ) continue;
                                    int qne = full_to_no_ext_[q];

                                    for (int hs = 0; hs < nirrep_; hs++) {
                                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                            long int s = ss + pitzer_offset_full[hs];

                                            bool s_external = ( !active_reference_orbitals_[s] && !inactive_reference_orbitals_[s] );
                                            if ( s_external ) continue;
                                            int sne = full_to_no_ext_[s];

                                            double eint = tei[INDEX(q,s)*nn1o2+INDEX(p,j)];

                                            long int pq = pne*nne+qne;
                                            long int ls = lne*nne+sne;

                                            long int qp = qne*nne+pne;
                                            long int sl = sne*nne+lne;

                                            dumabab -= eint * ( D2bb[pq*nne*nne+ls] + D2ab[qp*nne*nne+sl] );
                                            dumbaba -= eint * ( D2aa[pq*nne*nne+ls] + D2ab[pq*nne*nne+ls] );
                                        }
                                    }
                                }
                            }
                        }
                    }

                    long int id = ij * nh + kl;
                    newA[id] += dumabab;
                    //newA[id] += dumbaba;
    
                    // AED: I think second funky term is a transpose of the first
                    // A(ij,il)(2) = A(ji,li)(1)
    
                    long int ji = ibas_erpa_sym[h][j][i];
                    long int li = ibas_erpa_sym[h][l][i];
    
                    id = ji * nh + li;
                    newA[id] += dumabab;
                    //newA[id] += dumbaba;
                }

            }
        }

        // A2 exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)
        // i do not think there are any ba/ba contributions 

        // A2 exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)
        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_no_ext[hjl]; ik++) {

                    long int i = bas_no_ext[hjl][ik][0];
                    long int k = bas_no_ext[hjl][ik][1];

                    long int ine = full_to_no_ext_[i];
                    long int kne = full_to_no_ext_[k];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumabab = 0.0;
                    double dumbaba = 0.0;

                    for (int qs = 0; qs < gems_no_ext[hjl]; qs++ ) {
                        
                        long int q = bas_no_ext[hjl][qs][0];
                        long int s = bas_no_ext[hjl][qs][1];

                        long int qne = full_to_no_ext_[q];
                        long int sne = full_to_no_ext_[s];
                        
                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int qk = qne*nne+kne;
                        long int is = ine*nne+sne;

                        long int kq = kne*nne+qne;
                        long int si = sne*nne+ine;

                        dumabab -= eint * D2ab[kq*nne*nne+is];
                        dumbaba -= eint * D2ab[qk*nne*nne+si];

                    }

                    long int id = ij * nh + kl;
                    newA[id] += dumabab;
                    //newA[id] += dumbaba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    newA[id] += dumabab;
                    //newA[id] += dumbaba;

                }
            }
        }

        // A2: last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
        // i do not think there are any ba/ba contributions 

        // A2: last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_no_ext[hli]; li++) {

                long int l = bas_no_ext[hli][li][0];
                long int i = bas_no_ext[hli][li][1];

                long int lne = full_to_no_ext_[l];
                long int ine = full_to_no_ext_[i];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];


                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumbaba = 0.0;
                    double dumabab = 0.0;

                    for (int pq = 0; pq < gems_no_ext[hli]; pq++ ) {

                        long int p = bas_no_ext[hli][pq][0];
                        long int q = bas_no_ext[hli][pq][1];

                        long int pne = full_to_no_ext_[p];
                        long int qne = full_to_no_ext_[q];

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                       
                        dumabab += eint * D2ab[qne*nne*nne*nne+pne*nne*nne+ine*nne+lne];
                        dumbaba += eint * D2ab[pne*nne*nne*nne+qne*nne*nne+lne*nne+ine];
                    }

                    long int id = ij * nh + kl;
                    newA[id] += dumabab;
                    //newA[id] += dumbaba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    newA[id] += dumabab;
                    //newA[id] += dumbaba;
                }

            }

        }

        std::shared_ptr<Matrix> Amat ( new Matrix(totdim,totdim) );
        std::shared_ptr<Matrix> Bmat ( new Matrix(totdim,totdim) );
        for (int ij = 0; ij < totdim; ij++) {
            for (int kl = 0; kl < totdim; kl++) {
                Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                Bmat->pointer()[ij][kl] = newB[ij*totdim+kl];
            }
        }

        // symmetrize
        // A(kl,ji) = A(ij,kl)
        // A(kl,ji) = A(lk,ij)

        // project out null (and negative) space of A matrix
        std::shared_ptr<Matrix> eigvec ( new Matrix(totdim,totdim) );
        std::shared_ptr<Vector> eigval ( new Vector(totdim) );

        Amat->diagonalize(eigvec,eigval,descending);
        //eigval->print();

        bool prune = false;//true;
        long int newdim = 0;

        // get rid of small/negative eigenvalues of A?
        if ( prune ) {

            newdim = 0;
            int * map = (int*)malloc(totdim*sizeof(int));
            memset((void*)map,'\0',totdim*sizeof(int));
            int * skip = (int*)malloc(totdim*sizeof(int));
            memset((void*)skip,'\0',totdim*sizeof(int));
            for (long int i = 0; i < totdim; i++) {
                double val = eigval->pointer()[i];
                skip[i] = 1;
                if ( val > pruning_threshold ) {
                    map[newdim] = i;
                    skip[i] = 0;
                    newdim++;
                }
            }
            //if ( newdim % 2 == 1 ) { // eigenvalues should come in pairs
            //    newdim--;
            //}
            for (int ij = 0; ij < totdim; ij++) {
                for (int kl = 0; kl < totdim; kl++) {
                    Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                }
            }
            Amat->transform(eigvec);
            Bmat->transform(eigvec);
        }else {
            newdim = totdim;
        }

        for (long int i = 0; i < newdim; i++) {
            for (long int j = 0; j < newdim; j++) {
                newA[i*newdim+j] = Amat->pointer()[i][j];
                newB[i*newdim+j] = Bmat->pointer()[i][j];
            }
        }

        if ( newdim < totdim ) {
            outfile->Printf("\n");
            outfile->Printf("    <<< warning >>> ");
            outfile->Printf("\n");
            outfile->Printf("    reducing dimension of [ k*l, [H, j*i ] ] from %5i to %5i\n",totdim,newdim);
            outfile->Printf("\n");
            eigval->print();
        }

        int info = 0;
        if ( prune ) {
            info = SymmetricGeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }else {
            GeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }

        // back transform wave functions to original untruncated basis
        std::shared_ptr<Matrix>  Cmat (new Matrix(totdim,totdim));
        for (int i = 0; i < newdim; i++) {
            for (int j = 0; j < totdim; j++) {
                double dum = 0.0;
                for (int k = 0; k < newdim; k++) {
                    dum += eigvec->pointer()[j][k] * newB[i * newdim + k];
                }
                Cmat->pointer()[i][j] = dum;
            }
        }

        // oscillator strengths: 
        // symmetry of state is h, consider only cartesian component with same symmetry
       
        if ( info != 0 ) {
            outfile->Printf("\n");
            outfile->Printf("    Error: diagonalization failed.  info = %5i\n",info);
            outfile->Printf("\n");
            throw PsiException("Error: diagonalization failed.",__FILE__,__LINE__);
        }

        outfile->Printf("    state");
        outfile->Printf("          energy (Eh)");
        outfile->Printf("       ex energy (Eh)");
        outfile->Printf("       ex energy (eV)");
        outfile->Printf("       f, osc. strength\n");

        for (long int state = newdim-1; state >= 0; state--){
          
            if ( eig[state] > 0.0 ) {

                if ( prune ) eig[state] = 1.0 / eig[state];

                // avoid considering nonsensical states
                if ( fabs(eig[state] ) > 10000.0 ) continue;

                //printf("hey excited state: %20.12lf\n",eig[state]+e2+e1+enuc_);

                // normalize C:
                double nrm = 0.0;

                double * C = Cmat->pointer()[state];
                for (int ij = 0; ij < nh; ij++) {
                    int i  = bas_erpa_sym[h][ij][0];
                    int j  = bas_erpa_sym[h][ij][1];
                    int hi = symmetry[i];
                    int hj = symmetry[j];
                    for( int kk = 0; kk < amopi_[hi]; kk++){
                        int k   = kk + pitzer_offset[hi];
                        int kj  = ibas_erpa_sym[h][k][j];

                        // check if kj is a proper ERPA excitations
                        if ( kj == -999 ) continue;

                        nrm    += 0.5 * C[ij] * C[kj] * D1a[k*amo_+i];
                        nrm    += 0.5 * C[ij] * C[kj] * D1b[k*amo_+i];
                    }
                    for( int ll = 0; ll < amopi_[hj]; ll++){
                        int l  = ll + pitzer_offset[hj];
                        int il = ibas_erpa_sym[h][i][l];

                        // check if il is a proper ERPA excitations
                        if ( il == -999 ) continue;

                        nrm   -= 0.5 * C[ij] * C[il] * D1a[j*amo_+l];
                        nrm   -= 0.5 * C[ij] * C[il] * D1b[j*amo_+l];
                    }
                }

                if ( nrm < 0.0 ) continue;

                nrm = sqrt(nrm);
                for (int ij = 0; ij < nh; ij++) {
                    C[ij] /= nrm;
                }


                outfile->Printf(" SPINFLIP   %5i %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138);

                bool * skip = (bool*)malloc(nh*sizeof(bool));
                memset((void*)skip,'\0',nh*sizeof(bool));

                for (int k = 0; k < number_coefficients; k++) {

                    double max = -999.0e99;
                    int maxij = 0;
                    for (int ij = 0; ij < nh; ij++) {
                        if ( skip[ij] ) continue;
                        if ( fabs(Cmat->pointer()[state][ij]) > max ) {
                            max = fabs(Cmat->pointer()[state][ij]);
                            maxij = ij;
                        }
                    }
                    skip[maxij] = true;
                    int id = maxij;
                    int i = bas_erpa_sym[h][id][0];
                    int j = bas_erpa_sym[h][id][1];
                    int hi = symmetry[i];
                    int hj = symmetry[j];
                    //outfile->Printf("                %s %3i, symm %3i -> %3i, symm %3i:  %20.6lf\n",(id == maxij) ? "alpha" : "beta ",i, hi, j, hj, Cmat->pointer()[state][maxij]);
                    if ( fabs(Cmat->pointer()[state][maxij]) > coefficient_threshold ) {
                        outfile->Printf("                %3i -> %3i  %20.6lf\n",i+1, j+1, Cmat->pointer()[state][maxij]);
                    }
                }

                free(skip);

                outfile->Printf("\n");
            } 
        }
   
        free(newA);
        free(newB);
        free(cc);
        free(eig);
    }

    free(tei);
    free(h1);
    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);

}*/

void ERPASolver::SpinAdaptedExtendedRPA() {

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals
    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // build two-electron integrals

    // note from AED (9/11/17): since the ERPA code below was designed for use 
    // with 3-index integrals, I had to play some tricks to make it work with
    // 4-index integrals.  If you want to use 4-index integrals (which is a bad
    // idea), the code will perform an SVD of the ERI tensor to yield factorized 
    // integrals that are appropriate for the 3-index code.  The factorization 
    // isn't symmetric, so there are now two buffers for 3-index integrals: Qmo_ 
    // and Qmo2.  In the case of DF, Qmo2 = Qmo_.  I say using the 4-index code
    // is a bad idea because it requires more storage.

    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));

    int info;
    double * Qmo2;
    if ( is_df_ ) {
        F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
        Qmo2 = Qmo_;
    }else { 
        // unpack two-electron integrals (which were initially blocked by symmetry)
        long int offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        // now, compute SVD of tei_full_sym_ and store result in Qmo_, Qmo2
        double * U  = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * VT = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * S  = (double*)malloc(nn1o2*sizeof(double));
        memset((void*)U,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)VT,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)S,'\0',nn1o2*sizeof(double));

        Qmo_ = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        memset((void*)Qmo_,'\0',nn1o2*nn1o2*sizeof(double));

        int lwork      = 5 * nn1o2;
        double * work  = (double*)malloc(lwork*sizeof(double));

        //info = C_DGESVD('A','A',nn1o2,nn1o2,tei,nn1o2,S,VT,nn1o2,U,nn1o2,work,lwork);
        char jobu  = 'A';
        char jobvt = 'A';
        int info;
        int dim = (int)nn1o2;
        info = F_DGESVD(&jobvt,&jobu,&dim,&dim,tei,&dim,S,U,&dim,VT,&dim,work,&lwork,&info);
        if ( info != 0 ) {
            throw PsiException("something went wrong in the SVD of the 4-index ERI tensor.",__FILE__,__LINE__);
        }

        // repack two-electron integrals (which were destroyed by SVD)
        memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));
        offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        double err = 0.0;
        for (int Q = 0; Q < nn1o2; Q++) {
            for (int kl = 0; kl < nn1o2; kl++) {
                Qmo_[kl*nn1o2 + Q] = U[Q*nn1o2 + kl] * S[Q];
            }
        }
        free(U);
        free(S);
        Qmo2 = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        C_DCOPY(nn1o2*nn1o2,VT,1,Qmo2,1);
        free(VT);

        //for (int ij = 0; ij < nn1o2; ij++) {
        //    for (int kl = 0; kl < nn1o2; kl++) {
        //        double dum = 0.0;
        //        for (int Q = 0; Q < nn1o2; Q++) {
        //            //dum += U[ij*nn1o2 + Q] * VT[kl*nn1o2 + Q] * S[Q];
        //            //dum += U[Q*nn1o2 + ij] * VT[Q*nn1o2 + kl] * S[Q];
        //            //dum += VT[ij*nn1o2 + Q] * U[kl*nn1o2 + Q] * S[Q];
        //            //dum += VT[Q*nn1o2 + ij] * U[Q*nn1o2 + kl] * S[Q];

        //            //dum += U[ij*nn1o2 + Q] * VT[Q*nn1o2 + kl] * S[Q];
        //            // this one:
        //            //dum += VT[ij*nn1o2 + Q] * U[Q*nn1o2 + kl] * S[Q];
        //            //dum += U[Q*nn1o2 + ij] * VT[Q*nn1o2 + kl] * S[Q];

        //            dum += VT[ij*nn1o2 + Q] * Qmo_[kl*nn1o2 + Q];
        //        }
        //        double diff = fabs(tei[ij*nn1o2+kl] - dum);
        //        if ( diff > 1e-6 ) {
        //            printf("%5i %5i %20.12lf %20.12lf\n",ij,kl,tei[ij*nn1o2+kl],dum);
        //        }
        //        err += diff*diff;
        //    }
        //}
        //err = sqrt(err);
        //printf("%20.12lf\n",err);
        nQ_ = nn1o2;
        F_DGEMM('t','n',nn1o2,nn1o2,nQ_,1.0,Qmo_,nQ_,Qmo2,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
    }

    // read TDPM from disk
    ReadTPDM(D2aa,D2bb,D2ab,D1a,D1b,tei);

    //ReadOPDM(D1a,D1b);

    long int n = nmo_*nmo_;

    // check energy
    double e2_aa = 0.0;
    double e2_ab = 0.0;
    double e2_bb = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2_ab +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2_aa += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2_bb += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];
                }
            }
        }
    }
    double e2 = e2_aa + e2_ab + e2_bb;

    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }

    //printf("%20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",e1,e2_ab,e2_aa,e2_bb,e1+e2,e1+e2+enuc_);
    //printf("%5i\n",info);
    //exit(0);

    // build intermediates for funky 3-index sums:
    double * tempa  = (double*)malloc(nmo_*nmo_*nQ_*sizeof(double));
    double * tempb  = (double*)malloc(nmo_*nmo_*nQ_*sizeof(double));
    double * Ia = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * Ib = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)tempa,'\0',nmo_*nmo_*nQ_*sizeof(double));
    memset((void*)tempb,'\0',nmo_*nmo_*nQ_*sizeof(double));
    memset((void*)Ia,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Ib,'\0',nmo_*nmo_*sizeof(double));

    for (int l = 0; l < amo_; l++) {

        int hl = symmetry[l];

        for (int p = 0; p < amo_; p++) {

            int hp = symmetry[p];

            int hpl = hp ^ hl;

            for (int Q = 0; Q < nQ_; Q++) {

                double duma = 0.0;
                double dumb = 0.0;

                for (int qs = 0; qs < gems_ab[hpl]; qs++) {

                    long int q = bas_ab_sym[hpl][qs][0];
                    long int s = bas_ab_sym[hpl][qs][1];

                    long int pq = p*nmo_+q;
                    long int ls = l*nmo_+s;

                    long int qp = q*nmo_+p;
                    long int sl = s*nmo_+l;

                    duma += Qmo_[INDEX(q,s)*nQ_+Q] * ( D2aa[pq*nmo_*nmo_+ls] + D2ab[pq*nmo_*nmo_+ls] );
                    dumb += Qmo_[INDEX(q,s)*nQ_+Q] * ( D2bb[pq*nmo_*nmo_+ls] + D2ab[qp*nmo_*nmo_+sl] );

                }

                tempa[l*nmo_*nQ_ + p*nQ_ + Q] = duma;
                tempb[l*nmo_*nQ_ + p*nQ_ + Q] = dumb;
            }
        }
    }
    for (int j = 0; j < amo_; j++) {
        for (int l = 0; l < amo_; l++) {
            double duma = 0.0;
            double dumb = 0.0;
            for (int p = 0; p < amo_; p++) {
                duma += C_DDOT(nQ_,tempa + l*nmo_*nQ_ + p*nQ_,1,Qmo2 + INDEX(p,j)*nQ_,1);
                dumb += C_DDOT(nQ_,tempb + l*nmo_*nQ_ + p*nQ_,1,Qmo2 + INDEX(p,j)*nQ_,1);
            }
            Ia[j*nmo_+l] = duma;
            Ib[j*nmo_+l] = dumb;
        }
    }
    free(tempa);
    free(tempb);

    // correlated (ERPA) ground-state 2-RDM
    double * D2aa_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
    double * D2bb_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
    double * D2ab_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));

    memset((void*)D2aa_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));
    memset((void*)D2bb_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));
    memset((void*)D2ab_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));

    // Gamma(pqrs) = gamma(pr)gamma(qs) + sum_v gamma_0v(pr)gamma_0v(qs) - gamma(qr)delta(ps)

    // D2aa(prqs) = gamma_a(pr)gamma_a(qs) + sum_v gamma_a_0v(pr)gamma_a_0v(qs) - gamma_a(qr)delta(ps)
    // D2bb(prqs) = gamma_b(pr)gamma_b(qs) + sum_v gamma_b_0v(pr)gamma_b_0v(qs) - gamma_b(qr)delta(ps)
    // D2ab(prqs) = gamma_a(pr)gamma_b(qs) + sum_v gamma_a_0v(pr)gamma_b_0v(qs)

    double ac_e2_1rdm = 0.0;
    for (int p = 0; p < nso_; p++) {
        for (int q = 0; q < nso_; q++) {
            for (int r = 0; r < nso_; r++) {
                for (int s = 0; s < nso_; s++) {

                    double np_a = D1a[p*nso_+p];
                    double nq_a = D1a[q*nso_+q];
                    double nr_a = D1a[r*nso_+r];
                    double ns_a = D1a[s*nso_+s];

                    double np_b = D1b[p*nso_+p];
                    double nq_b = D1b[q*nso_+q];
                    double nr_b = D1b[r*nso_+r];
                    double ns_b = D1b[s*nso_+s];

                    // Eq(12) ... which apparently has no contributions from the same-spin blocks of D2
                    //D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += (np_a - 1.0) * nq_b * (r==q)*(p==s);
                    // wait, there are ONLY contributions from the same-spin blocks of D2
                    D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += (np_a - 1.0) * nq_a * (r==q)*(p==s);
                    D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += (np_b - 1.0) * nq_b * (r==q)*(p==s);

                    //D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += np_a * nq_b * (p==r)*(q==s);
                    //D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += np_a * nq_a * (p==r)*(q==s) - np_a * (q==r);
                    //D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += np_b * nq_b * (p==r)*(q==s) - np_b * (q==r);

                    //double ga_pr = D1a[p*nso_+r];
                    //double gb_qs = D1b[q*nso_+s];
                    //D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_pr * gb_qs;// * 0.5;
                    ////if ( p != q && r != s ) {
                    //    double ga_qs = D1a[q*nso_+s];
                    //    double ga_qr = D1a[q*nso_+r];

                    //    double gb_pr = D1b[p*nso_+r];
                    //    double gb_qr = D1b[q*nso_+r];
                    //    D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_pr * ga_qs - ga_qr * (p==s);// * 0.5;
                    //    D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gb_pr * gb_qs - gb_qr * (p==s);// * 0.5;
                    ////}
                }
            }
        }
    }


    for (int h = 0; h < nirrep_; h++) {

        long int nh = gems_erpa[h];
        //long int nh = gems_ab[h];
        if ( nh == 0 ) continue;

        long int totdim = nh;
        //if ( type == "AA" ) totdim = 2*nh;

        double * newB = (double*)malloc(totdim*totdim*sizeof(double));
        double * newA = (double*)malloc(totdim*totdim*sizeof(double));
        double * cc   = (double*)malloc(totdim*totdim*sizeof(double));
        double * eig  = (double*)malloc(totdim*sizeof(double));

        memset((void*)newA,'\0',totdim*totdim*sizeof(double));
        memset((void*)newB,'\0',totdim*totdim*sizeof(double));
        memset((void*)cc,'\0',totdim*totdim*sizeof(double));
        memset((void*)eig,'\0',totdim*sizeof(double));

        outfile->Printf("\n");
        outfile->Printf("    Symmetry: %5i\n",h);
        outfile->Printf("\n");

        // build B(ij,kl) = <| [ k*l, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {
            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];
            for (long int kl = 0; kl < nh; kl++) {
                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;
                double dumabab = 0.0;
                double dumbaba = 0.0;

                if ( j == l ) {
                    duma    += D1a[k*nmo_+i];
                    dumb    += D1b[k*nmo_+i];

                    //dumabab += D1a[k*nmo_+i];
                    //dumbaba += D1b[k*nmo_+i];
                }
                if ( i == k ) {
                    duma    -= D1a[j*nmo_+l];
                    dumb    -= D1b[j*nmo_+l];

                    //dumabab -= D1b[j*nmo_+l];
                    //dumbaba -= D1a[j*nmo_+l];
                }

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);
                    
                long int id = ij * nh + kl;
                newB[id] = 0.5 * (duma + dumb);
                //newB[aaaa_sym] = duma;
                //newB[bbbb_sym] = dumb;
            }
        }

        // build A1(ij,kl) = <| [ k*l, [H1, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;
                //double dumabab = 0.0;
                //double dumbaba = 0.0;

                duma    += D1a[k*nmo_+i] * h1[l*nmo_+j];
                duma    += D1a[j*nmo_+l] * h1[i*nmo_+k];

                //dumabab += D1a[k*nmo_+i] * h1[l*nmo_+j];
                //dumabab += D1b[j*nmo_+l] * h1[i*nmo_+k];

                //dumbaba += D1b[k*nmo_+i] * h1[l*nmo_+j];
                //dumbaba += D1a[j*nmo_+l] * h1[i*nmo_+k];

                dumb    += D1b[k*nmo_+i] * h1[l*nmo_+j];
                dumb    += D1b[j*nmo_+l] * h1[i*nmo_+k];

                if ( j == l ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                            dumb    -= D1b[k*nmo_+q] * h1[i*nmo_+q];

                            //dumabab -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                            //dumbaba -= D1b[k*nmo_+q] * h1[i*nmo_+q];
                        }
                    }
                }
                if ( i == k ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                            dumb    -= D1b[q*nmo_+l] * h1[q*nmo_+j];

                            //dumabab -= D1b[q*nmo_+l] * h1[q*nmo_+j];
                            //dumbaba -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                        }
                    }
                }

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);
                    
                long int id = ij * nh + kl;
                newA[id] += 0.5 * (duma + dumb);
                //newA[aaaa_sym] += duma;
                //newA[bbbb_sym] += dumb;
            }
        }

#if 0
        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma = 0.0;
                double dumb = 0.0;

                double dumab = 0.0;
                double dumba = 0.0;

                double dumabab = 0.0;
                double dumbaba = 0.0;

                // 2 coulomb-like terms:

                // 1: (ik|qs) D(jq;ls)
/*
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                //double eint = TEI(i,k,q,s,SymmetryPair(hi,hk));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,k),1,Qmo_+nQ_*INDEX(q,s),1);
                                double eint = tei[INDEX(i,k)*nn1o2+INDEX(q,s)];

                                long int jq = j*nmo_+q;
                                long int ls = l*nmo_+s;

                                long int qj = q*nmo_+j;
                                long int sl = s*nmo_+l;

                                duma    += eint * ( D2aa[jq*n+ls] + D2ab[jq*n+ls] );
                                dumb    += eint * ( D2bb[jq*n+ls] + D2ab[qj*n+sl] );

                                dumabab += eint * ( D2bb[jq*n+ls] + D2ab[qj*n+sl] );
                                dumbaba += eint * ( D2aa[jq*n+ls] + D2ab[jq*n+ls] );
                            }
                        }
                    }
                }
                // 2: (lj|qs) D(kq;is)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                //double eint = TEI(l,j,q,s,SymmetryPair(hl,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(l,j),1,Qmo_+nQ_*INDEX(q,s),1);
                                double eint = tei[INDEX(l,j)*nn1o2+INDEX(q,s)];

                                long int kq = k*nmo_+q;
                                long int is = i*nmo_+s;

                                long int qk = q*nmo_+k;
                                long int si = s*nmo_+i;

                                duma += eint * ( D2aa[kq*n+is] + D2ab[kq*n+is] );
                                dumb += eint * ( D2bb[kq*n+is] + D2ab[qk*n+si] );

                                dumabab += eint * ( D2aa[kq*n+is] + D2ab[kq*n+is] );
                                dumbaba += eint * ( D2bb[kq*n+is] + D2ab[qk*n+si] );
                            }
                        }
                    }
                }
*/

                // funky sums over 3 indices:
                // - dik (qs|pj) D(pq;ls) + djl (qs|ip) D(kq;ps)
/*
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    for (int hs = 0; hs < nirrep_; hs++) {
                                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                            long int s = ss + pitzer_offset_full[hs];
                                            //double eint = TEI(q,s,p,j,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,s),1,Qmo_+nQ_*INDEX(p,j),1);
                                            double eint = tei[INDEX(q,s)*nn1o2+INDEX(p,j)];

                                            long int pq = p*nmo_+q;
                                            long int ls = l*nmo_+s;

                                            long int qp = q*nmo_+p;
                                            long int sl = s*nmo_+l;

                                            duma -= eint * ( D2aa[pq*n+ls] + D2ab[pq*n+ls] );
                                            dumb -= eint * ( D2bb[pq*n+ls] + D2ab[qp*n+sl] );

                                            dumabab -= eint * ( D2bb[pq*n+ls] + D2ab[qp*n+sl] );
                                            dumbaba -= eint * ( D2aa[pq*n+ls] + D2ab[pq*n+ls] );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
*/
/*
                if ( j == l ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    for (int hs = 0; hs < nirrep_; hs++) {
                                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                            long int s = ss + pitzer_offset_full[hs];
                                            //double eint = TEI(q,s,i,p,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,s),1,Qmo_+nQ_*INDEX(i,p),1);
                                            double eint = tei[INDEX(q,s)*nn1o2+INDEX(i,p)];

                                            long int kq = k*nmo_+q;
                                            long int ps = p*nmo_+s;

                                            long int qk = q*nmo_+k;
                                            long int sp = s*nmo_+p;

                                            duma -= eint * ( D2aa[kq*n+ps] + D2ab[kq*n+ps] );
                                            dumb -= eint * ( D2bb[kq*n+ps] + D2ab[qk*n+sp] );

                                            dumabab -= eint * ( D2aa[kq*n+ps] + D2ab[kq*n+ps] );
                                            dumbaba -= eint * ( D2bb[kq*n+ps] + D2ab[qk*n+sp] );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
*/

/*
                // exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                long int kq = k*nmo_+q;
                                long int is = i*nmo_+s;

                                duma -= eint * D2aa[kq*n+is];
                                dumb -= eint * D2bb[kq*n+is];

                                eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                long int jq = j*nmo_+q;
                                long int ls = l*nmo_+s;

                                duma -= eint * D2aa[jq*n+ls];
                                dumb -= eint * D2bb[jq*n+ls];
                            }
                        }
                    }
                }
*/

/*
                // exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];

                                double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                long int qk = q*nmo_+k;
                                long int is = i*nmo_+s;

                                long int kq = k*nmo_+q;
                                long int si = s*nmo_+i;

                                dumab += eint * D2ab[qk*n+is];
                                dumba += eint * D2ab[kq*n+si];

                                dumabab -= eint * D2ab[kq*n+is];
                                dumbaba -= eint * D2ab[qk*n+si];

                                eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                long int jq = j*nmo_+q;
                                long int sl = s*nmo_+l;

                                long int qj = q*nmo_+j;
                                long int ls = l*nmo_+s;

                                dumab += eint * D2ab[jq*n+sl];
                                dumba += eint * D2ab[qj*n+ls];

                                dumabab -= eint * D2ab[qj*n+sl];
                                dumbaba -= eint * D2ab[jq*n+ls];
                            }
                        }
                    }
                }
*/

/*
                // last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
                for (int hp = 0; hp < nirrep_; hp++) {
                    for (long int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        for (int hq = 0; hq < nirrep_; hq++) {
                            for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                long int q = qq + pitzer_offset_full[hq];

                                //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                long int pq = p*nmo_+q;
                                long int li = l*nmo_+i;

                                duma += eint * D2aa[pq*n+li];
                                dumb += eint * D2bb[pq*n+li];

                                //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                eint = tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                long int kj = k*nmo_+j;

                                duma += eint * D2aa[kj*n+pq];
                                dumb += eint * D2bb[kj*n+pq];
                            }
                        }
                    }
                }
*/

/*
                // last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)
                for (int hp = 0; hp < nirrep_; hp++) {
                    for (long int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        for (int hq = 0; hq < nirrep_; hq++) {
                            for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                long int q = qq + pitzer_offset_full[hq];

                                //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                long int pq = p*nmo_+q;
                                long int il = i*nmo_+l;

                                long int qp = q*nmo_+p;
                                long int li = l*nmo_+i;

                                dumab -= eint * D2ab[pq*n+il];
                                dumba -= eint * D2ab[qp*n+li];

                                dumabab += eint * D2ab[qp*n+il];
                                dumbaba += eint * D2ab[pq*n+li];

                                //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                eint = tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                long int kj = k*nmo_+j;
                                long int jk = j*nmo_+k;

                                dumab -= eint * D2ab[jk*n+pq];
                                dumba -= eint * D2ab[kj*n+qp];

                                dumabab += eint * D2ab[kj*n+pq];
                                dumbaba += eint * D2ab[jk*n+qp];
                            }
                        }
                    }
                }

*/

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                long int aabb_sym = (ij + 0 * nh) * 2 * nh + (kl + 1 * nh);
                long int bbaa_sym = (ij + 1 * nh) * 2 * nh + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                newA[aaaa_sym] += duma;
                newA[aabb_sym] += dumab;
                newA[bbaa_sym] += dumba;
                newA[bbbb_sym] += dumb;

            }
        }
#endif

        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>

        // N^5 version of A2 coulomb-like terms (1) and (2)
        // 
        // 1: (ik|qs) D(jq;ls)
        // 2: (lj|qs) D(kq;is)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            if ( gems_ab[hjl] == 0 ) continue;

            double * tempa = (double*)malloc(nQ_*gems_ab[hjl]*sizeof(double));
            double * tempb = (double*)malloc(nQ_*gems_ab[hjl]*sizeof(double));

            memset((void*)tempa,'\0',nQ_*gems_ab[hjl]*sizeof(double));
            memset((void*)tempb,'\0',nQ_*gems_ab[hjl]*sizeof(double));

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                for (int Q = 0; Q < nQ_; Q++) {

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {

                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        long int jq = j*nmo_+q;
                        long int ls = l*nmo_+s;
                        
                        long int qj = q*nmo_+j;
                        long int sl = s*nmo_+l;
                                    
                        duma += Qmo_[INDEX(q,s)*nQ_+Q] * (D2aa[jq * nmo_*nmo_ + ls] + D2ab[jq * nmo_*nmo_ + ls]);
                        dumb += Qmo_[INDEX(q,s)*nQ_+Q] * (D2bb[jq * nmo_*nmo_ + ls] + D2ab[qj * nmo_*nmo_ + sl]);

                    }

                    tempa[Q*gems_ab[hjl] + jl] = duma;
                    tempb[Q*gems_ab[hjl] + jl] = dumb;
                }
            }

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {
                
                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;
                
                    int hi = symmetry[i];
                
                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    // check if ij/kl are proper ERPA excitations
                    if ( ij == -999 || kl == -999 ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int Q = 0; Q < nQ_; Q++) {
                        duma += tempa[Q*gems_ab[hjl] + jl] * Qmo2[INDEX(i,k)*nQ_+Q];
                        dumb += tempb[Q*gems_ab[hjl] + jl] * Qmo2[INDEX(i,k)*nQ_+Q];
                    }

                    long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * (duma + dumb);

                    //newA[aaaa_sym] += duma;
                    //newA[bbbb_sym] += dumb;

                    // AED: I think term (2) above is just a funny transpose of term (1)
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aaaa_sym = (lk + 0 * nh) * 2 * nh + (ji + 0 * nh);
                    bbbb_sym = (lk + 1 * nh) * 2 * nh + (ji + 1 * nh);

                    id = ji * nh + lk;
                    newA[id] += 0.5 * (duma + dumb);

                    //newA[aaaa_sym] += duma;
                    //newA[bbbb_sym] += dumb;
                }
            }
            free(tempa);
            free(tempb);
        }

        // A2 exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int kq = k*nmo_+q;
                        long int is = i*nmo_+s;

                        duma -= eint * D2aa[kq*n+is];
                        dumb -= eint * D2bb[kq*n+is];

                        //eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                        //long int jq = j*nmo_+q;
                        //long int ls = l*nmo_+s;

                        //duma -= eint * D2aa[jq*n+ls];
                        //dumb -= eint * D2bb[jq*n+ls];
                    }

                    long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * (duma + dumb);


                    //newA[aaaa_sym] += duma;
                    //newA[bbbb_sym] += dumb;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;


                    aaaa_sym = (lk + 0 * nh) * 2 * nh + (ji + 0 * nh);
                    bbbb_sym = (lk + 1 * nh) * 2 * nh + (ji + 1 * nh);

                    newA[id] += 0.5 * (duma + dumb);
                    //newA[aaaa_sym] += duma;
                    //newA[bbbb_sym] += dumb;

                }
            }
        }

        // A2 exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int qk = q*nmo_+k;
                        long int is = i*nmo_+s;

                        long int kq = k*nmo_+q;
                        long int si = s*nmo_+i;

                        dumab += eint * D2ab[qk*n+is];
                        dumba += eint * D2ab[kq*n+si];

                        //dumabab -= eint * D2ab[kq*n+is];
                        //dumbaba -= eint * D2ab[qk*n+si];

                    }

                    long int aabb_sym = (ij + 0 * nh) * 2 * nh + (kl + 1 * nh);
                    long int bbaa_sym = (ij + 1 * nh) * 2 * nh + (kl + 0 * nh);

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * (dumab + dumba);

                    //newA[aabb_sym] += dumab;
                    //newA[bbaa_sym] += dumba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    aabb_sym = (ji + 0 * nh) * 2 * nh + (lk + 1 * nh);
                    bbaa_sym = (ji + 1 * nh) * 2 * nh + (lk + 0 * nh);

                    newA[id] += 0.5 * (dumab + dumba);

                    //newA[aabb_sym] += dumab;
                    //newA[bbaa_sym] += dumba;

                }
            }
        }

        // N^5 version of A2 funky sums over 3 indices:
        // 
        // 1: -dik (qs|pj) D(pq;ls)
        // 2: -djl (qs|ip) D(kq;ps)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 
        // Also, the intermediates for this guy can be built
        // at N^4 cost and can be done outside of these loops (see above)

        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            int hj = symmetry[j];
            for (int ll = 0; ll < amopi_[hj]; ll++) {

                int l = ll + pitzer_offset[hj];

                if ( i == l ) continue;

                long int il = ibas_erpa_sym[h][i][l];

                // check if il is a proper ERPA excitation

                if ( il == -999 ) continue;

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (il + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (il + 1 * nh);

                long int id = ij * nh + il;

                newA[id] -= 0.5 * Ia[j*nmo_+l];
                newA[id] -= 0.5 * Ib[j*nmo_+l];

                //newA[aaaa_sym] -= Ia[j*nmo_+l];
                //newA[bbbb_sym] -= Ib[j*nmo_+l];

                // AED: I think second funky term is a transpose of the first
                // A(ij,il)(2) = A(ji,li)(1)

                // if ij/il are proper ERPA excitations, then ji/li should be as well

                long int ji = ibas_erpa_sym[h][j][i];
                long int li = ibas_erpa_sym[h][l][i];

                aaaa_sym = (ji + 0 * nh) * 2 * nh + (li + 0 * nh);
                bbbb_sym = (ji + 1 * nh) * 2 * nh + (li + 1 * nh);

                id = ji * nh + li;

                newA[id] -= 0.5 * Ia[j*nmo_+l];
                newA[id] -= 0.5 * Ib[j*nmo_+l];

                //newA[aaaa_sym] -= Ia[j*nmo_+l];
                //newA[bbbb_sym] -= Ib[j*nmo_+l];

            }
        }

        // A2: last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        duma += eint * D2aa[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];
                        dumb += eint * D2bb[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];

                    }

                    long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                    long int id = ij * nh + kl;

                    newA[id] += 0.5 * (duma + dumb);

                    //newA[aaaa_sym] += duma;
                    //newA[bbbb_sym] += dumb;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitaitons, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    id = ji * nh + lk;

                    aaaa_sym = (ji + 0 * nh) * 2 * nh + (lk + 0 * nh);
                    bbbb_sym = (ji + 1 * nh) * 2 * nh + (lk + 1 * nh);

                    newA[id] += 0.5 * (duma + dumb);
                    //newA[aaaa_sym] += duma;
                    //newA[bbbb_sym] += dumb;

                }

            }

        }

        // A2: last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)

        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        dumab -= eint * D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+i*nmo_+l];
                        dumba -= eint * D2ab[q*nmo_*nmo_*nmo_+p*nmo_*nmo_+l*nmo_+i];

                    }

                    long int aabb_sym = (ij + 0 * nh) * 2 * nh + (kl + 1 * nh);
                    long int bbaa_sym = (ij + 1 * nh) * 2 * nh + (kl + 0 * nh);

                    long int id = ij * nh + kl;
                    newA[id] += 0.5 * (dumab + dumba);

                    //newA[aabb_sym] += dumab;
                    //newA[bbaa_sym] += dumba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aabb_sym = (ji + 0 * nh) * 2 * nh + (lk + 1 * nh);
                    bbaa_sym = (ji + 1 * nh) * 2 * nh + (lk + 0 * nh);

                    id = ji * nh + lk;

                    newA[id] += 0.5 * (dumab + dumba);

                    //newA[aabb_sym] += dumab;
                    //newA[bbaa_sym] += dumba;

                }

            }

        }

        std::shared_ptr<Matrix> Amat ( new Matrix(totdim,totdim) );
        std::shared_ptr<Matrix> Bmat ( new Matrix(totdim,totdim) );
        for (int ij = 0; ij < totdim; ij++) {
            for (int kl = 0; kl < totdim; kl++) {
                Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                Bmat->pointer()[ij][kl] = newB[ij*totdim+kl];
            }
        }

        // symmetrize
        // A(kl,ji) = A(ij,kl)
        // A(kl,ji) = A(lk,ij)

        // project out null (and negative) space of A matrix
        std::shared_ptr<Matrix> eigvec ( new Matrix(totdim,totdim) );
        std::shared_ptr<Vector> eigval ( new Vector(totdim) );

        Amat->diagonalize(eigvec,eigval,descending);
        //eigval->print();

        bool prune = true;
        long int newdim = 0;

        // get rid of small/negative eigenvalues of A?
        if ( prune ) {

            newdim = 0;
            int * map = (int*)malloc(totdim*sizeof(int));
            memset((void*)map,'\0',totdim*sizeof(int));
            int * skip = (int*)malloc(totdim*sizeof(int));
            memset((void*)skip,'\0',totdim*sizeof(int));
            for (long int i = 0; i < totdim; i++) {
                double val = eigval->pointer()[i];
                skip[i] = 1;
                if ( val > pruning_threshold ) {
                    map[newdim] = i;
                    skip[i] = 0;
                    newdim++;
                }
            }
            //if ( newdim % 2 == 1 ) { // eigenvalues should come in pairs
            //    newdim--;
            //}
            for (int ij = 0; ij < totdim; ij++) {
                for (int kl = 0; kl < totdim; kl++) {
                    Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                }
            }
            Amat->transform(eigvec);
            Bmat->transform(eigvec);
        }else {
            newdim = totdim;
        }

        for (long int i = 0; i < newdim; i++) {
            for (long int j = 0; j < newdim; j++) {
                newA[i*newdim+j] = Amat->pointer()[i][j];
                newB[i*newdim+j] = Bmat->pointer()[i][j];
            }
        }

        outfile->Printf("\n");
        if ( newdim < totdim ) {
            outfile->Printf("    <<< warning >>> ");
            outfile->Printf("\n");
            outfile->Printf("    reducing dimension of [ k*l, [H, j*i ] ] from %5i to %5i\n",totdim,newdim);
            outfile->Printf("\n");
        }

        int info = 0;
        if ( prune ) {
            info = SymmetricGeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }else {
	    bool * energy_is_real = (bool*)malloc(totdim*sizeof(bool));
            GeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig,energy_is_real);
            free(energy_is_real);
        }

        // back transform wave functions to original untruncated basis
        std::shared_ptr<Matrix>  Cmat (new Matrix(totdim,totdim));
        for (int i = 0; i < newdim; i++) {
            for (int j = 0; j < totdim; j++) {
                double dum = 0.0;
                for (int k = 0; k < newdim; k++) {
                    dum += eigvec->pointer()[j][k] * newB[i * newdim + k];
                }
                Cmat->pointer()[i][j] = dum;
            }
        }

        // oscillator strengths: 
        // symmetry of state is h, consider only cartesian component with same symmetry
       
        if ( info != 0 ) {
            outfile->Printf("\n");
            outfile->Printf("    Error: diagonalization failed.  info = %5i\n",info);
            outfile->Printf("\n");
        }

        outfile->Printf("    state");
        outfile->Printf("          energy (Eh)");
        outfile->Printf("       ex energy (Eh)");
        outfile->Printf("       ex energy (eV)");
        outfile->Printf("       f, osc. strength\n");

        for (long int state = newdim-1; state >= 0; state--){
          
            if ( eig[state] > 0.0 ) {

                if ( prune ) eig[state] = 1.0 / eig[state];

                // avoid considering nonsensical states
                if ( fabs(eig[state] ) > 10000.0 ) continue;

                // normalize C:
                double nrm = 0.0;

                for (int ij = 0; ij < nh; ij++) {
                    int i  = bas_erpa_sym[h][ij][0];
                    int j  = bas_erpa_sym[h][ij][1];
                    int hi = symmetry[i];
                    int hj = symmetry[j];
                    for( int kk = 0; kk < amopi_[hi]; kk++){
                        int k   = kk + pitzer_offset[hi];
                        int kj  = ibas_erpa_sym[h][k][j];

                        // check if kj is a proper ERPA excitations
                        if ( kj == -999 ) continue;

                        nrm    += 0.5 * Cmat->pointer()[state][ij] * Cmat->pointer()[state][kj] * D1a[k*amo_+i];
                        nrm    += 0.5 * Cmat->pointer()[state][ij] * Cmat->pointer()[state][kj] * D1b[k*amo_+i];
                    }
                    for( int ll = 0; ll < amopi_[hj]; ll++){
                        int l  = ll + pitzer_offset[hj];
                        int il = ibas_erpa_sym[h][i][l];

                        // check if il is a proper ERPA excitations
                        if ( il == -999 ) continue;

                        nrm   -= 0.5 * Cmat->pointer()[state][ij] * Cmat->pointer()[state][il] * D1a[j*amo_+l];
                        nrm   -= 0.5 * Cmat->pointer()[state][ij] * Cmat->pointer()[state][il] * D1b[j*amo_+l];
                    }
                }


                if ( nrm < 0.0 ) continue;

                nrm = sqrt(nrm);
                for (int ij = 0; ij < nh; ij++) {
                    Cmat->pointer()[state][ij] /= nrm;
                }
                
                double dumx = 0.0;
                double dumy = 0.0;
                double dumz = 0.0;

                if ( h == hx_ ) {
                   
                    dumx = TransitionDipoleMoment_SpinAdapted(Cmat->pointer()[state], D1a, D1b, nh,h,"X");

                }
                if ( h == hy_ ) {
    
                    dumy = TransitionDipoleMoment_SpinAdapted(Cmat->pointer()[state], D1a, D1b, nh,h,"Y");

                }
                if ( h == hz_ ) {
    
                    dumz = TransitionDipoleMoment_SpinAdapted(Cmat->pointer()[state], D1a, D1b, nh,h,"Z");

                }

                double * TDa = (double*)malloc(nso_*nso_*sizeof(double));
                double * TDb = (double*)malloc(nso_*nso_*sizeof(double));

                memset((void*)TDa,'\0',nso_*nso_*sizeof(double));
                memset((void*)TDb,'\0',nso_*nso_*sizeof(double));

                // build transition density matrix:
                TransitionDensity_SpinAdapted(Cmat->pointer()[state], D1a, D1b, TDa,TDb,nh,h);

                // use transition density matrix to update ground-state 2-RDM
                for (int p = 0; p < nso_; p++) {
                    for (int q = 0; q < nso_; q++) {
                        for (int r = 0; r < nso_; r++) {
                            for (int s = 0; s < nso_; s++) {
                                double ga_rp = TDa[p*nso_+r]; // careful with transpose ...
                                double gb_qs = TDb[s*nso_+q];
                                // Eq(10)
                                D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_rp * gb_qs;// * 0.5;
                                //if ( p != q && r != s) {

                                    double ga_qs = TDa[s*nso_+q];
                                    D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_rp * ga_qs;// * 0.5;

                                    double gb_rp = TDb[p*nso_+r]; // careful with transpose ...
                                    D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gb_rp * gb_qs;// * 0.5;
                                //}

                                // Eq(17) ... only valid in the NO basis


                                // Eq(17) ... only valid in the NO basis

                                //double np = D1a[p*nso_+p];
                                //double nq = D1a[q*nso_+q];
                                //double nr = D1a[r*nso_+r];
                                //double nr = D1a[r*nso_+r];

                                //double ga_

                                //D2_17[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += 2.0 * (
                                //    - 0.5 * ( np*(1.0-nq) + nq*(1.0-np) ) * (p==r)*(q==s);
                            }
                        }
                    }
                }
                free(TDa);
                free(TDb);


                //outfile->Printf(" NEWWAY   %5i %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138);

                double val = 2./3. * eig[state] * (dumx*dumx+dumy*dumy+dumz*dumz);
                
                outfile->Printf(" NEWWAY   %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138, val);
                outfile->Printf("                Transition Moment: X: %10.6lf    Y: %10.6lf   Z: %10.6lf\n", dumx, dumy, dumz);
                
                bool * skip = (bool*)malloc(nh*sizeof(bool));
                memset((void*)skip,'\0',nh*sizeof(bool));
                
                for (int k = 0; k < number_coefficients; k++) {

                    double max = -999.0e99;
                    int maxij = 0;
                    for (int ij = 0; ij < nh; ij++) {
                        if ( skip[ij] ) continue;
                        if ( fabs(Cmat->pointer()[state][ij]) > max ) {
                            max = fabs(Cmat->pointer()[state][ij]);
                            maxij = ij;
                        }
                    }
                    skip[maxij] = true;
                    int id = maxij;
                    int i = bas_erpa_sym[h][id][0];
                    int j = bas_erpa_sym[h][id][1];
                    int hi = symmetry[i];
                    int hj = symmetry[j];
                    //outfile->Printf("                %s %3i, symm %3i -> %3i, symm %3i:  %20.6lf\n",(id == maxij) ? "alpha" : "beta ",i, hi, j, hj, Cmat->pointer()[state][maxij]);
                    if ( fabs(Cmat->pointer()[state][maxij]) > coefficient_threshold ) {
                        outfile->Printf("                %3i -> %3i  %20.6lf\n",i+1, j+1, Cmat->pointer()[state][maxij]);
                    }
                }

                free(skip);
                
                outfile->Printf("\n");                   
            } 
        }
   
        free(newA);
        free(newB);
        free(cc);
        free(eig);
    }

    // check energy
    double rpa_e2 = 0.0;
    double ac_e2 = 0.0;

    // recompute zeroth-order energy
    e2 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];

                    rpa_e2 +=       eint * D2ab_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                    rpa_e2 += 0.5 * eint * D2aa_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                    rpa_e2 += 0.5 * eint * D2bb_prime[(i*nmo_+j)*n+(k*nmo_+l)];

                    int nactive = 0;
                    if ( active_reference_orbitals_[i] ) nactive++;
                    if ( active_reference_orbitals_[j] ) nactive++;
                    if ( active_reference_orbitals_[k] ) nactive++;
                    if ( active_reference_orbitals_[l] ) nactive++;

                    int ninactive = 0;
                    if ( inactive_reference_orbitals_[i] ) ninactive++;
                    if ( inactive_reference_orbitals_[j] ) ninactive++;
                    if ( inactive_reference_orbitals_[k] ) ninactive++;
                    if ( inactive_reference_orbitals_[l] ) ninactive++;

                    if ( nactive != 4 && ninactive != 4 && ( nactive + ninactive ) != 0  ){

                        ac_e2 +=       eint * D2ab_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                        ac_e2 += 0.5 * eint * D2aa_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                        ac_e2 += 0.5 * eint * D2bb_prime[(i*nmo_+j)*n+(k*nmo_+l)];

                    }
                    e2 +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];
                }
            }
        }
    }
    e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }

    // check traces:
    double traa = 0.0;
    double trbb = 0.0;
    double trab = 0.0;
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = 0; j < nmo_; j++) {
            long int ij = i*nmo_+j;
            traa += D2ab_prime[(i*nmo_+j)*n+(i*nmo_+j)];// * 0.5;
            trbb += D2aa_prime[(i*nmo_+j)*n+(i*nmo_+j)];// * 0.5;
            trab += D2bb_prime[(i*nmo_+j)*n+(i*nmo_+j)];// * 0.5;
        }
    }

    //outfile->Printf("\n");                   
    //outfile->Printf("    ==> ERPA correlation information: <==\n");
    //outfile->Printf("\n");                   
    //outfile->Printf("        Tr(D2aa):            %20.12lf\n",traa);
    //outfile->Printf("        Tr(D2bb):            %20.12lf\n",trbb);
    //outfile->Printf("        Tr(D2ab):            %20.12lf\n",trab);
    //outfile->Printf("\n");                   
    //outfile->Printf("        * ERPA total energy: %20.12lf\n",e1+rpa_e2+enuc_);
    outfile->Printf("\n");                   
    outfile->Printf("    * AC1 total energy:  %20.12lf\n",e1+e2 + 0.5 * ac_e2 +enuc_);
    outfile->Printf("\n");                   

    free(Ia);
    free(Ib);
    free(tei);
    free(h1);
    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);
    free(D2aa_prime);
    free(D2bb_prime);
    free(D2ab_prime);
    if ( !is_df_ ) {
        free(Qmo2);
    }
}

void ERPASolver::NewExtendedRPA(std::string type) {

    if ( type != "AA" ) {
        throw PsiException("erpa currently only works for spin-conserving excitations",__FILE__,__LINE__);
    }

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals
    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // build two-electron integrals

    // note from AED (9/11/17): since the ERPA code below was designed for use 
    // with 3-index integrals, I had to play some tricks to make it work with
    // 4-index integrals.  If you want to use 4-index integrals (which is a bad
    // idea), the code will perform an SVD of the ERI tensor to yield factorized 
    // integrals that are appropriate for the 3-index code.  The factorization 
    // isn't symmetric, so there are now two buffers for 3-index integrals: Qmo_ 
    // and Qmo2.  In the case of DF, Qmo2 = Qmo_.  I say using the 4-index code
    // is a bad idea because it requires more storage.

    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));

    int info;
    double * Qmo2;
    if ( is_df_ ) {
        F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
        Qmo2 = Qmo_;
    }else { 
        // unpack two-electron integrals (which were initially blocked by symmetry)
        long int offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        // now, compute SVD of tei_full_sym_ and store result in Qmo_, Qmo2
        double * U  = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * VT = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        double * S  = (double*)malloc(nn1o2*sizeof(double));
        memset((void*)U,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)VT,'\0',nn1o2*nn1o2*sizeof(double));
        memset((void*)S,'\0',nn1o2*sizeof(double));

        Qmo_ = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        memset((void*)Qmo_,'\0',nn1o2*nn1o2*sizeof(double));

        int lwork      = 5 * nn1o2;
        double * work  = (double*)malloc(lwork*sizeof(double));

        //info = C_DGESVD('A','A',nn1o2,nn1o2,tei,nn1o2,S,VT,nn1o2,U,nn1o2,work,lwork);
        char jobu  = 'A';
        char jobvt = 'A';
        int info;
        int dim = (int)nn1o2;
        info = F_DGESVD(&jobvt,&jobu,&dim,&dim,tei,&dim,S,U,&dim,VT,&dim,work,&lwork,&info);
        if ( info != 0 ) {
            throw PsiException("something went wrong in the SVD of the 4-index ERI tensor.",__FILE__,__LINE__);
        }

        // repack two-electron integrals (which were destroyed by SVD)
        memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));
        offset = 0;
        for (int h = 0; h < nirrep_; h++) {
            for (long int ij = 0; ij < gems_full[h]; ij++) {
                long int i = bas_really_full_sym[h][ij][0];
                long int j = bas_really_full_sym[h][ij][1];
                for (long int kl = 0; kl < gems_full[h]; kl++) {
                    long int k = bas_really_full_sym[h][kl][0];
                    long int l = bas_really_full_sym[h][kl][1];
                    tei[INDEX(i,j)*nn1o2+INDEX(k,l)] = tei_full_sym_[offset + INDEX(ij,kl)];
                }
            }
            offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
        }

        double err = 0.0;
        for (int Q = 0; Q < nn1o2; Q++) {
            for (int kl = 0; kl < nn1o2; kl++) {
                Qmo_[kl*nn1o2 + Q] = U[Q*nn1o2 + kl] * S[Q];
            }
        }
        free(U);
        free(S);
        Qmo2 = (double*)malloc(nn1o2*nn1o2*sizeof(double));
        C_DCOPY(nn1o2*nn1o2,VT,1,Qmo2,1);
        free(VT);

        //for (int ij = 0; ij < nn1o2; ij++) {
        //    for (int kl = 0; kl < nn1o2; kl++) {
        //        double dum = 0.0;
        //        for (int Q = 0; Q < nn1o2; Q++) {
        //            //dum += U[ij*nn1o2 + Q] * VT[kl*nn1o2 + Q] * S[Q];
        //            //dum += U[Q*nn1o2 + ij] * VT[Q*nn1o2 + kl] * S[Q];
        //            //dum += VT[ij*nn1o2 + Q] * U[kl*nn1o2 + Q] * S[Q];
        //            //dum += VT[Q*nn1o2 + ij] * U[Q*nn1o2 + kl] * S[Q];

        //            //dum += U[ij*nn1o2 + Q] * VT[Q*nn1o2 + kl] * S[Q];
        //            // this one:
        //            //dum += VT[ij*nn1o2 + Q] * U[Q*nn1o2 + kl] * S[Q];
        //            //dum += U[Q*nn1o2 + ij] * VT[Q*nn1o2 + kl] * S[Q];

        //            dum += VT[ij*nn1o2 + Q] * Qmo_[kl*nn1o2 + Q];
        //        }
        //        double diff = fabs(tei[ij*nn1o2+kl] - dum);
        //        if ( diff > 1e-6 ) {
        //            printf("%5i %5i %20.12lf %20.12lf\n",ij,kl,tei[ij*nn1o2+kl],dum);
        //        }
        //        err += diff*diff;
        //    }
        //}
        //err = sqrt(err);
        //printf("%20.12lf\n",err);
        nQ_ = nn1o2;
        F_DGEMM('t','n',nn1o2,nn1o2,nQ_,1.0,Qmo_,nQ_,Qmo2,nQ_,0.0,tei,nmo_*(nmo_+1)/2);
    }

    // read TDPM from disk
    ReadTPDM(D2aa,D2bb,D2ab,D1a,D1b,tei);

    //ReadOPDM(D1a,D1b);

    long int n = nmo_*nmo_;

    // check energy
    double e2_aa = 0.0;
    double e2_ab = 0.0;
    double e2_bb = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2_ab +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2_aa += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2_bb += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];
                }
            }
        }
    }
    double e2 = e2_aa + e2_ab + e2_bb;

    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }

    //printf("%20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",e1,e2_ab,e2_aa,e2_bb,e1+e2,e1+e2+enuc_);
    //printf("%5i\n",info);
    //exit(0);

    // build intermediates for funky 3-index sums:
    double * tempa  = (double*)malloc(nmo_*nmo_*nQ_*sizeof(double));
    double * tempb  = (double*)malloc(nmo_*nmo_*nQ_*sizeof(double));
    double * Ia = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * Ib = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)tempa,'\0',nmo_*nmo_*nQ_*sizeof(double));
    memset((void*)tempb,'\0',nmo_*nmo_*nQ_*sizeof(double));
    memset((void*)Ia,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)Ib,'\0',nmo_*nmo_*sizeof(double));

    for (int l = 0; l < amo_; l++) {

        int hl = symmetry[l];

        for (int p = 0; p < amo_; p++) {

            int hp = symmetry[p];

            int hpl = hp ^ hl;

            for (int Q = 0; Q < nQ_; Q++) {

                double duma = 0.0;
                double dumb = 0.0;

                for (int qs = 0; qs < gems_ab[hpl]; qs++) {

                    long int q = bas_ab_sym[hpl][qs][0];
                    long int s = bas_ab_sym[hpl][qs][1];

                    long int pq = p*nmo_+q;
                    long int ls = l*nmo_+s;

                    long int qp = q*nmo_+p;
                    long int sl = s*nmo_+l;

                    duma += Qmo_[INDEX(q,s)*nQ_+Q] * ( D2aa[pq*nmo_*nmo_+ls] + D2ab[pq*nmo_*nmo_+ls] );
                    dumb += Qmo_[INDEX(q,s)*nQ_+Q] * ( D2bb[pq*nmo_*nmo_+ls] + D2ab[qp*nmo_*nmo_+sl] );

                }

                tempa[l*nmo_*nQ_ + p*nQ_ + Q] = duma;
                tempb[l*nmo_*nQ_ + p*nQ_ + Q] = dumb;
            }
        }
    }
    for (int j = 0; j < amo_; j++) {
        for (int l = 0; l < amo_; l++) {
            double duma = 0.0;
            double dumb = 0.0;
            for (int p = 0; p < amo_; p++) {
                duma += C_DDOT(nQ_,tempa + l*nmo_*nQ_ + p*nQ_,1,Qmo2 + INDEX(p,j)*nQ_,1);
                dumb += C_DDOT(nQ_,tempb + l*nmo_*nQ_ + p*nQ_,1,Qmo2 + INDEX(p,j)*nQ_,1);
            }
            Ia[j*nmo_+l] = duma;
            Ib[j*nmo_+l] = dumb;
        }
    }
    free(tempa);
    free(tempb);

    // correlated (ERPA) ground-state 2-RDM
    double * D2aa_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
    double * D2bb_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));
    double * D2ab_prime = (double*)malloc(nso_*nso_*nso_*nso_*sizeof(double));

    memset((void*)D2aa_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));
    memset((void*)D2bb_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));
    memset((void*)D2ab_prime,'\0',nso_*nso_*nso_*nso_*sizeof(double));

    // Gamma(pqrs) = gamma(pr)gamma(qs) + sum_v gamma_0v(pr)gamma_0v(qs) - gamma(qr)delta(ps)

    // D2aa(prqs) = gamma_a(pr)gamma_a(qs) + sum_v gamma_a_0v(pr)gamma_a_0v(qs) - gamma_a(qr)delta(ps)
    // D2bb(prqs) = gamma_b(pr)gamma_b(qs) + sum_v gamma_b_0v(pr)gamma_b_0v(qs) - gamma_b(qr)delta(ps)
    // D2ab(prqs) = gamma_a(pr)gamma_b(qs) + sum_v gamma_a_0v(pr)gamma_b_0v(qs)

    for (int p = 0; p < nso_; p++) {
        for (int q = 0; q < nso_; q++) {
            for (int r = 0; r < nso_; r++) {
                for (int s = 0; s < nso_; s++) {
                    double ga_pr = D1a[p*nso_+r];
                    double gb_qs = D1b[q*nso_+s];

                    double np_a = D1a[p*nso_+p];
                    double nq_a = D1a[q*nso_+q];
                    double nr_a = D1a[r*nso_+r];
                    double ns_a = D1a[s*nso_+s];

                    double np_b = D1b[p*nso_+p];
                    double nq_b = D1b[q*nso_+q];
                    double nr_b = D1b[r*nso_+r];
                    double ns_b = D1b[s*nso_+s];

                    // Eq(12) ... which apparently has no contributions from the same-spin blocks of D2
                    //D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += (np_a - 1.0) * nq_b * (r==q)*(p==s);
                    // wait, there are ONLY contributions from the same-spin blocks of D2
                    D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += (np_a - 1.0) * nq_a * (r==q)*(p==s);
                    D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += (np_b - 1.0) * nq_b * (r==q)*(p==s);

                    ////if ( p != q && r != s ) {
                    //    double ga_qs = D1a[q*nso_+s];
                    //    double ga_qr = D1a[q*nso_+r];

                    //    double gb_pr = D1b[p*nso_+r];
                    //    double gb_qr = D1b[q*nso_+r];
                    //    D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_pr * ga_qs - ga_qr * (p==s);// * 0.5;
                    //    D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gb_pr * gb_qs - gb_qr * (p==s);// * 0.5;
                    ////}
                }
            }
        }
    }


    for (int h = 0; h < nirrep_; h++) {

        long int nh = gems_erpa[h];
        //long int nh = gems_ab[h];
        if ( nh == 0 ) continue;

        long int totdim = nh;
        if ( type == "AA" ) totdim = 2*nh;

        double * newB = (double*)malloc(totdim*totdim*sizeof(double));
        double * newA = (double*)malloc(totdim*totdim*sizeof(double));
        double * cc   = (double*)malloc(totdim*totdim*sizeof(double));
        double * eig  = (double*)malloc(totdim*sizeof(double));

        memset((void*)newA,'\0',totdim*totdim*sizeof(double));
        memset((void*)newB,'\0',totdim*totdim*sizeof(double));
        memset((void*)cc,'\0',totdim*totdim*sizeof(double));
        memset((void*)eig,'\0',totdim*sizeof(double));

        outfile->Printf("\n");
        outfile->Printf("    Symmetry: %5i\n",h);
        outfile->Printf("\n");

        // build B(ij,kl) = <| [ k*l, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {
            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];
            for (long int kl = 0; kl < nh; kl++) {
                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;
                double dumabab = 0.0;
                double dumbaba = 0.0;

                if ( j == l ) {
                    duma    += D1a[k*nmo_+i];
                    dumb    += D1b[k*nmo_+i];

                    //dumabab += D1a[k*nmo_+i];
                    //dumbaba += D1b[k*nmo_+i];
                }
                if ( i == k ) {
                    duma    -= D1a[j*nmo_+l];
                    dumb    -= D1b[j*nmo_+l];

                    //dumabab -= D1b[j*nmo_+l];
                    //dumbaba -= D1a[j*nmo_+l];
                }

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);
                    
                newB[aaaa_sym] = duma;
                newB[bbbb_sym] = dumb;
            }
        }

        // build A1(ij,kl) = <| [ k*l, [H1, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma    = 0.0;
                double dumb    = 0.0;
                //double dumabab = 0.0;
                //double dumbaba = 0.0;

                duma    += D1a[k*nmo_+i] * h1[l*nmo_+j];
                duma    += D1a[j*nmo_+l] * h1[i*nmo_+k];

                //dumabab += D1a[k*nmo_+i] * h1[l*nmo_+j];
                //dumabab += D1b[j*nmo_+l] * h1[i*nmo_+k];

                //dumbaba += D1b[k*nmo_+i] * h1[l*nmo_+j];
                //dumbaba += D1a[j*nmo_+l] * h1[i*nmo_+k];

                dumb    += D1b[k*nmo_+i] * h1[l*nmo_+j];
                dumb    += D1b[j*nmo_+l] * h1[i*nmo_+k];

                if ( j == l ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                            dumb    -= D1b[k*nmo_+q] * h1[i*nmo_+q];

                            //dumabab -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                            //dumbaba -= D1b[k*nmo_+q] * h1[i*nmo_+q];
                        }
                    }
                }
                if ( i == k ) {
                    for (int hq = 0; hq < nirrep_; hq++) {
                        for (int qq = 0; qq < amopi_[hq]; qq++) {
                            long int q = qq + pitzer_offset_full[hq];
                            duma    -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                            dumb    -= D1b[q*nmo_+l] * h1[q*nmo_+j];

                            //dumabab -= D1b[q*nmo_+l] * h1[q*nmo_+j];
                            //dumbaba -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                        }
                    }
                }

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);
                    
                newA[aaaa_sym] += duma;
                newA[bbbb_sym] += dumb;
            }
        }

#if 0
        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>
        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            for (long int kl = 0; kl < nh; kl++) {

                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];

                double duma = 0.0;
                double dumb = 0.0;

                double dumab = 0.0;
                double dumba = 0.0;

                double dumabab = 0.0;
                double dumbaba = 0.0;

                // 2 coulomb-like terms:

                // 1: (ik|qs) D(jq;ls)
/*
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                //double eint = TEI(i,k,q,s,SymmetryPair(hi,hk));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,k),1,Qmo_+nQ_*INDEX(q,s),1);
                                double eint = tei[INDEX(i,k)*nn1o2+INDEX(q,s)];

                                long int jq = j*nmo_+q;
                                long int ls = l*nmo_+s;

                                long int qj = q*nmo_+j;
                                long int sl = s*nmo_+l;

                                duma    += eint * ( D2aa[jq*n+ls] + D2ab[jq*n+ls] );
                                dumb    += eint * ( D2bb[jq*n+ls] + D2ab[qj*n+sl] );

                                dumabab += eint * ( D2bb[jq*n+ls] + D2ab[qj*n+sl] );
                                dumbaba += eint * ( D2aa[jq*n+ls] + D2ab[jq*n+ls] );
                            }
                        }
                    }
                }
                // 2: (lj|qs) D(kq;is)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                //double eint = TEI(l,j,q,s,SymmetryPair(hl,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(l,j),1,Qmo_+nQ_*INDEX(q,s),1);
                                double eint = tei[INDEX(l,j)*nn1o2+INDEX(q,s)];

                                long int kq = k*nmo_+q;
                                long int is = i*nmo_+s;

                                long int qk = q*nmo_+k;
                                long int si = s*nmo_+i;

                                duma += eint * ( D2aa[kq*n+is] + D2ab[kq*n+is] );
                                dumb += eint * ( D2bb[kq*n+is] + D2ab[qk*n+si] );

                                dumabab += eint * ( D2aa[kq*n+is] + D2ab[kq*n+is] );
                                dumbaba += eint * ( D2bb[kq*n+is] + D2ab[qk*n+si] );
                            }
                        }
                    }
                }
*/

                // funky sums over 3 indices:
                // - dik (qs|pj) D(pq;ls) + djl (qs|ip) D(kq;ps)
/*
                if ( i == k ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    for (int hs = 0; hs < nirrep_; hs++) {
                                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                            long int s = ss + pitzer_offset_full[hs];
                                            //double eint = TEI(q,s,p,j,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,s),1,Qmo_+nQ_*INDEX(p,j),1);
                                            double eint = tei[INDEX(q,s)*nn1o2+INDEX(p,j)];

                                            long int pq = p*nmo_+q;
                                            long int ls = l*nmo_+s;

                                            long int qp = q*nmo_+p;
                                            long int sl = s*nmo_+l;

                                            duma -= eint * ( D2aa[pq*n+ls] + D2ab[pq*n+ls] );
                                            dumb -= eint * ( D2bb[pq*n+ls] + D2ab[qp*n+sl] );

                                            dumabab -= eint * ( D2bb[pq*n+ls] + D2ab[qp*n+sl] );
                                            dumbaba -= eint * ( D2aa[pq*n+ls] + D2ab[pq*n+ls] );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
*/
/*
                if ( j == l ) {
                    for (int hp = 0; hp < nirrep_; hp++) {
                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                            long int p = pp + pitzer_offset_full[hp];
                            for (int hq = 0; hq < nirrep_; hq++) {
                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                    long int q = qq + pitzer_offset_full[hq];
                                    for (int hs = 0; hs < nirrep_; hs++) {
                                        for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                            long int s = ss + pitzer_offset_full[hs];
                                            //double eint = TEI(q,s,i,p,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,s),1,Qmo_+nQ_*INDEX(i,p),1);
                                            double eint = tei[INDEX(q,s)*nn1o2+INDEX(i,p)];

                                            long int kq = k*nmo_+q;
                                            long int ps = p*nmo_+s;

                                            long int qk = q*nmo_+k;
                                            long int sp = s*nmo_+p;

                                            duma -= eint * ( D2aa[kq*n+ps] + D2ab[kq*n+ps] );
                                            dumb -= eint * ( D2bb[kq*n+ps] + D2ab[qk*n+sp] );

                                            dumabab -= eint * ( D2aa[kq*n+ps] + D2ab[kq*n+ps] );
                                            dumbaba -= eint * ( D2bb[kq*n+ps] + D2ab[qk*n+sp] );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
*/

/*
                // exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];
                                double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                long int kq = k*nmo_+q;
                                long int is = i*nmo_+s;

                                duma -= eint * D2aa[kq*n+is];
                                dumb -= eint * D2bb[kq*n+is];

                                eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                long int jq = j*nmo_+q;
                                long int ls = l*nmo_+s;

                                duma -= eint * D2aa[jq*n+ls];
                                dumb -= eint * D2bb[jq*n+ls];
                            }
                        }
                    }
                }
*/

/*
                // exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)
                for (int hq = 0; hq < nirrep_; hq++) {
                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                        long int q = qq + pitzer_offset_full[hq];
                        for (int hs = 0; hs < nirrep_; hs++) {
                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                long int s = ss + pitzer_offset_full[hs];

                                double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                long int qk = q*nmo_+k;
                                long int is = i*nmo_+s;

                                long int kq = k*nmo_+q;
                                long int si = s*nmo_+i;

                                dumab += eint * D2ab[qk*n+is];
                                dumba += eint * D2ab[kq*n+si];

                                dumabab -= eint * D2ab[kq*n+is];
                                dumbaba -= eint * D2ab[qk*n+si];

                                eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                long int jq = j*nmo_+q;
                                long int sl = s*nmo_+l;

                                long int qj = q*nmo_+j;
                                long int ls = l*nmo_+s;

                                dumab += eint * D2ab[jq*n+sl];
                                dumba += eint * D2ab[qj*n+ls];

                                dumabab -= eint * D2ab[qj*n+sl];
                                dumbaba -= eint * D2ab[jq*n+ls];
                            }
                        }
                    }
                }
*/

/*
                // last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
                for (int hp = 0; hp < nirrep_; hp++) {
                    for (long int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        for (int hq = 0; hq < nirrep_; hq++) {
                            for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                long int q = qq + pitzer_offset_full[hq];

                                //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                long int pq = p*nmo_+q;
                                long int li = l*nmo_+i;

                                duma += eint * D2aa[pq*n+li];
                                dumb += eint * D2bb[pq*n+li];

                                //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                eint = tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                long int kj = k*nmo_+j;

                                duma += eint * D2aa[kj*n+pq];
                                dumb += eint * D2bb[kj*n+pq];
                            }
                        }
                    }
                }
*/

/*
                // last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)
                for (int hp = 0; hp < nirrep_; hp++) {
                    for (long int pp = 0; pp < amopi_[hp]; pp++) {
                        long int p = pp + pitzer_offset_full[hp];
                        for (int hq = 0; hq < nirrep_; hq++) {
                            for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                long int q = qq + pitzer_offset_full[hq];

                                //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                long int pq = p*nmo_+q;
                                long int il = i*nmo_+l;

                                long int qp = q*nmo_+p;
                                long int li = l*nmo_+i;

                                dumab -= eint * D2ab[pq*n+il];
                                dumba -= eint * D2ab[qp*n+li];

                                dumabab += eint * D2ab[qp*n+il];
                                dumbaba += eint * D2ab[pq*n+li];

                                //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                eint = tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                long int kj = k*nmo_+j;
                                long int jk = j*nmo_+k;

                                dumab -= eint * D2ab[jk*n+pq];
                                dumba -= eint * D2ab[kj*n+qp];

                                dumabab += eint * D2ab[kj*n+pq];
                                dumbaba += eint * D2ab[jk*n+qp];
                            }
                        }
                    }
                }

*/

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                long int aabb_sym = (ij + 0 * nh) * 2 * nh + (kl + 1 * nh);
                long int bbaa_sym = (ij + 1 * nh) * 2 * nh + (kl + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                newA[aaaa_sym] += duma;
                newA[aabb_sym] += dumab;
                newA[bbaa_sym] += dumba;
                newA[bbbb_sym] += dumb;

            }
        }
#endif

        // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>

        // N^5 version of A2 coulomb-like terms (1) and (2)
        // 
        // 1: (ik|qs) D(jq;ls)
        // 2: (lj|qs) D(kq;is)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            if ( gems_ab[hjl] == 0 ) continue;

            double * tempa = (double*)malloc(nQ_*gems_ab[hjl]*sizeof(double));
            double * tempb = (double*)malloc(nQ_*gems_ab[hjl]*sizeof(double));

            memset((void*)tempa,'\0',nQ_*gems_ab[hjl]*sizeof(double));
            memset((void*)tempb,'\0',nQ_*gems_ab[hjl]*sizeof(double));

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                for (int Q = 0; Q < nQ_; Q++) {

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {

                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        long int jq = j*nmo_+q;
                        long int ls = l*nmo_+s;
                        
                        long int qj = q*nmo_+j;
                        long int sl = s*nmo_+l;
                                    
                        duma += Qmo_[INDEX(q,s)*nQ_+Q] * (D2aa[jq * nmo_*nmo_ + ls] + D2ab[jq * nmo_*nmo_ + ls]);
                        dumb += Qmo_[INDEX(q,s)*nQ_+Q] * (D2bb[jq * nmo_*nmo_ + ls] + D2ab[qj * nmo_*nmo_ + sl]);

                    }

                    tempa[Q*gems_ab[hjl] + jl] = duma;
                    tempb[Q*gems_ab[hjl] + jl] = dumb;
                }
            }

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {
                
                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;
                
                    int hi = symmetry[i];
                
                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    // check if ij/kl are proper ERPA excitations
                    if ( ij == -999 || kl == -999 ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int Q = 0; Q < nQ_; Q++) {
                        duma += tempa[Q*gems_ab[hjl] + jl] * Qmo2[INDEX(i,k)*nQ_+Q];
                        dumb += tempb[Q*gems_ab[hjl] + jl] * Qmo2[INDEX(i,k)*nQ_+Q];
                    }

                    long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                    // AED: I think term (2) above is just a funny transpose of term (1)
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aaaa_sym = (lk + 0 * nh) * 2 * nh + (ji + 0 * nh);
                    bbbb_sym = (lk + 1 * nh) * 2 * nh + (ji + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;
                }
            }
            free(tempa);
            free(tempb);
        }

        // A2 exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int kq = k*nmo_+q;
                        long int is = i*nmo_+s;

                        duma -= eint * D2aa[kq*n+is];
                        dumb -= eint * D2bb[kq*n+is];

                        //eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                        //long int jq = j*nmo_+q;
                        //long int ls = l*nmo_+s;

                        //duma -= eint * D2aa[jq*n+ls];
                        //dumb -= eint * D2bb[jq*n+ls];
                    }

                    long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(lk,ji)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aaaa_sym = (lk + 0 * nh) * 2 * nh + (ji + 0 * nh);
                    bbbb_sym = (lk + 1 * nh) * 2 * nh + (ji + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                }
            }
        }

        // A2 exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)

        for (long int hjl = 0; hjl < nirrep_; hjl++) {

            for (long int jl = 0; jl < gems_ab[hjl]; jl++) {

                long int j = bas_ab_sym[hjl][jl][0];
                long int l = bas_ab_sym[hjl][jl][1];

                int hj  = symmetry[j];

                for (long int ik = 0; ik < gems_ab[hjl]; ik++) {

                    long int i = bas_ab_sym[hjl][ik][0];
                    long int k = bas_ab_sym[hjl][ik][1];

                    if ( i == j || k == l ) continue;

                    int hi = symmetry[i];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int qs = 0; qs < gems_ab[hjl]; qs++ ) {
                        
                        long int q = bas_ab_sym[hjl][qs][0];
                        long int s = bas_ab_sym[hjl][qs][1];

                        double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                        long int qk = q*nmo_+k;
                        long int is = i*nmo_+s;

                        long int kq = k*nmo_+q;
                        long int si = s*nmo_+i;

                        dumab += eint * D2ab[qk*n+is];
                        dumba += eint * D2ab[kq*n+si];

                        //dumabab -= eint * D2ab[kq*n+is];
                        //dumbaba -= eint * D2ab[qk*n+si];

                    }

                    long int aabb_sym = (ij + 0 * nh) * 2 * nh + (kl + 1 * nh);
                    long int bbaa_sym = (ij + 1 * nh) * 2 * nh + (kl + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aabb_sym = (ji + 0 * nh) * 2 * nh + (lk + 1 * nh);
                    bbaa_sym = (ji + 1 * nh) * 2 * nh + (lk + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                }
            }
        }

        // N^5 version of A2 funky sums over 3 indices:
        // 
        // 1: -dik (qs|pj) D(pq;ls)
        // 2: -djl (qs|ip) D(kq;ps)
        // 
        // Since the 2nd term is a transpose of the 1st, it is not explicitly constructed.
        // 
        // Also, the intermediates for this guy can be built
        // at N^4 cost and can be done outside of these loops (see above)

        for (long int ij = 0; ij < nh; ij++) {

            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];

            int hj = symmetry[j];
            for (int ll = 0; ll < amopi_[hj]; ll++) {

                int l = ll + pitzer_offset[hj];

                if ( i == l ) continue;

                long int il = ibas_erpa_sym[h][i][l];

                // check if il is a proper ERPA excitation

                if ( il == -999 ) continue;

                long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (il + 0 * nh);
                long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (il + 1 * nh);

                newA[aaaa_sym] -= Ia[j*nmo_+l];
                newA[bbbb_sym] -= Ib[j*nmo_+l];

                // AED: I think second funky term is a transpose of the first
                // A(ij,il)(2) = A(ji,li)(1)

                // if ij/il are proper ERPA excitations, then ji/li should be as well

                long int ji = ibas_erpa_sym[h][j][i];
                long int li = ibas_erpa_sym[h][l][i];

                aaaa_sym = (ji + 0 * nh) * 2 * nh + (li + 0 * nh);
                bbbb_sym = (ji + 1 * nh) * 2 * nh + (li + 1 * nh);

                newA[aaaa_sym] -= Ia[j*nmo_+l];
                newA[bbbb_sym] -= Ib[j*nmo_+l];

            }
        }

        // A2: last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double duma = 0.0;
                    double dumb = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        duma += eint * D2aa[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];
                        dumb += eint * D2bb[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+l*nmo_+i];

                    }

                    long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitaitons, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aaaa_sym = (ji + 0 * nh) * 2 * nh + (lk + 0 * nh);
                    bbbb_sym = (ji + 1 * nh) * 2 * nh + (lk + 1 * nh);

                    newA[aaaa_sym] += duma;
                    newA[bbbb_sym] += dumb;

                }

            }

        }

        // A2: last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)

        for (long int hli = 0; hli < nirrep_; hli++) {

            for (long int li = 0; li < gems_ab[hli]; li++) {

                long int l = bas_ab_sym[hli][li][0];
                long int i = bas_ab_sym[hli][li][1];

                int hi  = symmetry[i];

                for (long int jk = 0; jk < gems_ab[hli]; jk++) {

                    long int j = bas_ab_sym[hli][jk][0];
                    long int k = bas_ab_sym[hli][jk][1];

                    if ( i == j || k == l ) continue;

                    int hj = symmetry[j];

                    int hij = hi ^ hj;

                    if ( hij != h ) continue;

                    // check if ij/kl are proper ERPA excitations

                    long int ij = ibas_erpa_sym[h][i][j];
                    long int kl = ibas_erpa_sym[h][k][l];

                    if ( ij == -999 || kl == -999 ) continue;

                    double dumab = 0.0;
                    double dumba = 0.0;

                    for (int pq = 0; pq < gems_ab[hli]; pq++ ) {

                        long int p = bas_ab_sym[hli][pq][0];
                        long int q = bas_ab_sym[hli][pq][1];

                        double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];
                        
                        dumab -= eint * D2ab[p*nmo_*nmo_*nmo_+q*nmo_*nmo_+i*nmo_+l];
                        dumba -= eint * D2ab[q*nmo_*nmo_*nmo_+p*nmo_*nmo_+l*nmo_+i];

                    }

                    long int aabb_sym = (ij + 0 * nh) * 2 * nh + (kl + 1 * nh);
                    long int bbaa_sym = (ij + 1 * nh) * 2 * nh + (kl + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                    // AED: I think second exchange is a transpose of the first
                    // A(ij,kl)(2) = A(ji,lk)(1)

                    // if ij/kl are proper ERPA excitations, then ji/lk should be as well

                    long int ji = ibas_erpa_sym[h][j][i];
                    long int lk = ibas_erpa_sym[h][l][k];

                    aabb_sym = (ji + 0 * nh) * 2 * nh + (lk + 1 * nh);
                    bbaa_sym = (ji + 1 * nh) * 2 * nh + (lk + 0 * nh);

                    newA[aabb_sym] += dumab;
                    newA[bbaa_sym] += dumba;

                }

            }

        }

        std::shared_ptr<Matrix> Amat ( new Matrix(totdim,totdim) );
        std::shared_ptr<Matrix> Bmat ( new Matrix(totdim,totdim) );
        for (int ij = 0; ij < totdim; ij++) {
            for (int kl = 0; kl < totdim; kl++) {
                Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                Bmat->pointer()[ij][kl] = newB[ij*totdim+kl];
            }
        }

        // symmetrize
        // A(kl,ji) = A(ij,kl)
        // A(kl,ji) = A(lk,ij)

        // project out null (and negative) space of A matrix
        std::shared_ptr<Matrix> eigvec ( new Matrix(totdim,totdim) );
        std::shared_ptr<Vector> eigval ( new Vector(totdim) );

        Amat->diagonalize(eigvec,eigval,descending);
        //eigval->print();

        bool prune = true;
        long int newdim = 0;

        // get rid of small/negative eigenvalues of A?
        if ( prune ) {

            newdim = 0;
            int * map = (int*)malloc(totdim*sizeof(int));
            memset((void*)map,'\0',totdim*sizeof(int));
            int * skip = (int*)malloc(totdim*sizeof(int));
            memset((void*)skip,'\0',totdim*sizeof(int));
            for (long int i = 0; i < totdim; i++) {
                double val = eigval->pointer()[i];
                skip[i] = 1;
                if ( val > pruning_threshold ) {
                    map[newdim] = i;
                    skip[i] = 0;
                    newdim++;
                }
            }
            //if ( newdim % 2 == 1 ) { // eigenvalues should come in pairs
            //    newdim--;
            //}
            for (int ij = 0; ij < totdim; ij++) {
                for (int kl = 0; kl < totdim; kl++) {
                    Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                }
            }
            Amat->transform(eigvec);
            Bmat->transform(eigvec);
        }else {
            newdim = totdim;
        }

        for (long int i = 0; i < newdim; i++) {
            for (long int j = 0; j < newdim; j++) {
                newA[i*newdim+j] = Amat->pointer()[i][j];
                newB[i*newdim+j] = Bmat->pointer()[i][j];
            }
        }

        outfile->Printf("\n");
        if ( newdim < totdim ) {
            outfile->Printf("    <<< warning >>> ");
            outfile->Printf("\n");
            outfile->Printf("    reducing dimension of [ k*l, [H, j*i ] ] from %5i to %5i\n",totdim,newdim);
            outfile->Printf("\n");
        }

        int info = 0;
        if ( prune ) {
            info = SymmetricGeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }else {
	    bool * energy_is_real = (bool*)malloc(totdim*sizeof(bool));
            GeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig,energy_is_real);
            free(energy_is_real);
        }

        // back transform wave functions to original untruncated basis
        std::shared_ptr<Matrix>  Cmat (new Matrix(totdim,totdim));
        for (int i = 0; i < newdim; i++) {
            for (int j = 0; j < totdim; j++) {
                double dum = 0.0;
                for (int k = 0; k < newdim; k++) {
                    dum += eigvec->pointer()[j][k] * newB[i * newdim + k];
                }
                Cmat->pointer()[i][j] = dum;
            }
        }

        // oscillator strengths: 
        // symmetry of state is h, consider only cartesian component with same symmetry
       
        if ( info != 0 ) {
            outfile->Printf("\n");
            outfile->Printf("    Error: diagonalization failed.  info = %5i\n",info);
            outfile->Printf("\n");
        }

        outfile->Printf("    state");
        outfile->Printf("          energy (Eh)");
        outfile->Printf("       ex energy (Eh)");
        outfile->Printf("       ex energy (eV)");
        outfile->Printf("       f, osc. strength\n");

        for (long int state = newdim-1; state >= 0; state--){
          
            if ( eig[state] > 0.0 ) {
                // normalize C:
                double nrm = 0.0;

                for (int ij = 0; ij < nh; ij++) {
                    int i  = bas_erpa_sym[h][ij][0];
                    int j  = bas_erpa_sym[h][ij][1];
                    int hi = symmetry[i];
                    int hj = symmetry[j];
                    for( int kk = 0; kk < amopi_[hi]; kk++){
                        int k   = kk + pitzer_offset[hi];
                        int kj  = ibas_erpa_sym[h][k][j];

                        // check if kj is a proper ERPA excitations
                        if ( kj == -999 ) continue;

                        nrm    += Cmat->pointer()[state][ij      ] * Cmat->pointer()[state][kj      ] * D1a[k*amo_+i];
                        nrm    += Cmat->pointer()[state][ij + nh ] * Cmat->pointer()[state][kj + nh ] * D1b[k*amo_+i];
                    }
                    for( int ll = 0; ll < amopi_[hj]; ll++){
                        int l  = ll + pitzer_offset[hj];
                        int il = ibas_erpa_sym[h][i][l];

                        // check if il is a proper ERPA excitations
                        if ( il == -999 ) continue;

                        nrm   -= Cmat->pointer()[state][ij      ] * Cmat->pointer()[state][il      ] * D1a[j*amo_+l];
                        nrm   -= Cmat->pointer()[state][ij + nh ] * Cmat->pointer()[state][il + nh ] * D1b[j*amo_+l];
                    }
                }


                nrm = sqrt(nrm);
                for (int ij = 0; ij < nh; ij++) {
                    Cmat->pointer()[state][ij     ] /= nrm;
                    Cmat->pointer()[state][ij + nh] /= nrm;
                }
                
                double dumx = 0.0;
                double dumy = 0.0;
                double dumz = 0.0;
    
                if ( h == hx_ ) {
                   
                    dumx = TransitionDipoleMoment(Cmat->pointer()[state], D1a, D1b, nh,h,"X");

                }
                if ( h == hy_ ) {
    
                    dumy = TransitionDipoleMoment(Cmat->pointer()[state], D1a, D1b, nh,h,"Y");

                }
                if ( h == hz_ ) {
    
                    dumz = TransitionDipoleMoment(Cmat->pointer()[state], D1a, D1b, nh,h,"Z");

                }

                double * TDa = (double*)malloc(nso_*nso_*sizeof(double));
                double * TDb = (double*)malloc(nso_*nso_*sizeof(double));

                memset((void*)TDa,'\0',nso_*nso_*sizeof(double));
                memset((void*)TDb,'\0',nso_*nso_*sizeof(double));

                // build transition density matrix:
                bool is_singlet = true;
                double max = -999.0;
                int ijmax = 0;
                for ( int ij = 0; ij < nh; ij++) {
                    if (fabs(Cmat->pointer()[state][ij]) > max ) {
                        max = fabs(Cmat->pointer()[state][ij]);
                        ijmax = ij;
                    }
                }
                //printf(" singlet? %20.12lf %20.12lf\n",Cmat->pointer()[state][ijmax],Cmat->pointer()[state][ijmax+nh]);
                //if (fabs(Cmat->pointer()[state][ijmax] - Cmat->pointer()[state][ijmax+nh]) < 1e-3)
                //if (fabs((Cmat->pointer()[state][ijmax] - Cmat->pointer()[state][ijmax+nh])/Cmat->pointer()[state][ijmax]) < 1e-3) {
                    if ( Cmat->pointer()[state][ijmax] * Cmat->pointer()[state][ijmax+nh] >= 0.0 ) {

                        //if (fabs(Cmat->pointer()[state][ijmax] - Cmat->pointer()[state][ijmax+nh]) > 1e-3) {
                        //    printf("%20.12lf %20.12lf\n",Cmat->pointer()[state][ijmax],Cmat->pointer()[state][ijmax+nh]);
                        //}

                        is_singlet = true;
                    //}
                    //printf("this one is a singlet\n");
                }else {
                    is_singlet = false;
                }

                //outfile->Printf("WARNING: AED removed restriction to include only singlets!\n");
                if ( is_singlet ) {
                    TransitionDensity(Cmat->pointer()[state], D1a, D1b, TDa,TDb,nh,h);
//if ( is_singlet ) printf("singlet:\n");
//else              printf("triplet:\n");
//for (int i = 0; i < nso_; i++) {
//    for (int j = 0; j < nso_; j++) {
//        printf("%5i %5i %20.12lf %20.12lf\n",i,j,TDa[i*nso_+j],TDb[i*nso_+j]);
//    }
//}

                    // use transition density matrix to update ground-state 2-RDM
                    for (int p = 0; p < nso_; p++) {
                        for (int q = 0; q < nso_; q++) {
                            for (int r = 0; r < nso_; r++) {
                                for (int s = 0; s < nso_; s++) {
                                    double ga_rp = TDa[p*nso_+r]; // careful with transpose ...
                                    double gb_qs = TDb[s*nso_+q];
                                    //double ga_rp = TDa[p*nso_+r]; // careful with transpose ...
                                    //double gb_qs = TDb[q*nso_+s];
                                    //double ga_rp = TDa[p*nso_+r]; // careful with transpose ...
                                    //double gb_qs = TDb[s*nso_+q];
                                    D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_rp * gb_qs;// * 0.5;
                                    //D2ab_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_rp * gb_qs * 0.5;
                                    //D2ab_prime[r*nso_*nso_*nso_+s*nso_*nso_+p*nso_+q] += ga_rp * gb_qs * 0.5;
                                    //if ( p != q && r != s) {

                                        double ga_qs = TDa[s*nso_+q];
                                        D2aa_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += ga_rp * ga_qs;// * 0.5;

                                        double gb_rp = TDb[p*nso_+r]; // careful with transpose ...
                                        D2bb_prime[p*nso_*nso_*nso_+q*nso_*nso_+r*nso_+s] += gb_rp * gb_qs;// * 0.5;
                                    //}
                                }
                            }
                        }
                    }
                }
                free(TDa);
                free(TDb);


                if ( prune ) eig[state] = 1.0 / eig[state];

                double val = 2./3. * eig[state] * (dumx*dumx+dumy*dumy+dumz*dumz);
                
                outfile->Printf(" OLDWAY   %5i %20.12lf %20.12lf %20.12lf %20.12lf\n",newdim-state,eig[state]+e2+e1+enuc_,eig[state],eig[state] * 27.21138, val);
                outfile->Printf("                Transition Moment: X: %10.6lf    Y: %10.6lf   Z: %10.6lf\n", dumx, dumy, dumz);
                
                bool * skip = (bool*)malloc(2*nh*sizeof(bool));
                memset((void*)skip,'\0',2*nh*sizeof(bool));
                
                for (int k = 0; k < number_coefficients; k++) {
                    double max = -999.0e99;
                    int maxij = 0;
                    for (int ij = 0; ij < 2*nh; ij++) {
                        if ( skip[ij] ) continue;
                        if ( fabs(Cmat->pointer()[state][ij]) > max ) {
                            max = fabs(Cmat->pointer()[state][ij]);
                            maxij = ij;
                        }
                    }
                    skip[maxij] = true;
                    int id = maxij;
                    if ( id >= nh ) id -= nh;
                    int i = bas_erpa_sym[h][id][0];
                    int j = bas_erpa_sym[h][id][1];
                    int hi = symmetry[i];
                    int hj = symmetry[j];
                    //outfile->Printf("                %s %3i, symm %3i -> %3i, symm %3i:  %20.6lf\n",(id == maxij) ? "alpha" : "beta ",i, hi, j, hj, Cmat->pointer()[state][maxij]);
                    if ( fabs(Cmat->pointer()[state][maxij]) > coefficient_threshold ) {
                        outfile->Printf("                %s %3i -> %3i  %20.6lf\n",(id == maxij) ? "alpha" : "beta ",i+1, j+1, Cmat->pointer()[state][maxij]);
                    }

                }

                free(skip);
                
                outfile->Printf("\n");                   
            } 
        }
   
        free(newA);
        free(newB);
        free(cc);
        free(eig);
    }

    // check energy
    double rpa_e2 = 0.0;
    double ac_e2 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    rpa_e2 +=       eint * D2ab_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                    rpa_e2 += 0.5 * eint * D2aa_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                    rpa_e2 += 0.5 * eint * D2bb_prime[(i*nmo_+j)*n+(k*nmo_+l)];

                    int nactive = 0;
                    if ( active_reference_orbitals_[i] ) nactive++;
                    if ( active_reference_orbitals_[j] ) nactive++;
                    if ( active_reference_orbitals_[k] ) nactive++;
                    if ( active_reference_orbitals_[l] ) nactive++;

                    int ninactive = 0;
                    if ( inactive_reference_orbitals_[i] ) ninactive++;
                    if ( inactive_reference_orbitals_[j] ) ninactive++;
                    if ( inactive_reference_orbitals_[k] ) ninactive++;
                    if ( inactive_reference_orbitals_[l] ) ninactive++;

                    //if ( nactive != 4 && ninactive != 4 && ( nactive + ninactive ) != 0  )
                    if ( (nactive + ninactive) != 4 && ( nactive + ninactive ) != 0  ) {

                        ac_e2 +=       eint * D2ab_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                        ac_e2 += 0.5 * eint * D2aa_prime[(i*nmo_+j)*n+(k*nmo_+l)];
                        ac_e2 += 0.5 * eint * D2bb_prime[(i*nmo_+j)*n+(k*nmo_+l)];

                    }
                }
            }
        }
    }
    e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }
    //rpa_e2 += e2;
    //rpa_e2 *= 0.5;

    // check traces:
    double traa = 0.0;
    double trbb = 0.0;
    double trab = 0.0;
    for (long int i = 0; i < nmo_; i++) {
        for (long int j = 0; j < nmo_; j++) {
            //double duma = 0.0;
            //double dumb = 0.0;
            //for (long int k = 0; k < nmo_; k++) {
            //    duma += D2ab_prime[(i*nmo_+k)*n+(j*nmo_+k)] * 0.5;
            //    dumb += D2ab_prime[(k*nmo_+i)*n+(k*nmo_+j)] * 0.5;
            //}
            //printf("%20.12lf %20.12lf %20.12lf\n",D1a[i*nmo_+j],duma/nbeta_,D1a[i*nmo_+j]-duma/nbeta_);
            long int ij = i*nmo_+j;
            traa += D2ab_prime[(i*nmo_+j)*n+(i*nmo_+j)];// * 0.5;
            trbb += D2aa_prime[(i*nmo_+j)*n+(i*nmo_+j)];// * 0.5;
            trab += D2bb_prime[(i*nmo_+j)*n+(i*nmo_+j)];// * 0.5;
        }
    }

    //outfile->Printf("\n");                   
    //outfile->Printf("    ==> ERPA correlation information: <==\n");
    //outfile->Printf("\n");                   
    //outfile->Printf("        Tr(D2aa):            %20.12lf\n",traa);
    //outfile->Printf("        Tr(D2bb):            %20.12lf\n",trbb);
    //outfile->Printf("        Tr(D2ab):            %20.12lf\n",trab);
    //outfile->Printf("\n");                   
    //outfile->Printf("        * ERPA total energy: %20.12lf\n",e1+rpa_e2+enuc_);
    outfile->Printf("\n");                   
    outfile->Printf("    * AC1 total energy:  %20.12lf\n",e1+e2 + 0.5 * ac_e2 +enuc_);
    outfile->Printf("\n");                   

    free(Ia);
    free(Ib);
    free(tei);
    free(h1);
    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);
    free(D2aa_prime);
    free(D2bb_prime);
    free(D2ab_prime);
    if ( !is_df_ ) {
        free(Qmo2);
    }
}

void ERPASolver::OldExtendedRPA(std::string type) {

    double * A   = (double*)malloc(16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * B   = (double*)malloc(16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * cc  = (double*)malloc(16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * eig = (double*)malloc(16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));

    double * h1   = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1a  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D1b  = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * D2aa = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2bb = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));
    double * D2ab = (double*)malloc(nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)A,'\0',16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)B,'\0',16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)cc,'\0',16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)eig,'\0',16 * nmo_*nmo_*nmo_*nmo_*sizeof(double));

    memset((void*)h1,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1a,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D1b,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)D2aa,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2bb,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));
    memset((void*)D2ab,'\0',nmo_*nmo_*nmo_*nmo_*sizeof(double));

    // get one- and two-electron integrals
    // one-electron integrals
    long int offset = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (long int i = 0; i < nmopi_[h] - frzvpi_[h]; i++) {
            int ii = i + pitzer_offset_full[h];
            for (long int j = i; j < nmopi_[h] - frzvpi_[h]; j++) {
                int jj = j + pitzer_offset_full[h];
                h1[ii*nmo_+jj] = oei_full_sym_[offset + INDEX(i,j)];
                h1[jj*nmo_+ii] = oei_full_sym_[offset + INDEX(i,j)];
            }
        }
        offset += ( nmopi_[h] - frzvpi_[h] ) * ( nmopi_[h] - frzvpi_[h] + 1 ) / 2;
    }

    // build two-electron integrals
    if ( !is_df_ ) {
        throw PsiException("erpa works only with scf_type df/cd for now.",__FILE__,__LINE__);
    }
    long int nn1o2 = nmo_*(nmo_+1)/2;
    double * tei = (double*)malloc(nn1o2*nn1o2*sizeof(double));
    memset((void*)tei,'\0',nn1o2*nn1o2*sizeof(double));
    F_DGEMM('t','n',nmo_*(nmo_+1)/2,nmo_*(nmo_+1)/2,nQ_,1.0,Qmo_,nQ_,Qmo_,nQ_,0.0,tei,nmo_*(nmo_+1)/2);

    // read TDPM from disk
    ReadTPDM(D2aa,D2bb,D2ab,D1a,D1b,tei);

    long int n = nmo_*nmo_;

    // check energy
    double e2 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            long int kl = k*nmo_+l;
            for (long int i = 0; i < nmo_; i++) {
                for (long int j = 0; j < nmo_; j++) {
                    long int ij = i*nmo_+j;
                    //double eint = TEI(i,k,j,l,SymmetryPair(symmetry_full[i],symmetry_full[k]));//
                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(j,l)];
                    e2 +=       eint * D2ab[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2aa[(i*nmo_+j)*n+(k*nmo_+l)];
                    e2 += 0.5 * eint * D2bb[(i*nmo_+j)*n+(k*nmo_+l)];
                }
            }
        }
    }
    double e1 = 0.0;
    for (long int k = 0; k < nmo_; k++) {
        for (long int l = 0; l < nmo_; l++) {
            e1 += D1a[k*nmo_+l] * h1[k*nmo_+l];
            e1 += D1b[k*nmo_+l] * h1[k*nmo_+l];
        }
    }

    //printf("%20.12lf %20.12lf %20.12lf\n",e1,e2,e1+e2);

    //Process::environment.globals["CURRENT ENERGY"] = e1 + e2 + enuc_;
    if ( ref_rdm_ == "CI" ) {
        double efci = Process::environment.globals["FCI TOTAL ENERGY"];
        if ( fabs (efci - (e1+e2+enuc_)) ) {
            printf("%20.12lf %20.12lf\n",efci,e1+e2+enuc_);
            printf("\n");
            printf("symmetry_energy_order [");
            for (int i = 0; i < nmo_-1; i++) {
                printf("%2i,",symmetry_energy_order[i]);
            }
            printf("%2i",symmetry_energy_order[nmo_-1]);
            printf("]\n");
            printf("\n");
            //throw PsiException("Energy from density on disk does not match FCI energy.  Check orbital symmetries.",__FILE__,__LINE__);
        }
    }

    // build B(ij,kl) = <| [ k*l, j*i ] |>
    for (int hi = 0; hi < nirrep_; hi++) {
        for (long int ii = 0; ii < amopi_[hi]; ii++) {
            long int i = ii + pitzer_offset_full[hi];
            for (int hj = 0; hj < nirrep_; hj++) {
                for (long int jj = 0; jj < amopi_[hj]; jj++) {
                    long int j = jj + pitzer_offset_full[hj];
                    long int ij = i*nmo_+j;
                    for (int hk = 0; hk < nirrep_; hk++) {
                        for (long int kk = 0; kk < amopi_[hk]; kk++) {
                            long int k = kk + pitzer_offset_full[hk];
                            for (int hl = 0; hl < nirrep_; hl++) {
                                for (long int ll = 0; ll < amopi_[hl]; ll++) {
                                    long int l = ll + pitzer_offset_full[hl];

                                    if ( SymmetryPair(hi,hj) != SymmetryPair(hk,hl) ) continue;

                                    double duma    = 0.0;
                                    double dumb    = 0.0;
                                    double dumabab = 0.0;
                                    double dumbaba = 0.0;

                                    if ( j == l ) {
                                        duma    += D1a[k*nmo_+i];
                                        dumb    += D1b[k*nmo_+i];

                                        dumabab += D1a[k*nmo_+i];
                                        dumbaba += D1b[k*nmo_+i];
                                    }
                                    if ( i == k ) {
                                        duma    -= D1a[j*nmo_+l];
                                        dumb    -= D1b[j*nmo_+l];

                                        dumabab -= D1b[j*nmo_+l];
                                        dumbaba -= D1a[j*nmo_+l];
                                    }

                                    // aa,aa
                                    long int aaaa = UIndex(i,j,k,l,0,0,nmo_);
                                    B[aaaa] += duma;

                                    // bb,bb
                                    long int bbbb = UIndex(i,j,k,l,3,3,nmo_);
                                    B[bbbb] += dumb;

                                    // ab,ab
                                    long int abab = UIndex(i,j,k,l,1,1,nmo_);
                                    B[abab] += dumabab;

                                    // ba,ba
                                    long int baba = UIndex(i,j,k,l,2,2,nmo_);
                                    B[baba] += dumbaba;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // build A1(ij,kl) = <| [ k*l, [H1, j*i ] |>
    for (int hi = 0; hi < nirrep_; hi++) {
        for (long int ii = 0; ii < amopi_[hi]; ii++) {
            long int i = ii + pitzer_offset_full[hi];
            for (int hj = 0; hj < nirrep_; hj++) {
                for (long int jj = 0; jj < amopi_[hj]; jj++) {
                    long int j = jj + pitzer_offset_full[hj];
                    long int ij = i*nmo_+j;
                    for (int hk = 0; hk < nirrep_; hk++) {
                        for (long int kk = 0; kk < amopi_[hk]; kk++) {
                            long int k = kk + pitzer_offset_full[hk];
                            for (int hl = 0; hl < nirrep_; hl++) {
                                for (long int ll = 0; ll < amopi_[hl]; ll++) {
                                    long int l = ll + pitzer_offset_full[hl];
                                    long int kl = k*nmo_+l;

                                    if ( SymmetryPair(hi,hj) != SymmetryPair(hk,hl) ) continue;

                                    double duma    = 0.0;
                                    double dumb    = 0.0;
                                    double dumabab = 0.0;
                                    double dumbaba = 0.0;

                                    duma    += D1a[k*nmo_+i] * h1[l*nmo_+j];
                                    duma    += D1a[j*nmo_+l] * h1[i*nmo_+k];

                                    dumabab += D1a[k*nmo_+i] * h1[l*nmo_+j];
                                    dumabab += D1b[j*nmo_+l] * h1[i*nmo_+k];

                                    dumbaba += D1b[k*nmo_+i] * h1[l*nmo_+j];
                                    dumbaba += D1a[j*nmo_+l] * h1[i*nmo_+k];

                                    dumb    += D1b[k*nmo_+i] * h1[l*nmo_+j];
                                    dumb    += D1b[j*nmo_+l] * h1[i*nmo_+k];

                                    if ( j == l ) {
                                        for (int hq = 0; hq < nirrep_; hq++) {
                                            for (int qq = 0; qq < amopi_[hq]; qq++) {
                                                long int q = qq + pitzer_offset_full[hq];
                                                duma    -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                                                dumb    -= D1b[k*nmo_+q] * h1[i*nmo_+q];

                                                dumabab -= D1a[k*nmo_+q] * h1[i*nmo_+q];
                                                dumbaba -= D1b[k*nmo_+q] * h1[i*nmo_+q];
                                            }
                                        }
                                    }
                                    if ( i == k ) {
                                        for (int hq = 0; hq < nirrep_; hq++) {
                                            for (int qq = 0; qq < amopi_[hq]; qq++) {
                                                long int q = qq + pitzer_offset_full[hq];
                                                duma    -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                                                dumb    -= D1b[q*nmo_+l] * h1[q*nmo_+j];

                                                dumabab -= D1b[q*nmo_+l] * h1[q*nmo_+j];
                                                dumbaba -= D1a[q*nmo_+l] * h1[q*nmo_+j];
                                            }
                                        }
                                    }
                                    // aa,aa
                                    long int aaaa = UIndex(i,j,k,l,0,0,nmo_);
                                    A[aaaa] += duma;

                                    // bb,bb
                                    long int bbbb = UIndex(i,j,k,l,3,3,nmo_);
                                    A[bbbb] += dumb;

                                    // ab,ab
                                    long int abab = UIndex(i,j,k,l,1,1,nmo_);
                                    A[abab] += dumabab;

                                    // ba,ba
                                    long int baba = UIndex(i,j,k,l,2,2,nmo_);
                                    A[baba] += dumbaba;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // build A2(ij,kl) = <| [ k*l, [H2, j*i ] |>
    for (int hk = 0; hk < nirrep_; hk++) {
        for (long int kk = 0; kk < amopi_[hk]; kk++) {
            long int k = kk + pitzer_offset_full[hk];
            for (int hl = 0; hl < nirrep_; hl++) {
                for (long int ll = 0; ll < amopi_[hl]; ll++) {
                    long int l = ll + pitzer_offset_full[hl];
                    long int kl = k*nmo_+l;
                    for (int hi = 0; hi < nirrep_; hi++) {
                        for (long int ii = 0; ii < amopi_[hi]; ii++) {
                            long int i = ii + pitzer_offset_full[hi];
                            for (int hj = 0; hj < nirrep_; hj++) {
                                for (long int jj = 0; jj < amopi_[hj]; jj++) {

                                    if ( SymmetryPair(hi,hj) != SymmetryPair(hk,hl) ) continue;

                                    long int j = jj + pitzer_offset_full[hj];
                                    long int ij = i*nmo_+j;

                                    double duma = 0.0;
                                    double dumb = 0.0;

                                    double dumab = 0.0;
                                    double dumba = 0.0;

                                    double dumabab = 0.0;
                                    double dumbaba = 0.0;

                                    // 2 coulomb-like terms:

                                    // 1: (ik|qs) D(jq;ls)
                                    for (int hq = 0; hq < nirrep_; hq++) {
                                        for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                            long int q = qq + pitzer_offset_full[hq];
                                            for (int hs = 0; hs < nirrep_; hs++) {
                                                for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                                    long int s = ss + pitzer_offset_full[hs];
                                                    //double eint = TEI(i,k,q,s,SymmetryPair(hi,hk));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,k),1,Qmo_+nQ_*INDEX(q,s),1);
                                                    double eint = tei[INDEX(i,k)*nn1o2+INDEX(q,s)];

                                                    long int jq = j*nmo_+q;
                                                    long int ls = l*nmo_+s;

                                                    long int qj = q*nmo_+j;
                                                    long int sl = s*nmo_+l;

                                                    duma    += eint * ( D2aa[jq*n+ls] + D2ab[jq*n+ls] );
                                                    dumb    += eint * ( D2bb[jq*n+ls] + D2ab[qj*n+sl] );

                                                    dumabab += eint * ( D2bb[jq*n+ls] + D2ab[qj*n+sl] );
                                                    dumbaba += eint * ( D2aa[jq*n+ls] + D2ab[jq*n+ls] );
                                                }
                                            }
                                        }
                                    }
                                    // 2: (lj|qs) D(kq;is)
                                    for (int hq = 0; hq < nirrep_; hq++) {
                                        for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                            long int q = qq + pitzer_offset_full[hq];
                                            for (int hs = 0; hs < nirrep_; hs++) {
                                                for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                                    long int s = ss + pitzer_offset_full[hs];
                                                    //double eint = TEI(l,j,q,s,SymmetryPair(hl,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(l,j),1,Qmo_+nQ_*INDEX(q,s),1);
                                                    double eint = tei[INDEX(l,j)*nn1o2+INDEX(q,s)];

                                                    long int kq = k*nmo_+q;
                                                    long int is = i*nmo_+s;

                                                    long int qk = q*nmo_+k;
                                                    long int si = s*nmo_+i;

                                                    duma += eint * ( D2aa[kq*n+is] + D2ab[kq*n+is] );
                                                    dumb += eint * ( D2bb[kq*n+is] + D2ab[qk*n+si] );

                                                    dumabab += eint * ( D2aa[kq*n+is] + D2ab[kq*n+is] );
                                                    dumbaba += eint * ( D2bb[kq*n+is] + D2ab[qk*n+si] );
                                                }
                                            }
                                        }
                                    }

                                    // funky sums over 3 indices:
                                    // dik [- (qs|pj) D(pq;ls) - (qs|ip) D(kq;ps)]

                                    if ( i == k ) {
                                        for (int hp = 0; hp < nirrep_; hp++) {
                                            for (long int pp = 0; pp < amopi_[hp]; pp++) {
                                                long int p = pp + pitzer_offset_full[hp];
                                                for (int hq = 0; hq < nirrep_; hq++) {
                                                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                                        long int q = qq + pitzer_offset_full[hq];
                                                        for (int hs = 0; hs < nirrep_; hs++) {
                                                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                                                long int s = ss + pitzer_offset_full[hs];
                                                                //double eint = TEI(q,s,p,j,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,s),1,Qmo_+nQ_*INDEX(p,j),1);
                                                                double eint = tei[INDEX(q,s)*nn1o2+INDEX(p,j)];

                                                                long int pq = p*nmo_+q;
                                                                long int ls = l*nmo_+s;

                                                                long int qp = q*nmo_+p;
                                                                long int sl = s*nmo_+l;

                                                                duma -= eint * ( D2aa[pq*n+ls] + D2ab[pq*n+ls] );
                                                                dumb -= eint * ( D2bb[pq*n+ls] + D2ab[qp*n+sl] );

                                                                dumabab -= eint * ( D2bb[pq*n+ls] + D2ab[qp*n+sl] );
                                                                dumbaba -= eint * ( D2aa[pq*n+ls] + D2ab[pq*n+ls] );
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if ( j == l ) {
                                        for (int hp = 0; hp < nirrep_; hp++) {
                                            for (long int pp = 0; pp < amopi_[hp]; pp++) {
                                                long int p = pp + pitzer_offset_full[hp];
                                                for (int hq = 0; hq < nirrep_; hq++) {
                                                    for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                                        long int q = qq + pitzer_offset_full[hq];
                                                        for (int hs = 0; hs < nirrep_; hs++) {
                                                            for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                                                long int s = ss + pitzer_offset_full[hs];
                                                                //double eint = TEI(q,s,i,p,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,s),1,Qmo_+nQ_*INDEX(i,p),1);
                                                                double eint = tei[INDEX(q,s)*nn1o2+INDEX(i,p)];

                                                                long int kq = k*nmo_+q;
                                                                long int ps = p*nmo_+s;

                                                                long int qk = q*nmo_+k;
                                                                long int sp = s*nmo_+p;

                                                                duma -= eint * ( D2aa[kq*n+ps] + D2ab[kq*n+ps] );
                                                                dumb -= eint * ( D2bb[kq*n+ps] + D2ab[qk*n+sp] );

                                                                dumabab -= eint * ( D2aa[kq*n+ps] + D2ab[kq*n+ps] );
                                                                dumbaba -= eint * ( D2bb[kq*n+ps] + D2ab[qk*n+sp] );
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    // exchange-like terms: - (qj|ls) D(kq;is) - (is|qk) D(jq;ls)
                                    for (int hq = 0; hq < nirrep_; hq++) {
                                        for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                            long int q = qq + pitzer_offset_full[hq];
                                            for (int hs = 0; hs < nirrep_; hs++) {
                                                for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                                    long int s = ss + pitzer_offset_full[hs];
                                                    //double eint = TEI(q,j,l,s,SymmetryPair(hq,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,j),1,Qmo_+nQ_*INDEX(l,s),1);
                                                    double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                                    long int kq = k*nmo_+q;
                                                    long int is = i*nmo_+s;

                                                    duma -= eint * D2aa[kq*n+is];
                                                    dumb -= eint * D2bb[kq*n+is];

                                                    //eint = TEI(i,s,q,k,SymmetryPair(hi,hs));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,s),1,Qmo_+nQ_*INDEX(q,k),1);
                                                    eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                                    long int jq = j*nmo_+q;
                                                    long int ls = l*nmo_+s;

                                                    duma -= eint * D2aa[jq*n+ls];
                                                    dumb -= eint * D2bb[jq*n+ls];
                                                }
                                            }
                                        }
                                    }

                                    // exchange-like terms (contribute to aa/bb blocks): - (qj|ls) D(qk;is) - (is|qk) D(jq;sl)
                                    for (int hq = 0; hq < nirrep_; hq++) {
                                        for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                            long int q = qq + pitzer_offset_full[hq];
                                            for (int hs = 0; hs < nirrep_; hs++) {
                                                for (long int ss = 0; ss < amopi_[hs]; ss++) {
                                                    long int s = ss + pitzer_offset_full[hs];

                                                    //double eint = TEI(q,j,l,s,SymmetryPair(hq,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(q,j),1,Qmo_+nQ_*INDEX(l,s),1);
                                                    double eint = tei[INDEX(q,j)*nn1o2+INDEX(l,s)];

                                                    long int qk = q*nmo_+k;
                                                    long int is = i*nmo_+s;

                                                    long int kq = k*nmo_+q;
                                                    long int si = s*nmo_+i;

                                                    dumab += eint * D2ab[qk*n+is];
                                                    dumba += eint * D2ab[kq*n+si];

                                                    dumabab -= eint * D2ab[kq*n+is];
                                                    dumbaba -= eint * D2ab[qk*n+si];

                                                    //eint = TEI(i,s,q,k,SymmetryPair(hi,hs));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,s),1,Qmo_+nQ_*INDEX(q,k),1);
                                                    eint = tei[INDEX(i,s)*nn1o2+INDEX(q,k)];

                                                    long int jq = j*nmo_+q;
                                                    long int sl = s*nmo_+l;

                                                    long int qj = q*nmo_+j;
                                                    long int ls = l*nmo_+s;

                                                    dumab += eint * D2ab[jq*n+sl];
                                                    dumba += eint * D2ab[qj*n+ls];

                                                    dumabab -= eint * D2ab[qj*n+sl];
                                                    dumbaba -= eint * D2ab[jq*n+ls];
                                                }
                                            }
                                        }
                                    }

                                    // last funny term: (pj|qk) D(pq;li) + (ip|lq) D(kj;pq)
                                    for (int hp = 0; hp < nirrep_; hp++) {
                                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                                            long int p = pp + pitzer_offset_full[hp];
                                            for (int hq = 0; hq < nirrep_; hq++) {
                                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                                    long int q = qq + pitzer_offset_full[hq];

                                                    //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                                    double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                                    long int pq = p*nmo_+q;
                                                    long int li = l*nmo_+i;

                                                    duma += eint * D2aa[pq*n+li];
                                                    dumb += eint * D2bb[pq*n+li];

                                                    //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                                    eint = tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                                    long int kj = k*nmo_+j;

                                                    duma += eint * D2aa[kj*n+pq];
                                                    dumb += eint * D2bb[kj*n+pq];
                                                }
                                            }
                                        }
                                    }
                                    // last funny term ( contributes to aa/bb ): (pj|qk) D(pq;il) + (ip|lq) D(jk;pq)
                                    for (int hp = 0; hp < nirrep_; hp++) {
                                        for (long int pp = 0; pp < amopi_[hp]; pp++) {
                                            long int p = pp + pitzer_offset_full[hp];
                                            for (int hq = 0; hq < nirrep_; hq++) {
                                                for (long int qq = 0; qq < amopi_[hq]; qq++) {
                                                    long int q = qq + pitzer_offset_full[hq];

                                                    //double eint = TEI(p,j,q,k,SymmetryPair(hp,hj));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(p,j),1,Qmo_+nQ_*INDEX(q,k),1);
                                                    double eint = tei[INDEX(p,j)*nn1o2+INDEX(q,k)];

                                                    long int pq = p*nmo_+q;
                                                    long int il = i*nmo_+l;

                                                    long int qp = q*nmo_+p;
                                                    long int li = l*nmo_+i;

                                                    dumab -= eint * D2ab[pq*n+il];
                                                    dumba -= eint * D2ab[qp*n+li];

                                                    dumabab += eint * D2ab[qp*n+il];
                                                    dumbaba += eint * D2ab[pq*n+li];

                                                    //eint = TEI(i,p,l,q,SymmetryPair(hi,hp));//C_DDOT(nQ_,Qmo_ + nQ_*INDEX(i,p),1,Qmo_+nQ_*INDEX(l,q),1);
                                                    eint = tei[INDEX(i,p)*nn1o2+INDEX(l,q)];

                                                    long int kj = k*nmo_+j;
                                                    long int jk = j*nmo_+k;

                                                    dumab -= eint * D2ab[jk*n+pq];
                                                    dumba -= eint * D2ab[kj*n+qp];

                                                    dumabab += eint * D2ab[kj*n+pq];
                                                    dumbaba += eint * D2ab[jk*n+qp];
                                                }
                                            }
                                        }
                                    }

                                    //A[kl*n+ij]         += duma + dumab;// + dumb + dumab + dumba;

                                    // aa,aa
                                    long int aaaa = UIndex(i,j,k,l,0,0,nmo_);
                                    A[aaaa] += duma;

                                    // bb,bb
                                    long int bbbb = UIndex(i,j,k,l,3,3,nmo_);
                                    A[bbbb] += dumb;

                                    // bb,aa
                                    long int bbaa = UIndex(i,j,k,l,3,0,nmo_);
                                    A[bbaa] += dumba;

                                    // aa,bb
                                    long int aabb = UIndex(i,j,k,l,0,3,nmo_);
                                    A[aabb] += dumab;

                                    // ab,ab
                                    long int abab = UIndex(i,j,k,l,1,1,nmo_);
                                    A[abab] += dumabab;

                                    // ba,ba
                                    long int baba = UIndex(i,j,k,l,2,2,nmo_);
                                    A[baba] += dumbaba;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    double * newB = (double*)malloc(16*n*n*sizeof(double));
    double * newA = (double*)malloc(16*n*n*sizeof(double));
    memset((void*)newA,'\0',16*n*n*sizeof(double));
    memset((void*)newB,'\0',16*n*n*sizeof(double));

    int target_state_sym = 7;

    for (int h = 0; h < nirrep_; h++) {

        long int nh = gems_erpa[h];
        //long int nh = gems_ab[h];
        if ( nh == 0 ) continue;

        long int totdim = nh;
        if ( type == "AA" ) totdim = 2*nh;

        memset((void*)newA,'\0',16*nh*nh*sizeof(double));
        memset((void*)newB,'\0',16*nh*nh*sizeof(double));

        outfile->Printf("\n");
        outfile->Printf("    Symmetry: %5i\n",h);
        outfile->Printf("\n");


        for (long int ij = 0; ij < nh; ij++) {
            long int i = bas_erpa_sym[h][ij][0];
            long int j = bas_erpa_sym[h][ij][1];
            //long int i = bas_ab_sym[h][ij][0];
            //long int j = bas_ab_sym[h][ij][1];
            for (long int kl = 0; kl < nh; kl++) {
                long int k = bas_erpa_sym[h][kl][0];
                long int l = bas_erpa_sym[h][kl][1];
                //long int k = bas_ab_sym[h][kl][0];
                //long int l = bas_ab_sym[h][kl][1];

                long int aaaa     = UIndex(i,j,k,l,0,0,nmo_);
                long int aabb     = UIndex(i,j,k,l,0,3,nmo_);
                long int bbaa     = UIndex(i,j,k,l,3,0,nmo_);
                long int abab     = UIndex(i,j,k,l,1,1,nmo_);
                long int baba     = UIndex(i,j,k,l,2,2,nmo_);
                long int bbbb     = UIndex(i,j,k,l,3,3,nmo_);

                //long int aaaa_sym = (ij + 0 * nh) * 4 * nh + (kl + 0 * nh);
                //long int aabb_sym = (ij + 0 * nh) * 4 * nh + (kl + 3 * nh);
                //long int abab_sym = (ij + 1 * nh) * 4 * nh + (kl + 1 * nh);
                //long int baba_sym = (ij + 2 * nh) * 4 * nh + (kl + 2 * nh);
                //long int bbaa_sym = (ij + 3 * nh) * 4 * nh + (kl + 0 * nh);
                //long int bbbb_sym = (ij + 3 * nh) * 4 * nh + (kl + 3 * nh);

// can limit types of excitations here

                if ( type == "AA" ) {

                    long int aaaa_sym = (ij + 0 * nh) * 2 * nh + (kl + 0 * nh);
                    long int aabb_sym = (ij + 0 * nh) * 2 * nh + (kl + 1 * nh);
                    long int bbaa_sym = (ij + 1 * nh) * 2 * nh + (kl + 0 * nh);
                    long int bbbb_sym = (ij + 1 * nh) * 2 * nh + (kl + 1 * nh);

                    newA[aaaa_sym] = A[aaaa];
                    newA[aabb_sym] = A[aabb];
                    newA[bbaa_sym] = A[bbaa];
                    newA[bbbb_sym] = A[bbbb];

                    newB[aaaa_sym] = B[aaaa];
                    newB[aabb_sym] = B[aabb];
                    newB[bbaa_sym] = B[bbaa];
                    newB[bbbb_sym] = B[bbbb];

                }else if (type == "AB" ) {

                    long int abab_sym = ij * nh + kl;

                    newA[abab_sym] = A[abab];

                    newB[abab_sym] = B[abab];

                }else if (type == "BA" ) {

                    long int baba_sym = ij * nh + kl;

                    newA[baba_sym] = A[baba];

                    newB[baba_sym] = B[baba];

                }else {
                    throw PsiException("invalid excitation type",__FILE__,__LINE__);
                }

            }
        }

        // project out null (and negative) space of A matrix
        std::shared_ptr<Matrix> Amat ( new Matrix(totdim,totdim) );
        std::shared_ptr<Matrix> Bmat ( new Matrix(totdim,totdim) );
        for (int ij = 0; ij < totdim; ij++) {
            for (int kl = 0; kl < totdim; kl++) {
                Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                Bmat->pointer()[ij][kl] = newB[ij*totdim+kl];
            }
        }

        // symmetrize
        // A(kl,ji) = A(ij,kl)
        // A(kl,ji) = A(lk,ij)

        std::shared_ptr<Matrix> Bsave ( new Matrix(Amat) );
        std::shared_ptr<Matrix> eigvec ( new Matrix(totdim,totdim) );
        std::shared_ptr<Vector> eigval ( new Vector(totdim) );

        Amat->diagonalize(eigvec,eigval,descending);
        //eigval->print();
        bool prune = true;
        long int newdim = 0;

        // get rid of small/negative eigenvalues of A?
        if ( prune ) {

            newdim = 0;
            int * map = (int*)malloc(totdim*sizeof(int));
            memset((void*)map,'\0',totdim*sizeof(int));
            int * skip = (int*)malloc(totdim*sizeof(int));
            memset((void*)skip,'\0',totdim*sizeof(int));
            for (long int i = 0; i < totdim; i++) {
                double val = eigval->pointer()[i];
                skip[i] = 1;
                if ( val > pruning_threshold ) {
                    map[newdim] = i;
                    skip[i] = 0;
                    newdim++;
                }
            }
            //if ( newdim % 2 == 1 ) { // eigenvalues should come in pairs
            //    newdim--;
            //}
            for (int ij = 0; ij < totdim; ij++) {
                for (int kl = 0; kl < totdim; kl++) {
                    Amat->pointer()[ij][kl] = newA[ij*totdim+kl];
                }
            }
            Amat->transform(eigvec);
            Bmat->transform(eigvec);
        }else {
            newdim = totdim;
        }

        for (long int i = 0; i < newdim; i++) {
            for (long int j = 0; j < newdim; j++) {
                newA[i*newdim+j] = Amat->pointer()[i][j];
                newB[i*newdim+j] = Bmat->pointer()[i][j];
            }
        }

        outfile->Printf("\n");
        if ( newdim < totdim ) {
            outfile->Printf("    <<< warning >>> ");
            outfile->Printf("\n");
            outfile->Printf("    reducing dimension of [ k*l, [H, j*i ] ] from %5i to %5i\n",totdim,newdim);
            outfile->Printf("\n");
        }

        int info = 0;
        if ( prune ) {
            info = SymmetricGeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig);
        }else {
	    bool * energy_is_real = (bool*)malloc(totdim*sizeof(bool));
            GeneralizedEigenvalueProblem(newdim,newA,newB,cc,eig,energy_is_real);
            free(energy_is_real);
        }

        outfile->Printf("    state");
        outfile->Printf("          energy (Eh)");
        outfile->Printf("       ex energy (Eh)");
        outfile->Printf("       ex energy (eV)\n");
        //outfile->Printf("            < n | 0 >\n");
        if ( info != 0 ) {
            outfile->Printf("\n");
            outfile->Printf("    Error: diagonalization failed.  info = %5i\n",info);
            outfile->Printf("\n");
        }else {
            for (long int i = newdim-1; i >= 0; i--) {
                if ( prune ) eig[i] = 1.0 / eig[i];
                if ( eig[i] > 0.0 ) {

                    outfile->Printf("    %5i %20.12lf %20.12lf %20.12lf\n",newdim-i,eig[i]+Process::environment.globals["CURRENT ENERGY"],eig[i],eig[i] * 27.21138);

                }
            }
        }
    }

    free(A);
    //free(B);
    free(newA);
    free(newB);
    free(cc);
    free(eig);

    free(D1a);
    free(D1b);
    free(D2aa);
    free(D2bb);
    free(D2ab);
    //free(D3aaa);
    //free(D3aab);
    //free(D3bba);
    //free(D3bbb);
    //free(D4aaaa);
    //free(D4aaab);
    //free(D4aabb);
    //free(D4bbba);
    //free(D4bbbb);

}

// void ERPASolver::TransitionDensity_SpinAdapted(double *C,  double * D1a, double * D1b, double * TDa, double * TDb, int nh, int h) {
//
//    // Tvo: 
//
//    // T_{pq} = <v| p*q |0>
//    // T_{pq} = <0| [ Ov , p*q ] |0>
//    // T_{pq} = D(j,q)c(p,j) - D(p,i)c(i,q)
//
//
//   // AED: 2/13/18 I'm not sure if the basis for the sum over pq should be the ERPA basis or the ab basis
//   for (int pq = 0; pq < nh; pq++) {
//   
//       int p = bas_erpa_sym[h][pq][0];
//       int q = bas_erpa_sym[h][pq][1];
//   
//       int hp = symmetry[p];
//       int hq = symmetry[q];
//
//       double duma = 0.0;
//       double dumb = 0.0;
//
//       for (int jj = 0; jj < nmopi_[hq]; jj++) {
//
//           int j  = jj + pitzer_offset[hq];
//
//           int pj = ibas_erpa_sym[h][p][j];
//
//           // check if pj is a proper ERPA excitation:
//           if ( pj == -999 ) continue;
//
//           duma += D1a[j*amo_+q] * C[pj];
//           dumb += D1b[j*amo_+q] * C[pj];
//       }
//
//       for (int ii = 0; ii < nmopi_[hp]; ii++) {
//
//           int i = ii + pitzer_offset[hp];
//
//           int iq = ibas_erpa_sym[h][i][q];
//
//           // check if iq is a proper ERPA excitation:
//           if ( iq == -999 ) continue;
//
//           duma -= D1a[p*amo_+i] * C[iq];
//           dumb -= D1b[p*amo_+i] * C[iq];
//       }
//
//       TDa[p*amo_+q] = 1.0/sqrt(2.0) * duma;
//       TDb[p*amo_+q] = 1.0/sqrt(2.0) * dumb;
//   }
//  }
void ERPASolver::TransitionDensity_SpinAdaptedl(double *C,  double * D1a, double * D1b, double * TDla, double * TDlb, int nh, int h) {
 //Build left transition density 
 // T_{pq} = C(j,i)<0| [j*i, p*q ] |0>
 //Tl[pq]=<v|p*q|0>= C(j,p)*D(j,q)-C(q,i)*D(p,i)

    for (int pq = 0; pq < nh; pq++) {
    
        int p = bas_erpa_sym[h][pq][0];
        int q = bas_erpa_sym[h][pq][1];
    
        int hp = symmetry[p];
        int hq = symmetry[q];

        double duma = 0.0;
        double dumb = 0.0;

        for (int jj = 0; jj < nmopi_[hq]; jj++) {

            int j  = jj + pitzer_offset[hq];

            int jp = ibas_erpa_sym[h][j][p];

            // check if pj is a proper ERPA excitation:
            if ( jp == -999 ) continue;

            duma += D1a[j*amo_+q] * C[jp];
            dumb += D1b[j*amo_+q] * C[jp];
        }

        for (int ii = 0; ii < nmopi_[hp]; ii++) {

            int i = ii + pitzer_offset[hp];

            int qi = ibas_erpa_sym[h][q][i];

            // check if iq is a proper ERPA excitation:
            if ( qi == -999 ) continue;

            duma -= D1a[p*amo_+i] * C[qi];
            dumb -= D1b[p*amo_+i] * C[qi];
        }

        TDla[p*amo_+q] = 1.0/sqrt(2.0) * duma;
       // outfile->Printf("\ndumal %20.12lf\n", duma);
        TDlb[p*amo_+q] = 1.0/sqrt(2.0) * dumb;
    }
  }
void ERPASolver::TransitionDensity_SpinAdapted(double *C,  double * D1a, double * D1b, double * TDa, double * TDb, int nh, int h) {
  //Right transtion density
  // T_{pq} = <0| [ p*q, j*i ] |0>*C(i,j)
  // T_{pq} = -D(j,q)c(p,j) + D(p,i)c(i,q)


    for (int pq = 0; pq < nh; pq++) {
    
        int p = bas_erpa_sym[h][pq][0];
        int q = bas_erpa_sym[h][pq][1];
    
        int hp = symmetry[p];
        int hq = symmetry[q];

        double duma = 0.0;
        double dumb = 0.0;

        for (int jj = 0; jj < nmopi_[hq]; jj++) {

            int j  = jj + pitzer_offset[hq];

            int pj = ibas_erpa_sym[h][p][j];

            // check if pj is a proper ERPA excitation:
            if ( pj == -999 ) continue;

            duma -= D1a[j*amo_+q] * C[pj];
            dumb -= D1b[j*amo_+q] * C[pj];
        }

        for (int ii = 0; ii < nmopi_[hp]; ii++) {

            int i = ii + pitzer_offset[hp];

            int iq = ibas_erpa_sym[h][i][q];

            // check if iq is a proper ERPA excitation:
            if ( iq == -999 ) continue;

            duma += D1a[p*amo_+i] * C[iq];
            dumb += D1b[p*amo_+i] * C[iq];
        }

        TDa[p*amo_+q] = 1.0/sqrt(2.0) * duma;
       // outfile->Printf("\nduma %20.12lf\n", duma);
        TDb[p*amo_+q] = 1.0/sqrt(2.0) * dumb;
    }
}
                                                      


void ERPASolver::TransitionDensity(double * C, double * D1a, double * D1b, double * TDa, double * TDb, int nh, int h) {

    // Tvo: 

    // T_{pq} = <v| p*q |0>
    // T_{pq} = <0| [ Ov , p*q ] |0>
    // T_{pq} = D(j,q)c(p,j) - D(p,i)c(i,q)


    // AED: 2/13/18 I'm not sure if the basis for the sum over pq should be the ERPA basis or the ab basis
    for (int pq = 0; pq < nh; pq++) {
    
        int p = bas_erpa_sym[h][pq][0];
        int q = bas_erpa_sym[h][pq][1];
    
        int hp = symmetry[p];
        int hq = symmetry[q];

        double duma = 0.0;
        double dumb = 0.0;

        for (int jj = 0; jj < nmopi_[hq]; jj++) {

            int j  = jj + pitzer_offset[hq];

            int pj = ibas_erpa_sym[h][p][j];

            // check if pj is a proper ERPA excitation:
            if ( pj == -999 ) continue;

            duma += D1a[j*amo_+q] * C[pj     ];
            dumb += D1b[j*amo_+q] * C[pj + nh];
        }

        for (int ii = 0; ii < nmopi_[hp]; ii++) {

            int i = ii + pitzer_offset[hp];

            int iq = ibas_erpa_sym[h][i][q];

            // check if iq is a proper ERPA excitation:
            if ( iq == -999 ) continue;

            duma -= D1a[p*amo_+i] * C[iq     ];
            dumb -= D1b[p*amo_+i] * C[iq + nh];
        }

        TDa[p*amo_+q] = duma;
        TDb[p*amo_+q] = dumb;
    }

}

double ERPASolver::TransitionDipoleMoment_SpinAdapted(double * C, double * D1a, double * D1b, int nh, int h,std::string component) {

    double dum = 0.0;
    
    // D1 term
    for (int ij = 0; ij < nh; ij++) {
    
        int i = bas_erpa_sym[h][ij][0];
        int j = bas_erpa_sym[h][ij][1];
    
        int hi = symmetry[i];
        int hj = symmetry[j];

        double **mua;
        double **mub;

        if ( component == "X" ) {

            mua = mux_moa_->pointer(hi);
            mub = mux_mob_->pointer(hi);

        }else if ( component == "Y" ) {

            mua = muy_moa_->pointer(hi);
            mub = muy_mob_->pointer(hi);

        }else if ( component == "Z" ) {

            mua = muz_moa_->pointer(hi);
            mub = muz_mob_->pointer(hi);

        }else {

            throw PsiException("ERROR: cartesian component not recognized",__FILE__,__LINE__);

        }

        int jj = j - pitzer_offset[hj];
        int ii = i - pitzer_offset[hi];
        for( int KK = 0; KK < amopi_[hi]; KK++){
            int K   = KK + pitzer_offset[hi];
            int Kj  = ibas_erpa_sym[h][K][j];
            dum    += C[ij] * mua[KK][jj] * D1a[K*amo_+i];
            dum    += C[ij] * mub[KK][jj] * D1b[K*amo_+i];
        }
        for( int LL = 0; LL < amopi_[hj]; LL++){
            int L  = LL + pitzer_offset[hj];
            int iL = ibas_erpa_sym[h][i][L];
            dum   -= C[ij] * mua[ii][LL] * D1a[j*amo_+L];
            dum   -= C[ij] * mub[ii][LL] * D1b[j*amo_+L];
        }
    }


    return 1.0/sqrt(2.0) * dum;
}

double ERPASolver::TransitionDipoleMoment(double * C, double * D1a, double * D1b, int nh, int h,std::string component) {

    double dum = 0.0;
    
    // D1 term
    for (int ij = 0; ij < nh; ij++) {
    
        int i = bas_erpa_sym[h][ij][0];
        int j = bas_erpa_sym[h][ij][1];
    
        int hi = symmetry[i];
        int hj = symmetry[j];

        double **mua;
        double **mub;

        if ( component == "X" ) {

            mua = mux_moa_->pointer(hi);
            mub = mux_mob_->pointer(hi);

        }else if ( component == "Y" ) {

            mua = muy_moa_->pointer(hi);
            mub = muy_mob_->pointer(hi);

        }else if ( component == "Z" ) {

            mua = muz_moa_->pointer(hi);
            mub = muz_mob_->pointer(hi);

        }else {

            throw PsiException("ERROR: cartesian component not recognized",__FILE__,__LINE__);

        }

        int jj = j - pitzer_offset[hj];
        int ii = i - pitzer_offset[hi];
        for( int KK = 0; KK < amopi_[hi]; KK++){
            int K   = KK + pitzer_offset[hi];
            int Kj  = ibas_erpa_sym[h][K][j];
            dum    += C[ij      ] * mua[KK][jj] * D1a[K*amo_+i];
            dum    += C[ij + nh ] * mub[KK][jj] * D1b[K*amo_+i];
        }
        for( int LL = 0; LL < amopi_[hj]; LL++){
            int L  = LL + pitzer_offset[hj];
            int iL = ibas_erpa_sym[h][i][L];
            dum   -= C[ij      ] * mua[ii][LL] * D1a[j*amo_+L];
            dum   -= C[ij + nh ] * mub[ii][LL] * D1b[j*amo_+L];
        }
    }


    return dum;
}

/*void ERPASolver::ReconstructD3(double * D3aaa, double * D3aab, double * D3bba, double * D3bbb,
                                double * D2aa, double * D2bb, double * D2ab, double * D1a, double * D1b) {

    long int n1 = nmo_;
    long int n2 = n1 * nmo_;
    long int n3 = n2 * nmo_;
    long int n4 = n3 * nmo_;
    long int n5 = n4 * nmo_;

    // (3 D2 - D1^2)^D1
    
    // 3 D2 - D1^2
    double * taa = (double*)malloc(n4*sizeof(double));
    double * tbb = (double*)malloc(n4*sizeof(double));
    double * tab = (double*)malloc(n4*sizeof(double));

    memset((void*)taa,'\0',n4*sizeof(double));
    memset((void*)tbb,'\0',n4*sizeof(double));
    memset((void*)tab,'\0',n4*sizeof(double));

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {
                                               
                    taa[i*n3 + j*n2 + k*n1 + l] = 3.0 * D2aa[i*n3 + j*n2 + k*n1 + l] - 2.0 * ( D1a[i*n1 + k] * D1a[j*n1 + l] - D1a[j*n1 + k] * D1a[i*n1 + l] );
                                                
                    tbb[i*n3 + j*n2 + k*n1 + l] = 3.0 * D2bb[i*n3 + j*n2 + k*n1 + l] - 2.0 * ( D1b[i*n1 + k] * D1b[j*n1 + l] - D1b[j*n1 + k] * D1b[i*n1 + l] );
                                               
                    tab[i*n3 + j*n2 + k*n1 + l] = 3.0 * D2ab[i*n3 + j*n2 + k*n1 + l] - 2.0 * D1a[i*n1 + k] * D1b[j*n1 + l];
                                               
                }
            }
        }
    }

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            for (int k = 0; k < nmo_; k++) {
                for (int l = 0; l < nmo_; l++) {
                    for (int m = 0; m < nmo_; m++) {
                        for (int n = 0; n < nmo_; n++) {

                            double dumaaa = + taa[i*n3 + j*n2 + l*n1 + m] * D1a[k*n1 + n]  // +ijk,lmn
                                            - taa[i*n3 + k*n2 + l*n1 + m] * D1a[j*n1 + n]  // -ikj,lmn
                                            - taa[k*n3 + j*n2 + l*n1 + m] * D1a[i*n1 + n]  // -kji,lmn

                                            - taa[i*n3 + j*n2 + l*n1 + n] * D1a[k*n1 + m]  // -ijk,lnm
                                            + taa[i*n3 + k*n2 + l*n1 + n] * D1a[j*n1 + m]  // +ikj,lnm
                                            + taa[k*n3 + j*n2 + l*n1 + n] * D1a[i*n1 + m]  // +kji,lnm

                                            - taa[i*n3 + j*n2 + n*n1 + m] * D1a[k*n1 + l]  // -ijk,nml
                                            + taa[i*n3 + k*n2 + n*n1 + m] * D1a[j*n1 + l]  // +ikj,nml
                                            + taa[k*n3 + j*n2 + n*n1 + m] * D1a[i*n1 + l]; // +kji,nml

                            D3aaa[i*n5 + j*n4 + k*n3 + l*n2 + m*n1 + n] = dumaaa;

                            double dumbbb = + tbb[i*n3 + j*n2 + l*n1 + m] * D1b[k*n1 + n]  // +ijk,lmn
                                            - tbb[i*n3 + k*n2 + l*n1 + m] * D1b[j*n1 + n]  // -ikj,lmn
                                            - tbb[k*n3 + j*n2 + l*n1 + m] * D1b[i*n1 + n]  // -kji,lmn

                                            - tbb[i*n3 + j*n2 + l*n1 + n] * D1b[k*n1 + m]  // -ijk,lnm
                                            + tbb[i*n3 + k*n2 + l*n1 + n] * D1b[j*n1 + m]  // +ikj,lnm
                                            + tbb[k*n3 + j*n2 + l*n1 + n] * D1b[i*n1 + m]  // +kji,lnm

                                            - tbb[i*n3 + j*n2 + n*n1 + m] * D1b[k*n1 + l]  // -ijk,nml
                                            + tbb[i*n3 + k*n2 + n*n1 + m] * D1b[j*n1 + l]  // +ikj,nml
                                            + tbb[k*n3 + j*n2 + n*n1 + m] * D1b[i*n1 + l]; // +kji,nml

                            D3bbb[i*n5 + j*n4 + k*n3 + l*n2 + m*n1 + n] = dumbbb;

                            double dumaab = + taa[i*n3 + j*n2 + l*n1 + m] * D1b[k*n1 + n]  // +ijk,lmn

                                            + tab[i*n3 + k*n2 + l*n1 + n] * D1a[j*n1 + m]  // +ikj,lnm
                                            - tab[j*n3 + k*n2 + l*n1 + n] * D1a[i*n1 + m]  // -jki,lnm

                                            - tab[i*n3 + k*n2 + m*n1 + n] * D1a[j*n1 + l]  // +ikj,nml
                                            + tab[j*n3 + k*n2 + m*n1 + n] * D1a[i*n1 + l]; // +kji,nml

                            D3aab[i*n5 + j*n4 + k*n3 + l*n2 + m*n1 + n] = dumaab;

                            double dumbba = + tbb[i*n3 + j*n2 + l*n1 + m] * D1a[k*n1 + n]  // +ijk,lmn

                                            + tab[k*n3 + i*n2 + n*n1 + l] * D1b[j*n1 + m]  // +ikj,lnm
                                            - tab[k*n3 + j*n2 + n*n1 + l] * D1b[i*n1 + m]  // -jki,lnm

                                            - tab[k*n3 + i*n2 + n*n1 + m] * D1b[j*n1 + l]  // +ikj,nml
                                            + tab[k*n3 + j*n2 + n*n1 + m] * D1b[i*n1 + l]; // +kji,nml

                            D3bba[i*n5 + j*n4 + k*n3 + l*n2 + m*n1 + n] = dumbba;
                        }
                    }
                }
            }
        }
    }

    free(taa);
    free(tbb);
    free(tab);

}*/



}} // End namespaces


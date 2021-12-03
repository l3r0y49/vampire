//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) A.Christison 2021. All rights reserved.
//
//   Email: ac2071@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"

namespace molecular_dynamics{

   namespace internal{
         
      void md_main(){
         
         //read_in_params()
         
         //set up data structures
         
         //generate coords
         crystal_generate();
         //derive potentals
         define_potental_tables();
         //inital printout
         inital_printout();
         //equaliberation run
         evolve_sample(N_equi_steps);
         //production run
         evolve_sample(N_prod_steps);
         //end run
         terminate();
      }       
   }

} // end of molecular_dynamics namespace


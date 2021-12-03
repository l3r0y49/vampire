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

// moleculardynamics module headers
#include "internal.hpp"

namespace moleculardynamics{

   //----------------------------------------------------------------------------
   // Function to initialize molecular_dynamics module
   //----------------------------------------------------------------------------
   void initialize(){
         
         //call routine to read in system params and define necessary arrays
         //set up data structures
         
         //read_in_params()
      
         // call md driver
         main_md();
         
      return;

   }

} // end of molecular_dynamics namespace


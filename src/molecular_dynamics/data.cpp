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
#include <cmath>
// Vampire headers
#include "molecular_dynamics.hpp"

// moleculardynamics module headers
#include "internal.hpp"

namespace molecular_dynamics{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside molecular_dynamics module
      //------------------------------------------------------------------------
      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties
      
      vel_acc=false;
      N=0;
      
      //resize box_size now that dimensions has a value
      box_size.resize(dimensions);
      
      temperature_sum = 0.0;
      energy_kinetic_sum = 0.0;
      energy_potental_sum = 0.0;
      pressure_sum = 0.0;
      

      
   } // end of internal namespace

} // end of molecular_dynamics namespace


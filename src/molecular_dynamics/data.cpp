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
namespace mdi=molecular_dynamics::internal;

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
      
      bool vel_acc=false;
      int N=0;
      
      //resize box_size now that dimensions has a value
      std::vector<double> box_size(dimensions);
      
      double temperature_sum = 0.0;
      double energy_kinetic_sum = 0.0;
      double energy_potental_sum = 0.0;
      double pressure_sum = 0.0;

   } // end of internal namespace

} // end of molecular_dynamics namespace


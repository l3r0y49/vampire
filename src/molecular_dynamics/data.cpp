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
      
      const int dimensions=3; //switches between 3D & 2D
      bool vel_acc=false;
      int n=0;
      
      //resize box_size now that dimensions has a value
      box_size.resize(dimensions)
      
      const double r_cuttoff = 2.5;   
      const double phi_cuttoff = 4.0/pow(r_cuttoff,12) - 4.0/pow(r_cuttoff,6); 
      const int table_size = 2001;
      const double r_min = 0.5;
      const double r_sq_min = pow(r_min,2);
      const double delta_r_sq = (pow(r_cuttoff,2)-r_sq_min) /(table_size - 1);
      const double inv_delat_r_sq = 1.0/delta_r_sq;
      
      double temperature_sum = 0.0;
      double energy_kinetic_sum = 0.0;
      double energy_potental_sum = 0.0;
      double pressure_sum = 0.0;
      
      const int max_pairs_per_atom = 100;
      
   } // end of internal namespace

} // end of molecular_dynamics namespace


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

// molecular_dynamics module headers
#include "internal.hpp"

namespace molecular_dynamics{
   
   namespace internal{
      //----------------------------------------------------------------------------
      // Function to generate tables for Lennard-Jones potental
      //----------------------------------------------------------------------------
      void define_potental_tables(){
         //temp routine variables
         double r_sq;
         double rm_2;
         double rm_6;
         double rm_12;
         
         for(i=0,i<table_size,i++){ //may need internal::, whatch for i ob1 error
            
            r_sq = r_sq_min + i * delta_r_sq;
            rm_2 = 1.0/r_sq;     //1/r^2
            rm_6 = pow(rm_2,3)   //1/r^6
            rm_12 = pow(rm_6,3)  //1/r^12
            
            phi_tab(i+1) = 4.0*(rm_12-rm_6)-phi_cuttoff  //4(1/r^12 - 1/r^6)-phi(Rc)
            
            //dphi = -(1/r)(dv/dr)
            d_phi_tab(i+1) = 24.0*rm_2*(2.0*rm_12-rm_6)  //24(1/r^14 - 1/r^8)
         }
         return;
      }
      
   }
} // end of molecular_dynamics namespace


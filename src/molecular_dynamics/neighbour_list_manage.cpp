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
      //----------------------------------------------------------------------------
      // Function to refold particles that left the box back into the box due to the periodic boundary conditions
      //----------------------------------------------------------------------------
      void refold_positions(){
         //loop over array element wise
         //c++ arrays go across then down
         int i;
         int j;
      
         for(i=0;i<positions.size();i++){
         
            for(j=0;j<positions.size();j++){
            
               //refold if out of box
               if(positions[j][i]>0.5){
                  positions[j][i]=positions[j][i]-1.0
               }
            
               //refold if out of box
               if(positions[j][i]<-0.5){
                  positions[j][i]=positions[j][i]+1.0
               }
            }
         }

         return;
      }
      }

} // end of molecular_dynamics namespace


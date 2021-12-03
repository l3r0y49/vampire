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

         // call md driver
         main_md();
         
      return;

   }
   
   void inital_printout(){
      
      printf("# Number of steps: %i time step: %f total time %f \n",Nsteps,deltat,double(Nsteps)*deltat);
      printf("# Number of atoms: %i \n",N);
      printf("# Box size: %f %f %f volume: %f \n",box_size[0],box_size[1],box_size[1],volume);
      if(change_rho){
         printf("# Density: %f (changed)\n",density);
      }else{
         printf("# Density: %f (unchanged)\n",density);
      }
      if(t_constat){
         printf("# Constant T run with T =: %f \n",t_constat);
      }else{
         printf("# Free evolution run\n");
      }
      printf("# skin: %f maximum neighbor list length: %i \n",skin,max_list_length);
   }

} // end of molecular_dynamics namespace


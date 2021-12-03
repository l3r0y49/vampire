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
   
   void initalize_values(){
      
   }
   
   void read_input(){
      //sample values
      //hardwired values to be read in from file once file I/O implemented
      /*
      **implement reminder in file I/0**
      cut_off_LJ = 2.5;
       */
      atomic_lattice=1.5496;
      nx=4;
      ny=4;
      nz=4;
      dispalc=0.1;

      //simulation values
      //hardwired values to be read in from file once file I/O implemented
       N_equi_steps=10000;
       N_prod_steps=pow(2.0,13); //exponent of 13-18 known to work from MD module
       delta_t=0.004;
       rho_requested=0.7;
       t_requested=1.5;
       skin=0.1;
       
       if(rho_requested>0.0){
         change_rho = true;
       }else{
          change_rho = false;
       }
       if(t_requested>0.0){
         t_constat = true;
       }else{
          t_constat = false;
       }
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
      
      printf("#\n");
      printf("# Step   Temperature     Kinetic      Potential   Total Energy    Pressure\n");
      printf("# -----  ------------  ------------  ------------  ------------  ------------\n");
   }

} // end of molecular_dynamics namespace


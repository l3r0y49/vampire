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
      std::vector<double>pos_at_real(dimensions), mass_center(dimensions);
      double scale,x_sum,y_sum,z_sum;
      int i;
      
      if(N <=0 ){
         printf("FATAL ERROR: N is %i \n",N);
      }
      
      //compute volume and density, do not change in run
      volume = box_size[0]*box_size[1]*box_size[2];
      density = N/volume;
      
      //if user changes density, do it here
      if(change_rho){
         scale = pow((density/rho_requested),(1.0/dimensions));
         box_size = scale*box_size;
         density = N / volume;
      }
      
      //can now allocate arrays with atomic info
      dispalcement.resize(dimensions,std::vector<double>(N));
      positions.resize(dimensions,std::vector<double>(N));
      velocities.resize(dimensions,std::vector<double>(N));
      accelerations.resize(dimensions,std::vector<double>(N));
      energy_potental.resize(N);
      energy_kinetic.resize(N);
      
      //and neighbor list
      max_list_length = max_pairs_per_atom*N;
      list.resize(max_list_length);
      advance.resize(N);
      marker_1.resize(N);
      marker_2.resize(N);
      dispalcement_list.resize(dimensions,std::vector<double>(N));
      
      //positions written to array in gen_coords file
      //don't need to manage velocities and accelerations as they are all now kept internal
      
      //compute center of mass coordianates - porobably a more afficent way to do all this
      for(i=0;i<N;i++){
         mass_center[0]=positions[i][0];
         mass_center[1]=positions[i][1];
         mass_center[2]=positions[i][2];
      }
      mass_center[0]/=N;
      mass_center[1]/=N;
      mass_center[2]/=N;
      //translate atoms to center of mass is at origin
      for(i=0;i<N;i++){
         positions[i][0]-=mass_center[0];
         positions[i][1]-=mass_center[1];
         positions[i][2]-=mass_center[2];
      }
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


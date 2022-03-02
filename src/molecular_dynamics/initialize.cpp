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
#include <cstdio>
// Vampire headers
#include "molecular_dynamics.hpp"

// moleculardynamics module headers
#include "internal.hpp"
namespace mdi=molecular_dynamics::internal;

namespace moleculardynamics{
   
   namespace internal{

      //----------------------------------------------------------------------------
      // Function to initialize molecular_dynamics module
      //----------------------------------------------------------------------------
      void initialize(){
      
         // call md driver
         main_md();
         return;
      }
   
   void initalize_values(){
      std::vector<double>pos_at_real(mdi::dimensions), mass_center(mdi::dimensions);
      double scale,x_sum,y_sum,z_sum;
      int i;
      
      if(mdi::N <=0 ){
         printf("FATAL ERROR: N is %i \n",mdi::N);
      }
      
      //compute volume and density, do not change in run
      mdi::volume = mdi::box_size[0]*mdi::box_size[1]*mdi::box_size[2];
      mdi::density = mdi::N/mdi::volume;
      
      //if user changes density, do it here
      if(mdi::change_rho){
         scale = pow((mdi::density/mdi::rho_requested),(1.0/mdi::dimensions));
         mdi::box_size = scale*mdi::box_size;
         mdi::density = mdi::N / mdi::volume;
      }
      
      //can now allocate arrays with atomic info
      mdi::dispalcement.resize(dimensions,std::vector<double>(mdi::N));
      mdi::positions.resize(dimensions,std::vector<double>(mdi::N));
      mdi::velocities.resize(dimensions,std::vector<double>(mdi::N));
      mdi::accelerations.resize(dimensions,std::vector<double>(mdi::N));
      mdi::energy_potental.resize(mdi::N);
      mdi::energy_kinetic.resize(mdi::N);
      
      //and neighbor list
      mdi::max_list_length = mdi::max_pairs_per_atom*mdi::N;
      mdi::list.resize(mdi::max_list_length);
      mdi::advance.resize(mdi::N);
      mdi::marker_1.resize(mdi::N);
      mdi::marker_2.resize(mdi::N);
      mdi::dispalcement_list.resize(dimensions,std::vector<double>(mdi::N));
      
      //positions written to array in gen_coords file
      //don't need to manage velocities and accelerations as they are all now kept internal
      
      //compute center of mass coordianates - porobably a more afficent way to do all this
      for(i=0;i<mdi::N;i++){
         mass_center[0]=mdi::positions[i][0];
         mass_center[1]=mdi::positions[i][1];
         mass_center[2]=mdi::positions[i][2];
      }
      mass_center[0]/=double(mdi::N);
      mass_center[1]/=double(mdi::N);
      mass_center[2]/=double(mdi::N);
      //translate atoms to center of mass is at origin
      for(i=0;i<mdi::N;i++){
         mdi::positions[i][0]-=mass_center[0];
         mdi::positions[i][1]-=mass_center[1];
         mdi::positions[i][2]-=mass_center[2];
      }
   }
   
   void read_input(){
      //sample values
      //hardwired values to be read in from file once file I/O implemented
      /*
      **implement reminder in file I/0**
      cut_off_LJ = 2.5;
       */
      mdi::atomic_lattice=1.5496;
      mdi::nx=4;
      mdi::ny=4;
      mdi::nz=4;
      mdi::dispalc=0.1;
      
      mdi::box_size[0]=mdi::nx;
      mdi::box_size[1]=mdi::ny;
      mdi::box_size[2]=mdi::nz;

      //simulation values
      //hardwired values to be read in from file once file I/O implemented
       mdi::N_equi_steps=10000;
       mdi::N_prod_steps=pow(2.0,13); //exponent of 13-18 known to work from MD module
       mdi::delta_t=0.004;
       mdi::rho_requested=0.7;
       mdi::t_requested=1.5;
       mdi::skin=0.1;
       
       if(mdi::rho_requested>0.0){
         mdi::change_rho = true;
       }else{
          mdi::change_rho = false;
       }
       if(mdi::t_requested>0.0){
         mdi::t_constat = true;
       }else{
          mdi::t_constat = false;
       }
   }
   
   void inital_printout(){
      
      printf("# Number of steps: %i time step: %f total time %f \n",mdi::N_prod_steps,mdi::delta_t,double(mdi::N_prod_steps)*mdi::delta_t
);
      printf("# Number of atoms: %i \n",mdi::N);
      printf("# Box size: %f %f %f volume: %f \n",mdi::box_size[0],mdi::box_size[1],mdi::box_size[1],mdi::volume);
      if(mdi::change_rho){
         printf("# Density: %f (changed)\n",mdi::density);
      }else{
         printf("# Density: %f (unchanged)\n",mdi::density);
      }
      if(mdi::t_constat){
         printf("# Constant T run with T =: %f \n",mdi::t_constat);
      }else{
         printf("# Free evolution run\n");
      }
      printf("# skin: %f maximum neighbor list length: %i \n",mdi::skin,mdi::max_list_length
);
      
      printf("#\n");
      printf("# Step   Temperature     Kinetic      Potential   Total Energy    Pressure\n");
      printf("# -----  ------------  ------------  ------------  ------------  ------------\n");
   }
   }

} // end of molecular_dynamics namespace


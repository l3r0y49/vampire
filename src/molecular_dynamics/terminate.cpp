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
namespace mdi=molecular_dynamics::internal;

namespace molecular_dynamics{

   namespace internal{
         // output routines and deallocate vectors if necessary
      
         void terminate(){
            
            //output routines
            print_statistics();
            write_sample(); 
            //- may want to call at step intervals during production stage to get visulisation of atom movements over time
            
            //clear vectors to clean memory
            mdi::box_size.clear();
            mdi::dispalcement.clear();
            mdi::positions.clear();
            mdi::velocities.clear();
            mdi::accelerations.clear();
            mdi::energy_potental.clear();
            mdi::energy_kinetic.clear();
            mdi::advance.clear();
            mdi::marker_1.clear();
            mdi::marker_2.clear();
            mdi::list.clear();
            mdi::dispalcement_list.clear();
            
            return;
         }
         
      void print_statistics(){
         
         if(mdi::N_prod_steps <= 0){ //protect from N=0 situation
            return;
         }
         //print all stats for production phase
         printf("# Means %f %f %f %f %f \n",mdi::temperature_sum/mdi::N_prod_steps, mdi::energy_kinetic_sum/mdi::N_prod_steps, mdi::energy_potental_sum/mdi::N_prod_steps, (mdi::energy_kinetic_sum+mdi::energy_potental_sum)/mdi::N_prod_steps, mdi::pressure_sum/mdi::N_prod_steps);
         
         return;
      }
      
      void write_sample(){
         int i;
         vel_acc=true;
         //print out all real space positions
         
         //use file I/O later to do so
         
         for(i=0;i<mdi::N;i++){
            printf("%f %f %f\n",mdi::positions[i][0]*mdi::box_size[0],mdi::positions[i][1]*mdi::box_size[0],mdi::positions[i][2]*mdi::box_size[0]);
         }
         printf("\n");
         for(i=0;i<mdi::N;i++){
            printf("%f %f %f\n",mdi::velocities[i][0]*mdi::box_size[0],mdi::velocities[i][1]*mdi::box_size[0],mdi::velocities[i][2]*mdi::box_size[0]);
         }
         printf("\n");
         for(i=0;i<mdi::N;i++){
            printf("%f %f %f\n",mdi::accelerations[i][0]*mdi::box_size[0],mdi::accelerations[i][1]*mdi::box_size[0],mdi::accelerations[i][2]*mdi::box_size[0]);
         }
         printf("\n");
         
         
         return;
      }
         
   }

} // end of molecular_dynamics namespace


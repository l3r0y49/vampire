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

namespace molecular_dynamics{

   namespace internal{
         // output routines and deallocate vectors if necessary
      
         void terminate(){
            //output routines
            
            print_statistcs();
            
            //write_sample(); 
            //- may want to call at step intervals during production stage to get visulisation of atom movements over time
            
            //clear vectors to clean memory
            box_size.clear();
            dispalcement.clear();
            positions.clear();
            velocities.clear();
            accelerations.clear();
            energy_potental.clear();
            energy_kinetic.clear();
            advance.clear();
            marker_1.clear();
            marker_2.clear();
            list.clear();
            dispalcement_list.clear();
            
            return;
         }
         
      void print_statistics(){
         
         if(N_steps <= 0){
            return;
         }
         printf("# Means %f %f %f %f %f \n",temperature_sum/Nsteps, ene_kin_sum/Nsteps, ene_pot_sum/ Nsteps, (ene_kin_sum+ ene_pot_sum)/Nsteps, pressure_sum / Nsteps)
         
         return;
      }
      
      void write_sample(){
         int i;
         vel_acc=true;
         
         for(i=0;i<N;i++){
            printf("%f %f %f\n",positions[i][0]*box_size,positions[i][1]*box_size,positions[i][2]*box_size);
         }
         printf("\n")
         for(i=0;i<N;i++){
            printf("%f %f %f\n",velocities[i][0]*box_size,velocities[i][1]*box_size,velocities[i][2]*box_size);
         }
         printf("\n")
         for(i=0;i<N;i++){
            printf("%f %f %f\n",accelerations[i][0]*box_size,accelerations[i][1]*box_size,accelerations[i][2]*box_size);
         }
         printf("\n")
         
         
         return;
      }
         
   }

} // end of molecular_dynamics namespace


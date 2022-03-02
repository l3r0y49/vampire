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
#include <numeric>
#include <vector>

// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"
namespace mdi=molecular_dynamics::internal;

namespace molecular_dynamics{
   
      namespace internal{

         //----------------------------------------------------------------------------
         // Function to assess the two largest displacements of the system, if their sum is larger than skin returns true
         //----------------------------------------------------------------------------
         moved_too_much(skin){
            double skin;
            double displ;
            double displ1=0.0;
            double displ2=0.0;
            ///std::vector <double> temp_dispalcement_list(3);
            int i;
      
            for(i=0;i<mdi::N;i++){
               //set size to 0 and then to required to 
               ///temp_dispalcement_list.resize(0);
               ///temp_dispalcement_list.resize(3);
               //temp_dispalcement_list.resize(dispalcement_list[0].size());
               
               ///populate_1d_with_column_int(temp_dispalcement_list,dispalcement_list,i);
               
//                displ=std::sqrt(std::inner_product(temp_dispalcement_list.begin(), temp_dispalcement_list.end(), temp_dispalcement_list.begin(),0));
               displ=std::sqrt(std::inner_product(dispalcement_list[i].begin(), dispalcement_list[i].end(), dispalcement_list[i].begin(),0));
               
               if(displ>=displ1){
                  displ2=displ1;
                  displ1=displ;
                  
               }else if(displ>=displ2){
                  displ2=displ;
               }
            }
            
            if(displ1+displ2>mdi::skin){
               return true;
            }else{
               return false;
            }
         }
         
         //----------------------------------------------------------------------------
         // Main routine that controls the time evolution of the system
         //----------------------------------------------------------------------------
         
         void evolve_sample(int N_steps){
            int step,i;
            double ene_kin_aver,ene_pot_aver,ene_tot_aver,temperature,pressure,chi;
            bool list_update_requested = true;
            std::vector<double> dispalcement_sqr
            
            compute_temperature(ene_kin_aver,temperature);
            
            //"Velocity Verlet" Integrator (see e.g. Allen and Tildesley book, p. 81)
            // velocity scaling applied when t_constat is enabled
            
            for(step=1;step<=N_steps;step++){
               
               refold_positions();
               //genarate square of displacements
               
               //errors================================================== need for loop
               std::transform(mdi::dispalcement.begin(),mdi::dispalcement.end(),dispalcement_sqr,multiplies<double>());
               mdi::dispalcement =mdi::dispalcement*mdi::velocities + 0.5*(mdi::dispalcement)*mdi::accelerations;   //dr = r(t+dt) -r (t)
               
               mdi::positions += mdi::dispalcement;  // r(t+dt)
               
               //============================================================= - for loop required
               
               if(mdi::t_constat && (temperature>0)){  //velocity rescale for constant temp.
                  compute_temperature(ene_kin_aver,temperature);
                  chi = sqrt(mdi::t_requested/temperature);
                  mdi::velocities = chi*mdi::velocities + 0.5*mdi::dispalcement*mdi::accelerations; //v(t+dt/2) (scaled)
               } else {
                  mdi::velocities = mdi::velocities + 0.5*mdi::dispalcement*mdi::accelerations; //v(t+dt/2)
               }
               if(list_update_requested){       //update required
                  update_list(mdi::r_cuttoff+mdi::skin);  //do update
                  list_update_requested=false;  //updated
               }
               compute_forces(); //a(t+dt)
               //=================for loop required ============
               mdi::velocities=mdi::velocities + 0.5*mdi::dispalcement*mdi::accelerations;    //v(t+dt)
               compute_temperature(ene_kin_aver,temperature);              // at t+dt, energy_kinetic
               ene_kin_aver = std::accumulate(mdi::energy_potental.begin(), mdi::energy_potental.end(),0)/mdi::N;
               ene_tot_aver = ene_kin_aver + ene_pot_aver;
               
               //pressure calculation, see the Allen and Tildesley book, section 2.4
               
               pressure = mdi::density*temperature + mdi::virial/mdi::volume;
               
               //update displacement list
               for(i=0;i<mdi::N;i++){
                  mdi::dispalcement_list[i] = mdi::dispalcement_list[j] + mdi::box_size*mdi::displacement[j];
               }
               //deterioration test, if moved too much relative to skin
               //list update scheduled for next step
               list_update_requested=moved_too_much(skin);
               
               printf("%i %f %f %f %f \n",step,temperature,ene_kin_aver,ene_pot_aver,pressure); //**make wtite out to file once I/O routines are set up**
               
               //accumlate stats
               mdi::temperature_sum += temperature;
               mdi::energy_kinetic_sum += ene_kin_aver;
               mdi::energy_potental_sum += ene_pot_aver;
               mdi::pressure_sum += pressure;
            }
            
            return;
         }
   
         //----------------------------------------------------------------------------
         // Routine to populate an 1d vector of doubles with the given column of a 2d vector of doubles (of the same size)
         //----------------------------------------------------------------------------
   
//          void populate_1d_with_column_doubles(std::vector<int>& one_d_vector,std::vector<std::vector<int> >& two_d_vector,int index){
//       
//             int i;
//             //copy the column of the 2d vector into the row of the 1d vector
//             for(i=0;i<two_d_vector[index].size();i++){
//                one_d_vector[i] = two_d_vector[i][index];
//             }
//       
//             return;
//          }
//       }
   
   } // end of molecular_dynamics namespace
}

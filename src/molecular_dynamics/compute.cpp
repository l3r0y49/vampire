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
#include <vector>
#include <numeric>   // std::inner_product
// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"
namespace mdi=molecular_dynamics::internal;

namespace molecular_dynamics{
   
   namespace internal{
      //----------------------------------------------------------------------------
      // Function to generate tables for Lennard-Jones potental
      //----------------------------------------------------------------------------
      void define_potental_tables(){
         //temp routine variables
         double r_sq,rm_2,rm_6,rm_12;
         int i;
         
         for(i=0;i<table_size;i++){ //may need internal::, watch for i ob1 error
            
            r_sq = r_sq_min + i * delta_r_sq;
            rm_2 = 1.0/r_sq;     //1/r^2
            rm_6 = pow(rm_2,3);   //1/r^6
            rm_12 = pow(rm_6,3);  //1/r^12
            
            mdi::phi_tab[i+1] = 4.0*(rm_12-rm_6)-phi_cuttoff;  //4(1/r^12 - 1/r^6)-phi(Rc)
            
            //d_phi = -(1/r)(dv/dr)
            mdi::d_phi_tab[i+1] = 24.0*rm_2*(2.0*rm_12-rm_6);  //24(1/r^14 - 1/r^8)
         }
         return;
      }
      //----------------------------------------------------------------------------
      // Function to calculate forces between atoms
      //----------------------------------------------------------------------------
      void compute_forces(){
         std::vector<double> sij(dimensions),rij(dimensions);
//          std::vector<double> temp_1d_pos_i(3),temp_1d_pos_j(3);
         double r_sqij,phi,d_phi,rk,weight;
         int i,j,k,l,m;
         
         //reset force, potental energy and virial terms
         mdi::accelerations.resize(0,std::vector<double>(0));
         mdi::accelerations.resize(dimensions,std::vector<double>(N));
         mdi::energy_potental.resize(0);
         mdi::energy_potental.resize(N);
         mdi::virial=0.0;
         
         for(i=0;i<mdi::N;i++){
            //useful part of neighbour list
            for(l=mdi::marker_1[i];l<=mdi::marker_2[i];l++){
               j = mdi::list[l];
               
               //reset 1d position vector containers values without deallocating memory
//                temp_1d_pos_1.resize(0);
//                temp_1d_pos_2.resize(0);
//                temp_1d_pos_1.resize(3);
//                temp_1d_pos_2.resize(3);
               //extract i and j atom positions
               
//                populate_1d_with_column_doubles(temp_1d_pos_i,positions,i);
//                populate_1d_with_column_doubles(temp_1d_pos_j,positions,j);
               //distance between i and j
//                sij = temp_1d_pos_i - temp_1d_pos_j;
               sij = mdi::positions[i] - mdi::positions[j];
               
               //apply boundary conditions where needed
               for(m=0;m<sij.size();m++){
                  if(sij[m]>0.5){
                     sij[m]-=1.0 ;
                  }else if(sij[m]<0.5){
                     sij[m]+=1.0;
                  }
               }
               
               //==============errors=================== - need to mutiply 2d vectors -- for loop??
               rij=mdi::box_size[0]*sij;  //real space units
               r_sqij=std::inner_product(rij.begin(),rij.end(),rij.begin(),0);  //square distance
               //========================================
               
               if(r_sqij < pow(mdi::r_cuttoff,2)){   //are particles interacting?
                  rk = (r_sqij - r_sq_min)* inv_delat_r_sq +1.0; // "continuous index in tables"
                  k= int(rk);    //descrete index
                  if(k<1){
                     k = 1;      //to protect
                  }
                  weight = rk - double(k);            //fractional part [0,1]
                  phi = weight* mdi::phi_tab[k+1] + (1.0-weight)* mdi::phi_tab[k]; //linear interpolation
                  d_phi = weight* mdi::d_phi_tab[k+1] + (1.0-weight)* mdi::d_phi_tab[k];
                  mdi::energy_potental[i] += 0.5*phi;     //accumulate energy
                  mdi::energy_potental[j] += 0.5*phi;     //shared bewteen i&j
                  mdi::virial -= d_phi*r_sqij;             //accum. virial sum r(dv/dr)
                  mdi::acc[i] += d_phi*sij;                //accum. forces (Fij = -Fji)
                  mdi::acc[j] -= d_phi*sij; 
               }
            }
         }
      mdi::virial = -mdi::virial/mdi::dimensions;                    //virial term
      }
      
      void compute_temperature(double energy_kin_aver,double temperature){
         std::vector<double> real_vel(mdi::dimensions);
         int i;
         
         for(i=0;i<mdi::N;i++){ //ob1 error?
            
            //for rool required need sum of all velcitites
            
            //vector operation
            real_vel = mdi::box_size*mdi::velocities[i];

            mdi::energy_kinetic[i] = 0.5*std::inner_product(real_vel.begin(),real_vel.end(),real_vel.begin(),0);
         }
         energy_kin_aver = std::accumulate(mdi::energy_kinetic.begin(),mdi::energy_kinetic.end(),0)/mdi::N;
         temperature = 2.0*energy_kin_aver/mdi::dimensions;
      }
      
   }
} // end of molecular_dynamics namespace


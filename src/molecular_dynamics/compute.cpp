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
         double r_sq,rm_2,rm_6,rm_12;
         int i;
         
         for(i=0;i<table_size;i++){ //may need internal::, watch for i ob1 error
            
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
      //----------------------------------------------------------------------------
      // Function to calculate forces between atoms
      //----------------------------------------------------------------------------
      void compute_forces(){
         std::vector<double> sij(dimensions),rij(dimensions);
//          std::vector<double> temp_1d_pos_i(3),temp_1d_pos_j(3);
         double r_sqij,phi,d_phi,rk,weight;
         int i,j,k,l,m;
         
         //reset force, potental energy and virial terms
         accelerations.resize(0,std::vector<double>(0));
         accelerations.resize(dimensions,std::vector<double>(n));
         energy_potental.resize(0);
         energy_potental.resize(n);
         virial=0.0;
         
         for(i=0;i<n;i++){
            //useful part of neighbour list
            for(l=marker_1(i);l<=marker_2(i);l++){
               j = list(l);
               
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
               sij = positions[i] - positions[j];
               
               //apply boundary conditions where needed
               for(m=0;m<sij;m++){
                  if(sij(m)>0.5){
                     sij(m)-=1.0 ;
                  }else if(sij(m)<0.5){
                     sij(m)+=1.0;
                  }
               }
               
               rij=box_size*sij;  //real space units
               r_sqij=std::inner_product(rij,rij);  //square distance
               
               if(r_sqij < pow(rcutoff,2)){   //are particles interacting?
                  rk = (r_sqij - r_sq_min)* inv_delat_r_sq +1.0; // "continuous index in tables"
                  k= int(rk);    //descrete index
                  if(k<1){
                     k = 1;      //to protect
                  }
                  weight = rk - double(k);            //fractional part [0,1]
                  phi = weight* phi_tab(k+1) + (1.0-weight)* phi_tab(k); //linear interpolation
                  dphi = weight* d_phi_tab(k+1) + (1.0-weight)* d_phi_tab(k);
                  energy_potental(i) += 0.5*phi;     //accumulate energy
                  energy_potental(j) += 0.5*phi;     //shared bewteen i&j
                  virial -= dphi*r_sqij;             //accum. virial sum r(dv/dr)
                  acc[i] += d_phi*sij                //accum. forces (Fij = -Fji)
                  acc[j] -= d_phi*sij 
               }
            }
         }
      viral = -viral/dimensions                    //virial term
      }
      
   }
} // end of molecular_dynamics namespace


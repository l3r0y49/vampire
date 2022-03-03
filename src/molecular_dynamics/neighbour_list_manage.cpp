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
#include <cstdio>
#include <cmath>
#include <algorithm>
// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"
namespace mdi=molecular_dynamics::internal;

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
      
         for(i=0;i<mdi::positions.size();i++){
            for(j=0;j<mdi::positions.size();j++){
               //refold if out of box
               if(mdi::positions[j][i]>0.5){
                  mdi::positions[j][i]-=1.0;
               }
               //refold if out of box
               if(mdi::positions[j][i]<-0.5){
                  mdi::positions[j][i]+=1.0;
               }
            }
         }

         return;
      }
      //----------------------------------------------------------------------------
      // Update the neighbour list, include all atom paris that are
      // within range of each other
      
      // The Fincham-Ralston method is used (D. Fincham and B.J. Ralston,
      // Comp.Phys.Comm. 23, 127 (1981) ), which allows vectorization of the
      // first of the two inner loop on vector computers and therefore gives
      // good performance.
      //----------------------------------------------------------------------------
      void update_list(double range){
         double range_sq;
         std::vector<double> sij(dimensions),rij(dimensions);
         double r_sqij;
         int i,j,l,m;
         
         printf("Neighbour list update \n");
         range_sq = pow(range,2);
         
         //Fincham_Ralston loop for list updates
         l=0;
         for(i=0;i<N;i++){
            for(j=i+1;j<N;j++){
             sij = mdi::positions[i] - mdi::positions[j];
             
             //apply boundary conditions where needed
               for(m=0;m<sij.size();m++){
                  if(sij[m]>0.5){
                     sij[m]-=1.0 ;
                  }else if(sij[m]<0.5){
                        sij[m]+=1.0;
                  }
               }
               
               rij=std::transform( mdi::box_size.begin(), mdi::box_size.end(),sij.begin(),std::multiplies<double>() );
               
               
               r_sqij=std::inner_product(rij.begin(),rij.end(),rij.begin(),0);  //square distance
               
               if(r_sqij<range_sq){ //is j an neighbour of i?
                  mdi::advance[j]=1.0;    //yes
               }else{
                  mdi::advance[j]=0.0;     //no
               }
            }
            mdi::marker_1[i]=l;          //start list for i
            for(j=i+1;j<mdi::N;j++){
               if(l>mdi::max_list_length){
                  printf("update_list: FATAL: list too small for skin \n");
                  printf("%i value of parameter max_pairs_per_atom needs increasing\n",max_pairs_per_atom);
                  std::exit;
               }
               mdi::list[l]=j;     //j included
               l+=mdi::advance[j]; //only if advance(j) is 1 
            }
            mdi::marker_2[i]=l-1;  //end of list for l
         }
      mdi::list_length= l-1;    //final lenght of list
               
      printf("%i index in list \n",mdi::list_length);
      
      mdi::dispalcement_list.resize(0,std::vector<double>(0));
      mdi::dispalcement_list.resize(dimensions,std::vector<double>(N));
      
      return;
      }
   }

} // end of molecular_dynamics namespace


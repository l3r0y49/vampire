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
// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"

namespace molecular_dynamics{
   
      namespace internal{

         //----------------------------------------------------------------------------
         // Function to assess the two largest displacements of the system, if their sum is larger than skin returns true
         //----------------------------------------------------------------------------
         bool moved_too_much(skin){
            double skin;
            double displ;
            double displ1=0.0;
            double displ2=0.0;
            std::vector <double> temp_dispalcement_list;
            int i;
      
            for(i=0;i<N;i++){
               temp_dispalcement_list.clear();
               temp_dispalcement_list.resize(dispalcement_list[0].size());
               
               populate_1d_with_column_int(temp_dispalcement_list,dispalcement_list,i);
               
               displ=std::sqrt(std::inner_product(temp_dispalcement_list.begin(), temp_dispalcement_list.end(), temp_dispalcement_list.begin(),0));
               
               if(displ>=displ1){
                  displ2=displ1;
                  displ1=displ;
                  
               }else if(displ>=displ2){
                  displ2=displ;
               }
            }
            
            if(displ1+displ2>skin){
               return true;
            }else{
               return false;
            }
         }
   
         //----------------------------------------------------------------------------
         // Routine to populate an 1d vector of doubles with the given column of a 2d vector of doubles (of the same size)
         //----------------------------------------------------------------------------
   
         void populate_1d_with_column_doubles(std::vector<int>& one_d_vector,std::vector<std::vector<int> >& two_d_vector,int index){
      
            int i;
            //copy the column of the 2d vector into the row of the 1d vector
            for(i=0;i<two_d_vector[index].size();i++){
               one_d_vector[i] = two_d_vector[i][index];
            }
      
            return;
         }
      }
   
} // end of molecular_dynamics namespace


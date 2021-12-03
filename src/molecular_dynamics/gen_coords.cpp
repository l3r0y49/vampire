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
#include <random>
// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"

namespace molecular_dynamics{

   namespace internal{
         
      void generate_crystal(){
         //double atomic_lattice,int nx,int ny,int nz,double dispalc
         const int crtout = 6;
         std::vector<double>rands(3);
         int i,j,k,l,n_count;
         double x,y,z;
         const int nbase = 4; //atoms in fcc cell
         std::vector<std::vector<double> > rcell = {{0.0,0.0,0.0}, //coords of atom
                                                    {0.5,0.5,0.0},
                                                    {0.0,0.0,0.5},
                                                    {0.5,0.0,0.5}}
         // no file I/O needed at this stage
         
         //define no. of particles and box size 
         N=4*nx*ny*nz                  //internal N
         box_size[1]=nx*atomic_lattice //internal box_size
         box_size[2]=ny*atomic_lattice //internal box_size
         box_size[3]=nx*atomic_lattice //internal box_size
         
         N_count=0; //number of atom positions written into array

         for(k=0,k<nz,k++){  //OB1??
            for(j=0,j<ny,j++){ //OB1??
               for(i=0,i<nx,i++){ //OB1??
                  for(l=0,l<nbase,l++){ //OB1??
                     
                     rands = {random(),random(),random()},
                     x = atomic_lattice*(i + rcell(1,L)) + 2.0*dispalc(rands[1]-0.5)
                     y = atomic_lattice*(j + rcell(2,L)) + 2.0*dispalc(rands[2]-0.5)
                     z = atomic_lattice*(k + rcell(3,L)) + 2.0*dispalc(rands[3]-0.5)
            
                     //write direct to displacemnet array
                     displacement[N_count] = {x,y,z};
                     N_count++;
                  }
               }
            }
         }
      }
      
      //random number template function
      /*credit - https://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range, 03/12/21*/
      template<typename T>
      T random() {
         std::random_device                  rand_dev;
         std::mt19937                        generator(rand_dev());
         std::uniform_int_distribution<T>    distr(0.0, 1.0);
         return distr(generator);
      }
         
   }

} // end of molecular_dynamics namespace


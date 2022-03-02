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
#include <cmath>
#include <ctime>
// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"
namespace mdi=molecular_dynamics::internal;

namespace molecular_dynamics{

   namespace internal{
      
      class rng{
         
         std::mt19937 mt;  //mersenne twister
         std::uniform_real_distribution<double> dist;
         
      public:
         
         //seed rng
         
         void seed(unsigned int random_seed){
            dist = std::uniform_real_distribution<double>(0.0,1.0);
            std::mt19937::result_type mt_seed = random_seed;
            mt.seed(mt_seed); //seed gen
         }
         
         //wrapper function
         double grnd(){
            return dist(mt);
            
         }
      };
         
      void generate_crystal(){
         //double atomic_lattice,int nx,int ny,int nz,double dispalc - use mdi:: namespace variables
         const int crtout = 6;
//          std::vector<double>rands(3);
         int i,j,k,l,N_count;
         double x,y,z;
         const int nbase = 4; //atoms in fcc cell
         std::vector<std::vector<double> > rcell = {{0.0,0.0,0.0}, //coords of atom
                                                    {0.5,0.5,0.0},
                                                    {0.0,0.0,0.5},
                                                    {0.5,0.0,0.5}};
         // no file I/O needed at this stage
         
         std::vector<rng> random_generators(3);
         
         //define no. of particles and box size 
         mdi::N=4*mdi::nx*mdi::ny*mdi::nz;                  //internal N
         mdi::box_size[0]=mdi::nx*mdi::atomic_lattice; //internal box_size
         mdi::box_size[1]=mdi::ny*mdi::atomic_lattice; //internal box_size
         mdi::box_size[2]=mdi::nx*mdi::atomic_lattice; //internal box_size
         
         N_count=0; //number of atom positions written into array

         for(k=0;k<mdi::nz;k++){  //OB1??
            for(j=0;j<mdi::ny;j++){ //OB1??
               for(i=0;i<mdi::nx;i++){ //OB1??
                  for(l=0;l<nbase;l++){ //OB1??
                     
                     //seed each rng
                     random_generators[0].seed(l*i);
                     random_generators[1].seed(l*j);
                     random_generators[2].seed(l*k);
                     
                     //generate array with 3 rn
//                      rands = {random_generators[0].grnd(),random_generators[1].grnd(),random_generators[2].grnd()};
                     
                     x = mdi::atomic_lattice*(i + rcell[1,L]) + 2.0*mdi::dispalc[random_generators[0].grnd()-0.5];
                     y = mdi::atomic_lattice*(j + rcell[2,L]) + 2.0*mdi::dispalc[random_generators[1].grnd()-0.5];
                     z = mdi::atomic_lattice*(k + rcell[3,L]) + 2.0*mdi::dispalc[random_generators[2].grnd()-0.5];
            
                     //write direct to positions array
                     mdi::positions[N_count] = {x/mdi::box_size[0],y/mdi::box_size[2],z/mdi::box_size[2]};
                     //reset elements to 0 in other atomistic arrays
                     mdi::energy_potental[N_count]=0.0;
                     mdi::energy_kinetic[N_count]=0.0;
                     mdi::velocities[N_count]={0.0,0.0,0.0};
                     mdi::accelerations[N_count]={0.0,0.0,0.0};
                     mdi::N_count++;
                  }
               }
            }
         }
      }
      
      //random number template function
      /*credit - https://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range, 03/12/21*/
//       template<typename T>
//       T random() {
//          std::random_device                  rand_dev;
//          std::mt19937                        generator(rand_dev());
//          std::uniform_int_distribution<T>    distr(0.0, 1.0);
//          return distr(generator);
//       }
         
   }

} // end of molecular_dynamics namespace


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

#ifndef MOLECULAR_DYNAMICS_INTERNAL_H_
#define MOLECULAR_DYNAMICS_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the molecular_dynamics module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>
#include <numeric>
#include <math.h> 
#include <random>
// Vampire headers
#include "molecular_dynamics.hpp"

// molecular_dynamics module headers
#include "internal.hpp"

namespace molecular_dynamics{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------
       
       //particle properties
       
       //select 2D or 3D simulation
       const int dimensions=3; //switches between 3D & 2D
       //velocities and acclerations known? i.e starting equaliberration or production stages
       extern bool vel_acc;
       extern int N;
       extern std::vector<double> box_size;
       std::vector<std::vector <double> > dispalcement;
       std::vector<std::vector <double> > positions;
       std::vector<std::vector <double> > velocities;
       std::vector<std::vector <double> > accelerations;
       std::vector<double> energy_potental;
       std::vector<double> energy_kinetic;
       
       //double vector array structure 
       
       /*
        array[0] = {a,b,c}
             x    y    z
        0   [a]  [b]  [c]
        1   [d]  [e]  [f]
        2   [g]  [h]  [i]
        .
        .
       N-2  [k]  [l]  [m]
       N-1  [n]  [o]  [p]
        */
       
       extern double volume;
       extern double density;
       extern double virial;   //virial term to compute pressure
       
       //!!!!!!!!!!may require some file name delcrations to store coords at various stages of run, ie in between equil and prod runs!!!!!!!!!!
       //md simulation control properties
       int N_equi_steps;        //steps in the equiliberation stage
       int N_prod_steps;        //steps in the production stage
       double delta_t;          //time steps
       double rho_requested;    //desired density
       double t_requested;      //desired temperature
       bool change_rho;         //true when user is changing density
       bool t_constat;          //becomes true when t_requested >=0
       double skin;             //additional range in neigbour list
       
       //Parameters and tables for the Lennard-Jones energy_potental
       
      const double r_cuttoff = 2.5;   
      const double phi_cuttoff = 4.0/pow(r_cuttoff,12) - 4.0/pow(r_cuttoff,6); 
      const int table_size = 2001;
      const double r_min = 0.5;
      const double r_sq_min = pow(r_min,2);
      const double delta_r_sq = (pow(r_cuttoff,2)-r_sq_min) /(table_size - 1);
      const double inv_delat_r_sq = 1.0/delta_r_sq;
       
       std::vector<double> phi_tab;
       std::vector<double> d_phi_tab;   //phi(r),1/r dphi/dr
       
       //Statistical quantities accumulated during run
       
       extern double temperature_sum;
       extern double energy_kinetic_sum;
       extern double energy_potental_sum;
       extern double pressure_sum;
       
       //Neighbour list properties
        const int max_pairs_per_atom = 100;
       
       std::vector<int> advance;
       std::vector<int> marker_1;
       std::vector<int> marker_2;
       std::vector<int> list;
       extern int max_list_length;
       int list_length;
       std::vector<std::vector <double> > dispalcement_list;
       
       //sample required values
       double atomic_lattice;
       int nx;
       int ny;
       int nz;
       double dispalc;
       const double cut_off_LJ = 2.5;
       const double cut_off_Al = 5.55805441821810;
      
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
       
      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             double test;

             // constructor
             mp_t (const unsigned int max_materials = 100):
                test(0.0) // constructor initialisation of test variable
             {
                // constructor body for initialising more complex data/arrays
             }; // end of constructor

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties
      
      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------
      
      //neighbour_list_manage
      void refold_positions();
      void update_list(double range);
      
      //compute
      void define_potental_tables();
      void compute_forces();
      void compute_temperature(double energy_kin_aver,double temperature);
      
      //evolve
      bool moved_too_much(double skin);
      void evolve_sample(int N_steps);
//       void populate_1d_with_column_doubles(std::vector<int>& one_d_vector,std::vector<std::vector<int> >& two_d_vector,int index);
      
      //main
      void md_main();
      
      //generate_coords
      //double atomic_lattice,int nx,int ny,int nz,double dispalc
      void generate_crystal();
//       template<typename T> T random();
      
      
      //terminate
      void terminate();
      void print_statistics();
      void write_sample();
      
      //initalise
      void inital_printout();
      void read_input();
      void initalize_values();

   } // end of internal namespace

} // end of molecular_dynamics namespace

#endif //MOLECULAR_DYNAMICS_INTERNAL_H_

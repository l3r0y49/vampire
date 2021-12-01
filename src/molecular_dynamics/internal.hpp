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
       const int dimensions;
       //velocities and acclerations known? i.e starting equaliberration or production stages
       bool vel_acc;
       int n;
       std::vector<double> box_size;
       std::vector<std::vector <double> > dispalcement;
       std::vector<std::vector <double> > positions;
       std::vector<std::vector <double> > velocities;
       std::vector<std::vector <double> > accelerations;
       std::vector<double> energy_potental;
       std::vector<double> energy_kinetic;
       
       double volume;
       double density;
       double virial;   //virial term to compute pressure
       
       //!!!!!!!!!!may require some file name delcrations to store coords at various stages of run, ie in between equil and prod runs!!!!!!!!!!
       //md simulation control properties
       int N_equi_steps;        //steps in the equiliberation stage
       int N_prod_steps;        //steps in the production stage
       double delta_t;          //time steps
       double rho_requested;    //desired density
       double t_requested;      //desired temperature
       bool change_roh;         //true when user is changing density
       bool t_constat;          //becomes true when t_requested >=0
       double skin;             //additional range in neigbour list
       
       //Parameters and tables for the Lennard-Jones energy_potental
       
       const double r_cuttoff;   //cuttof distance
       const double phi_cuttoff; //potental at cuttoff
       const int table_size;
       const double r_min;
       const double r_sq_min;
       const double delta_r_sq;
       const double inv_delat_r_sq;
       
       std::vector<double> phi_tab;
       std::vector<double> d_phi_tab;   //phi(r),1/r dphi/dr
       
       //Statistical quantities accumulated during run
       
       double temperature_sum;
       double energy_kinetic_sum;
       double energy_potental_sum;
       double pressure_sum;
       
       //Neighbour list properties
       const int max_pairs_per_atom;
       std::vector<double> advance;
       std::vector<double> marker_1;
       std::vector<double> marker_2;
       std::vector<double> list;
       int max_list_length;
       int list_length;
       std::vector<std::vector <double> > dispalcement_list;
       
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
      
      //compute
      void define_potental_tables();
      
      //evolve
      bool moved_too_much(skin);
      void populate_1d_with_column_doubles(std::vector<int>& one_d_vector,std::vector<std::vector<int> >& two_d_vector,int index);

   } // end of internal namespace

} // end of molecular_dynamics namespace

#endif //MOLECULAR_DYNAMICS_INTERNAL_H_

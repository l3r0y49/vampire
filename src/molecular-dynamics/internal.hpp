//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Alex_Christison 2021. All rights reserved.
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

   } // end of internal namespace

} // end of molecular_dynamics namespace

#endif //MOLECULAR_DYNAMICS_INTERNAL_H_

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

#ifndef MOLECULARDYNAMICS_INTERNAL_H_
#define MOLECULARDYNAMICS_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the moleculardynamics module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "moleculardynamics.hpp"

// moleculardynamics module headers
#include "internal.hpp"

namespace moleculardynamics{

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

} // end of moleculardynamics namespace

#endif //MOLECULARDYNAMICS_INTERNAL_H_

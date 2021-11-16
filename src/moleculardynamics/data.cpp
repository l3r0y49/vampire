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

// Vampire headers
#include "moleculardynamics.hpp"

// moleculardynamics module headers
#include "internal.hpp"

namespace moleculardynamics{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside moleculardynamics module
      //------------------------------------------------------------------------

      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

   } // end of internal namespace

} // end of moleculardynamics namespace


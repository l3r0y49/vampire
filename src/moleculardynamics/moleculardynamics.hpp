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

#ifndef MOLECULARDYNAMICS_H_
#define MOLECULARDYNAMICS_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "moleculardynamics.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for moleculardynamics module
//--------------------------------------------------------------------------------
namespace moleculardynamics{

   //-----------------------------------------------------------------------------
   // Function to initialise moleculardynamics module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for moleculardynamics module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of moleculardynamics namespace

#endif //MOLECULARDYNAMICS_H_

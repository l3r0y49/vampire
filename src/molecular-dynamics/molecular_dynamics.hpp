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

#ifndef MOLECULAR_DYNAMICS_H_
#define MOLECULAR_DYNAMICS_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "molecular_dynamics.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for molecular_dynamics module
//--------------------------------------------------------------------------------
namespace molecular_dynamics{

   //-----------------------------------------------------------------------------
   // Function to initialise molecular_dynamics module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for molecular_dynamics module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of molecular_dynamics namespace

#endif //MOLECULAR_DYNAMICS_H_

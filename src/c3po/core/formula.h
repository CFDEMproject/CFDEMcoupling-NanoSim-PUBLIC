/*-----------------------------------------------------------------------------*\
                  ___   _____   _____   _____     ___   
                 / ___\/\  __`\/\  __`\/\  __`\  / __`\ 
                /\ \__/\ \ \_\ \ \ \_\ \ \ \_\ \/\ \_\ \
                \ \____\\ \  __/\ \  __/\ \  __/\ \____/
                 \/____/ \ \ \/  \ \ \/  \ \ \/  \/___/ 
                          \ \_\   \ \_\   \ \_\         
                           \/_/    \/_/    \/_/         

         A Compilation for Fluid-Particle Data Post PrOcessing

Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
               2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria
---------------------------------------------------------------------------------
License
    CPPPO is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with CPPPO. If not, see <http://www.gnu.org/licenses/lgpl.html>.

	This code is designed for on-the-fly post processing of fluid-particle
	data (e.g., of velocity, pressure, concentration, temperature field).

	Parts of the code were developed in the frame of the NanoSim project funded
	by the European Commission through FP7 Grant agreement no. 604656.
\*-----------------------------------------------------------------------------
    Description
    tool for interpretation and execution of user-defined formulas
\*-----------------------------------------------------------------------------*/
#ifndef C3PO_FORMULA_H
#define C3PO_FORMULA_H

#include "stdio.h"
#include <string>
#include <vector>
#include "mpi.h"

namespace C3PO_NS
{
 class Formula
 {
  public:
  
  Formula(const char* formula);
  ~Formula();

  void            interpretFormula(int, int) const; //interprets the strings in the formula
  std::string     getFormula() const; //just for display
  std::string     getNumerator() const; //just for display
  std::string     getDenominator() const; //just for display
  void            evaluate(std::vector< std::vector<double>  > *  inputFields_, std::vector< double > * outputField_ );
  
  private:
  
  mutable std::string raw_formula_;
  mutable std::vector<char> formula_;
  mutable std::vector<char> numerator_; 
  mutable std::vector<char> denominator_;   //by default empty
  
  void ParseFormula() const;
  void CheckFormula() const;
  
  void throw_error(std::string function,std::string msg_) const;
  
 };
}
#endif

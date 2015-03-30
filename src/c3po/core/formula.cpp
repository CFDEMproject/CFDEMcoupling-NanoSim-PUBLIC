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

#include "formula.h"
#include "stdlib.h"
#include <cstdio>

using namespace C3PO_NS;

Formula::Formula(const char* formula)
:
raw_formula_(formula)
{
 ParseFormula();
 CheckFormula();
}
/*-----------------------------------------------------------------------------*/
Formula::~Formula()
{
}
/*-----------------------------------------------------------------------------*/
void Formula::ParseFormula() const
{
 for (unsigned int it=0; it<raw_formula_.size();it++)
  {
        
        if (raw_formula_[it] != ' ' )
        {  
          formula_.push_back(raw_formula_[it]);
         
        }
  }
}
/*-----------------------------------------------------------------------------*/
void Formula::CheckFormula() const
{
 bool operator_before=true;

 for (unsigned int it=0; it<formula_.size();it++)
  {
   if(operator_before )    
   {
 
    operator_before=false;
   
    if((formula_[it]-'0')<0 || (formula_[it]-'0')>9)
    { 
   
      printf("\nFormula::CheckFormula -> ERROR: '%c' is not a valid entry for CPPPO formulas OR the formula is not correct!\n Valid entries for position %i are:\n  a number between 0 and 9\n",formula_[it],it);
    
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    
    }
   }
   else if(!operator_before)
   {
    if(formula_[it]!='*' && formula_[it]!='/' && formula_[it]!='+' && formula_[it]!='-'  )
    {
     printf("\nFormula::CheckFormula -> ERROR: '%c' is not a valid entry for CPPPO formulas OR the formula is not correct!\n Valid entries for position %i are:\n '+'\n '-'\n '*'\n '/'\n",formula_[it],it);
    
     MPI_Barrier(MPI_COMM_WORLD);
     MPI_Finalize();
     exit(1);
    }
   
    operator_before=true;
    
   }   
  }
  

}

/*-----------------------------------------------------------------------------*/
std::string Formula::getFormula()
{
 std::string out_;
 
 for(unsigned int i=0;i<formula_.size();i++)
  out_.push_back(formula_[i]);
  
 return out_;
 
}
/*-----------------------------------------------------------------------------*/
void Formula::evaluate(std::vector< std::vector<double>  > * inputFields_, std::vector< double > * outputField_ )
{
 int size_=(*inputFields_)[0].size();

 //start calculation
 for(int i=0;i<size_;i++)
 {
  double value_=(*inputFields_)[(formula_[0]-'0')][i];
  
  for(int n=1;n<int(formula_.size());n++)
  {
   if(formula_[n]=='+') value_+=(*inputFields_)[(formula_[n+1]-'0')][i];
   if(formula_[n]=='-') value_-=(*inputFields_)[(formula_[n+1]-'0')][i];
   if(formula_[n]=='*') value_= value_ * (*inputFields_)[(formula_[n+1]-'0')][i];
   if(formula_[n]=='/') value_= value_ / ( (*inputFields_)[(formula_[n+1]-'0')][i] + 1e-15);
   n++;
  }
 
 outputField_->push_back(value_); 
 }
}


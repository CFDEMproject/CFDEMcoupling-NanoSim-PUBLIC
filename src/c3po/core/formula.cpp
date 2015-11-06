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
#include "mpi.h"

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
void Formula::throw_error(std::string function ,std::string msg_) const
{
 std::cout << "\nFormula::"<< function << " -> ERROR: " << msg_;
 MPI_Barrier(MPI_COMM_WORLD);
 //MPI_Finalize();
 exit(1);
}
/*-----------------------------------------------------------------------------*/
void Formula::ParseFormula() const
{
  bool foundNumerator_   = false;
  bool foundDenominator_ = false;
  for (unsigned int it=0; it<raw_formula_.size();it++)
  {
        if (raw_formula_[it] == '(' && !foundDenominator_ )
            foundNumerator_ = true;

        if (raw_formula_[it] == '%')
        {
            if(foundDenominator_ == true) throw_error("ParseFormula","Only one numerator/denominator separator is allowed!");
            foundNumerator_   = false;
            foundDenominator_ = true;
        }

        if ( raw_formula_[it] != ' ' )
        {  
          formula_.push_back(raw_formula_[it]);

          if(     foundNumerator_ 
               && !(raw_formula_[it]=='(') 
               && !(raw_formula_[it]==')')  
            )
            numerator_.push_back(raw_formula_[it]);

          if(      foundDenominator_ 
                && !(raw_formula_[it]=='(') 
                && !(raw_formula_[it]==')')
                && !(raw_formula_[it]=='%')  
            )
            denominator_.push_back(raw_formula_[it]);

        }
  }

  //If user did not specify brackets, fill numerator
  if(numerator_.size()==0)
      numerator_.assign(formula_.begin(),formula_.end());


   std::cout << "formula: " << (getFormula()).c_str() << " \n";

}
/*-----------------------------------------------------------------------------*/
void Formula::interpretFormula(int numVecs, int numScalars) const
{
    //Replace words starting with 'vec' and 'scalar' with ids
    for(unsigned int i=(numerator_.size()-1);i>0;i--) //start at the back
    {
        if( numerator_[i]=='c' )
        {
         if(i>1)
         {
          if( numerator_[i-1]=='e' && numerator_[i-2]=='v')
          {
           if(i==numerator_.size()-1) throw_error("interpretFormula","vector field number missing!");
           
           numerator_[i+1]=numerator_[i+1]-1; //for vectors just decrease the index
           
           if(numerator_[i+1]-'0'>9) throw_error("interpretFormula","vector field number missing!");
           
            
           numerator_.erase(numerator_.begin()+i-2,numerator_.begin()+i+1);
          }
         }
         else
         {
          throw_error("interpretFormula","vector field not correctly defined in your formula numerator!");
         }
        } 
        
       
        if( numerator_[i]=='r' )
        {
         if(i>4)
         {
          if( numerator_[i-1]=='a' && numerator_[i-2]=='l'  && numerator_[i-3]=='a' && numerator_[i-4]=='c' && numerator_[i-5]=='s' )
          {
            if(i==numerator_.size()-1) throw_error("interpretFormula","scalar field number missing!");
            numerator_[i+1]=numerator_[i+1]-1+numVecs;
            if(numerator_[i+1]-'0'>9) throw_error("interpretFormula","scalar field number missing!"); 
            numerator_.erase(numerator_.begin()+i-5,numerator_.begin()+i+1);
          }
         }
         else
         {
          throw_error("interpretFormula","scalar field not correctly defined in your formula numerator!");      
         }
        }
        
        if(   ( numerator_[i]=='v' && numerator_[i+2]!='c')  
            ||( numerator_[i]=='e' && numerator_[i+1]!='c') 
            ||( numerator_[i]=='s' && numerator_[i+5]!='r') 
            ||( numerator_[i]=='a' && numerator_[i+3]!='r' && numerator_[i+1]!='r')  
            ||( numerator_[i]=='l' && numerator_[i+2]!='r') 
            ||( numerator_[i]=='c' && numerator_[i+4]!='r' && numerator_[i-1]!='e') 
          )
         throw_error("interpretFormula","fields are not correctly defined in your formula numerator! Use 'vec' or 'scalar' to identify them.");
         
    }
    std::cout << "interpreted numerator: " << (getNumerator()).c_str() << " \n";
    

    //Replace words starting with 'vec' and 'scalar' with ids
    if(denominator_.size()>0)
    for(unsigned int i=(denominator_.size()-1);i>0;i--) //start at the back
    {
      
       if( denominator_[i]=='c' )
       {
        if(i>1)
        {
        if( denominator_[i-1]=='e' && denominator_[i-2]=='v' )
         {
            if(i==denominator_.size()-1) throw_error("interpretFormula","vector field number missing!");
            
            denominator_[i+1]=denominator_[i+1]-1; //for vectors just decrease the index
            
             if(denominator_[i+1]-'0'>9) throw_error("interpretFormula","vector field number missing!");
            
            denominator_.erase(denominator_.begin()+i-2,denominator_.begin()+i+1);
         }
         else
         {
          throw_error("interpretFormula","scalar field not correctly defined in your formula denominator!");      
         }
        }
       }
       
       
       if( denominator_[i]=='r' )
       {
        if(i>4)
        {
         if( denominator_[i-1]=='a' && denominator_[i-2]=='l'  && denominator_[i-3]=='a' && denominator_[i-4]=='c' && denominator_[i-5]=='s' )
         {
            if(i==denominator_.size()-1) throw_error("interpretFormula","scalar field number missing!");
            denominator_[i+1]=denominator_[i+1]-1+numVecs; //for 
             if(denominator_[i+1]-'0'>9) throw_error("interpretFormula","scalar field number missing!");
            denominator_.erase(denominator_.begin()+i-5,denominator_.begin()+i+1);
         }
        }
        else
        {
          throw_error("interpretFormula","scalar field not correctly defined in your formula denominator!");      
        }
        
       }
        
        
        if(   ( denominator_[i]=='v' && denominator_[i+2]!='c')  
            ||( denominator_[i]=='e' && denominator_[i+1]!='c') 
            ||( denominator_[i]=='s' && denominator_[i+5]!='r') 
            ||( denominator_[i]=='a' && denominator_[i+3]!='r' && denominator_[i+1]!='r')  
            ||( denominator_[i]=='l' && denominator_[i+2]!='r') 
            ||( denominator_[i]=='c' && denominator_[i+4]!='r' && denominator_[i-1]!='e') 
          )
         throw_error("interpretFormula","fields are not correctly defined in your formula denominator! Use 'vec' or 'scalar' to identify them.");
    }

    std::cout << "interpreted denominator: " << (getDenominator()).c_str() << " \n";
}

/*-----------------------------------------------------------------------------*/
void Formula::CheckFormula() const
{
 bool operator_before=true;

 for (unsigned int it=0; it<formula_.size();it++)
  {
   if(formula_[it]=='(' || formula_[it]==')' 
      || formula_[it]=='v'|| formula_[it]=='e'|| formula_[it]=='c'
      || formula_[it]=='s'|| formula_[it]=='a'|| formula_[it]=='l'|| formula_[it]=='r'
     )
       continue;

   if(formula_[it]=='%')
   {
       operator_before=true;
       continue;
   }

   if(operator_before )    
   {
    operator_before=false;
   }
   else
   {
    if(formula_[it]!='*' && formula_[it]!='/' && formula_[it]!='+' && formula_[it]!='-'  )
    {
    
     char buf[512];
     sprintf(buf,"\nFormula::CheckFormula -> ERROR: '%c' is not a valid entry for CPPPO formulas OR the formula is not correct (i.e., does not contain 'scalar' or 'vec' quantities!\n Valid entries for position %i are:\n '+'\n '-'\n '*'\n '/'\n",formula_[it],it);
     std::string msg_(buf);

     throw_error("CheckFormula",msg_); 
    
    }
   
    operator_before=true;
    
   }   
  }
  

}

/*-----------------------------------------------------------------------------*/
std::string Formula::getFormula() const
{
 std::string out_;
 for(unsigned int i=0;i<formula_.size();i++)
  out_.push_back(formula_[i]);
  
 return out_;
}


/*-----------------------------------------------------------------------------*/
std::string Formula::getNumerator() const
{
 std::string out_;
 for(unsigned int i=0;i<numerator_.size();i++)
  out_.push_back(numerator_[i]);
  
 return out_;
}


/*-----------------------------------------------------------------------------*/
std::string Formula::getDenominator() const
{
 std::string out_;
 for(unsigned int i=0;i<denominator_.size();i++)
  out_.push_back(denominator_[i]);
  
 return out_;
}
/*-----------------------------------------------------------------------------*/
void Formula::evaluate(std::vector< std::vector<double>  > * inputFields_, std::vector< double > * outputField_ )
{
 int size_=(*inputFields_)[0].size();

 //start calculation
 for(int i=0;i<size_;i++)
 {
  //Compute numerator
  double valueNumerator_=(*inputFields_)[(numerator_[0]-'0')][i];
  for(int n=1;n<int(numerator_.size());n++)
  {
   if(numerator_[n]=='+') valueNumerator_+= (*inputFields_)[(numerator_[n+1]-'0')][i];
   if(numerator_[n]=='-') valueNumerator_-= (*inputFields_)[(numerator_[n+1]-'0')][i];
   if(numerator_[n]=='*') valueNumerator_*= (*inputFields_)[(numerator_[n+1]-'0')][i];
   if(numerator_[n]=='/') valueNumerator_/= (*inputFields_)[(numerator_[n+1]-'0')][i] + 1e-15;
   n++; //Advance another value because of operator
  }

  //Compute denominator
  double valueDenominator_ = 1.0;

  if(denominator_.size()>0)
  {
    valueDenominator_=(*inputFields_)[(denominator_[0]-'0')][i];
    for(int n=1;n<int(denominator_.size());n++)
    {
     if(denominator_[n]=='+') valueDenominator_+= (*inputFields_)[(denominator_[n+1]-'0')][i];
     if(denominator_[n]=='-') valueDenominator_-= (*inputFields_)[(denominator_[n+1]-'0')][i];
     if(denominator_[n]=='*') valueDenominator_*= (*inputFields_)[(denominator_[n+1]-'0')][i];
     if(denominator_[n]=='/') valueDenominator_/= (*inputFields_)[(denominator_[n+1]-'0')][i] + 1e-15;
     n++; //Advance another value because of operator
    }
  }
  
  outputField_->push_back(valueNumerator_ / valueDenominator_); 
 }
}


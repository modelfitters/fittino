/* $Id: ProspinoCalculator.cpp 2775 2016-11-15 12:41:43Z sarrazin $ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        ProspinoCalculator.cpp                                           *
*                                                                              *
* Description Wrapper around MadGraph                                          *
* Authors                                                                      *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include "TMath.h"
#include "ProspinoCalculator.h"
#include "ConfigurationException.h"
#include "PhysicsModel.h"
#include "ModelParameter.h"
#include "Factory.h"
#include "SimplePrediction.h"
#include <iostream>
#include <fstream>
#include <CalculatorException.h>
#include <boost/regex.hpp>
#include "CalculatorBase.h"
#include <cstdlib>

Fittino::ProspinoCalculator::ProspinoCalculator( const ModelBase* model, const boost::property_tree::ptree& ptree ) 
  : SLHACalculatorBase( model, &ptree ),
    
    
    
    _ProsOldSLHAFile ( ptree.get<std::string>("ProOldSLHAFile") ),	
    _SLHAFile( ptree.get<std::string>("SLHAFile") ),
    _executable ( ptree.get<std::string>("Executable") ) 

{

    _ProOutput = ptree.get<std::string>("ProOutput");

AddOutput("Pro_XSLO_05m_gg");
AddOutput("Pro_XSLO_05m_bb1");
AddOutput("Pro_XSLO_05m_bb2");
AddOutput("Pro_XSLO_05m_tb1");
AddOutput("Pro_XSLO_05m_tb2");
AddOutput("Pro_XSLO_05m_sg");
AddOutput("Pro_XSLO_05m_sb");
AddOutput("Pro_XSLO_05m_ss");

AddOutput("Pro_XSLO_1m_gg");
AddOutput("Pro_XSLO_1m_bb1");
AddOutput("Pro_XSLO_1m_tb1");
AddOutput("Pro_XSLO_1m_bb2");
AddOutput("Pro_XSLO_1m_tb2");
AddOutput("Pro_XSLO_1m_sg");
AddOutput("Pro_XSLO_1m_sb");
AddOutput("Pro_XSLO_1m_ss");

AddOutput("Pro_XSLO_2m_gg");
AddOutput("Pro_XSLO_2m_bb1");
AddOutput("Pro_XSLO_2m_tb1");
AddOutput("Pro_XSLO_2m_bb2");
AddOutput("Pro_XSLO_2m_tb2");
AddOutput("Pro_XSLO_2m_sg");
AddOutput("Pro_XSLO_2m_sb");
AddOutput("Pro_XSLO_2m_ss");


AddOutput("Pro_XSNLO_05m_gg");
AddOutput("Pro_XSNLO_05m_bb1");
AddOutput("Pro_XSNLO_05m_bb2");
AddOutput("Pro_XSNLO_05m_tb1");
AddOutput("Pro_XSNLO_05m_tb2");
AddOutput("Pro_XSNLO_05m_sg");
AddOutput("Pro_XSNLO_05m_sb");
AddOutput("Pro_XSNLO_05m_ss");

AddOutput("Pro_XSNLO_1m_gg");
AddOutput("Pro_XSNLO_1m_bb1");
AddOutput("Pro_XSNLO_1m_tb1");
AddOutput("Pro_XSNLO_1m_bb2");
AddOutput("Pro_XSNLO_1m_tb2");
AddOutput("Pro_XSNLO_1m_sg");
AddOutput("Pro_XSNLO_1m_sb");
AddOutput("Pro_XSNLO_1m_ss");

AddOutput("Pro_XSNLO_2m_gg");
AddOutput("Pro_XSNLO_2m_bb1");
AddOutput("Pro_XSNLO_2m_tb1");
AddOutput("Pro_XSNLO_2m_bb2");
AddOutput("Pro_XSNLO_2m_tb2");
AddOutput("Pro_XSNLO_2m_sg");
AddOutput("Pro_XSNLO_2m_sb");
AddOutput("Pro_XSNLO_2m_ss");

AddOutput("Pro_ErrNLO_gg");
AddOutput("Pro_ErrNLO_bb1");
AddOutput("Pro_ErrNLO_tb1");
AddOutput("Pro_ErrNLO_bb2");
AddOutput("Pro_ErrNLO_tb2");
AddOutput("Pro_ErrNLO_sg");
AddOutput("Pro_ErrNLO_sb");
AddOutput("Pro_ErrNLO_ss");



}

Fittino::ProspinoCalculator::~ProspinoCalculator() {
  
}

void Fittino::ProspinoCalculator::CalculatePredictions() {
int n;
// removes the old SLHA files
n=0;
std:: string remove;
remove = "rm "+ _ProsOldSLHAFile;
const char *cstr = remove.c_str();
// gives the new created SPheno SLHA file, the name Prospino needs
std:: string copypaste;
copypaste = "cp " + _SLHAFile + " " + _ProsOldSLHAFile ;
const char *cstrr = copypaste.c_str();

// calls prospino
std:: string execute;
execute = "cd " + _executable + " ; ./prospino_2.run" ; 
const char *cstrrr = execute.c_str();


std:: system( cstr );
std:: system ( cstrr );


std:: system( cstrrr );
// reads in the output file

ReadFile();

}
void Fittino::ProspinoCalculator::ReadFile() {
 
        std::string dummy, temp1, temp5, temp6;
	std::string line;
        double temp2 = 0, temp3=0;
	int temp4 =0;
	std::vector<std::string> process;
	std::vector<double> XSLO;
	std::vector<double> XSNLO;
	std::vector<int> SQFam;
	std::ifstream infile(_ProOutput);
	while(std::getline(infile, line)){
         if(line.empty()) break  ;
	    std::istringstream ss(line); 	
	    ss>>temp1>>temp4>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy>>temp2>>dummy>>temp3>>dummy>>dummy>>dummy>>dummy;
	    process.push_back(temp1);
	    XSLO.push_back(temp2);
	    XSNLO.push_back(temp3);
	    SQFam.push_back(temp4);
	    std::cout<<temp4<<std::endl; 

}
       infile.close();

   // for every process the LOSX at the three scale is read in
       std::map<std::string, double> proXSLO05m;
   for (int i=0 ; i<process.size() ; i+=3){
     if(process[i]=="bb"|| process[i]=="tb"){
	 temp5 = std::to_string(SQFam[i]); 	
         temp6 = process[i] + temp5;	
          proXSLO05m[temp6]=XSLO[i];
	}
      else{
	proXSLO05m[process[i]] = XSLO[i];
	} 
    
   }

       std::map<std::string, double> proXSLO1m;
   for (int i=1 ; i<process.size() ; i+=3){
     if(process[i]=="bb"|| process[i]=="tb"){
	 temp5 = std::to_string(SQFam[i]); 	
         temp6 = process[i] + temp5;	
          proXSLO1m[temp6]=XSLO[i];
		}
      else{
 	proXSLO1m[process[i]] = XSLO[i];
	}	 
    
}
       std::map<std::string, double> proXSLO2m;
   for (int i=2 ; i<process.size() ; i+=3){
     if(process[i]=="bb"|| process[i]=="tb"){
	 temp5 = std::to_string(SQFam[i]); 	
         temp6 = process[i] + temp5;	
          proXSLO2m[temp6]=XSLO[i];
		}
      else{
 	proXSLO2m[process[i]] = XSLO[i];
      }	 
    
}

  for( const auto& elem : proXSLO05m ) {

   
    SetOutput( "Pro_XSLO_05m_" + elem.first, elem.second  );
  }
  for( const auto& elem : proXSLO1m ) {

   
    SetOutput( "Pro_XSLO_1m_" + elem.first, elem.second  );
  }
  for( const auto& elem : proXSLO2m ) {

   
    SetOutput( "Pro_XSLO_2m_" + elem.first, elem.second  );
  }



// for every process the NLOXS at the three scale is read in

       std::map<std::string, double> proXSNLO05m;
   for (int i=0 ; i<process.size() ; i+=3){
     if(process[i]=="bb"|| process[i]=="tb"){
	 temp5 = std::to_string(SQFam[i]); 	
         temp6 = process[i] + temp5;	
          proXSNLO05m[temp6]=XSNLO[i];
	}
      else{
	proXSNLO05m[process[i]] = XSNLO[i];
	} 
    
   }

       std::map<std::string, double> proXSNLO1m;
   for (int i=1 ; i<process.size() ; i+=3){
     if(process[i]=="bb"|| process[i]=="tb"){
	 temp5 = std::to_string(SQFam[i]); 	
         temp6 = process[i] + temp5;	
          proXSNLO1m[temp6]=XSNLO[i];
		}
      else{
 	proXSNLO1m[process[i]] = XSNLO[i];
	}	 
    
}
       std::map<std::string, double> proXSNLO2m;
       std::map<std::string, double> proErr1NLO;
       std::map<std::string, double> proErr2NLO;

// NLOXS error is calcualted as the abs(difference of the scale differences in the NLOXS)

   for (int i=2 ; i<process.size() ; i+=3){
     if(process[i]=="bb"|| process[i]=="tb"){
	 temp5 = std::to_string(SQFam[i]); 	
         temp6 = process[i] + temp5;	
         proXSNLO2m[temp6]=XSNLO[i];
       	 proErr1NLO[temp6] = proXSNLO05m[temp6] - proXSNLO1m[temp6];
 	 proErr2NLO[temp6] = proXSNLO2m[temp6] - proXSNLO1m[temp6];
		}
      else{
 	proXSNLO2m[process[i]] = XSNLO[i];
 	proErr1NLO[process[i]] = proXSNLO05m[process[i]] - proXSNLO1m[process[i]];
 	proErr2NLO[process[i]] = proXSNLO2m[process[i]] - proXSNLO1m[process[i]];
      }	 
    
}

    std::map<std::string, double> proErrNLO;
  for( int i = 0; i<process.size(); i+=3 ) {
    
     if(process[i]=="bb"|| process[i]=="tb"){
	 temp5 = std::to_string(SQFam[i]); 	
         temp6 = process[i] + temp5;	
      if( std::abs(proErr1NLO[temp6]) > std::abs(proErr2NLO[temp6])){
        proErrNLO[temp6] = std::abs(proErr1NLO[temp6]); 
       }
   
      else{
        proErrNLO[temp6] = std::abs(proErr2NLO[temp6]); 
     }
   }
     else{
      if( std::abs(proErr1NLO[process[i]]) > std::abs(proErr2NLO[process[i]])){
        proErrNLO[process[i]] = std::abs(proErr1NLO[process[i]]); 
       }
   
      else{
        proErrNLO[process[i]] = std::abs(proErr2NLO[process[i]]); 
     }
   }
  }

  for( const auto& elem : proXSNLO05m ) {

 

// outputs are set 
  
    SetOutput( "Pro_XSNLO_05m_" + elem.first, elem.second  );
  }
  for( const auto& elem : proXSNLO1m ) {

   
    SetOutput( "Pro_XSNLO_1m_" + elem.first, elem.second  );
  }
  for( const auto& elem : proXSNLO2m ) {

   
    SetOutput( "Pro_XSNLO_2m_" + elem.first, elem.second  );
  }

  for( const auto& elem : proErrNLO ) {

   
    SetOutput( "Pro_ErrNLO_" + elem.first, elem.second  );
  }

}


void Fittino::ProspinoCalculator::SetupMeasuredValues() {

}

void Fittino::ProspinoCalculator::Initialize() {
  
}

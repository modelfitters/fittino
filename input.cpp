/***************************************************************************
                                  input.cpp
                             -------------------    
    Class that parses a given input file, stores the read values and
    makes them available.
                             -------------------
    $Id$
    $Name$
                             -------------------
    begin                : November 30, 2003
    authors              : Philip Bechtle, Peter Wienemann
    email                : philip.bechtle@desy.de, peter.wienemann@desy.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
#include <input.h>
#include <yy.h>
#include <units.h>
#include <values.h>
#include <misc.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMath.h>

vector<doubleVec_t> CrossSectionProduction;

Input::Input(const char* inputfile)
{
  // Set defaults
  yyCalculator = SPHENO;
  yyRelicDensityCalculator = NORELICDENSITYCALCULATOR;
  yyLEOCalculator = NOLEOCALCULATOR;
  yyISR       = true;
  FillNameMap();
  
  // Open input file and parse it
  yyin = fopen(inputfile, "r");
  if (!yyin) {
      cerr<<"Input file "<<inputfile<<" not found"<<endl;
      exit(EXIT_FAILURE);
  }
  yyInputFileLineNo = 1;
  yyparse();
  fclose(yyin);

  if (yyParseError) {
    cerr<<"ERROR: Syntax error in Fittino input file "<<inputfile<<endl;
    cerr<<"       Check above messages for details"<<endl;
    exit(EXIT_FAILURE);
  }
  if (yySetErrorFlag) {
    cerr<<"ERROR: Syntax error in Fittino input file "<<inputfile<<endl;
    cerr<<"       Check above messages for details"<<endl;
    exit(EXIT_FAILURE);
  }

  yyInputFileLineNo = 1;

  /* check whether defaults are given */
  CheckDefault("massZ"     , ID_Z  ,  91.1876    , 0.0021     );
  //CheckDefault("massW"     , ID_W  ,  80.423     , 0.039      );
  CheckDefault("massCharm" , ID_c  ,   1.27      , 0.2        );
  CheckDefault("massBottom", ID_b  ,   4.5       , 0.25       );
  CheckDefault("massTop"   , ID_t  , 174.3       , 5.1        );
  CheckDefault("massTau"   , ID_tau,   1.77699   , 0.00029    );
  CheckDefault("alphas"    , 0     ,   0.1187    , 0.002      );
  // G_F is skipped because its small uncertainty causes numerical problems with covariance matrix inversion.
  CheckDefault("G_F"       , 0     ,   1.16639e-5, 0.000001e-5);
  CheckDefault("alphaem"   , 0     , 127.934     , 0.027      );
  //  CheckDefault("QEWSB"     , 0     , 1000.       , 50.        );

  // cout << "mass = " << mass << endl; 

  //  cout << "after yyparse" << endl;
  //  DumpMeasuredVector();

  /* Fill IDs, bounds, positions */
  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
    if ( yyMeasuredVec[i].bound_low == 0. && yyMeasuredVec[i].bound_up == 0. ) {
      yyMeasuredVec[i].bound_low = -1.e+6; // MINDOUBLE;
      yyMeasuredVec[i].bound_up  = 1.e+6;  // MAXDOUBLE;	
      if (!strncmp(yyMeasuredVec[i].name.c_str(), "mass", 4)) {
	yyMeasuredVec[i].bound_low = 0.;
	yyMeasuredVec[i].bound_up  = 10.*yyMeasuredVec[i].value;
      }
    }
   }    

  /* fill IDs */
  Fill_IDs();

  //filling positions in vector to be fitted
  for(unsigned int i=0; i<yyMeasuredVec.size(); i++) {
     //if (yyMeasuredVec[i].nofit == true) {
     //   yyMeasuredVec[i].position = -1;
     //}
     //else {
	yyMeasuredVec[i].position = i;
     //}
  }

  //  cout << "after Fill_IDs" << endl;
  // DumpMeasuredVector();

  /* fill diagonal elements of correlation matrix, calculate covariance matrix and its inverse */
  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
     yyMeasuredCorrelationMatrix.add(i, i, 1);
  }

  TagMap* yyTagMap = new TagMap(&yyMeasuredVec);

  if (!yyTagMap->empty()) {
     CorrelationMatrix* correlationMatrix = new CorrelationMatrix(&yyMeasuredVec, yyTagMap);
     correlationMatrix->Calculate();

     for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
	for (unsigned int j=i; j<yyMeasuredVec.size(); j++) {
	   yyMeasuredCorrelationMatrix.add(i,j,correlationMatrix->GetCorrelation(i,j));
	}
     }
  }

  //cout << "filling measured covariances "<< endl;
  yyMeasuredCorrelationMatrix.CalculateCovarianceMatrix();
  //cout << endl << "correlation and covariance matrix before introducing direct input covariance:" << endl;
  yyMeasuredCorrelationMatrix.Print();
  //yyMeasuredCorrelationMatrix.PrintCovariance();
  //yyMeasuredCorrelationMatrix.PrintInverseCovariance();

  // cout << "calling AddFullCovarianceMatrix" << endl;
  // TODO: add DirectInputCovariance Matrix here
  if (yyDirectInputCovarianceMatrix.TestCovarianceMatrixExistence()) {
    yyMeasuredCorrelationMatrix.AddFullCovarianceMatrix(yyDirectInputCovarianceMatrix);

    // test
    // cout << endl << "correlation and covariance matrix after adding direct input covariance matrix to yyMeasuredCorrelationMatrix:" << endl;
    //yyMeasuredCorrelationMatrix.Print();
    //yyMeasuredCorrelationMatrix.PrintCovariance();
    
    // TODO: calculate correlation matrix again
    yyMeasuredCorrelationMatrix.CalculateInverseCovarianceMatrix();
    yyMeasuredCorrelationMatrix.TransformCovarianceMatrixIntoCorrelationMatrix();
  }

  // test
  //cout << endl << "correlation and covariance matrix after calculating correlation and inverse covariance matrix::" << endl;
  //yyTagMap->Print();
  //yyMeasuredCorrelationMatrix.Print();
  //yyMeasuredCorrelationMatrix.PrintCovariance();
  //yyMeasuredCorrelationMatrix.PrintInverseCovariance();

  fMeasuredVec = yyMeasuredVec;
  fMeasuredCorrelationMatrix = yyMeasuredCorrelationMatrix;

  //  cout << "before Fill Cross SectionProd" << endl;
  //DumpMeasuredVector();

  FillCrossSectionProduction();

  // Check that each alias number occurs only once per observable type
  for (unsigned int itype = 1; itype<10; itype++) { // loop over ObservableType (excluding masses = 0)
     for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
	if (yyMeasuredVec[i].type == (int)itype) {
	   for (unsigned int j=0; j<yyMeasuredVec.size(); j++) {
	      if ((int)i != (int)j && yyMeasuredVec[j].type == (int)itype &&
		    yyMeasuredVec[i].alias == yyMeasuredVec[j].alias) {
		 cerr<<"ERROR: Alias "<<yyMeasuredVec[i].alias<<" is used for "
		    <<yyMeasuredVec[i].name<<" and "<<yyMeasuredVec[j].name<<endl;
		 exit(EXIT_FAILURE);
	      }
	   }
	}
     }
  }

  //  cout << "after Fill Cross SectionProd" << endl;
  //DumpMeasuredVector();

  // list vector
  cout << "------------------------------------------------------------" << endl;
  cout << "Listing input observables" << endl;
  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
     if (yyMeasuredVec[i].nofit == false) {
	cout << yyMeasuredVec[i].name << " " << yyMeasuredVec[i].value << " +- " << yyMeasuredVec[i].error << endl;
     }
  }


  //cout << "before End" << endl;
  //DumpMeasuredVector();
}


Input::~Input()
{
   fMeasuredVec.clear();
}

void Input::CheckDefault(string name, int ID, double value, double error)
{
   MeasuredValue tmpVec;
   int found;
   found = 0;
   for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
      if (!yyMeasuredVec[i].name.compare(name)) {
	 found = 1;
	 break;
      }
   }    
   if (found == 0) {
      tmpVec.nofit = false;
      tmpVec.name = name;
      tmpVec.value = value;
      tmpVec.error = error;
      tmpVec.id = ID;
      tmpVec.alias = 0;
      tmpVec.theovalue = 0.;
      tmpVec.bound_low = -1.e+6;
      tmpVec.bound_up = 1.e+6;
      if (!strncmp(name.c_str(), "mass", 4)) {
	 tmpVec.type = mass;
	 //      cout << "assigning type mass " << mass << " to " << name << endl; 
	 tmpVec.bound_low = 0.;
	 tmpVec.bound_up = 10.*tmpVec.value;	
      }
      else tmpVec.type = other;
      yyMeasuredVec.push_back(tmpVec);
      cout << "Adding default value: "<<name<<" = "<<value<<" +- "<<error<<endl;
   }

   return;
}


void Input::Fill_IDs()
{

   string tmpname;

   for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
      if (!strncmp(yyMeasuredVec[i].name.c_str(), "mass", 4)) {
	 tmpname = yyMeasuredVec[i].name;
	 tmpname.erase(0,4);
	 yyMeasuredVec[i].id = yyParticleIDs[tmpname];
	 //cout << " FILLIDS: " << tmpname << " id = " << yyMeasuredVec[i].id << endl;
      }
   }

   return;

}

void Input::FillNameMap()
{
   /*  
   // THIS PART IS DISABLED: NMSSM STUFF PARTLY OVERWRITTEN BY NON-NMSSM STUFF 
   yyParticleNames[ID_d       ] = "Down";
   yyParticleNames[ID_u       ] = "Up"             ;   
   yyParticleNames[ID_s       ] = "Strange"	    ;
   yyParticleNames[ID_c       ] = "Charm"	    ;
   yyParticleNames[ID_b       ] = "Bottom"	    ;
   yyParticleNames[ID_t       ] = "Top"	    ;

   yyParticleNames[ID_e       ] = "Electron"	    ;
   yyParticleNames[ID_nue     ] = "Nue"	    ;
   yyParticleNames[ID_mu      ] = "Muon"	    ;
   yyParticleNames[ID_numu    ] = "Numu"	    ;
   yyParticleNames[ID_tau     ] = "Tau"	    ;
   yyParticleNames[ID_nutau   ] = "Nutau"	    ;

   yyParticleNames[ID_g       ] = "Gluon"	    ;
   yyParticleNames[ID_gamma   ] = "Gamma"	    ;
   yyParticleNames[ID_Z       ] = "Z"		    ;
   yyParticleNames[ID_W       ] = "W"		    ;

   // ATTENTION: NMSSM STUFF PARTLY OVERWRITTEN BY NON-NMSSM STUFF 
   yyParticleNames[ID_H1      ] = "H1"		    ; // NMSSM
   yyParticleNames[ID_H2      ] = "H2"		    ; // NMSSM
   yyParticleNames[ID_H3      ] = "H3"		    ; // NMSSM
   yyParticleNames[ID_A1      ] = "A1"		    ; // NMSSM
   yyParticleNames[ID_A2      ] = "A2"		    ; // NMSSM

   yyParticleNames[ID_h       ] = "h0"		    ;
   yyParticleNames[ID_H       ] = "H0"		    ;
   yyParticleNames[ID_A       ] = "A0"		    ;
   yyParticleNames[ID_Hplus   ] = "Hplus"	    ;

   yyParticleNames[ID_G       ] = "Graviton"	    ;

   yyParticleNames[ID_sdL     ] = "SdownL"         ;
   yyParticleNames[ID_suL     ] = "SupL"           ;
   yyParticleNames[ID_ssL     ] = "SstrangeL"      ;
   yyParticleNames[ID_scL     ] = "ScharmL"        ;
   yyParticleNames[ID_sb1     ] = "Sbottom1"       ;
   yyParticleNames[ID_st1     ] = "Stop1"          ;

   yyParticleNames[ID_seL     ] = "SelectronL"     ;
   yyParticleNames[ID_snueL   ] = "SnueL"          ;
   yyParticleNames[ID_smuL    ] = "SmuL"           ;
   yyParticleNames[ID_snumuL  ] = "SnumuL"         ;
   yyParticleNames[ID_stau1   ] = "Stau1"          ;
   yyParticleNames[ID_snutauL ] = "Snutau1"        ;

   yyParticleNames[ID_sdR     ] = "SdownR"    	    ;
   yyParticleNames[ID_suR     ] = "SupR"           ;
   yyParticleNames[ID_ssR     ] = "SstrangeR"      ;
   yyParticleNames[ID_scR     ] = "ScharmR"        ;
   yyParticleNames[ID_sb2     ] = "Sbottom2"       ;
   yyParticleNames[ID_st2     ] = "Stop2"          ;

   yyParticleNames[ID_seR     ] = "SelectronR"     ;
   yyParticleNames[ID_snueR   ] = "SnueR"          ;
   yyParticleNames[ID_smuR    ] = "SmuR"           ;
   yyParticleNames[ID_snumuR  ] = "SnumuR"         ;
   yyParticleNames[ID_stau2   ] = "Stau2"          ;
   yyParticleNames[ID_snutauR ] = "Snutau2"        ;

   yyParticleNames[ID_gluino  ] = "Gluino"  	    ;
   yyParticleNames[ID_chi01   ] = "Neutralino1"    ;
   yyParticleNames[ID_chi02   ] = "Neutralino2"    ;
   yyParticleNames[ID_chip1   ] = "Chargino1"      ;
   yyParticleNames[ID_chi03   ] = "Neutralino3"    ;
   yyParticleNames[ID_chi04   ] = "Neutralino4"    ;
   yyParticleNames[ID_chi05   ] = "Neutralino5"    ; // for NMSSM
   yyParticleNames[ID_chip2   ] = "Chargino2" 	    ;
   yyParticleNames[ID_gravitino ] = "Gravitino"  ;
   */
      //---------------------------------------------------
      //---------------------------------------------------

      yyParticleIDs["Down"       ] = ID_d          ;
   yyParticleIDs["Up"         ] = ID_u          ;   
   yyParticleIDs["Strange"	 ] = ID_s          ;
   yyParticleIDs["Charm"	 ] = ID_c          ;
   yyParticleIDs["Bottom"	 ] = ID_b          ;
   yyParticleIDs["Top"	 ] = ID_t          ;

   yyParticleIDs["Electron"	 ] = ID_e          ;
   yyParticleIDs["Nue"	 ] = ID_nue        ;
   yyParticleIDs["Muon"	 ] = ID_mu         ;
   yyParticleIDs["Numu"	 ] = ID_numu       ;
   yyParticleIDs["Tau"	 ] = ID_tau        ;
   yyParticleIDs["Nutau"	 ] = ID_nutau      ;

   yyParticleIDs["Gluon"	 ] = ID_g          ;
   yyParticleIDs["Gamma"	 ] = ID_gamma      ;
   yyParticleIDs["Z"		 ] = ID_Z          ;
   yyParticleIDs["W"		 ] = ID_W          ;

   yyParticleIDs["h0"	 ] = ID_h          ;
   yyParticleIDs["H0"	 ] = ID_H          ;
   yyParticleIDs["A0"	 ] = ID_A          ;
   yyParticleIDs["Hplus"	 ] = ID_Hplus      ;

   yyParticleIDs["H1"	 ] = ID_H1          ; // NMSSM
   yyParticleIDs["H2"	 ] = ID_H2          ; // NMSSM
   yyParticleIDs["H3"	 ] = ID_H3          ; // NMSSM
   yyParticleIDs["A1"	 ] = ID_A1          ; // NMSSM
   yyParticleIDs["A2"	 ] = ID_A2          ; // NMSSM

   yyParticleIDs["Graviton"	 ] = ID_G          ;

   yyParticleIDs["SdownL"       ] = ID_sdL       ;
   yyParticleIDs["SupL"         ] = ID_suL       ;
   yyParticleIDs["SstrangeL"    ] = ID_ssL       ;
   yyParticleIDs["ScharmL"      ] = ID_scL       ;
   yyParticleIDs["Sbottom1"     ] = ID_sb1       ;
   yyParticleIDs["Stop1"        ] = ID_st1       ;

   yyParticleIDs["SelectronL"   ] = ID_seL       ;
   yyParticleIDs["SnueL"        ] = ID_snueL     ;
   yyParticleIDs["SmuL"         ] = ID_smuL      ;
   yyParticleIDs["SnumuL"       ] = ID_snumuL    ;
   yyParticleIDs["Stau1"        ] = ID_stau1     ;
   yyParticleIDs["Snutau1"      ] = ID_snutauL   ;

   yyParticleIDs["SdownR"       ] = ID_sdR       ;
   yyParticleIDs["SupR"         ] = ID_suR       ;
   yyParticleIDs["SstrangeR"    ] = ID_ssR       ;
   yyParticleIDs["ScharmR"      ] = ID_scR       ;
   yyParticleIDs["Sbottom2"     ] = ID_sb2       ;
   yyParticleIDs["Stop2"        ] = ID_st2       ;

   yyParticleIDs["SelectronR"   ] = ID_seR       ;
   yyParticleIDs["SnueR"        ] = ID_snueR     ;
   yyParticleIDs["SmuR"         ] = ID_smuR      ;
   yyParticleIDs["SnumuR"       ] = ID_snumuR    ;
   yyParticleIDs["Stau2"        ] = ID_stau2     ;
   yyParticleIDs["Snutau2"      ] = ID_snutauR   ;

   yyParticleIDs["Gluino"  	   ] = ID_gluino    ;
   yyParticleIDs["Neutralino1"  ] = ID_chi01     ;
   yyParticleIDs["Neutralino2"  ] = ID_chi02     ;
   yyParticleIDs["Chargino1"    ] = ID_chip1     ;
   yyParticleIDs["Neutralino3"  ] = ID_chi03     ;
   yyParticleIDs["Neutralino4"  ] = ID_chi04     ;
   yyParticleIDs["Neutralino5"  ] = ID_chi05     ; // NMSSM
   yyParticleIDs["Chargino2"    ] = ID_chip2     ;
   yyParticleIDs["Gravitino"    ] = ID_gravitino ;



   return;
}

int Input::ReturnParticleID(string name)
{
   return yyParticleIDs[name];
}


string Input::ReturnParticleName(int ID)
{
   cerr << "ReturnParticleName has been disabled due to NMSSM name clash" << endl;
   cerr << "(see warning in implementation of Input::FillNameMap)" << endl;
   exit(EXIT_FAILURE);
   //  return yyParticleNames[ID];
}


void Input::DumpMeasuredVector()
{
   for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
      cout<<yyMeasuredVec[i].name<<" = "<<yyMeasuredVec[i].value<<" +- "
	 <<yyMeasuredVec[i].error<<endl;
      cout<<"    theovalue = "<<yyMeasuredVec[i].theovalue<<endl;
      cout<<"    id        = "<<yyMeasuredVec[i].id<<endl;
      cout<<"    bound_low = "<<yyMeasuredVec[i].bound_low<<endl;
      cout<<"    bound_up  = "<<yyMeasuredVec[i].bound_up<<endl;
      cout<<"    type      = "<<yyMeasuredVec[i].type<<endl;
      cout<<"    alias     = "<<yyMeasuredVec[i].alias<<endl;
      cout<<"    sqrts     = "<<yyMeasuredVec[i].sqrts<<endl;
      cout<<"    pol1      = "<<yyMeasuredVec[i].polarisation1<<endl;
      cout<<"    pol2      = "<<yyMeasuredVec[i].polarisation2<<endl;
      cout<<"    daughter  = ";
      for (unsigned int j=0; j<yyMeasuredVec[i].daughters.size(); j++) {
	 cout<<yyMeasuredVec[i].daughters[j]<<" ";
      }
      cout<<endl;
      cout<<"    products  = ";
      for (unsigned int j=0; j<yyMeasuredVec[i].products.size(); j++) {
	 cout<<yyMeasuredVec[i].products[j]<<" ";
      }

      cout<< endl;
   }
}

void Input::FillCrossSectionProduction()
{
   vector<double> tmpVec;
   unsigned int k;
   bool found;

   CrossSectionProduction.clear();
   for (unsigned int i = 0; i < yyMeasuredVec.size(); i++) {
      if (yyMeasuredVec[i].type == xsection) {
	 k = CrossSectionProduction.size();
	 found = false;
	 for ( unsigned int j = 0; j < k; j++ ) {
	    if ((CrossSectionProduction[j][0]==yyMeasuredVec[i].sqrts)
		  && (CrossSectionProduction[j][1]==yyMeasuredVec[i].polarisation1)
		  && (CrossSectionProduction[j][2]==yyMeasuredVec[i].polarisation2) ) {
	       found = true;
	       break;
	    }
	 }
	 if (!found) {
	    tmpVec.clear();
	    tmpVec.push_back(yyMeasuredVec[i].sqrts);
	    tmpVec.push_back(yyMeasuredVec[i].polarisation1);
	    tmpVec.push_back(yyMeasuredVec[i].polarisation2);
	    CrossSectionProduction.push_back(tmpVec);
	 }
      }
   }
   if (CrossSectionProduction.size()==0) {
      tmpVec.clear();
      tmpVec.push_back(1000.);
      tmpVec.push_back(0.0);
      tmpVec.push_back(0.0);
      CrossSectionProduction.push_back(tmpVec);    
   }

   return;

}

void Input::ScatterObservablesBefore() {

   cout << yyDashedLine << endl;
   cout << "Scattering observables around their true value" << endl;

   // fill local vector trueValueObservableVector with true observable values from input file
   TVectorD trueValueObservableVector(yyMeasuredVec.size());
   for (unsigned int j=0; j<yyMeasuredVec.size(); j++) {
      trueValueObservableVector(j) = yyMeasuredVec[j].value;
   }

   // perform scattering
   gRandom->SetSeed(yySeedForObservableScattering); // Use fixed seed here
   TVectorD pseudoMeasuredObservableVector(yyMeasuredVec.size());
   pseudoMeasuredObservableVector = getCorrelatedRandomNumbers( trueValueObservableVector, this->GetMeasuredCorrelationMatrix().GetCovarianceMatrix() );
   for (unsigned int j = 0; j < yyMeasuredVec.size(); j++) {
      yyMeasuredVec[j].value = pseudoMeasuredObservableVector(j);
   }

   // test printouts for initial value scattering
   cout << "new mean values for" << endl;
   for (unsigned int j = 0; j < yyMeasuredVec.size(); j++) {
      if (yyMeasuredVec[j].nofit == false) {
         cout << yyMeasuredVec[j].name << ": " << yyMeasuredVec[j].value << " (deviation from true value " << trueValueObservableVector(j) << " within " << TMath::Ceil((TMath::Abs(yyMeasuredVec[j].value-trueValueObservableVector(j))/yyMeasuredVec[j].error)) <<" sigma)" <<endl;
      }
   }

}
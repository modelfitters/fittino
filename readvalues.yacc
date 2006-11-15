%{
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <map>
#include <leshouches.h>
#include <yy.h>
#include <misc.h>

char yyBlockName[255];
double yyScale;
// SMINPUTS
double yyG_F;		     
double yyoneoveralpha_em_mz; 
double yyalpha_s_mz;	     
double yyMassZ;		     
double yyMass_b;	     
double yyMass_t;	     
double yyMass_tau;              
//MODEL
int    yyModel;
// MINPAR
double yy_Qmax;
double yytanb;
double yy_m0;
double yy_m12;
double yy_signmu;
double yy_A;
// EXTPAR
double yyMinput;
double yyM1;
double yyM2;
double yyM3;
double yy_At;
double yy_Ab;
double yy_Atau;
double yy_m2Hd;
double yy_m2Hu;
double yy_mu;
double yy_MassA;
double yy_mse1L;
double yy_mse2L;
double yy_mse3L;
double yy_mse1R;
double yy_mse2R;
double yy_mse3R;
double yy_msq1L;
double yy_msq2L;
double yy_msq3L;
double yy_msu1R;
double yy_msu2R;
double yy_msu3R;
double yy_msd1R;
double yy_msd2R;
double yy_msd3R;
double yy_N51;
double yy_N52;
double yy_N53;

// SPINFO
char yy_spectrum_calc_name[256];
char yy_spectrum_calc_version[256];
char yy_warnings[1024];
char yy_errors[1024];
char yy_error_codes[1024];
typedef struct {
    char tmpname[256];
} tmpstr_t;
tmpstr_t tmpstr;
vector<tmpstr_t> tmpStrings;
vector<int> tmpNumbers;


// MASS
map<int,double> yyMass;


// DECAY
decay_element_t tmp_branch;
map<int,decay_element_t> branching_ratios;


// XS
xs_element_t tmp_xs;
map<doubleVec_t,xs_element_t> cross_sections;



doubleVec_t tmpVec, tmp2Vec;
vector<doubleVec_t> tmpParams;

int found, skip;
 
vector<MeasuredValue> yyMeasuredVec; // contains measured observables
vector<MeasuredValue> yyThrownVec;   // contains smeared observables
 
vector<MeasuredValue> yyFittedVec; // contains mu, M1, M2 ...
vector<MeasuredValue> yyUniversalityVec;
vector<MeasuredValue> yyFixedVec;  // contains mu, M1, M2 ...
 
CorrelationMatrix     yyMeasuredCorrelationMatrix(&yyMeasuredVec);
CorrelationMatrix     yyFittedCorrelationMatrix(&yyFittedVec);

vector <parameter_t> yyFixedPar;
vector <parameter_t> yyFittedPar;
vector<ScanParameter> yyScanPar; // contains parameters to scan
vector <parameter_t> indchisq_vec;
 
bool          yyUseLoopCorrections = true;
bool          yyCalcPullDist;
bool          yyScanParameters = false;
bool          yyPerformFit = true;
bool          yyISR;
bool          yyCalculatorError;
bool          yyUseMinos;
bool          yyGetContours = false;
bool          yyUseHesse;
bool          yyUseSimAnnBefore = false;
bool          yyUseSimAnnWhile = false;
bool          yyUseGivenStartValues;
bool          yyFitAllDirectly;
bool          yyCalcIndChisqContr;
bool          yyBoundsOnX = true;
bool          yySepFitTanbX = true;
bool          yySepFitTanbMu = false;
bool          yySepFitmA = false;
bool          yyScanX = true;
bool          yyVerbose = true;
bool          yyAdaptiveSimAnn = false;
bool          yyNoBoundsAtAll = false;
bool          yySimAnnUncertainty = false;
bool          yySimAnnUncertaintyRunDown = false;
bool          yyRandomDirUncertainties = false;
bool          yyPerformSingleFits = false;
bool          yyUseMicrOmegas = false;

unsigned int yyCalculator;
string       yyCalculatorPath = "";

unsigned int yyFitModel = MSSM;
 
map<int,string> yyParticleNames;
map<string,int> yyParticleIDs;

int           yyParseError = 0;
unsigned int  yyNumberOfMinimizations = 1;
double        yyErrDef = 1.;
int           yyMaxCallsSimAnn = 300000;
double        yyTempRedSimAnn = 0.4;
double        yyInitTempSimAnn = -1.;
int           yyNumberPulls = 0;

int           yyNumberOfDirections = 10000;

int           yyRandomGeneratorSeed = -1;

double        yyXscanlow = -6000.;
double        yyXscanhigh = 2000.;

//##############################################################
double        yybsg    = -10000.;   //  1   BR(b -> s gamma) 
double        yybsmm   = -10000.;   //  2   BR(b -> s mu+ mu-)
double        yyB_smm  = -10000.;   //  3   BR(Bs -> mu+ mu-)
double        yyB_utn  = -10000.;   //  4   BR(B_u -> tau nu)
double        yydMB_d  = -10000.;   //  5   |Delta(M_Bd)| [ps^-1]
double        yydMB_s  = -10000.;   //  6   |Delta(M_Bs)| [ps^-1]
double        yygmin2e = -10000.;   // 10   Delta(g-2)_electron
double        yygmin2m = -10000.;   // 11   Delta(g-2)_muon
double        yygmin2t = -10000.;   // 12   Delta(g-2)_tau
double        yydrho   = -10000.;   // 30   Delta(rho_parameter)
//==============================================================
double        yyOmega  = -10000.;   //      relic density
//##############################################################

double        yyMaxCalculatorTime = 20.;

std::vector<InputFileLine> yyInputFile;
struct InputFileLine yyInputFileLine;
int           yyInputFileLineNo = 1;

%}

%union {
    int      integer;
    double   real;
    char     name[255];
};

%token <name> T_KEY T_WORD
%token <real> T_NUMBER
%token <integer> T_ENERGYUNIT T_SWITCHSTATE T_CROSSSECTIONUNIT
%token T_ERRORSIGN T_BRA T_KET T_COMMA T_GOESTO T_ALIAS T_NOFIT
%token T_BLOCK T_SCALE T_DECAY T_NEWLINE T_BR T_LEO T_XS T_CALCULATOR T_XSBR T_BRRATIO
%token <name> T_COMPARATOR T_UNIVERSALITY T_PATH T_NEWLINE
 
%type <name> sentence
%type <real> value err

%%

input:
            /* empty */
            | input T_NEWLINE
              {
	        if (yyInputFileLine.prevalue.size()) yyInputFileLine.postvalue += "\t";
		yyInputFileLine.postvalue += $2;
		yyInputFile.push_back(yyInputFileLine);
		yyInputFileLine.prevalue.erase();
//		yyInputFileLine.error = -1;
		yyInputFileLine.error.clear();
		yyInputFileLine.postvalue.erase();
		yyInputFileLineNo++;
              }
            | input block
            | input decay
            | input xs
              {
//                  printf("xs input\n");
              }
	    | input T_KEY value
	      {
		char c[1000];
		yyInputFileLine.prevalue  = $2;
                yyInputFileLine.prevalue += "\t";
		sprintf(c, "%f", $3);
                yyInputFileLine.prevalue += c;

	        if (!strcmp($2,"NumberOfMinimizations")) {
		  // cout << "FOUND NUMBER OF MINIMIZATIONS "<<$3<<endl;
		  if ($3>0) {
                    yyNumberOfMinimizations = (int)$3;
                  }
		}
	        else if (!strcmp($2,"ErrDef")) {
		  // cout << "FOUND ErrDef "<<$3<<endl;
		  if ($3>0) {
                    yyErrDef = $3;
                  }
                } 	        
	        else if (!strcmp($2,"MaxCallsSimAnn")) {
		  // cout << "FOUND MaxCallsSimAnn "<<$3<<endl;
		  if ($3>0) {
                    yyMaxCallsSimAnn = (int)$3;
                  }
                } 	        
	        else if (!strcmp($2,"TempRedSimAnn")) {
		  // cout << "FOUND TempRedSimAnn "<<$3<<endl;
		  if ($3>0) {
                    yyTempRedSimAnn = $3;
                  }
                } 	// yyInitTempSimAnn        
	        else if (!strcmp($2,"InitTempSimAnn")) {
		  // cout << "FOUND TempRedSimAnn "<<$3<<endl;
		  if ($3>0) {
                    yyInitTempSimAnn = $3;
                  }
                } 	
		else if (!strcmp($2,"RandomGeneratorSeed")) {
		  // cout << "FOUND RandomGeneratorSeed "<<$3<<endl;
		  yyRandomGeneratorSeed = (int)$3;
                }
		else if (!strcmp($2,"NumberOfDirections")) {
		  // cout << "FOUND NumberOfDirections "<<$3<<endl;
		  yyNumberOfDirections = (int)$3;
                }
                // yyMaxCalculatorTime
		else if (!strcmp($2,"MaxCalculatorTime")) {
		  yyMaxCalculatorTime = (double)$3;
                }
                else if (!strcmp($2,"NumberPulls")) {
		  // cout << "FOUND NumberPulls "<<$3<<endl;
		  if ($3>0) {
                    yyNumberPulls = (int)$3;
                  }
       		}
                else {
		  found = 0;
		  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
		    if (!yyMeasuredVec[i].name.compare($2)) {
		      found = 1;
		      yyMeasuredVec[i].value = $3;
		      yyMeasuredVec[i].error = -1;
		      break;
		    }
		  }
		  if (found == 0) {
		      MeasuredValue tmpValue;
                      tmpValue.nofit = false;
		      tmpValue.name = $2;
		      tmpValue.type = mass;
		      tmpValue.value = $3;
		      tmpValue.error = -1;
		      tmpValue.theovalue  = 0;
		      tmpValue.alias = 0;
		      tmpValue.id = 0;
		      tmpValue.bound_low = 0.;
		      tmpValue.bound_up = 0.;
		      if (!strncmp($2, "mass", 4))tmpValue.type = mass;
		      else tmpValue.type = other;
		      yyMeasuredVec.push_back(tmpValue);
		      if (!strncmp($2, "cos", 3)) {
			cout << "COS: " << tmpValue.name << " " << tmpValue.type << 
			  " " << other << " " << Pedge << endl;
		      }
//		      cout << "filling "<< $2 <<endl;
		  }
                }
	      }
	    | input T_KEY value err
	      {
		  yyInputFileLine.prevalue = $2;
		  yyInputFileLine.value = $3;
//		  yyInputFileLine.error = $4;

		  found = 0;
		  skip = 0;
		  if ($4*$4 < 2.2204e-16) {
                      cout<<"WARNING: Measurement of "<<$2<<" will not be used in fit. Uncertainty too small."<<endl;
		      cout<<"         Too small uncertainties cause numerical problems with covariance matrix inversion."<<endl;
		      cout<<"         Press any key to continue."<<endl;
                      getchar();
		      skip = 1;
                  }
		  if (!skip) {
		      for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
			if (!yyMeasuredVec[i].name.compare($2)) {
			  found = 1;
			  yyMeasuredVec[i].value = $3;
			  yyMeasuredVec[i].error = $4;
			  break;
			}
		      }
		      if (found == 0) {
			  MeasuredValue tmpValue;
			  tmpValue.nofit = false;
			  tmpValue.name = $2;
			  tmpValue.value = $3;
			  tmpValue.error = $4;
			  tmpValue.type = mass;
			  tmpValue.theovalue  = 0;
			  tmpValue.bound_low = 0.;
			  tmpValue.bound_up = 0.;
			  tmpValue.alias = 0;
			  tmpValue.id = 0;
			  if (!strncmp($2, "mass", 4))tmpValue.type = mass;
			  else tmpValue.type = other;
			  yyMeasuredVec.push_back(tmpValue);
			  if (!strncmp($2, "cos", 3)) {
			    cout << "COS: " << tmpValue.name << " " << tmpValue.type << endl;
			  }
		      }
                  }
	      }
	    | input T_UNIVERSALITY sentence
	      {
		 yyInputFileLine.prevalue  = $2;
		 yyInputFileLine.prevalue += "\t";
	         yyInputFileLine.prevalue += $3;

		 MeasuredValue tmpValue;
        	 string str;
	         str.erase();
		 str = $3;
//		 cout << "decompositing string "<<str<< endl;
		 char tmpstr5[255];
		 unsigned int pos = 0, newpos = 0;
//		 int anti = 1;
//	         int i = 0;
	         MeasuredValue tmpval;
	         string        firstname;
		 while ((newpos = str.find(" ", pos)) != string::npos) {
		      str.copy(tmpstr5, newpos - pos, pos);
	              tmpstr5[newpos-pos] = '\0';
                      if (!pos) {
                          firstname.erase();
                          firstname.append(tmpstr5);
                      }
	              else {
	                  tmpval.name.erase();
	                  tmpval.name.append(tmpstr5);
                          tmpval.universality = firstname;
                          yyUniversalityVec.push_back(tmpval);
                      }
		      pos = newpos + 1;
                 }		 					
              }	
	    | input T_KEY T_WORD T_WORD value
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += "\t";
		  yyInputFileLine.prevalue += $3;
		  yyInputFileLine.prevalue += "\t";
		  yyInputFileLine.prevalue += $4;
		  yyInputFileLine.prevalue += "\t";
	 	  sprintf(c, "%f", $5);
		  yyInputFileLine.prevalue += c;

		  if (!strcmp($2, "correlationCoefficient")) {
	              int i = 0;
		      int j = 0;
                      int aliasnumber;
                      char* charnumber;
		      if (!strncmp($3,"sigma",5)) {
                        charnumber = strchr($3,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == xsection) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    i = k;
//                            cout << "found correlation with xs " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
		      else if (!strncmp($3,"br",2)) {
                        charnumber = strchr($3,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    i = k;
//                            cout << "found correlation with br " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
		      else if (!strncmp($3,"edge",2)) {
                        charnumber = strchr($3,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == Pedge) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    i = k;
//                            cout << "found correlation with edge " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
		      else if (!strncmp($3,"width",2)) {
                        charnumber = strchr($3,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == Pwidth) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    i = k;
//                            cout << "found correlation with width " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
                      else {
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if (!yyMeasuredVec[k].name.compare($3)) {
			    i = k;
//                            cout << "found correlation with " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }			
                      }
		      if (!strncmp($4,"sigma",5)) {
                        charnumber = strchr($4,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == xsection) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    j = k;
//                            cout << "found correlation with xs " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
		      else if (!strncmp($4,"br",2)) {
                        charnumber = strchr($4,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    j = k;
//                            cout << "found correlation with br " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
		      else if (!strncmp($4,"edge",2)) {
                        charnumber = strchr($4,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == Pedge) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    j = k;
//                            cout << "found correlation with edge " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
		      else if (!strncmp($4,"width",2)) {
                        charnumber = strchr($4,'_');
                        if (charnumber == 0) yyerror ("Underscore not found");
		        aliasnumber = atoi((charnumber+1));
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if ((yyMeasuredVec[k].type == Pwidth) && (yyMeasuredVec[k].alias == aliasnumber)) {
			    j = k;
//                            cout << "found correlation with width " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }
                      }
                      else {
                        for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
                          if (!yyMeasuredVec[k].name.compare($4)) {
			    j = k;
//                            cout << "found correlation with " <<  yyMeasuredVec[k].name << endl;
			    break;
			  }
                        }			
                      }

//	              cout << "filling correlation coefficient for " << $3 << " " << $4 << " " << $5 << endl;
		      yyMeasuredCorrelationMatrix.add(i, j, $5);
		  }
	      }
	    | input T_KEY value value
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += "\t";
		  sprintf(c, "%f", $3);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += " ";
		  sprintf(c, "%f", $4);
		  yyInputFileLine.prevalue += c;

		  if (!strcmp($2, "XScanRange")) {
		    yyXscanlow = $3;
//		    cout << "yyXscanlow = " << yyXscanlow << endl;
		    yyXscanhigh = $4;
//		    cout << "yyXscanhigh = " << yyXscanhigh << endl;
		  }
	      }
	    | input T_KEY T_SWITCHSTATE
	      {
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += "\t";
		  if ($3 == on) yyInputFileLine.prevalue += "on";
	  	  else yyInputFileLine.prevalue += "off";

		  if (!strcmp($2, "LoopCorrections")) {
		      if ($3 == on) yyUseLoopCorrections = true;
		      else yyUseLoopCorrections = false;

		  }
		  if (!strcmp($2, "CalcPullDist")) {
		      if ($3 == on) yyCalcPullDist = true;
		      else yyCalcPullDist = false;
		  }
		  if (!strcmp($2, "ScanParameters")) {
		      if ($3 == on) yyScanParameters = true;
		      else yyScanParameters = false;
		  }
		  if (!strcmp($2, "PerformFit")) {
		      if ($3 == on) yyPerformFit = true;
		      else yyPerformFit = false;
		  }
		  if (!strcmp($2, "CalcIndChisqContr")) {
		      if ($3 == on) yyCalcIndChisqContr = true;
		      else yyCalcIndChisqContr = false;
		  }
		  if (!strcmp($2, "ISR")) {
		      if ($3 == on) yyISR = true;
		      else yyISR = false;
		  }
		  if (!strcmp($2, "UseMinos")) {
		      if ($3 == on) yyUseMinos = true;
		      else yyUseMinos = false;
		  } 
		  //yyGetContours
		  if (!strcmp($2, "GetContours")) {
		      if ($3 == on) yyGetContours = true;
		      else yyGetContours = false;
		  } 		  
		  if (!strcmp($2, "UseHesse")) {
		      if ($3 == on) yyUseHesse = true;
		      else yyUseHesse = false;
		  }
		  if (!strcmp($2, "UseSimAnnBefore")) {
		    if ($3 == on) yyUseSimAnnBefore = true;
		      else yyUseSimAnnBefore = false;
		  }
		  if (!strcmp($2, "UseSimAnnWhile")) {
		      if ($3 == on) yyUseSimAnnWhile = true;
		      else yyUseSimAnnWhile = false;
		  }
		  if (!strcmp($2, "UseGivenStartValues")) {
		      if ($3 == on) yyUseGivenStartValues = true;
		      else yyUseGivenStartValues = false;
		  }
		  if (!strcmp($2, "BoundsOnX")) {
		      if ($3 == off) yyBoundsOnX = false;
		      else yyBoundsOnX = true;
		  }
		  if (!strcmp($2, "FitAllDirectly")) {
		      if ($3 == on) yyFitAllDirectly = true;
		      else yyFitAllDirectly = false;
		  }
		  if (!strcmp($2, "SepFitTanbX")) {
		      if ($3 == on) yySepFitTanbX = true;
		      else yySepFitTanbX = false;
		  }
		  if (!strcmp($2, "SepFitTanbMu")) {
		      if ($3 == on) yySepFitTanbMu = true;
		      else yySepFitTanbMu = false;
		  } 
		  // yySepFitmA
		  if (!strcmp($2, "SepFitmA")) {
		      if ($3 == on) yySepFitmA = true;
		      else yySepFitmA = false;
		  } 		  
		  if (!strcmp($2, "ScanX")) {
		      if ($3 == on) yyScanX = true;
		      else yyScanX = false;
		  }
		  if (!strcmp($2, "Verbose")) {
		      if ($3 == on) yyVerbose = true;
		      else yyVerbose = false;
		  } // yyAdaptiveSimAnn
		  if (!strcmp($2, "AdaptiveSimAnn")) {
		      if ($3 == on) yyAdaptiveSimAnn = true;
		      else yyAdaptiveSimAnn = false;
		  } // yyAdaptiveSimAnn yyNoBoundsAtAll
		  if (!strcmp($2, "NoBoundsAtAll")) {
		      if ($3 == on) yyNoBoundsAtAll = true;
		      else yyNoBoundsAtAll = false;
		  } // yySimAnnUncertainty
		  if (!strcmp($2, "SimAnnUncertainty")) {
		      if ($3 == on) yySimAnnUncertainty = true;
		      else yySimAnnUncertainty = false;
		  }
		  if (!strcmp($2, "SimAnnUncertaintyRunDown")) {
		      if ($3 == on) yySimAnnUncertaintyRunDown = true;
		      else yySimAnnUncertaintyRunDown = false;
		  }
		  if (!strcmp($2, "RandomDirUncertainties")) {
		      if ($3 == on) yyRandomDirUncertainties = true;
		      else yyRandomDirUncertainties = false;
		  } // yyPerformSingleFits
 		  if (!strcmp($2, "PerformSingleFits")) {
		      if ($3 == on) yyPerformSingleFits = true;
		      else yyPerformSingleFits = false;
		  } 
 		  if (!strcmp($2, "UseMicrOmegas")) {
		      if ($3 == on) yyUseMicrOmegas = true;
		      else yyUseMicrOmegas = false;
		  } 	      }
	    | input T_CALCULATOR T_WORD
	      {
                   yyInputFileLine.prevalue  = "Calculator";
		   yyInputFileLine.prevalue += "\t";
		   yyInputFileLine.prevalue += $3;

		   if (!strcmp($3, "SPHENO")) {
		      yyCalculator = SPHENO;
		   }
		   if (!strcmp($3, "SUSPECT")) {
		      yyCalculator = SUSPECT;
		   }
	      }
	    | input T_CALCULATOR T_WORD T_PATH
	      {
                   yyInputFileLine.prevalue  = "Calculator";
		   yyInputFileLine.prevalue += "\t";
		   yyInputFileLine.prevalue += $3;
		   yyInputFileLine.prevalue += " ";
                   yyInputFileLine.prevalue += $4;

		   if (!strcmp($3, "SPHENO")) {
		      yyCalculator = SPHENO;
		   }
		   if (!strcmp($3, "SUSPECT")) {
		      yyCalculator = SUSPECT;
		   }
	           yyCalculatorPath = $4;
	      }
	    | input T_KEY T_WORD
	      {
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += "\t";
		  yyInputFileLine.prevalue += $3;

		  if (!strcmp($2, "fitParameter")) {
		     parameter_t tmpparam;
		     tmpparam.name = $3;
		     yyFittedPar.push_back(tmpparam);
		  }  
		  else if (!strcmp($2,"fitModel")) {
		    cout<<"Fitting "<<$3<<" model"<<endl;
		    if (!strcmp($3, "MSSM"))
		      yyFitModel = MSSM;
		    else if (!strcmp($3, "mSUGRA"))
		      yyFitModel = mSUGRA;
		    else if (!strcmp($3, "GMSB"))
		      yyFitModel = GMSB;
		    else if (!strcmp($3, "AMSB"))
		      yyFitModel = AMSB;
		    else {
		      cerr<<"Unknown fit model: "<<$3<<endl;
		      exit(EXIT_FAILURE);
		    }
		  }
	      }
	    | input T_KEY T_WORD value
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += "\t";
		  yyInputFileLine.prevalue += $3;
		  yyInputFileLine.prevalue += "\t";
		  sprintf(c, "%f", $4);
		  yyInputFileLine.prevalue += c;

		  if (!strcmp($2, "fixParameter")) {
	            parameter_t tmpparam;
                    tmpparam.name = $3;
                    tmpparam.value = $4;
                    yyFixedPar.push_back(tmpparam);
		  }
		  else if (!strcmp($2, "fitParameter")) {
	            parameter_t tmpparam;
                    tmpparam.name = $3;
                    tmpparam.value = $4;
                    yyFittedPar.push_back(tmpparam);
		  }
                  else {
                     yyerror ("Syntax error in Fixed or Fitted Par 1");
                  }
	      }
	    | input T_KEY T_WORD value value T_NUMBER
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += "\t";
		  yyInputFileLine.prevalue += $3;
		  yyInputFileLine.prevalue += "\t";
		  sprintf(c, "%f", $4);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += "\t";
		  sprintf(c, "%f", $5);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += "\t";
		  sprintf(c, "%f", $6);
		  yyInputFileLine.prevalue += c;

		  if (!strcmp($2, "scanParameter")) {
		    if (yyScanPar.size() == 2) {
		       cout<<"At most two-dimensinal parameter scans can be performed"<<endl;
		       exit(EXIT_FAILURE);
                    }
	            ScanParameter tmpscanpar;
                    tmpscanpar.name = $3;
                    tmpscanpar.min = $4;
                    tmpscanpar.max = $5;
                    tmpscanpar.numberOfSteps = (unsigned int)$6;
		    if (tmpscanpar.numberOfSteps < 2) {
		       yyerror("Number of scan steps must be at least 2");
                    }
                    yyScanPar.push_back(tmpscanpar);
		  }
                  else {
                     yyerror ("Syntax error in scanParameter");
                  }
	      }
	    | input T_KEY T_WORD T_COMPARATOR value
	      {
		  char c[1000];
                  yyInputFileLine.prevalue  = $2;
                  yyInputFileLine.prevalue += "\t";
                  yyInputFileLine.prevalue += $3;
                  yyInputFileLine.prevalue += "\t";
                  yyInputFileLine.prevalue += $4;
		  sprintf(c, "%f", $5);
		  yyInputFileLine.prevalue += c;

		  if (!strcmp($2, "limit")) {
		      MeasuredValue tmpValue;
                      tmpValue.nofit = false;
		      tmpValue.theovalue  = 0;
		      tmpValue.name = $3;
		      tmpValue.error = -1.;
                      tmpValue.nofit = false;
		      if (strcmp($4, ">")) {
			  tmpValue.value = 1.2*$5;
			  tmpValue.bound_low = $5;
			  tmpValue.bound_up = 1.E+6;
		      } else {
			  tmpValue.value = 0.8*$5;
			  tmpValue.bound_low = -1.E+6;
			  tmpValue.bound_up = $5;
		      }
		      yyMeasuredVec.push_back(tmpValue);
		  }
	      }
	    | input T_NOFIT T_WORD value err
	      {
                  yyInputFileLine.prevalue  = "nofit";
                  yyInputFileLine.prevalue += "\t";
                  yyInputFileLine.prevalue += $3;
	          yyInputFileLine.value = $4;
//	          yyInputFileLine.error = $5;

		  MeasuredValue tmpValue;
		  tmpValue.theovalue  = 0;
		  tmpValue.name = $3;
		  tmpValue.value = $4;
		  tmpValue.error = $5;
		  tmpValue.id = 0;
	          tmpValue.bound_low = -1.E+6;
	          tmpValue.bound_up = 1.E+6;
                  tmpValue.nofit = true;
		  if (!strncmp($3, "mass", 4)) tmpValue.type = mass;
		  else if (!strncmp($3, "cos", 3)) tmpValue.type = other;
		  yyMeasuredVec.push_back(tmpValue);
              }
	    | input T_KEY T_WORD value T_ERRORSIGN value
	      {
		  if (!strcmp($2, "fitParameter")) {
		    char c[1000];
		    yyInputFileLine.prevalue  = $2;
		    yyInputFileLine.prevalue += "\t";
		    yyInputFileLine.prevalue += $3;
		    yyInputFileLine.prevalue += "\t";
		    sprintf(c, "%f", $4);
		    yyInputFileLine.prevalue += c;
		    yyInputFileLine.prevalue += "\t+-\t";
		    sprintf(c, "%f", $6);
		    yyInputFileLine.prevalue += c;
	            parameter_t tmpparam;
                    tmpparam.name = $3;
                    tmpparam.value = $4;
                    tmpparam.error = $6;
                    yyFittedPar.push_back(tmpparam);
		  }
/*
		  else if (!strcmp($2, "fixParameter")) {
		    char c[1000];
		    yyInputFileLine.prevalue  = $2;
		    yyInputFileLine.prevalue += "\t";
		    yyInputFileLine.prevalue += $3;
		    sprintf(c, "%f", $4);
		    yyInputFileLine.prevalue += c;
		    yyInputFileLine.prevalue += "\t";
		    sprintf(c, "%f", $6);
		    yyInputFileLine.prevalue += c;
	            parameter_t tmpparam;
                    tmpparam.name = $3;
                    tmpparam.value = $4;
                    yyFixedPar.push_back(tmpparam);
		  }
*/
                  else {
		     cout << "name = " << $3 << endl;
                     yyerror ("Syntax error in Fixed or Fitted Par 2");
                  }

	      }
	    | input T_KEY T_NUMBER sentence value err T_ALIAS T_NUMBER
	      {
		char c[1000];
		yyInputFileLine.prevalue  = $2;
		yyInputFileLine.prevalue += " ";
		sprintf(c, "%d", (int)$3);
		yyInputFileLine.prevalue += c;
		yyInputFileLine.prevalue += " ";
		yyInputFileLine.prevalue += $4;
                yyInputFileLine.value = $5;
//                yyInputFileLine.error = $6;
		yyInputFileLine.postvalue = "\talias ";
		sprintf(c, "%d", (int)$8);
		yyInputFileLine.postvalue += c;

		if (!strcmp($2, "edge")) {
		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.theovalue  = 0;
		  tmpValue.name = "edge type ";
		  char srtt[20];
		  sprintf(srtt,"%d",(int)$3);
		  tmpValue.name.append(srtt);
		  tmpValue.name.append(" alias ");
		  sprintf(srtt,"%d:",(int)$8);
		  tmpValue.name.append(srtt);
		  tmpValue.type  = Pedge;
		  tmpValue.id    = (int)$3;
		  if (tmpValue.id<1 || tmpValue.id>5) {
                    yyerror("edge type not implemented");
                  } 
		  string str;
	          str.erase();
		  str = $4;
//		  cout << "decompositing string "<<str<< endl;
		  char tmpstr2[255];
		  unsigned int pos = 0, newpos = 0;
//		  int anti = 1;
	          int i = 0;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		      str.copy(tmpstr2, newpos - pos, pos);
//	              cout << "tmpstr2 after str.copy "<<tmpstr2<< endl;
	              tmpstr2[newpos-pos] = '\0';
		      if (tmpstr2[newpos-pos-1] == '~') {
			  tmpstr2[newpos-pos-1] = '\0';
		      }
                      if (!strncmp(tmpstr2,"mass",4)) {
//			cout << "substracting mass" << endl;
                        while (tmpstr2[i] != '\0') {
                          tmpstr2[i] = tmpstr2[i+4];
			  i++;
                        }
                        i = 0;
                      }
		      else {
			  tmpstr2[newpos-pos] = '\0';
		      }
		      pos = newpos + 1;
//		      cout << "Adding Element: " << tmpstr2 << " with ID " <<  yyParticleIDs[tmpstr2] << endl;
		      tmpValue.daughters.push_back(yyParticleIDs[tmpstr2]);
		  }
		  tmpValue.value = $5;
		  tmpValue.error = $6;
		  tmpValue.alias = (int)$8;
		  tmpValue.bound_up = 1e+6;
		  tmpValue.bound_low = 0.;
		  yyMeasuredVec.push_back(tmpValue);
		  strcpy($4,"");
		}
		
              }
	     | input T_BR T_BRA T_WORD T_GOESTO sentence T_KET value err T_ALIAS T_NUMBER
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = "BR(";
		  yyInputFileLine.prevalue += $4;
		  yyInputFileLine.prevalue += " -> ";
		  yyInputFileLine.prevalue += $6;
		  yyInputFileLine.prevalue += ")";
		  yyInputFileLine.value = $8;
//		  yyInputFileLine.error = $9;
		  yyInputFileLine.postvalue = "\talias ";
		  sprintf(c, "%d", (int)$11);
		  yyInputFileLine.postvalue += c;

		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.type  = br;
		  tmpValue.theovalue  = 0;
		  tmpValue.name  = $4;
		  tmpValue.id    = yyParticleIDs[$4];
		  tmpValue.name.append(" -> ");
		  tmpValue.name.append($6);
		  string str;
	          str.erase();
		  str = $6;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		      str.copy(tmpstr3, newpos - pos, pos);
		      if (tmpstr3[newpos-pos-1] == '~') {
			  tmpstr3[newpos-pos-1] = '\0';
			  anti = -1;
		      }
		      else {
			  tmpstr3[newpos-pos] = '\0';
			  anti = 1;
		      }
		      pos = newpos + 1;
		      tmpValue.daughters.push_back(anti * yyParticleIDs[tmpstr3]);
		  }
		  tmpValue.value = $8;
		  tmpValue.error = $9;
		  tmpValue.alias = (int)$11;
		  tmpValue.bound_up = 1e+6;
		  tmpValue.bound_low = 0.;
		  yyMeasuredVec.push_back(tmpValue);
//		  cout << "Added branching ratio " << tmpValue.name << endl;
	          strcpy($6,"");
	      }
	     | input T_LEO T_BRA T_WORD T_KET value err T_ALIAS T_NUMBER
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = "LEObs ( ";
		  yyInputFileLine.prevalue += $4;
		  yyInputFileLine.prevalue += " ) ";
		  yyInputFileLine.value = $6;
//		  yyInputFileLine.error = $7;
		  yyInputFileLine.postvalue = "\talias ";
		  sprintf(c, "%d", (int)$9);
		  yyInputFileLine.postvalue += c;

		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.type  = LEObs;
		  tmpValue.theovalue  = 0;
		  tmpValue.name  = $4;
                  //=================================
		  if (!strcmp($4, "bsg")) {
		    tmpValue.id    = bsg;
		  } else if (!strcmp($4, "bsmm")) {
		    tmpValue.id    = bsmm;
		  } else if (!strcmp($4, "B_smm")) {
		    tmpValue.id    = B_smm;
		  } else if (!strcmp($4, "B_utn")) {
		    tmpValue.id    = B_utn;
		  } else if (!strcmp($4, "dMB_d")) {
		    tmpValue.id    = dMB_d;
		  } else if (!strcmp($4, "dMB_s")) {
		    tmpValue.id    = dMB_s;
		  } else if (!strcmp($4, "gmin2e")) {
		    tmpValue.id    = gmin2e;
		  } else if (!strcmp($4, "gmin2m")) {
		    tmpValue.id    = gmin2m;
		  } else if (!strcmp($4, "gmin2t")) {
		    tmpValue.id    = gmin2t;
		  } else if (!strcmp($4, "drho")) {
		    tmpValue.id    = drho;
		  } else if (!strcmp($4, "omega")) {
		    tmpValue.id    = omega;
		  }
                  //=================================
		  tmpValue.value = $6;
		  tmpValue.error = $7;
		  tmpValue.alias = (int)$9;
		  tmpValue.bound_up = 1e+6;
		  tmpValue.bound_low = -1e+6;
		  yyMeasuredVec.push_back(tmpValue);
		  //		  cout << "Added le obs " << tmpValue.name << " = " << tmpValue.value << endl;
	      }
	     | input T_LEO T_BRA T_WORD T_KET T_COMPARATOR value T_ALIAS T_NUMBER
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = "bound on LEObs ( ";
		  yyInputFileLine.prevalue += $4;
		  yyInputFileLine.prevalue += " ) ";
		  yyInputFileLine.value = $7;
		  //		  yyInputFileLine.error = 0.01*$7;
		  yyInputFileLine.postvalue = "\talias ";
		  sprintf(c, "%d", (int)$9);
		  yyInputFileLine.postvalue += c;
		  
		  MeasuredValue tmpValue;
		  tmpValue.nofit = false;
		  tmpValue.type  = LEObs;
		  tmpValue.theovalue  = 0;
		  tmpValue.name  = $4;
                  //==================================
		  if (!strcmp($4, "bsg")) {
		    tmpValue.id    = bsg;
		  } else if (!strcmp($4, "bsmm")) {
		    tmpValue.id    = bsmm;
		  } else if (!strcmp($4, "B_smm")) {
		    tmpValue.id    = B_smm;
		  } else if (!strcmp($4, "B_utn")) {
		    tmpValue.id    = B_utn;
		  } else if (!strcmp($4, "dMB_d")) {
		    tmpValue.id    = dMB_d;
		  } else if (!strcmp($4, "dMB_s")) {
		    tmpValue.id    = dMB_s;
		  } else if (!strcmp($4, "gmin2e")) {
		    tmpValue.id    = gmin2e;
		  } else if (!strcmp($4, "gmin2m")) {
		    tmpValue.id    = gmin2m;
		  } else if (!strcmp($4, "gmin2t")) {
		    tmpValue.id    = gmin2t;
		  } else if (!strcmp($4, "drho")) {
		    tmpValue.id    = drho;
		  } else if (!strcmp($4, "omega")) {
		    tmpValue.id    = omega;
		  }
                  //==================================
		  cout << "T_COMPARATOR " << $4 << " " << $6 << " " << $7 << endl;
		  if (!strcmp($6, ">")) {
		    tmpValue.value = 1.2*$7;
		    tmpValue.error = -1.;
		    tmpValue.bound_up = 1e+6;
		    tmpValue.bound_low = $7;
		  }
		  else if (!strcmp($6, "<")) {
		    tmpValue.value = 0.8*$7;
		    tmpValue.error = -1.;
		    tmpValue.bound_up = $7;
		    tmpValue.bound_low = -1e+6;
		  } else {
		    yyerror ("Syntax Error in LEObs");
		  }
		  tmpValue.alias = (int)$9;
		  yyMeasuredVec.push_back(tmpValue);
		  //		  cout << "Added le obs " << tmpValue.name << " = " << tmpValue.value << endl;
	      }
	     | input T_BR T_BRA T_WORD T_GOESTO sentence T_KET T_ALIAS T_NUMBER
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = "BR ( ";
		  yyInputFileLine.prevalue += $4;
		  yyInputFileLine.prevalue += " -> ";
		  yyInputFileLine.prevalue += $6;
		  yyInputFileLine.prevalue += " ) ";
		  yyInputFileLine.prevalue += "\talias ";
		  sprintf(c, "%d", (int)$9);
		  yyInputFileLine.prevalue += c;

		  MeasuredValue tmpValue;
                  tmpValue.nofit = true;
		  tmpValue.type  = br;
		  tmpValue.theovalue  = -1;
		  tmpValue.name  = $4;
		  tmpValue.id    = yyParticleIDs[$4];
		  tmpValue.name.append(" -> ");
		  tmpValue.name.append($6);
		  string str;
	          str.erase();
		  str = $6;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		      str.copy(tmpstr3, newpos - pos, pos);
		      if (tmpstr3[newpos-pos-1] == '~') {
			  tmpstr3[newpos-pos-1] = '\0';
			  anti = -1;
		      }
		      else {
			  tmpstr3[newpos-pos] = '\0';
			  anti = 1;
		      }
		      pos = newpos + 1;
		      tmpValue.daughters.push_back(anti * yyParticleIDs[tmpstr3]);
		  }
		  tmpValue.value = -1;
		  tmpValue.error = -1;
		  tmpValue.alias = (int)$9;
		  tmpValue.bound_up = 1e+6;
		  tmpValue.bound_low = 0.;
		  yyMeasuredVec.push_back(tmpValue);
//		  cout << "Added branching ratio alias for " << tmpValue.name << endl;
	          strcpy($6,"");
	      }
	     | input T_KEY T_BRA sentence T_GOESTO sentence T_COMMA value T_COMMA value T_COMMA value T_KET value err T_ALIAS T_NUMBER
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += " ( ";
		  yyInputFileLine.prevalue += $4;
		  yyInputFileLine.prevalue += " -> ";
		  yyInputFileLine.prevalue += $6;
		  yyInputFileLine.prevalue += ", ";
		  sprintf(c, "%f GeV", $8);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += ", ";
		  sprintf(c, "%f", $10);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += ", ";
		  sprintf(c, "%f", $12);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += " ) ";
		  yyInputFileLine.value = $14;
//		  yyInputFileLine.error = $15;
		  yyInputFileLine.postvalue  = "\talias ";
		  sprintf(c, "%d", (int)$17);
		  yyInputFileLine.postvalue += c;

		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.name.erase();
		  tmpValue.type  = xsection;
		  tmpValue.theovalue  = 0;
		  tmpValue.name  = $4;
//	          cout << "in xs: $4 = " << $4 << endl;
		  tmpValue.id    = 11;
		  tmpValue.name.append(" -> ");
		  tmpValue.name.append($6);
		  string str;
	          str.erase();
	          str = $6;
		  char tmpstr4[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		      str.copy(tmpstr4, newpos - pos, pos);
		      if (tmpstr4[newpos-pos-1] == '~') {
			  tmpstr4[newpos-pos-1] = '\0';
			  anti = -1;
		      }
		      else {
			  tmpstr4[newpos-pos] = '\0';
			  anti = 1;
		      }
		      pos = newpos + 1;
		      tmpValue.products.push_back(anti * yyParticleIDs[tmpstr4]);
		  }
		  tmpValue.sqrts = $8;
		  tmpValue.polarisation1 = $10;
		  tmpValue.polarisation2 = $12;
		  tmpValue.value = $14;
		  tmpValue.error = $15;
		  tmpValue.alias = (int)$17;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
		  yyMeasuredVec.push_back(tmpValue);
//		  cout << "Added cross section " << tmpValue.name << endl;
		  tmpStrings.clear();
	          tmpNumbers.clear();
		  strcpy($4,"                     ");
	          strcpy($4,"");
		  strcpy($6,"                     ");
	          strcpy($6,"");
	      }
	     | input T_KEY T_BRA sentence T_GOESTO sentence T_COMMA value T_COMMA value T_COMMA value T_KET T_ALIAS T_NUMBER
	      {
		  char c[1000];
		  yyInputFileLine.prevalue  = $2;
		  yyInputFileLine.prevalue += " ( ";
		  yyInputFileLine.prevalue += $4;
		  yyInputFileLine.prevalue += " -> ";
		  yyInputFileLine.prevalue += $6;
		  yyInputFileLine.prevalue += ", ";
		  sprintf(c, "%f GeV", $8);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += ", ";
		  sprintf(c, "%f", $10);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += ", ";
		  sprintf(c, "%f", $12);
		  yyInputFileLine.prevalue += c;
		  yyInputFileLine.prevalue += " ) ";
		  yyInputFileLine.prevalue += "\talias ";
		  sprintf(c, "%d", (int)$15);
		  yyInputFileLine.prevalue += c;

		  MeasuredValue tmpValue;
                  tmpValue.nofit = true;
		  tmpValue.name.erase();
		  tmpValue.type  = xsection;
		  tmpValue.theovalue  = -1;
		  tmpValue.name  = $4;
//	          cout << "in xs: $4 = " << $4 << endl;
		  tmpValue.id    = 11;
		  tmpValue.name.append(" -> ");
		  tmpValue.name.append($6);
		  string str;
	          str.erase();
	          str = $6;
		  char tmpstr4[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		      str.copy(tmpstr4, newpos - pos, pos);
		      if (tmpstr4[newpos-pos-1] == '~') {
			  tmpstr4[newpos-pos-1] = '\0';
			  anti = -1;
		      }
		      else {
			  tmpstr4[newpos-pos] = '\0';
			  anti = 1;
		      }
		      pos = newpos + 1;
		      tmpValue.products.push_back(anti * yyParticleIDs[tmpstr4]);
		  }
		  tmpValue.sqrts = $8;
		  tmpValue.polarisation1 = $10;
		  tmpValue.polarisation2 = $12;
		  tmpValue.value = -1;
		  tmpValue.error = -1;
		  tmpValue.alias = (int)$15;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
		  yyMeasuredVec.push_back(tmpValue);
//		  cout << "Added cross section alias for " << tmpValue.name << endl;
		  tmpStrings.clear();
	          tmpNumbers.clear();
		  strcpy($4,"                     ");
	          strcpy($4,"");
		  strcpy($6,"                     ");
	          strcpy($6,"");
	      }
	     | input T_KEY T_BRA sentence T_KET value err T_ALIAS T_NUMBER
	      {
		char c[1000];
		yyInputFileLine.prevalue  = $2;
		yyInputFileLine.prevalue += " ( ";
		yyInputFileLine.prevalue += $4;
		yyInputFileLine.prevalue += " ) ";
		yyInputFileLine.value = $6;
//		yyInputFileLine.error = $7;
		yyInputFileLine.postvalue  = "\talias ";
		sprintf(c, "%d", (int)$9);
		yyInputFileLine.postvalue += c;

		if (!strcmp($2,"xsbr")) {
		  int i = 0;
//		  int j = 0;
		  int aliasnumber;
		  char* charnumber;
		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.name.erase();
		  tmpValue.type  = xsbr;
		  tmpValue.theovalue  = 0;
		  tmpValue.value = $6;
		  tmpValue.error = $7;
		  tmpValue.name  = $4;
		  tmpValue.alias = (int)$9;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
//		  cout << "xsbr found" << endl;
		  // break sentence in pieces:
		  string str;
	          str.erase();
		  str = $4;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  bool found_br = false;
		  int anti = 1;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		    str.copy(tmpstr3, newpos - pos, pos);
		    if (tmpstr3[newpos-pos-1] == '~') {
		      tmpstr3[newpos-pos-1] = '\0';
		      anti = -1;
		    }
		    else {
		      tmpstr3[newpos-pos] = '\0';
		      anti = 1;
		    }
		    pos = newpos + 1;
		    // break up:
		    if (!strncmp(tmpstr3,"sigma",5)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == xsection) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found xs " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "sigma " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove cross section from the list
//		      cout << "removing xs from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.sqrts = yyMeasuredVec[i].sqrts;
		      tmpValue.products.push_back(i);
		    }
		    else if (!strncmp(tmpstr3,"brsum",5)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == brsum) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found brsum " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "brsum " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else if (!strncmp(tmpstr3,"brprod",6)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == brprod) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found brprod " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "brprod " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else if (!strncmp(tmpstr3,"br",2)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found br " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "br " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		  }
		  yyMeasuredVec.push_back(tmpValue);
		}
		else if (!strcmp($2,"brratio")) {
		  int i = 0;
//		  int j = 0;
		  int aliasnumber;
		  char* charnumber;
		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.name.erase();
		  tmpValue.type  = brratio;
		  tmpValue.theovalue  = 0;
		  tmpValue.value = $6;
		  tmpValue.error = $7;
		  tmpValue.name  = $4;
		  tmpValue.alias = (int)$9;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
//		  cout << "brratio found" << endl;
		  // break sentence in pieces:
		  string str;
	          str.erase();
		  str = $4;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  bool found_br = false;
		  int countbrs = 0;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		    countbrs++;
		    if (countbrs>=3) {
		      cout << "syntax error in brratio: too many br, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");		      
		    }
		    str.copy(tmpstr3, newpos - pos, pos);
		    if (tmpstr3[newpos-pos-1] == '~') {
		      tmpstr3[newpos-pos-1] = '\0';
		      anti = -1;
		    }
		    else {
		      tmpstr3[newpos-pos] = '\0';
		      anti = 1;
		    }
		    pos = newpos + 1;
		    // break up:
		    if (!strncmp(tmpstr3,"brsum",5)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == brsum) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found brsum " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "brsum " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else if (!strncmp(tmpstr3,"brprod",6)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == brprod) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found brprod " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "brprod " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else if (!strncmp(tmpstr3,"br",2)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found br " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "br " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else {
		      cout << "syntax error in brratio, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");
		    }
		  }
		  if (countbrs<2) {
		      cout << "syntax error in brratio: not enough br, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");		      		    
		  }
		  yyMeasuredVec.push_back(tmpValue);
		}
		else if (!strcmp($2,"brsum")) {
		  int i = 0;
//		  int j = 0;
		  int aliasnumber;
		  char* charnumber;
		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.name.erase();
		  tmpValue.type  = brsum;
		  tmpValue.theovalue  = 0;
		  tmpValue.value = $6;
		  tmpValue.error = $7;
		  tmpValue.name  = $4;
		  tmpValue.alias = (int)$9;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
//		  cout << "brratio found" << endl;
		  // break sentence in pieces:
		  string str;
	          str.erase();
		  str = $4;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  bool found_br = false;
		  int countbrs = 0;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		    countbrs++;
		    str.copy(tmpstr3, newpos - pos, pos);
		    if (tmpstr3[newpos-pos-1] == '~') {
		      tmpstr3[newpos-pos-1] = '\0';
		      anti = -1;
		    }
		    else {
		      tmpstr3[newpos-pos] = '\0';
		      anti = 1;
		    }
		    pos = newpos + 1;
		    // break up:
		    if (!strncmp(tmpstr3,"br",2)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found br " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "br " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else {
		      cout << "syntax error in brsum, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");
		    }
		  }
		  if (countbrs<2) {
		      cout << "syntax error in brsum: not enough br, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");		      		    
		  }
		  yyMeasuredVec.push_back(tmpValue);
		}
		else if (!strcmp($2,"brprod")) {
		  int i = 0;
//		  int j = 0;
		  int aliasnumber;
		  char* charnumber;
		  MeasuredValue tmpValue;
                  tmpValue.nofit = false;
		  tmpValue.name.erase();
		  tmpValue.type  = brprod;
		  tmpValue.theovalue  = 0;
		  tmpValue.value = $6;
		  tmpValue.error = $7;
		  tmpValue.name  = $4;
		  tmpValue.alias = (int)$9;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
//		  cout << "brprod found" << endl;
		  // break sentence in pieces:
		  string str;
	          str.erase();
		  str = $4;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  bool found_br = false;
		  int countbrs = 0;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		    countbrs++;
		    str.copy(tmpstr3, newpos - pos, pos);
		    if (tmpstr3[newpos-pos-1] == '~') {
		      tmpstr3[newpos-pos-1] = '\0';
		      anti = -1;
		    }
		    else {
		      tmpstr3[newpos-pos] = '\0';
		      anti = 1;
		    }
		    pos = newpos + 1;
		    // break up:
		    if (!strncmp(tmpstr3,"br",2)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found br " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "br " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else {
		      cout << "syntax error in brprod, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");
		    }
		  }
		  if (countbrs<2) {
		      cout << "syntax error in brprod: not enough br, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");		      		    
		  }
		  yyMeasuredVec.push_back(tmpValue);
		}
	      }
	     | input T_KEY T_BRA sentence T_KET T_ALIAS T_NUMBER
	      {
		char c[1000];
		yyInputFileLine.prevalue  = $2;
		yyInputFileLine.prevalue += " ( ";
		yyInputFileLine.prevalue += $4;
		yyInputFileLine.prevalue += " ) ";
		yyInputFileLine.prevalue  = "\talias ";
		sprintf(c, "%d", (int)$7);
		yyInputFileLine.prevalue += c;

		if (!strcmp($2,"brsum")) {
		  int i = 0;
//		  int j = 0;
		  int aliasnumber;
		  char* charnumber;
		  MeasuredValue tmpValue;
                  tmpValue.nofit = true;
		  tmpValue.name.erase();
		  tmpValue.type  = brsum;
		  tmpValue.theovalue  = 0;
		  tmpValue.value = -1;
		  tmpValue.error = -1;
		  tmpValue.name  = $4;
		  tmpValue.alias = (int)$7;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
                  //cout << "alias for brratio found" << endl;
		  // break sentence in pieces:
		  string str;
	          str.erase();
		  str = $4;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  bool found_br = false;
		  int countbrs = 0;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		    countbrs++;
		    str.copy(tmpstr3, newpos - pos, pos);
		    if (tmpstr3[newpos-pos-1] == '~') {
		      tmpstr3[newpos-pos-1] = '\0';
		      anti = -1;
		    }
		    else {
		      tmpstr3[newpos-pos] = '\0';
		      anti = 1;
		    }
		    pos = newpos + 1;
		    // break up:
		    if (!strncmp(tmpstr3,"br",2)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found br " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "br " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else {
		      cout << "syntax error in brsum, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");
		    }
		  }
		  if (countbrs<2) {
		      cout << "syntax error in brsum: not enough br, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");		      		    
		  }
		  yyMeasuredVec.push_back(tmpValue);
		}
		else if (!strcmp($2,"brprod")) {
		  int i = 0;
//		  int j = 0;
		  int aliasnumber;
		  char* charnumber;
		  MeasuredValue tmpValue;
                  tmpValue.nofit = true;
		  tmpValue.name.erase();
		  tmpValue.type  = brprod;
		  tmpValue.theovalue  = 0;
		  tmpValue.value = -1;
		  tmpValue.error = -1;
		  tmpValue.name  = $4;
		  tmpValue.alias = (int)$7;
		  tmpValue.bound_up = 1e+9;
		  tmpValue.bound_low = 0.;
                  //cout << "alias for brratio found" << endl;
		  // break sentence in pieces:
		  string str;
	          str.erase();
		  str = $4;
		  char tmpstr3[255];
		  unsigned int pos = 0, newpos = 0;
		  int anti = 1;
		  bool found_br = false;
		  int countbrs = 0;
		  while ((newpos = str.find(" ", pos)) != string::npos) {
		    countbrs++;
		    str.copy(tmpstr3, newpos - pos, pos);
		    if (tmpstr3[newpos-pos-1] == '~') {
		      tmpstr3[newpos-pos-1] = '\0';
		      anti = -1;
		    }
		    else {
		      tmpstr3[newpos-pos] = '\0';
		      anti = 1;
		    }
		    pos = newpos + 1;
		    // break up:
		    if (!strncmp(tmpstr3,"br",2)) {
		      charnumber = strchr(tmpstr3,'_');
		      if (charnumber == 0) yyerror ("Underscore not found");
		      aliasnumber = atoi((charnumber+1));
		      found_br = false;
		      for (unsigned int k = 0; k<yyMeasuredVec.size();k++) {
			if ((yyMeasuredVec[k].type == br) && (yyMeasuredVec[k].alias == aliasnumber)) {
			  i = k;
			  found_br = true;
//			  cout << "found br " <<  yyMeasuredVec[k].name << endl;
			  break;
			}
		      }
		      if (!found_br) {
			cout << "br " << aliasnumber << " not found, alias: " << tmpValue.alias <<  endl;
			yyerror (" ");			
		      }
		      // remove br from the list
//		      cout << "removing br from the list" << endl;
		      yyMeasuredVec[i].nofit = true;
		      tmpValue.daughters.push_back(i);
		    }
		    else {
		      cout << "syntax error in brprod, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");
		    }
		  }
		  if (countbrs<2) {
		      cout << "syntax error in brprod: not enough br, alias: " << tmpValue.alias <<  endl;
		      yyerror (" ");		      		    
		  }
		  yyMeasuredVec.push_back(tmpValue);
		}
		else {
		  cout<<"No value and error given for "<<$2<<" ( "<<$4<<" )"<<endl;
		  yyerror("");
                }
	      }
            ;

//========================================================================
//========================================================================
//========================================================================
//========================================================================
block:      T_BLOCK T_WORD T_NEWLINE parameters
              {
//                  printf("Reading block %s...\n", $2);

//========================================================================
                  if (!strcmp($2, "MODSEL")) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
            
                          switch (int(tmpParams[i][0])) {
                          case 1:
                              yyModel = (int)tmpParams[i][1];
                              break;
			  default:
                              yyerror("Parameter in MODSEL not known");
                          }
                      }
                  }
//========================================================================
                  if (!strcmp($2, "SMINPUTS")) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
            
                          switch (int(tmpParams[i][0])) {
                          case 1:
                              yyoneoveralpha_em_mz =  tmpParams[i][1];
                              break;
                          case 2:
                              yyG_F =  tmpParams[i][1];
                              break;
                          case 3:
                              yyalpha_s_mz =  tmpParams[i][1];
                              break;
                          case 4:
                              yyMassZ = tmpParams[i][1];
			      break;
                          case 5:
                              yyMass_b = tmpParams[i][1];
			      break;
                          case 6:
                              yyMass_t = tmpParams[i][1];
			      break;
                          case 7:
                              yyMass_tau = tmpParams[i][1];
			      break;
			  default:
                              yyerror("Parameter in SMINPUTS not known");
                          }
                      }
                  }
//========================================================================
                  if (!strcmp($2, "MINPAR")) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
            		  if (int(tmpParams[i][0])==3) {
                              yytanb =  tmpParams[i][1];
                              continue;
                          }
                          if (yyModel==1) {              // mSUGRA
                              switch (int(tmpParams[i][0])) {
                              case 0:
                                  yy_Qmax =  tmpParams[i][1];
                                  break;
                              case 1:
                                  yy_m0 =  tmpParams[i][1];
                                  break;
                              case 2:
                                  yy_m12 =  tmpParams[i][1];
                                  break;
                              case 4:
                                  yy_signmu =  tmpParams[i][1];
                                  break;
                              case 5:
                                  yy_A =  tmpParams[i][1];
                                  break;
			      default:
                                  yyerror("Parameter in MINPAR not known");
                              }
                          }
                      }
                  }
//========================================================================
                  if (!strcmp($2, "SPINFO")) {
//		      cout << "tmpStrings.size()" << tmpStrings.size() << endl;
                      for (unsigned int i=0; i<tmpStrings.size(); i++) {
            
                          switch (int(tmpNumbers[i])) {
                          case 1:
                              strcpy(yy_spectrum_calc_name,tmpStrings[i].tmpname);
                              break;
                          case 2:
                              strcpy(yy_spectrum_calc_version,tmpStrings[i].tmpname);
                              break;
                          case 3:
                              strcat(yy_warnings,tmpStrings[i].tmpname);
                              break;
                          case 4:
                              strcat(yy_errors,tmpStrings[i].tmpname);
                              break;
                          case 5:
                              strcpy(yy_error_codes,tmpStrings[i].tmpname);
                              break;
			  default:
//			      cout << "tmpNumbers["<<i<<"]: "<< int(tmpNumbers[i]) << endl;
			    return 1;
                              yyerror("Parameter in SPINFO not known");
                          }
                      }
		      tmpNumbers.clear();
                      tmpStrings.clear();
//		      cout << "after clear: tmpStrings.size()" << tmpStrings.size() << endl;
                  }
//========================================================================
                  if (!strcmp($2, "SPhenoINFO")) {
		    // cout << "reading SPhenoINFO " << endl;
		  }
//========================================================================
                  if (!strcmp($2, "EXTPAR")) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
            
                          switch (int(tmpParams[i][0])) {
                          case 0:
                              yyMinput =  tmpParams[i][1];
                              break;
                          case 1:
                              yyM1 =  tmpParams[i][1];
                              break;
                          case 2:
                              yyM2 =  tmpParams[i][1];
                              break;
                          case 3:
                              yyM3 =  tmpParams[i][1];
                              break;
                          case 11:
                              yy_At =  tmpParams[i][1];
                              break;
                          case 12:
                              yy_Ab =  tmpParams[i][1];
                              break;
                          case 13:
                              yy_Atau =  tmpParams[i][1];
                              break;
                          case 21:
                              yy_m2Hd =  tmpParams[i][1];
                              break;
                          case 22:
                              yy_m2Hu =  tmpParams[i][1];
                              break;
                          case 23:
                              yy_mu =  tmpParams[i][1];
                              break;
                          case 24:
                              yy_MassA =  tmpParams[i][1];
                              break;
                          case 31:
                              yy_mse1L =  tmpParams[i][1];
                              break;
                          case 32:
                              yy_mse2L =  tmpParams[i][1];
                              break;
                          case 33:
                              yy_mse3L =  tmpParams[i][1];
                              break;
                          case 34:
                              yy_mse1R =  tmpParams[i][1];
                              break;
                          case 35:
                              yy_mse2R =  tmpParams[i][1];
                              break;
                          case 36:
                              yy_mse3R =  tmpParams[i][1];
                              break;
                          case 41:
                              yy_msq1L =  tmpParams[i][1];
                              break;
                          case 42:
                              yy_msq2L =  tmpParams[i][1];
                              break;
                          case 43:
                              yy_msq3L =  tmpParams[i][1];
                              break;
                          case 44:
                              yy_msu1R =  tmpParams[i][1];
                              break;
                          case 45:
                              yy_msu2R =  tmpParams[i][1];
                              break;
                          case 46:
                              yy_msu3R =  tmpParams[i][1];
                              break;
                          case 47:
                              yy_msd1R =  tmpParams[i][1];
                              break;
                          case 48:
                              yy_msd2R =  tmpParams[i][1];
                              break;
                          case 49:
                              yy_msd3R =  tmpParams[i][1];
                              break;
                          case 51:
                              yy_N51 =  tmpParams[i][1];
                              break;
                          case 52:
                              yy_N52 =  tmpParams[i][1];
                              break;
                          case 53:
                              yy_N53 =  tmpParams[i][1];
                              break;
			  
			  default:
                              yyerror("Parameter in EXTPAR not known");
                          }
                      }
                  }
//========================================================================
                  if (!strcmp($2, "MASS")) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
     			  yyMass[(unsigned int)tmpParams[i][0]]=tmpParams[i][1];     
                      }
                  }

		  // SPhenoLowEnergy 
//========================================================================
                  //==========================================
                  if (!strcmp($2, "SPhenoLowEnergy")) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
			if ((unsigned int)tmpParams[i][0]==1) {
     			  yybsg=tmpParams[i][1];
			}
			else if ((unsigned int)tmpParams[i][0]==2) {
     			  yybsmm=tmpParams[i][1]*1E6;
			}
			else if ((unsigned int)tmpParams[i][0]==3) {
     			  yyB_smm=tmpParams[i][1]*1E9;
			}
			else if ((unsigned int)tmpParams[i][0]==4) {
     			  yyB_utn=tmpParams[i][1];
			}
			else if ((unsigned int)tmpParams[i][0]==5) {
     			  yydMB_d=tmpParams[i][1];
			}
			else if ((unsigned int)tmpParams[i][0]==6) {
     			  yydMB_s=tmpParams[i][1];
			}
			else if ((unsigned int)tmpParams[i][0]==10) {
     			  yygmin2e=tmpParams[i][1]*1E12;
			}
			else if ((unsigned int)tmpParams[i][0]==11) {
     			  yygmin2m=tmpParams[i][1]*1E9;
			}
			else if ((unsigned int)tmpParams[i][0]==12) {
     			  yygmin2t=tmpParams[i][1]*1E9;
			}
			else if ((unsigned int)tmpParams[i][0]==30) {
     			  yydrho=tmpParams[i][1];
			}
                      }
                  }
                  //==========================================
		  // MicrOmegas 
//========================================================================
                  if (!strcmp($2, "MicrOmegas")) {
		    for (unsigned int i=0; i<tmpParams.size(); i++) {
		      if ((unsigned int)tmpParams[i][0]==1) {
			yyOmega=tmpParams[i][1];
		      }
		    }
                  }
//========================================================================
                  if (!strcmp($2, "SPhenoCrossSections")) {
		    // cout << "starting to read XS..." << endl;
                  }

                  tmpParams.clear();
              }

//========================================================================
//========================================================================
            | T_BLOCK T_WORD T_SCALE T_NUMBER T_NEWLINE parameters
              {
                  yyScale = $4;
//                  printf("Reading block %s at scale %f...\n", $2, $4);

                  tmpParams.clear();
              }
            ;

//========================================================================
//========================================================================
//========================================================================
//========================================================================
decay:      T_DECAY T_NUMBER T_NUMBER T_NEWLINE parameters
              {
//                  printf("Reading decay %i...\n", (int)$2);
                  tmp_branch.decays.clear();
//                  cout << "branch cleared, tmpParams.size() = "<<tmpParams.size() << endl;
		  if (tmpParams.size()>0) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
//                         cout << "loop "<<i<<endl;
                         tmpVec.clear();
                         tmpVec.push_back(tmpParams[i][0]);
                         tmpVec.push_back(tmpParams[i][1]);
//                         cout << "filled BR and no of daughters" << endl;
        	         for (unsigned int j=2; j<2+tmpParams[i][1]; j++) {
//                             cout << "tmpParams["<<i<<"]["<<j<<"] = "<<tmpParams[i][j]<<endl;
    		             tmpVec.push_back(tmpParams[i][j]);
                         }
                         tmp_branch.decays.push_back(tmpVec);
                      }
    		  tmp_branch.TWidth = $3;
    		  branching_ratios[(int)$2] = tmp_branch;
		  } else {
    		      tmp_branch.TWidth = 0.;	
		      tmpVec.clear();
                      tmpVec.push_back(0.);
                      tmpVec.push_back(0.);
                      tmp_branch.decays.push_back(tmpVec);
       		      branching_ratios[(int)$2] = tmp_branch;
                  }
		  tmpParams.clear();
              }
            ;	

//========================================================================
//========================================================================
//========================================================================
//========================================================================
xs:         T_XS T_NUMBER T_NUMBER T_NUMBER T_NUMBER T_NUMBER T_NUMBER T_NEWLINE parameters
              {
		// printf("Reading xs %i...\n", (int)$2);

                  tmp_xs.xs.clear();
//                  cout << "branch cleared, tmpParams.size() = "<<tmpParams.size() << endl;
		  if (tmpParams.size()>0) {
                      for (unsigned int i=0; i<tmpParams.size(); i++) {
//                         cout << "loop "<<i<<endl;
                         tmpVec.clear();
                         tmpVec.push_back(tmpParams[i][0]);
                         tmpVec.push_back(tmpParams[i][1]);
//                         cout << "filled XS and no of primaries" << endl;
	                 for (unsigned int j=2; j<2+tmpParams[i][1]; j++) {
//                             cout << "tmpParams["<<i<<"]["<<j<<"] = "<<tmpParams[i][j]<<endl;
		             tmpVec.push_back(tmpParams[i][j]);
                         }
                         tmp_xs.xs.push_back(tmpVec);
                      }
		      tmp_xs.initials[0] = (int)$2;
		      tmp_xs.initials[1] = (int)$3;
		      tmp_xs.ecm = $4;
		      tmp_xs.polarisation[0] = $5;
		      tmp_xs.polarisation[1] = $6;
		      tmp_xs.isr = (int)$7;
		      tmp2Vec.clear();
		      tmp2Vec.push_back($2);
		      tmp2Vec.push_back($4);
		      tmp2Vec.push_back($5);
		      tmp2Vec.push_back($6);
		      cross_sections[tmp2Vec] = tmp_xs;
//                      cout << "filling xs for " << $2 << " " << $5 << " " << $6 << endl;
//		      for (unsigned int k=0; k<tmp_xs.xs.size();k++ ) {
//		       cout <<  " k = " << k << " products " << tmp_xs.xs[k][2] << 
//                              " " << tmp_xs.xs[k][3] << endl;
//                      }
		  } else {
		      tmp_xs.initials[0] = (int)$2;
		      tmp_xs.initials[1] = (int)$3;
		      tmp_xs.ecm = $4;
		      tmp_xs.polarisation[0] = $5;
		      tmp_xs.polarisation[1] = $6;
		      tmp_xs.isr = (int)$7;
		      tmpVec.clear();
                      tmpVec.push_back(0.);
                      tmpVec.push_back(0.);
                      tmp_xs.xs.push_back(tmpVec);
		      tmp2Vec.clear();
		      tmp2Vec.push_back($2);
		      tmp2Vec.push_back($5);
		      tmp2Vec.push_back($6);
		      cross_sections[tmp2Vec] = tmp_xs;
                  }

		  tmpParams.clear();
              }
            ;	

//========================================================================
//========================================================================
//========================================================================
//========================================================================
parameters: /* empty */
            | parameters T_NEWLINE
            | parameters T_NUMBER T_NUMBER T_NEWLINE
              {
                  tmpVec.clear();
                  tmpVec.push_back($2);
                  tmpVec.push_back($3);                  
                  tmpParams.push_back(tmpVec);
              }
            | parameters T_NUMBER T_NUMBER T_NUMBER T_NEWLINE
              {
                  tmpVec.clear();
                  tmpVec.push_back($2);
                  tmpVec.push_back($3);                  
                  tmpVec.push_back($4);
                  tmpParams.push_back(tmpVec);
              }
            | parameters T_NUMBER T_NUMBER T_NUMBER T_NUMBER T_NEWLINE
              {
                  tmpVec.clear();
                  tmpVec.push_back($2);
                  tmpVec.push_back($3);                  
                  tmpVec.push_back($4);
                  tmpVec.push_back($5);
                  tmpParams.push_back(tmpVec);
              }
            | parameters T_NUMBER T_NUMBER T_NUMBER T_NUMBER T_NUMBER T_NEWLINE
              {
                  tmpVec.clear();
                  tmpVec.push_back($2);
                  tmpVec.push_back($3);                  
                  tmpVec.push_back($4);
                  tmpVec.push_back($5);
                  tmpVec.push_back($6);
                  tmpParams.push_back(tmpVec);
              }
            | parameters T_NUMBER sentence T_NEWLINE
              {
//	         cout << "parsing string for tmpStrings" << endl;
//                 cout << "$1 = " << $1 << "$2 = " << $2 << endl;
		 strcpy(tmpstr.tmpname,$3);
                 tmpStrings.push_back(tmpstr);
                 tmpNumbers.push_back((int)$2);
                 strcpy($3,"");
              }
           | parameters T_NUMBER T_NEWLINE
              {
	         
	      }
            ;

//===================================================================
//===================================================================
//===================================================================
//===================================================================
sentence:   T_WORD
              {
//                 cout << "found lonely word" << endl;
                 strcpy($$,$1);			 
                 strcat($$," ");			 
              }
            | sentence T_WORD 
              {
//                  cout << "appending string "<< $2 << " to string " << $$ << endl;
                  strcat($$,$2);
		  strcpy($2,"");
                  strcat($$," ");
//                  cout << "having appended "<< $2 << " to string " << $$ << endl;
              }
            ; 

err:       T_ERRORSIGN value
             {
                 $$ = $2;
		 yyInputFileLine.error.push_back($2);
             }
           | err T_ERRORSIGN value
             {
                 $$ = TMath::Sqrt($1 * $1 + $3 * $3);
		 yyInputFileLine.error.push_back($3);
             }
           ;
 
value:     T_NUMBER                        { $$ = $1; }
           | T_NUMBER T_ENERGYUNIT         { $$ = $1 * TMath::Power(10, $2); }
           | T_NUMBER T_CROSSSECTIONUNIT   { $$ = $1 * TMath::Power(10, $2); }
           ;

%%

void yyerror(char* s) {
    fprintf(stderr, "Error while reading input file: %s\n", s);
    exit(EXIT_FAILURE);
}


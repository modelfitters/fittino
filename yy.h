/***************************************************************************
                                     yy.h
                             -------------------    
    Declarations used for data transfer between yacc and Input class.
                             -------------------
    $Id$
    $Name$
                             -------------------
    begin                : September 14, 2003
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

#ifndef FITTINO_YY_H
#define FITTINO_YY_H

#include <stdio.h>
#include <misc.h>


extern FILE* yyin;

//extern "C" {
extern int yylex(void);
extern void yyerror(char*);
extern int yyparse(void);
extern int yywrap(void);
//}

extern vector<MeasuredValue> yyMeasuredVec;
extern vector<MeasuredValue> yyThrownVec;

extern vector<MeasuredValue> yyFittedVec;
extern vector<MeasuredValue> yyUniversalityVec;
extern vector<MeasuredValue> yyFixedVec;  // contains mu, M1, M2 ...

extern CorrelationMatrix     yyMeasuredCorrelationMatrix;
extern CorrelationMatrix     yyDirectInputCovarianceMatrix;
extern CorrelationMatrix     yyFittedCorrelationMatrix;
		        
extern bool yyUseLoopCorrections;
extern bool yyCalcPullDist;
extern bool yyScanParameters;
extern bool yyPerformFit;
extern bool yyISR;
extern bool yyCalculatorError;
extern bool yyUseMinos;
extern bool yyGetContours;
extern bool yyUseHesse;
extern bool yyUseSimAnnBefore;
extern bool yyUseSimAnnWhile;
extern bool yyUseMarkovChains;
extern bool yyMarkovChainReadjustWidth;
extern bool yyUseGivenStartValues;
extern bool yyFitAllDirectly;
extern bool yyCalcIndChisqContr;
extern bool yyBoundsOnX;
extern bool yySepFitTanbX;
extern bool yySepFitTanbMu;
extern bool yySepFitmA;
extern bool yyScanX;
extern bool yyVerbose;
extern bool yyAdaptiveSimAnn;
extern bool yyNoBoundsAtAll;
extern bool yySimAnnUncertainty;
extern bool yySimAnnUncertaintyRunDown;
extern bool yyGetTempFromFirstChiSqr;
extern bool yyRandomDirUncertainties;
extern bool yyPerformSingleFits;
extern bool yyUseHiggsLimits;
extern bool yyQuarkFlavourViolation;
extern bool yyRandomParameters;
extern bool yyUseObservableScatteringBefore;
extern bool yyToyFitShakeFittedVecBeforeStart;
extern bool          yySetErrorFlag ;
extern bool yyUseSimplexMinOnly;
extern bool yyUseSimplexMin;
extern bool yyRequireNeut1LSP;
extern bool yyPreliminaryScan;

extern unsigned int yyNumberOfMinimizations;
extern unsigned int yySeedForObservableScattering;

typedef struct {
  string name;
  double value;
  double error;
} parameter_t;

extern vector <parameter_t> yyFixedPar;
extern vector <parameter_t> yyFittedPar;
extern vector <parameter_t> yyRandomPar;
extern vector<ScanParameter> yyScanPar;

extern vector<parameter_t> yyCorrelatedErr;

extern string yyCalculatorPath;
extern string yyRelicDensityCalculatorPath;
extern string yyLEOCalculatorPath;
extern string yyDashedLine;
extern string yyMarkovInterfaceFilePath;


extern int    yyParseError;
extern int    yyNaN;
extern int    yyInfinity;
extern double        yyErrDef;
extern int           yyNumberPulls;
extern int           yyMaxCallsSimAnn;
extern int           yySimAnnEpochLength;
extern int           yyRandomGeneratorSeed;
extern int           yyNumberOfDirections;
extern int           yyMaxMarkovChain;
extern int           yyMarkovChainReadjustWidthPeriod;
extern double        yyMarkovChainChi2Scale;
extern double        yyTempRedSimAnn;
extern double        yyInitTempSimAnn;
extern double        yyXscanlow;
extern double        yyXscanhigh;
extern double        yyMaxCalculatorTime;
extern double        yyLimitSlope;

extern int           yyMinuitStrategy;
extern double        yyMachinePrecision;

extern unsigned int  yyFitModel;

//==========================================================
// Low Energy Observables: (SPhenoLowEnergy)
extern double        yybsg;     //  1  BR(b -> s gamma) 
extern double        yybsmm;    //  2  BR(b -> s mu+ mu-)
extern double        yyB_smm;   //  3  BR(Bs -> mu+ mu-)
extern double        yyB_utn;   //  4  BR(B_u -> tau nu)
extern double        yydMB_d;   //  5  |Delta(M_Bd)| [ps^-1]
extern double        yydMB_s;   //  6  |Delta(M_Bs)| [ps^-1]
extern double        yygmin2e;  // 10  Delta(g-2)_electron
extern double        yygmin2m;  // 11  Delta(g-2)_muon
extern double        yygmin2t;  // 12  Delta(g-2)_tau
extern double        yydrho;    // 30  Delta(rho_parameter)
//---------------------------------------------------------
extern double        yyOmega;   //     relic density
extern double        yygSquaredAh;  // hAZ coupling correction
extern double        yygSquaredZh;  // hZZ coupling correction

//==========================================================
extern double        yyYtau; // tau Yukawa coupling
extern double        yyg; // SU2 coupling constant
extern double        yygprime; // U1 coupling constant
extern double        yyN11; // neutralino mixing matrix element
extern double        yyN12; // neutralino mixing matrix element
extern double        yyN13; // neutralino mixing matrix element
extern double        yyThetaStau; // stau mixing angle

//==========================================================
// Low Energy Observables: (PREDICT)
extern double yyBsg_npf;         // R(B->s gamma)
extern double yydm_s_npf;        // R(Delta m_s)
extern double yyB_smm_npf;       // B(Bs->mumu)
extern double yyBtn_npf;         // R(B->tau nu)
extern double yyB_sXsll_npf;     // R(Bs->Xsll)
extern double yyKlnu_npf;        // R(K->l nu)
extern double yygmin2m_npf;      // D(g-2)
extern double yyMassW_npf;       // m(W)
extern double yysin_th_eff_npf;  // sin(th_eff(Qfb))
extern double yyGammaZ_npf;      // Gamma(Z)
extern double yyR_l_npf;         // R_l (l=e,mu)
extern double yyR_b_npf;         // R_b
extern double yyR_c_npf;         // R_c
extern double yyA_fbb_npf;       // A_fb(b)
extern double yyA_fbc_npf;       // A_fb(c)
extern double yyA_b_npf;         // A_b
extern double yyA_c_npf;         // A_c
extern double yyA_l_npf;         // A_l(SLD)
extern double yyMassh0_npf;      // m(h0)
extern double yyOmega_npf;       // Omega_h
extern double yyA_tau_npf;       // A_l(P_tau)
extern double yyA_fbl_npf;       // A_fb(l)
extern double yysigma_had0_npf;  // sigma_had^0
extern double yydm_d_npf;        // R(Delta m_d)
extern double yydm_k_npf;        // R(Delta m_k)
extern double yyKppinn_npf;      // R(K->pi nu nu)
extern double yyB_dll_npf;       // B(Bd->ll)
extern double yyDmsDmd_npf;      // R(Dms)/R(Dmd)
extern double yyD_0_npf;         // D_0(K*gamma)
extern double yybsg_npf;         // B(b->sg)

extern std::vector<InputFileLine> yyInputFile;
extern int yyInputFileLineNo;

#endif

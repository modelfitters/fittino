/***************************************************************************
                                    misc.h
                             -------------------    
    Miscellaneous definitions of general use.
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

#ifndef FITTINO_MISC_H
#define FITTINO_MISC_H

#include <string>
#include <iostream>
//#include <stream.h>
#include <vector>
#include <map>
#include <utility>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include<leshouches.h>

using namespace std;

#define sqr(x) ((x)*(x))


// Paricle IDs
#define ID_d      1
#define ID_u      2
#define ID_s      3
#define ID_c      4
#define ID_b      5
#define ID_t      6

#define ID_e     11
#define ID_nue   12
#define ID_mu    13
#define ID_numu  14
#define ID_tau   15
#define ID_nutau 16

#define ID_g     21
#define ID_gamma 22
#define ID_Z     23
#define ID_W     24

// for NMSSM
#define ID_H1    25
#define ID_H2    35
#define ID_H3    45
#define ID_A1    36
#define ID_A2    46
// else
#define ID_h     25
#define ID_H     35
#define ID_A     36
#define ID_Hplus 37
#define ID_G     39

#define ID_sdL     1000001
#define ID_suL     1000002
#define ID_ssL     1000003
#define ID_scL     1000004
#define ID_sb1     1000005
#define ID_st1     1000006

#define ID_seL     1000011
#define ID_snueL   1000012
#define ID_smuL    1000013
#define ID_snumuL  1000014
#define ID_stau1   1000015
#define ID_snutauL 1000016

#define ID_sdR     2000001
#define ID_suR     2000002
#define ID_ssR     2000003
#define ID_scR     2000004
#define ID_sb2     2000005
#define ID_st2     2000006
		   
#define ID_seR     2000011
#define ID_snueR   2000012
#define ID_smuR    2000013
#define ID_snumuR  2000014
#define ID_stau2   2000015
#define ID_snutauR 2000016

#define ID_gluino  1000021
#define ID_chi01   1000022
#define ID_chi02   1000023
#define ID_chip1   1000024
#define ID_chi03   1000025
#define ID_chi04   1000035
#define ID_chi05   1000045 // for NMSSM
#define ID_chip2   1000037
#define ID_gravitino   1000039


enum SwitchState    { off, on };
enum Calculator     { SPHENO, SUSPECT, SPHENO_FH };
enum RelicDensityCalculator { NORELICDENSITYCALCULATOR, MICROMEGAS };
enum LEOCalculator  { NOLEOCALCULATOR, NPFITTER};
enum ObservableType { mass, Pwidth, xsection, br, xsectionbr, Pedge, xsbr, brratio, brsum, LEObs, brprod, tauFromStau1Polarisation, SMPrecision, other, weighted };
//   LEObservables  { bsg,                                           gmin2,          drho, omega };
enum LEObservables  { bsg, bsmm, B_smm, B_utn, dMB_d, dMB_s, gmin2e, gmin2m, gmin2t, drho, omega, Bsg_npf, dm_s_npf, B_smm_npf, Btn_npf, B_sXsll_npf, Klnu_npf, gmin2m_npf, MassW_npf, sin_th_eff_npf, GammaZ_npf, R_l_npf, R_b_npf, R_c_npf, A_fbb_npf, A_fbc_npf, A_b_npf, A_c_npf, A_l_npf, Massh0_npf, Omega_npf, A_tau_npf, A_fbl_npf, sigma_had0_npf, dm_d_npf, dm_k_npf, Kppinn_npf, B_dll_npf, DmsDmd_npf, D_0_npf, bsg_npf};
enum SMPrecisionObservables { alphas, alphaem, gf };
enum FitModel       { MSSM, mSUGRA, GMSB, AMSB, NMSSM, XMSUGRA };
enum higgsLimitProcesses 
  { Zh_Zbb, 
    Zh_Ztt,
    Zh_ZAA_Zbbbb,
    Zh_ZAA_Ztttt,
    Ah_bbbb,
    Ah_bbtt,
    Ah_ttbb,
    Ah_tttt,
    Ah_AAA_bbbbbb,
    Ah_AAA_tttttt,
    higgsLimitProcessesLast
  };

extern bool         yyUseLoopCorrections;
extern bool         yyCalcPullDist;
extern bool         yyScanParameters;
extern unsigned int yyCalculator;
extern unsigned int yyRelicDensityCalculator;
extern unsigned int yyLEOCalculator;
extern map<int,string> yyParticleNames;
extern map<string,int> yyParticleIDs;

extern vector<doubleVec_t> CrossSectionProduction;

struct MeasuredValue {
  double          value;
  double          error;  // total error (added in quadrature)
  double          staterr; // statistical error
  map<string, double> correrror;
  int             position;
  double          theovalue;
  bool            theoset;
  bool            first;
  bool            hasbeenset;
  operator        double() const { return value; }
  string          name;
  int             id;
  bool            bound;
  double          bound_low;
  double          bound_up;
  ObservableType  type;
  int             alias;
  vector<int>     daughters;
  vector<int>     products;
  double          sqrts;
  double          polarisation1;
  double          polarisation2;
  bool            nofit;
  double          positive_error;
  double          negative_error;
  string          universality;
  string          calculator;
  bool            temp_nofit;
  double          scaling;
  bool            setScaling;
};

struct ScanParameter {
  string          name;
  double          min;
  double          max;
  unsigned int    numberOfSteps;
};

struct InputFileLine {
  string          prevalue;
  string          postvalue;
  double          value;
  vector<double>  error;
};

extern MeasuredValue fitterM1;
extern MeasuredValue fitterM2;
extern MeasuredValue fitterAbsM3;
extern MeasuredValue fitterMassA0;
extern MeasuredValue fitterTanBeta;
extern MeasuredValue fitterMu;
extern MeasuredValue fitterMassSelR;
extern MeasuredValue fitterMassSelL;
extern MeasuredValue fitterMassSupR;
extern MeasuredValue fitterMassSupL;
extern MeasuredValue fitterA;
extern MeasuredValue fitterMassTop;

// observables received from suspect:
extern double        sRecMassChargino1;
extern double        sRecMassChargino2;
extern double        sRecMassNeutralino1;
extern double        sRecMassNeutralino2;
extern double        sRecMassNeutralino3;
extern double        sRecMassNeutralino4;
extern double        sRecMassGluino;
extern double        sRecMassh0;
extern double        sRecMassH0;
extern double        sRecMassA0;
extern double        sRecMassHpm;
extern double        sRecMassSelR;
extern double        sRecMassSelL;
extern double        sRecMassSupR;
extern double        sRecMassSupL;
extern double        sRecCos2PhiL;
extern double        sRecCos2PhiR;
extern double        sRecMassTop;

class TagMap : public map<string, vector<MeasuredValue*> > {
   private:
      vector<MeasuredValue>* fMeasuredVector;
   public:
      TagMap(vector<MeasuredValue>* measuredvector);
      void Print();
};

class CorrelationMatrix {
   protected:
      const vector<MeasuredValue>* fMeasuredVector;
      TagMap*                      fTagMap;
      TMatrixDSym*                 fInverseCovarianceMatrix;
      TMatrixDSym*                 fCorrelationMatrix;
      TMatrixDSym*                 fCovarianceMatrix;
   public:
      CorrelationMatrix(); 
      CorrelationMatrix(const vector<MeasuredValue>* measuredvector);
      CorrelationMatrix(const vector<MeasuredValue>* measuredvector, TagMap* tagmap);
      CorrelationMatrix(const CorrelationMatrix& another);

      virtual ~CorrelationMatrix();

      const double& operator()(int rown, int coln) const;
      CorrelationMatrix& operator=(const CorrelationMatrix& source);

      const TMatrixDSym& GetCovarianceMatrix() const;
      double GetCorrelation(const string& name1, const string& name2) const;
      double GetCorrelation(unsigned int index1, unsigned int index2) const;
      double GetCovariance(const string& name1, const string& name2) const;
      double GetCovariance(unsigned int index1, unsigned int index2) const;
      double GetInverseCovariance(const string& name1, const string& name2) const;
      double GetInverseCovariance(unsigned int index1, unsigned int index2) const;
      void add(const string& name1, const string& name2, double entry);
      void add(unsigned int index1, unsigned int index2, double entry);
      void AddCovariance(const string& name1, const string& name2, double entry);
      void AddCovariance(unsigned int index1, unsigned int index2, double entry);
      void Calculate();
      void CalculateCovarianceMatrix();
      void Print();
      void PrintCovariance();
      void PrintInverseCovariance();
      void AddFullCovarianceMatrix(const CorrelationMatrix& addedCovarianceMatrix);
      void TransformCovarianceMatrixIntoCorrelationMatrix();
      void CalculateInverseCovarianceMatrix();
      bool TestCovarianceMatrixExistence();
};

TVectorD getCorrelatedRandomNumbers(const TVectorD& mean, const TMatrixDSym& covarianceMatrix);

#endif

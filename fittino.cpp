/***************************************************************************
                                 fittino.cpp
                             -------------------    
    Main class which does Fittino's calculations.
                             -------------------
    $Id$
    $Name$
                             -------------------

    begin                : September 10, 2003
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

#include <fittino.h>
#include <indchisq.h>
#include <input.h>
#include <smearedinput.h>
#include <units.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <yy.h>
#include <fstream>
#include <stdio.h>
#include <leshouches.h>

#include <TRandom.h>
#include <TMath.h>
#include <TFitter.h>
#include <TMinuit.h>
#include <TGraph.h>

#include <TFile.h>
#include <TH1.h>

MeasuredValue fitterM1;
MeasuredValue fitterM2;
MeasuredValue fitterAbsM3;
MeasuredValue fitterMassA0;
MeasuredValue fitterTanBeta;
MeasuredValue fitterMu;
MeasuredValue fitterMassSelR;
MeasuredValue fitterMassSelL;
MeasuredValue fitterMassSupR;
MeasuredValue fitterMassSupL;
MeasuredValue fitterA;
MeasuredValue fitterMassTop;

// observables received from suspect:
double        sRecMassChargino1;
double        sRecMassChargino2;
double        sRecMassNeutralino1;
double        sRecMassNeutralino2;
double        sRecMassNeutralino3;
double        sRecMassNeutralino4;
double        sRecMassGluino;
double        sRecMassh0;
double        sRecMassH0;
double        sRecMassA0;
double        sRecMassHpm;
double        sRecMassSelR;
double        sRecMassSelL;
double        sRecMassSupR;
double        sRecMassSupL;
double        sRecCos2PhiL;
double        sRecCos2PhiR;
double        sRecMassTop;


extern "C" { void susygen_call_test_(int *, double [],double *);}
extern "C" { void mssmsusygen_init_(int *, int *, float *, float *,float *,float *,float *,float *);}
extern "C" { void mssmsusygen_();}
extern "C" { void call_suspect_free_(int *,int *,double *,double *,double *,
				     double [],double *,double *,double *,double *,double *,double *);}

//extern "C" { void get_suspect_(double [],double [],double *,double [],double [2][2],double [2][2],double [4][4],double []);}
void callSuspect();
int callSPheno();
int checkcall(char programpath[1024], unsigned int runtime);
void WriteLesHouches(double* x);
void fitterFCN(int &npar, double *gin, double &f, double *x, int iflag); 
void chi2Function(int& npar, double* gin, double& f, double* x, int iflag);
MeasuredValue* ReturnMeasuredValue (string name);
MeasuredValue* ReturnFittedValue (string name);
MeasuredValue* ReturnFixedValue (string name);
int ReturnFittedPosition (string name);
bool FindInFitted (string name);
bool FindInFixed (string name);
bool FindInUniversality(string name);
MeasuredValue* ReturnUniversality (string name);
void  ReadLesHouches();
void  ParseLesHouches();
double give_xs (doubleVec_t initial, int channel, int element );
double give_br ( int id, int decay, int element );
int FindInFittedPar (string name);
double gchisq;
int    gstat;

Fittino::Fittino(const Input* input)
{
    fInput = input;
    fFittedCovarianceMatrix = 0;
}


Fittino::~Fittino()
{

}


void Fittino::calculateM1(const SmearedInput* input)
{
    double ka[3],kb[3],kc[3],kd[3];
    double rad[8], mka[4],mkb[4],mkc[4];
    double mst[4][2];
    double mnt[4];

    double sin2ThetaW = input->GetValue("sin2ThetaW");
    double cos2ThetaW = 1. - sin2ThetaW;
    double beta       = TMath::ATan(fTanBeta);

    mnt[0] = input->GetValue("massNeutralino1");
    mnt[1] = input->GetValue("massNeutralino2");
    //    mnt[2] = input->GetValue("massNeutralino3");
    //    mnt[3] = input->GetValue("massNeutralino4");

    for (int jsss=0;jsss<3;jsss++) {
      ka[jsss]=0.;
      kb[jsss]=0.;
      kc[jsss]=0.;
      kd[jsss]=0.;
    }

    double mz = input->GetValue("massZ");

    ka[2]=sqr(fM2)+2.*sqr(fAbsMu)+2.*sqr(mz);

    kb[0]=sqr(fM2)+2.*sqr(fAbsMu)+2.*sqr(mz)
          *cos2ThetaW;
    kb[1]=-2.*TMath::Sin(2.*beta)*fAbsMu*sqr(mz)
          *sin2ThetaW;
    kb[2]= 2.*sqr(fAbsMu)*sqr(fM2)+sqr(sqr(fAbsMu)+sqr(mz))
      -2.*sqr(mz)*TMath::Sin(2.*beta)*cos2ThetaW*fAbsMu*fM2
      +2.*sin2ThetaW*sqr(mz)*sqr(fM2);

    kd[0]= TMath::Power(fAbsMu,4)*sqr(fM2)+TMath::Power(mz,4)
	   *sqr(fAbsMu)*sqr(TMath::Sin(2.*beta))*sqr(cos2ThetaW)
           -2.*sqr(mz)*TMath::Power(fAbsMu,3)*fM2
	   *TMath::Sin(2.*beta)*cos2ThetaW;
    kd[1]= 2.*TMath::Power(mz,4)*sqr(fAbsMu)*fM2
	   *sin2ThetaW
	   *cos2ThetaW*sqr(TMath::Sin(2.*beta))
           -2.*sqr(mz)*TMath::Power(fAbsMu,3)*sqr(fM2)
           *TMath::Sin(2.*beta)*sin2ThetaW;
    kd[2]= TMath::Power(mz,4)*sqr(fAbsMu)
           *sqr(TMath::Sin(2.*beta))
	   *sqr(fM2)*sqr(sin2ThetaW);

    // --------------------------------------------------
    //   kc in version vom 22.5.01 nach Vergleich mit math
    kc[0]=TMath::Power(fAbsMu,4)+2.*sqr(fAbsMu)*sqr(fM2)+2.
          *sqr(mz)*sqr(fAbsMu)*cos2ThetaW
          +TMath::Power(mz,4)*sqr(cos2ThetaW)
          -2.*sqr(mz)*fAbsMu*TMath::Sin(2.*beta)*fM2*cos2ThetaW;
    kc[1]= 2.*TMath::Power(mz,4)*fM2*sin2ThetaW
           *cos2ThetaW-2.*sqr(mz)*TMath::Sin(2.*beta)
           *TMath::Power(fAbsMu,3)*sin2ThetaW-2.
           *sqr(mz)*TMath::Sin(2.*beta)*fAbsMu
           *sqr(fM2)*sin2ThetaW;
    kc[2]= TMath::Power(fAbsMu,4)*sqr(fM2)+TMath::Power(mz,4)
           *sqr(fAbsMu)*sqr(TMath::Sin(2.*beta))
           +2.*sqr(mz)*sqr(fAbsMu)*sqr(fM2)
           *sin2ThetaW
           -2.*sqr(mz)*TMath::Sin(2.*beta)*fM2
           *TMath::Power(fAbsMu,3)*cos2ThetaW
           +TMath::Power(mz,4)*sqr(fM2)
           *sqr(sin2ThetaW);

    //   If all four neutralino masses useful:
    //   then loop over is = 1,4
    //   Otherwise:
    //   Only use LSP Mass:

    int is = 0;

    rad[is+4]=0.;
    mka[is]=0.;
    mkb[is]=0.;
    mkc[is]=0.;
    mka[is]=-TMath::Power(mnt[is],6)+kb[0]*TMath::Power(mnt[is],4)
      -kc[0]*sqr(mnt[is])+kd[0];
    mkb[is]=kb[1]*TMath::Power(mnt[is],4)-kc[1]*sqr(mnt[is])+kd[1];
    mkc[is]=TMath::Power(mnt[is],8)-ka[2]*TMath::Power(mnt[is],6)
            +kb[2]*TMath::Power(mnt[is],4)
            -kc[2]*sqr(mnt[is])+kd[2];

    rad[is+4]=sqr(mkb[is])-4.*mka[is]*mkc[is];
    if (rad[is+4]<0.) {
      //    write(*,'(a23)') '# rad(mst)<0';
    }

    mst[is][0]=(-mkb[is]+TMath::Sqrt(sqr(mkb[is])-4.*mka[is]*mkc[is]))
               /2./mka[is];
    mst[is][1]=(-mkb[is]-TMath::Sqrt(sqr(mkb[is])-4.*mka[is]*mkc[is]))
               /2./mka[is];

    fM1.value = mst[0][0];
}


void Fittino::calculateM2(const SmearedInput* input)
{
    double massChargino1 = input->GetValue("massChargino1");
    double massChargino2 = input->GetValue("massChargino2");
    double massW         = input->GetValue("massW");
    double cos2phil      = input ->GetValue("cos2PhiL");
    double cos2phir      = input ->GetValue("cos2PhiR");

    double sigma = ((sqr(massChargino2)
		   +sqr(massChargino1))
                   /(2.*sqr(massW)))-1.;
    double delta = (sqr(massChargino2)
                   -sqr(massChargino1))
                   /(4.*sqr(massW));

    fM2.value = massW * TMath::Sqrt(sigma - delta
                * (cos2phil+cos2phir));
}


void Fittino::calculateAbsM3(const SmearedInput* input)
{
    fAbsM3.value = input->GetValue("massGluino");
}  


void Fittino::calculateTanBeta(const SmearedInput* input)
{
    double massChargino1 = input->GetValue("massChargino1");
    double massChargino2 = input->GetValue("massChargino2");
    double massW         = input->GetValue("massW");
    double cos2phil      = input ->GetValue("cos2PhiL");
    double cos2phir      = input ->GetValue("cos2PhiR");

    double delta = (sqr(massChargino2)
		   -sqr(massChargino1))
                   /(4.*sqr(massW));
    double tan2Beta = (1.+delta*(cos2phir
                      -cos2phil))/
		      (1.-delta*(cos2phir
                     -cos2phil));

    (tan2Beta >= 0.) ? fTanBeta.value = TMath::Sqrt(tan2Beta)
                     : fTanBeta.value = 0.;
}  


void Fittino::calculateAbsMu(const SmearedInput* input)
{
    double massChargino1 = input->GetValue("massChargino1");
    double massChargino2 = input->GetValue("massChargino2");
    double massW         = input->GetValue("massW");
    double cos2phil      = input ->GetValue("cos2PhiL");
    double cos2phir      = input ->GetValue("cos2PhiR");

    double sigma = (sqr(massChargino2)
                   +sqr(massChargino1))
                   /(2.*sqr(massW))-1.;
    double delta = (sqr(massChargino2)
                   -sqr(massChargino1))
                   /(4.*sqr(massW));

    fAbsMu.value = massW * TMath::Sqrt(sigma
                   + delta*(cos2phil+cos2phir));
}


void Fittino::calculateSignMu(const SmearedInput* input)
{
    double massChargino1 = input->GetValue("massChargino1");
    double massChargino2 = input->GetValue("massChargino2");
    double massW         = input->GetValue("massW");

    double beta  = TMath::ATan(fTanBeta);
    double delta = (sqr(massChargino2)
		    -sqr(massChargino1))
	           /(4.*sqr(massW));

    fSignMu.value = -(sqr(delta)-sqr(sqr(fM2)-sqr(fAbsMu))-
                    4.*sqr(massW)*(sqr(fM2)+sqr(fAbsMu)) - 
                    4.*TMath::Power(massW,4)
		    *(sqr(TMath::Cos(2.*beta))))
                    /(8.*sqr(massW)*fM2*
		    fAbsMu*TMath::Sin(2.*beta));
}

void Fittino::calculateMSelL(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mElectron = 0.000511;
  double mSelectron1 = input->GetValue("massSelectronL");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(0.5-sin2ThetaW) - sqr(mElectron) + sqr(mSelectron1);
  fMSelL.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMSelR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mElectron = 0.000511;
  double mSelectron1 = input->GetValue("massSelectronR");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-sin2ThetaW) - sqr(mElectron) + sqr(mSelectron1);
  fMSelR.value = TMath::Sqrt(Msq);

  return;
}

void Fittino::calculateMSmuL(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mMuon = 0.105658;
  double mSmuon1 = input->GetValue("massSmuL");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(0.5-sin2ThetaW) - sqr(mMuon) + sqr(mSmuon1);
  fMSmuL.value = TMath::Sqrt(Msq);
  //  cout << fMSmuL.value << endl;
  return;
}
void Fittino::calculateMSmuR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mMuon = 0.105658;
  double mSmuon1 = input->GetValue("massSmuR");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-sin2ThetaW) - sqr(mMuon) + sqr(mSmuon1);
  fMSmuR.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMStauL(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mTau = 1.77;
  double mStau1 = input->GetValue("massStau2");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(0.5-sin2ThetaW) - sqr(mTau) + sqr(mStau1);
  fMStauL.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMStauR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mTau = 1.77;
  double mStau1 = input->GetValue("massStau1");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-sin2ThetaW) - sqr(mTau) + sqr(mStau1);
  fMStauR.value = TMath::Sqrt(Msq);

  return;
}

void Fittino::calculateMSupL(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  //double mElectron = 0.000511;
  double mSup1 = input->GetValue("massSupL");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(0.5-(2./3.)*sin2ThetaW) + sqr(mSup1);
  fMSupL.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMScharmL(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  //double mElectron = 0.000511;
  double mScharm1 = input->GetValue("massScharmL");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(0.5-(2./3.)*sin2ThetaW) + sqr(mScharm1);
  fMScharmL.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMStopL(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mTop = 174.3;
  double mStop1 = input->GetValue("massStop2");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(0.5-(2./3.)*sin2ThetaW) -sqr(mTop) + sqr(mStop1);
  fMStopL.value = TMath::Sqrt(Msq);

  return;
}

void Fittino::calculateMSupR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  //double mTop = 174.3;
  double mSup1 = input->GetValue("massSupR");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-(2./3.)*sin2ThetaW) + sqr(mSup1);
  fMSupR.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMScharmR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  //double mTop = 174.3;
  double mScharm1 = input->GetValue("massScharmR");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-(2./3.)*sin2ThetaW) + sqr(mScharm1);
  fMScharmR.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMStopR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mTop = 174.3;
  double mStop1 = input->GetValue("massStop1");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-(2./3.)*sin2ThetaW) -sqr(mTop) + sqr(mStop1);
  fMStopR.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMSdownR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  //double mTop = 174.3;
  double mSdown1 = input->GetValue("massSdownR");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-(1./3.)*sin2ThetaW) + sqr(mSdown1);
  fMSdownR.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMSstrangeR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  //double mTop = 174.3;
  double mSstrange1 = input->GetValue("massSstrangeR");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-(1./3.)*sin2ThetaW) + sqr(mSstrange1);
  fMSstrangeR.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateMSbottomR(const SmearedInput* input)
{
  double massZ = input->GetValue("massZ");
  double beta  = TMath::ATan(fTanBeta);
  double sin2ThetaW = input->GetValue("sin2ThetaW");
  double mBottom = 4.2;
  double mSbottom1 = input->GetValue("massSbottom2");
  double Msq;

  Msq = sqr(massZ)*TMath::Cos(2.*beta)*(-(1./3.)*sin2ThetaW) -sqr(mBottom)+ sqr(mSbottom1);
  fMSbottomR.value = TMath::Sqrt(Msq);

  return;
}
void Fittino::calculateXtau(const SmearedInput* input)
{
  double mu = fMu.value;
  
  fXtau.value = mu * fTanBeta;

  return;
}
void Fittino::calculateXtop(const SmearedInput* input)
{
  double mu = fMu.value;
  
  fXtop.value = mu / fTanBeta;

  return;
}

void Fittino::calculateXbottom(const SmearedInput* input)
{
  double mu = fMu.value;
  
  fXbottom.value = mu * fTanBeta;

  return;
}

//==================================================================
//
//
//
//
//
//
//==================================================================

void Fittino::calculateTreeLevelValues(int nthrows)
{
    double sumM1  = 0, sumM2    = 0, sumM3       = 0;
    double sumMu  = 0, sumBeta  = 0, sumTanBeta  = 0.;
    double sum2M1 = 0, sum2M2   = 0, sum2M3      = 0;
    double sum2Mu = 0, sum2Beta = 0, sum2TanBeta = 0;
    double sumMSelL = 0, sum2MSelL = 0; 
    double sumMSmuL = 0, sum2MSmuL = 0; 
    double sumMStauL = 0, sum2MStauL = 0; 
    double sumMSelR = 0, sum2MSelR = 0; 
    double sumMSmuR = 0, sum2MSmuR = 0; 
    double sumMStauR = 0, sum2MStauR = 0; 
    double sumMSupL = 0, sum2MSupL = 0; 
    double sumMScharmL = 0, sum2MScharmL = 0; 
    double sumMStopL = 0, sum2MStopL = 0; 
    double sumMSupR = 0, sum2MSupR = 0; 
    double sumMScharmR = 0, sum2MScharmR = 0; 
    double sumMStopR = 0, sum2MStopR = 0; 
    double sumMSdownR = 0, sum2MSdownR = 0; 
    double sumMSstrangeR = 0, sum2MSstrangeR = 0; 
    double sumMSbottomR = 0, sum2MSbottomR = 0; 
    double sumXtau = 0, sum2Xtau = 0; 
    double sumXtop = 0, sum2Xtop = 0; 
    double sumXbottom = 0, sum2Xbottom = 0; 

    int    nvalid = 0;
    MeasuredValue tmpValue;
    MeasuredValue fnew;
    bool par_already_found;


    SmearedInput* smearedinput = new SmearedInput(fInput);

    for (int i = 0; i < nthrows; i++) {
	smearedinput->throwDice();

        // NOTICE: order of calculation is important
	calculateTanBeta(smearedinput);
	calculateM2(smearedinput);
	calculateAbsMu(smearedinput);
	calculateSignMu(smearedinput);
	fMu.value = TMath::Sign(fAbsMu, fSignMu);
	calculateM1(smearedinput);
	calculateAbsM3(smearedinput);

 	calculateMSelL(smearedinput);
 	calculateMSmuL(smearedinput);
 	calculateMStauL(smearedinput);
 	calculateMSelR(smearedinput);
 	calculateMSmuR(smearedinput);
 	calculateMStauR(smearedinput);

 	calculateMSupL(smearedinput);
 	calculateMScharmL(smearedinput);
 	calculateMStopL(smearedinput);
 	calculateMSupR(smearedinput);
 	calculateMScharmR(smearedinput);
 	calculateMStopR(smearedinput);
 	calculateMSdownR(smearedinput);
 	calculateMSstrangeR(smearedinput);
 	calculateMSbottomR(smearedinput);

 	calculateXtau(smearedinput);
 	calculateXtop(smearedinput);
 	calculateXbottom(smearedinput);
	
	if (fTanBeta > 0) {
	    nvalid++;
	    sumM1 += fM1;
	    sumM2 += fM2;
	    sumM3 += fAbsM3;
	    sumBeta += TMath::ATan(fTanBeta);
	    sumTanBeta += fTanBeta;
	    sumMu += fMu;

	    sumMSelL += fMSelL;
	    sumMSmuL += fMSmuL;
	    sumMStauL += fMStauL;
	    sumMSelR += fMSelR;
	    sumMSmuR += fMSmuR;
	    sumMStauR += fMStauR;
	    sumMSupL += fMSupL;
	    sumMScharmL += fMScharmL;
	    sumMStopL += fMStopL;
	    sumMSupR += fMSupR;
	    sumMScharmR += fMScharmR;
	    sumMStopR += fMStopR;
	    sumMSdownR += fMSdownR;
	    sumMSstrangeR += fMSstrangeR;
	    sumMSbottomR += fMSbottomR;
	    sumXtau += fXtau;
	    sumXtop += fXtop;
	    sumXbottom += fXbottom;

	    sum2M1 += sqr(fM1);
	    sum2M2 += sqr(fM2);
	    sum2M3 += sqr(fAbsM3);
	    sum2Beta += sqr(TMath::ATan(fTanBeta));
	    sum2TanBeta += sqr(fTanBeta);
	    sum2Mu += sqr(fMu);

	    sum2MSelL += sqr(fMSelL);
	    sum2MSmuL += sqr(fMSmuL);
	    sum2MStauL += sqr(fMStauL);
	    sum2MSelR += sqr(fMSelR);
	    sum2MSmuR += sqr(fMSmuR);
	    sum2MStauR += sqr(fMStauR);
	    sum2MSupL      += sqr(fMSupL);
	    sum2MScharmL   += sqr(fMScharmL);
	    sum2MStopL     += sqr(fMStopL);
	    sum2MSupR      += sqr(fMSupR);
	    sum2MScharmR   += sqr(fMScharmR);
	    sum2MStopR     += sqr(fMStopR);
	    sum2MSdownR    += sqr(fMSdownR);
	    sum2MSstrangeR += sqr(fMSstrangeR);
	    sum2MSbottomR  += sqr(fMSbottomR);
	    sum2Xtau       += sqr(fXtau);
	    sum2Xtop       += sqr(fXtop);
	    sum2Xbottom    += sqr(fXbottom);
	}
    }

    delete smearedinput;


    // look at given observables, decide which parameters can be fitted, ask 
    // user for list of parameters to be fitted

    fM1.name  = "M1";
    fM1.value = sumM1 / double(nvalid);
    fM1.error = TMath::Sqrt(sum2M1 / double(nvalid) - sqr(fM1.value));
    if (yyUseGivenStartValues && (FindInFittedPar("M1") >= 0)) fM1.value = yyFittedPar[FindInFittedPar("M1")].value ;
    fM1.bound_low = 0.;
    fM1.bound_up = 10000.;

    fM2.name  = "M2";
    fM2.value = sumM2 / double(nvalid);
    fM2.error = TMath::Sqrt(sum2M2 / double(nvalid) - sqr(fM2.value));
    if (yyUseGivenStartValues && (FindInFittedPar("M2") >= 0)) fM2.value = yyFittedPar[FindInFittedPar("M2")].value ;
    fM2.bound_low = 0.;
    fM2.bound_up = 10000.;

    fAbsM3.name  = "M3";
    fAbsM3.value = sumM3 / double(nvalid);
    fAbsM3.error = TMath::Sqrt(sum2M3 / double(nvalid) - sqr(fAbsM3.value));
    if (yyUseGivenStartValues && (FindInFittedPar("M3") >= 0)) fAbsM3.value = yyFittedPar[FindInFittedPar("M3")].value ;
    fAbsM3.bound_low = 0.;
    fAbsM3.bound_up = 10000.;

    fTanBeta.name  = "TanBeta";
    fTanBeta.value = sumTanBeta / double(nvalid);
    fTanBeta.error = TMath::Sqrt(sum2TanBeta / double(nvalid)
				 - sqr(fTanBeta.value));
    if (yyUseGivenStartValues && (FindInFittedPar("TanBeta") >= 0)) fTanBeta.value = yyFittedPar[FindInFittedPar("TanBeta")].value ;
    fTanBeta.bound_low = 0.;
    fTanBeta.bound_up = 1000.;

    fMu.name = "Mu";
    fMu.value = sumMu / double(nvalid);
    fMu.error = TMath::Sqrt(sum2Mu / double(nvalid) - sqr(fMu.value));
    if (yyUseGivenStartValues && (FindInFittedPar("Mu") >= 0)) fMu.value = yyFittedPar[FindInFittedPar("Mu")].value ;
    fMu.bound_low = 0.;
    fMu.bound_up = 10000.;

    fMSelL.name = "MSelectronL";
    fMSelL.value = sumMSelL/ double(nvalid);
    fMSelL.error = TMath::Sqrt(sum2MSelL / double(nvalid) - sqr(fMSelL.value));
    if ((TMath::Abs(fMSelL.error)>TMath::Abs(fMSelL.value/2.)) || 
	(TMath::Abs(fMSelL.error)<TMath::Abs(fMSelL.value/10.))) {
      fMSelL.error = TMath::Abs(fMSelL.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSelectronL") >= 0)) fMSelL.value = yyFittedPar[FindInFittedPar("MSelectronL")].value ;
    fMSelL.bound_low = 0.;
    fMSelL.bound_up = 10000.;

    fMSmuL.name = "MSmuL";
    fMSmuL.value = sumMSmuL/ double(nvalid);
    fMSmuL.error = TMath::Sqrt(sum2MSmuL / double(nvalid) - sqr(fMSmuL.value));
    if ((TMath::Abs(fMSmuL.error)>TMath::Abs(fMSmuL.value/2.)) || 
	(TMath::Abs(fMSmuL.error)<TMath::Abs(fMSmuL.value/10.))) {
      fMSmuL.error = TMath::Abs(fMSmuL.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSmuL") >= 0)) fMSmuL.value = yyFittedPar[FindInFittedPar("MSmuL")].value ;
    fMSmuL.bound_low = 0.;
    fMSmuL.bound_up = 10000.;

    fMStauL.name = "MStauL";
    fMStauL.value = sumMStauL/ double(nvalid);
    fMStauL.error = TMath::Sqrt(sum2MStauL / double(nvalid) - sqr(fMStauL.value));
    if ((TMath::Abs(fMStauL.error)>TMath::Abs(fMStauL.value/2.)) || 
	(TMath::Abs(fMStauL.error)<TMath::Abs(fMStauL.value/10.))) {
      fMStauL.error = TMath::Abs(fMStauL.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MStauL") >= 0)) fMStauL.value = yyFittedPar[FindInFittedPar("MStauL")].value ;
    fMStauL.bound_low = 0.;
    fMStauL.bound_up = 10000.;
    
    fMSelR.name = "MSelectronR";
    fMSelR.value = sumMSelR/ double(nvalid);
    fMSelR.error = TMath::Sqrt(sum2MSelR / double(nvalid) - sqr(fMSelR.value));
    if ((TMath::Abs(fMSelR.error)>TMath::Abs(fMSelR.value/2.)) || 
	(TMath::Abs(fMSelR.error)<TMath::Abs(fMSelR.value/10.))) {
      fMSelR.error = TMath::Abs(fMSelR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSelectronR") >= 0)) fMSelR.value = yyFittedPar[FindInFittedPar("MSelectronR")].value ;
    fMSelR.bound_low = 0.;
    fMSelR.bound_up = 10000.;

    fMSmuR.name = "MSmuR";
    fMSmuR.value = sumMSmuR/ double(nvalid);
    fMSmuR.error = TMath::Sqrt(sum2MSmuR / double(nvalid) - sqr(fMSmuR.value));
    if ((TMath::Abs(fMSmuR.error)>TMath::Abs(fMSmuR.value/2.)) || 
	(TMath::Abs(fMSmuR.error)<TMath::Abs(fMSmuR.value/10.))) {
      fMSmuR.error = TMath::Abs(fMSmuR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSmuR") >= 0)) fMSmuR.value = yyFittedPar[FindInFittedPar("MSmuR")].value ;
    fMSmuR.bound_low = 0.;
    fMSmuR.bound_up = 10000.;

    fMStauR.name = "MStauR";
    fMStauR.value = sumMStauR/ double(nvalid);
    fMStauR.error = TMath::Sqrt(sum2MStauR / double(nvalid) - sqr(fMStauR.value));
    if ((TMath::Abs(fMStauR.error)>TMath::Abs(fMStauR.value/2.)) || 
	(TMath::Abs(fMStauR.error)<TMath::Abs(fMStauR.value/10.))) {
      fMStauR.error = TMath::Abs(fMStauR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MStauR") >= 0)) fMStauR.value = yyFittedPar[FindInFittedPar("MStauR")].value ;
    fMStauR.bound_low = 0.;
    fMStauR.bound_up = 10000.;
    
    fMSupL.name = "MSupL";
    fMSupL.value = sumMSupL/ double(nvalid);
    fMSupL.error = TMath::Sqrt(sum2MSupL / double(nvalid) - sqr(fMSupL.value));
    if ((TMath::Abs(fMSupL.error)>TMath::Abs(fMSupL.value/2.)) || 
	(TMath::Abs(fMSupL.error)<TMath::Abs(fMSupL.value/10.))) {
      fMSupL.error = TMath::Abs(fMSupL.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSupL") >= 0)) fMSupL.value = yyFittedPar[FindInFittedPar("MSupL")].value ;
    fMSupL.bound_low = 0.;
    fMSupL.bound_up = 10000.;

    fMScharmL.name = "MScharmL";
    fMScharmL.value = sumMScharmL/ double(nvalid);
    fMScharmL.error = TMath::Sqrt(sum2MScharmL / double(nvalid) - sqr(fMScharmL.value));
    if ((TMath::Abs(fMScharmL.error)>TMath::Abs(fMScharmL.value/2.)) || 
	(TMath::Abs(fMScharmL.error)<TMath::Abs(fMScharmL.value/10.))) {
      fMScharmL.error = TMath::Abs(fMScharmL.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MScharmL") >= 0)) fMScharmL.value = yyFittedPar[FindInFittedPar("MScharmL")].value ;
    fMScharmL.bound_low = 0.;
    fMScharmL.bound_up = 10000.;

    fMStopL.name = "MStopL";
    fMStopL.value = sumMStopL/ double(nvalid);
    fMStopL.error = TMath::Sqrt(sum2MStopL / double(nvalid) - sqr(fMStopL.value));
    if ((TMath::Abs(fMStopL.error)>TMath::Abs(fMStopL.value/2.)) || 
	(TMath::Abs(fMStopL.error)<TMath::Abs(fMStopL.value/10.))) {
      fMStopL.error = TMath::Abs(fMStopL.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MStopL") >= 0)) fMStopL.value = yyFittedPar[FindInFittedPar("MStopL")].value ;
    fMStopL.bound_low = 0.;
    fMStopL.bound_up = 10000.;

    fMSupR.name = "MSupR";
    fMSupR.value = sumMSupR/ double(nvalid);
    fMSupR.error = TMath::Sqrt(sum2MSupR / double(nvalid) - sqr(fMSupR.value));
    if ((TMath::Abs(fMSupR.error)>TMath::Abs(fMSupR.value/2.)) || 
	(TMath::Abs(fMSupR.error)<TMath::Abs(fMSupR.value/10.))) {
      fMSupR.error = TMath::Abs(fMSupR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSupR") >= 0)) fMSupR.value = yyFittedPar[FindInFittedPar("MSupR")].value ;
    fMSupR.bound_low = 0.;
    fMSupR.bound_up = 10000.;

    fMScharmR.name = "MScharmR";
    fMScharmR.value = sumMScharmR/ double(nvalid);
    fMScharmR.error = TMath::Sqrt(sum2MScharmR / double(nvalid) - sqr(fMScharmR.value));
    if ((TMath::Abs(fMScharmR.error)>TMath::Abs(fMScharmR.value/2.)) || 
	(TMath::Abs(fMScharmR.error)<TMath::Abs(fMScharmR.value/10.))) {
      fMScharmR.error = TMath::Abs(fMScharmR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MScharmR") >= 0)) fMScharmR.value = yyFittedPar[FindInFittedPar("MScharmR")].value ;
    fMScharmR.bound_low = 0.;
    fMScharmR.bound_up = 10000.;

    fMStopR.name = "MStopR";
    fMStopR.value = sumMStopR/ double(nvalid);
    fMStopR.error = TMath::Sqrt(sum2MStopR / double(nvalid) - sqr(fMStopR.value));
    if ((TMath::Abs(fMStopR.error)>TMath::Abs(fMStopR.value/2.)) || 
	(TMath::Abs(fMStopR.error)<TMath::Abs(fMStopR.value/10.))) {
      fMStopR.error = TMath::Abs(fMStopR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MStopR") >= 0)) fMStopR.value = yyFittedPar[FindInFittedPar("MStopR")].value ;
    fMStopR.bound_low = 0.;
    fMStopR.bound_up = 10000.;

    fMSdownR.name = "MSdownR";
    fMSdownR.value = sumMSdownR/ double(nvalid);
    fMSdownR.error = TMath::Sqrt(sum2MSdownR / double(nvalid) - sqr(fMSdownR.value));
    if ((TMath::Abs(fMSdownR.error)>TMath::Abs(fMSdownR.value/2.)) || 
	(TMath::Abs(fMSdownR.error)<TMath::Abs(fMSdownR.value/10.))) {
      fMSdownR.error = TMath::Abs(fMSdownR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSdownR") >= 0)) fMSdownR.value = yyFittedPar[FindInFittedPar("MSdownR")].value ;
    fMSdownR.bound_low = 0.;
    fMSdownR.bound_up = 10000.;

    fMSstrangeR.name = "MSstrangeR";
    fMSstrangeR.value = sumMSstrangeR/ double(nvalid);
    fMSstrangeR.error = TMath::Sqrt(sum2MSstrangeR / double(nvalid) - sqr(fMSstrangeR.value));
    if ((TMath::Abs(fMSstrangeR.error)>TMath::Abs(fMSstrangeR.value/2.)) || 
	(TMath::Abs(fMSstrangeR.error)<TMath::Abs(fMSstrangeR.value/10.))) {
      fMSstrangeR.error = TMath::Abs(fMSstrangeR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSstrangeR") >= 0)) fMSstrangeR.value = yyFittedPar[FindInFittedPar("MSstrangeR")].value ;
    fMSstrangeR.bound_low = 0.;
    fMSstrangeR.bound_up = 10000.;

    fMSbottomR.name = "MSbottomR";
    fMSbottomR.value = sumMSbottomR/ double(nvalid);
    fMSbottomR.error = TMath::Sqrt(sum2MSbottomR / double(nvalid) - sqr(fMSbottomR.value));
    if ((TMath::Abs(fMSbottomR.error)>TMath::Abs(fMSbottomR.value/2.)) || 
	(TMath::Abs(fMSbottomR.error)<TMath::Abs(fMSbottomR.value/10.))) {
      fMSbottomR.error = TMath::Abs(fMSbottomR.value/10.);
    }
    if (yyUseGivenStartValues && (FindInFittedPar("MSbottomR") >= 0)) fMSbottomR.value = yyFittedPar[FindInFittedPar("MSbottomR")].value ;
    fMSbottomR.bound_low = 0.;
    fMSbottomR.bound_up = 10000.;

    fXtau.name = "Xtau";
    fXtau.value = sumXtau/ double(nvalid);
    fXtau.value = -fMu.value*fTanBeta.value;
    fXtau.error = TMath::Sqrt(sum2Xtau / double(nvalid) - sqr(fXtau.value));
    //    if ( fXtau.value > 1000. ) fXtau.value = 1000.;
    //    if ( fXtau.value < -1000. ) fXtau.value = -1000.;
    if ((TMath::Abs(fXtau.error)>TMath::Abs(fXtau.value/2.)) || 
	(TMath::Abs(fXtau.error)<TMath::Abs(fXtau.value/10.))) {
      fXtau.error = TMath::Abs(fXtau.value/10.);
    }
    fXtau.error = 100.;
    if (yyUseGivenStartValues && (FindInFittedPar("Xtau") >= 0)) fXtau.value = yyFittedPar[FindInFittedPar("Xtau")].value ;
    if (yyBoundsOnX) {
      fXtau.bound_low = yyXscanlow;
      fXtau.bound_up = yyXscanhigh;
    } else {
      fXtau.bound_low = -10000.;
      fXtau.bound_up = 10000.;
    }

    fXtop.name = "Xtop";
    fXtop.value = sumXtop/ double(nvalid);
    fXtop.value = 0.;
    fXtop.value = -fMu.value/fTanBeta.value;
    fXtop.error = TMath::Sqrt(sum2Xtop / double(nvalid) - sqr(fXtop.value));
    if ((TMath::Abs(fXtop.error)>TMath::Abs(fXtop.value/2.)) || 
	(TMath::Abs(fXtop.error)<TMath::Abs(fXtop.value/10.))) {
      fXtop.error = TMath::Abs(fXtop.value/10.);
    }
    fXtop.error = 100.;
    if (yyUseGivenStartValues && (FindInFittedPar("Xtop") >= 0)) fXtop.value = yyFittedPar[FindInFittedPar("Xtop")].value ;
    if (yyBoundsOnX) {
      fXtop.bound_low = yyXscanlow;
      fXtop.bound_up = yyXscanhigh;
    } else {
      fXtop.bound_low = -10000.;
      fXtop.bound_up = 10000.;
    }
    if (yyScanX && (fXtop.value > yyXscanhigh)) { 
      fXtop.value = yyXscanhigh;
    }
    if (yyScanX && (fXtop.value < yyXscanlow)) { 
      fXtop.value = yyXscanlow;
    }

    fXbottom.name = "Xbottom";
    fXbottom.value = sumXbottom/ double(nvalid);
    fXbottom.value = 0.;
    fXbottom.value = -fMu.value*fTanBeta.value;
    fXbottom.error = TMath::Sqrt(sum2Xbottom / double(nvalid) - sqr(fXbottom.value));
    if ((TMath::Abs(fXbottom.error)>TMath::Abs(fXbottom.value/2.)) || 
	(TMath::Abs(fXbottom.error)<TMath::Abs(fXbottom.value/10.))) {
      fXbottom.error = TMath::Abs(fXbottom.value/10.);
    }
    fXbottom.error = 100.;
    if (yyUseGivenStartValues && (FindInFittedPar("Xbottom") >= 0)) fXbottom.value = yyFittedPar[FindInFittedPar("Xbottom")].value ;
    if (yyBoundsOnX) {
      fXbottom.bound_low = yyXscanlow;
      fXbottom.bound_up = yyXscanhigh;
    } else {
      fXbottom.bound_low = -10000.;
      fXbottom.bound_up = 10000.;
    }
    if (yyScanX && (fXbottom.value > yyXscanhigh)) { 
      fXbottom.value = yyXscanhigh;
    }
    if (yyScanX && (fXbottom.value < yyXscanlow)) { 
      fXbottom.value = yyXscanlow;
    }


//    cout << "MSelectronL = " << fMSelL.value << " +- " << fMSelL.error << endl;
//    cout << "MSmuL = " << fMSmuL.value << " +- " << fMSmuL.error << endl;
//    cout << "MStauL = " << fMStauL.value << " +- " << fMStauL.error << endl;
//    cout << "MSelectronR = " << fMSelR.value << " +- " << fMSelR.error << endl;
//    cout << "MSmuR = " << fMSmuR.value << " +- " << fMSmuR.error << endl;
//    cout << "MStauR = " << fMStauR.value << " +- " << fMStauR.error << endl;

    //-----------------------------------------------------------------
    // fill yyFittedVec
    // first fill in the measured parameters
    for (unsigned int  i=0; i < yyFittedPar.size(); i++ ) {
      par_already_found = false;
      for (unsigned int j = 0; j < fInput->GetMeasuredVector().size(); j++ ) {
	if (!yyFittedPar[i].name.compare(fInput->GetMeasuredVector()[j].name)) {
	  yyFittedVec.push_back((fInput->GetMeasuredVector())[j]);
	  if (yyUseGivenStartValues) {
	    unsigned int ilength;
	    ilength = yyFittedVec.size();
	    yyFittedVec[ilength-1].value = yyFittedPar[i].value;
	    //    	    cout << " parameter " <<  yyFittedPar[i].name << " " << yyFittedVec[ilength-1].value << " " 
	    //	 << " to value " << yyFittedPar[i].value << endl;	   
	  }
	  par_already_found = true;
	  break;
	}
      }
      if (!par_already_found) {
	if (!yyFittedPar[i].name.compare("TanBeta")) {
	  yyFittedVec.push_back(fTanBeta);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("Mu")){
	  yyFittedVec.push_back(fMu);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("M1")) {
	  yyFittedVec.push_back(fM1);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("M2")) {
	  yyFittedVec.push_back(fM2);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("M3")) {
	  yyFittedVec.push_back(fAbsM3);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSelectronL")) {
	  yyFittedVec.push_back(fMSelL);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSelectronR")) {
	  yyFittedVec.push_back(fMSelR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSmuL")) {
	  yyFittedVec.push_back(fMSmuL);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSmuR")) {
	  yyFittedVec.push_back(fMSmuR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MStauL")) {
	  yyFittedVec.push_back(fMStauL);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MStauR")) {
	  yyFittedVec.push_back(fMStauR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSupL")) {
	  yyFittedVec.push_back(fMSupL);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSupR")) {
	  yyFittedVec.push_back(fMSupR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MScharmL")) {
	  yyFittedVec.push_back(fMScharmL);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MScharmR")) {
	  yyFittedVec.push_back(fMScharmR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MStopL")) {
	  yyFittedVec.push_back(fMStopL);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MStopR")) {
	  yyFittedVec.push_back(fMStopR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSdownR")) {
	  yyFittedVec.push_back(fMSdownR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSstrangeR")) {
	  yyFittedVec.push_back(fMSstrangeR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("MSbottomR")) {
	  yyFittedVec.push_back(fMSbottomR);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("Xtau")) {
	  yyFittedVec.push_back(fXtau);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("Xtop")) {
	  yyFittedVec.push_back(fXtop);
	  par_already_found = true;
	}
	else if (!yyFittedPar[i].name.compare("Xbottom")) {
	  yyFittedVec.push_back(fXbottom);
	  par_already_found = true;
	}
	
      }
      if (!par_already_found) {
	fnew.name  = yyFittedPar[i].name;
	fnew.value = yyFittedPar[i].value;
	fnew.error = yyFittedPar[i].value*0.05;
	fnew.bound_low = -1.E+6;
	fnew.bound_up = 1.E+6;
	yyFittedVec.push_back(fnew);
      }     
    }

//    for (unsigned int  i=0; i < fInput->GetMeasuredVector().size(); i++ ) {
//	if ( !fInput->GetMeasuredVector()[i].name.compare("massTop") || //
//	     !fInput->GetMeasuredVector()[i].name.compare("massA0")  || //
//	     !fInput->GetMeasuredVector()[i].name.compare("massCharm")  || //
//	     !fInput->GetMeasuredVector()[i].name.compare("massBottom") ) {
//	  yyFittedVec.push_back((fInput->GetMeasuredVector())[i]);
//	}
//    }
//    // then fill the estimated parameters
//    yyFittedVec.push_back(fTanBeta);
//    yyFittedVec.push_back(fMu);
//    yyFittedVec.push_back(fM1);
//    yyFittedVec.push_back(fM2);
//    yyFittedVec.push_back(fAbsM3);


    cout << "fill fixed parameters " << endl;

    FillFixedParameters();

    cout << "=======================================================" << endl;
    cout << "nthrows " << nthrows << " nvalid " << nvalid << endl;
    cout << "Estimated Parameters:" << endl;
    for (unsigned int  i=0; i < yyFittedVec.size(); i++ ) {
      cout << i << " " << yyFittedVec[i].name << " = " << yyFittedVec[i].value << " +- " << yyFittedVec[i].error 
	   << " bounds " << yyFittedVec[i].bound_low << " -- " <<  yyFittedVec[i].bound_up << endl;
    }
    cout << "=======================================================" << endl;

    

}



void Fittino::calculateLoopLevelValues() 
{
  Double_t amin, edm, errdef;
  Int_t    nvpar, nparx;
  double vlow, vhigh;
  double arguments[4];
  TString parname;
  double eminus, eplus, eparab, globcc;
  TGraph* TGraphDummy;
  string TGraphTitle;
  const int    npoints = 20;
  double xarray[npoints];
  double yarray[npoints];
  double closedxarray[npoints+1];
  double closedyarray[npoints+1];
  int ierr = 0;
  char GraphName[255];
  vector <double> saved_uncertainties;

  char ntuplename[256];
  char ntupletext[256];
  char ntuplevars[2048];

  //-------------------------------------------------------------------------
  // switch all temp_nofits on
  for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
    yyMeasuredVec[i].temp_nofit = false;
    yyMeasuredVec[i].first = true;
    yyMeasuredVec[i].hasbeenset = false;
  }  


  //-------------------------------------------------------------------------
  // eventually call simulated annealing
  system ("rm SimAnnNtupFile.root");
  if (yyUseSimAnnBefore) {
    unsigned int i = 0;
    // open the ntuple file
    TFile *SimAnnNtupFile = new TFile("SimAnnNtupFile.root","UPDATE");
    sprintf ( ntuplename, "ntuple%i", i );
    sprintf ( ntupletext, "path of the simulated annealing No. %i", i );
    sprintf ( ntuplevars, "n:t:f:acc:xopt" );
    for (unsigned int j=0; j < yyFittedVec.size(); j++ ) {
      sprintf ( ntuplevars, "%s:%s", ntuplevars, yyFittedVec[j].name.c_str() );
    }
    TNtuple *simannntuple = new TNtuple(ntuplename,ntupletext,ntuplevars);
    
    simulated_annealing(0,simannntuple);
    SimAnnNtupFile->Write();
    SimAnnNtupFile->Close();
  }

  //-------------------------------------------------------------------------
  // Set Up TMinuit
  //  TFitter* fitter = new TFitter(yyFittedVec.size());
  TMinuit* fitter = new TMinuit (yyFittedVec.size());
  fitter->SetFCN(fitterFCN);
  
  arguments[0] = 0;
  fitter->mnexcm("SET PRI", arguments, 1,ierr);
  arguments[0] = 1.;
  fitter->mnexcm("SET ERR", arguments, 1,ierr);
  
  for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    cout << "adding parameter " << yyFittedVec[i].name << endl;
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    if ( yyFittedVec[i].error <= 0. ) {
      fitter->mnparm(i, yyFittedVec[i].name.c_str(),
      			   yyFittedVec[i].value, TMath::Abs(yyFittedVec[i].value/10.), 
      			   yyFittedVec[i].bound_low, yyFittedVec[i].bound_up,ierr);
      saved_uncertainties.push_back(TMath::Abs(yyFittedVec[i].value/10.));
    } else {
      fitter->mnparm(i, yyFittedVec[i].name.c_str(),
			   yyFittedVec[i].value, yyFittedVec[i].error, 
			   yyFittedVec[i].bound_low, yyFittedVec[i].bound_up,ierr);
      saved_uncertainties.push_back(yyFittedVec[i].error);
    }
  }
  

  //==================================================================================
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //==================================================================================
  // fit everything directly or first do some subsample fits?
  if (!yyFitAllDirectly) {

    //-------------------------------------------------------------------------
    // first fit just mA to mA
    if (yySepFitmA) {
      for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
	if ( !( !yyMeasuredVec[i].name.compare("massA0") || !yyMeasuredVec[i].name.compare("massh0") ) ) {
	  yyMeasuredVec[i].temp_nofit = true;
	}
      }        
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if ( yyFittedVec[i].name.compare("massA0") ) {
	  arguments[0] = i+1;
	  cout << "fixing parameter "<< yyFittedVec[i].name << endl; 
	  fitter->mnexcm("FIX", arguments, 1,ierr);
	}
      }
      
      arguments[0] = 2;
      fitter->mnexcm("SET STRATEGY", arguments, 1, ierr);
      arguments[0] = 2000;
      arguments[1] = 0.1;
      fitter->mnexcm("MINIMIZE", arguments, 2,ierr);
      
      for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
	yyMeasuredVec[i].temp_nofit = false;
      }  
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if ( yyFittedVec[i].name.compare("massA0") ) {
	  arguments[0] = i+1;
	  fitter->mnexcm("RELEASE", arguments, 1,ierr);
	}
      }
    }
    //-------------------------------------------------------------------------
    // first fit just tanb and mu
    if (yySepFitTanbMu) {
      for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
	if ( !( !yyMeasuredVec[i].name.compare("massChargino1") || !yyMeasuredVec[i].name.compare("massChargino2") ||
		!yyMeasuredVec[i].name.compare("massNeutralino1") || !yyMeasuredVec[i].name.compare("massNeutralino2") || 
		!yyMeasuredVec[i].name.compare("massNeutralino3") || !yyMeasuredVec[i].name.compare("massNeutralino4") || 
		((yyMeasuredVec[i].type == xsection) && ((TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Neutralino1"]) || 
							 (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Neutralino2"]) || 
							 (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Chargino1"]) || 
							 (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Chargino2"]) ) ) ) ) {
	  yyMeasuredVec[i].temp_nofit = true;
	}
      }        
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if ( !( !yyFittedVec[i].name.compare("TanBeta") || !yyFittedVec[i].name.compare("Mu") ) ) {
	  arguments[0] = i+1;
	  cout << "fixing parameter "<< yyFittedVec[i].name << endl; 
	  fitter->mnexcm("FIX", arguments, 1,ierr);
	}
      }
      
      arguments[0] = 2;
      fitter->mnexcm("SET STRATEGY", arguments, 1, ierr);
      arguments[0] = 2000;
      arguments[1] = 0.1;
      fitter->mnexcm("MINIMIZE", arguments, 2,ierr);
      
      for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
	yyMeasuredVec[i].temp_nofit = false;
      }  
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if ( !( !yyFittedVec[i].name.compare("TanBeta") || !yyFittedVec[i].name.compare("Mu") ) ) {
	  arguments[0] = i+1;
	  fitter->mnexcm("RELEASE", arguments, 1,ierr);
	}
      }
      
    }


//  //-------------------------------------------------------------------------
//  // Perform the fit
//  // first fix everything but TanBeta, Mu, M1, M2, M3, massA0
//  for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
//    if ( !( !yyFittedVec[i].name.compare("TanBeta") || !yyFittedVec[i].name.compare("Mu") 
//	    || !yyFittedVec[i].name.compare("M1") || !yyFittedVec[i].name.compare("M2")
//	    || !yyFittedVec[i].name.compare("M3") || !yyFittedVec[i].name.compare("massA0") ) ) {
//      arguments[0] = i+1;
//      cout << "fixing parameter "<< yyFittedVec[i].name << endl; 
//      fitter->mnexcm("FIX", arguments, 1,ierr);
//    }
//  }
//  for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
//    if ( !( !yyMeasuredVec[i].name.compare("massChargino1") || !yyMeasuredVec[i].name.compare("massChargino2") ||
//	    !yyMeasuredVec[i].name.compare("massNeutralino1") || !yyMeasuredVec[i].name.compare("massNeutralino2") || 
//	    !yyMeasuredVec[i].name.compare("massGluino") || !yyMeasuredVec[i].name.compare("massA0") ||  
//	    !yyMeasuredVec[i].name.compare("massh0") || 
//	    ((yyMeasuredVec[i].type == xsection) && ((TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Neutralino1"]) || 
//						     (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Neutralino2"]) || 
//						     (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Chargino1"]) || 
//						     (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Chargino2"]) ) ) ) ) {
//      yyMeasuredVec[i].temp_nofit = true;
//    }
//  }  
//  
//  arguments[0] = 2000;
//  arguments[1] = 0.1;
//  fitter->mnexcm("MINIMIZE", arguments, 2,ierr);
//  for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
//    if ( !( !yyFittedVec[i].name.compare("TanBeta") || !yyFittedVec[i].name.compare("Mu") 
//	    || !yyFittedVec[i].name.compare("M1") || !yyFittedVec[i].name.compare("M2")
//	    || !yyFittedVec[i].name.compare("M3") || !yyFittedVec[i].name.compare("massA0") ) ) {
//      arguments[0] = i+1;
//      fitter->mnexcm("RELEASE", arguments, 1,ierr);
//    }
//  }
//
//  for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
//    yyMeasuredVec[i].temp_nofit = false;
//  }  
  
    for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
      if ( !( !yyMeasuredVec[i].name.compare("massSelectronL") ||  !yyMeasuredVec[i].name.compare("massSelectronR") || 
	      !yyMeasuredVec[i].name.compare("massSmuL") ||  !yyMeasuredVec[i].name.compare("massSmuR") || 
	      !yyMeasuredVec[i].name.compare("massStau1") ||  !yyMeasuredVec[i].name.compare("massStau2") ||
	      !yyMeasuredVec[i].name.compare("massSnueL") || !yyMeasuredVec[i].name.compare("massSnumuL") || 
	      !yyMeasuredVec[i].name.compare("massSnutauL") || 
	      ((yyMeasuredVec[i].type == xsection) && ((TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SelectronL"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SelectronR"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SmuL"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SmuR"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Stau1"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Stau2"]) ) ) ) ) {
	yyMeasuredVec[i].temp_nofit = true;
      }
    } 
    
//    //-------------------------------------------------------------------------
//    // scan Xtau
//    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
//	if (   !yyFittedVec[i].name.compare("Xtau") ) {
//	  arguments[0] = i+1;
//	  cout << "scanning parameter no " << i+1 << " " << yyFittedVec[i].name << endl;
//	  arguments[1] = 100;
//	  arguments[2] = yyFittedVec[i].bound_low;
//	  arguments[3] = yyFittedVec[i].bound_up;
//	  fitter->mnexcm("SCAN", arguments, 4,ierr);
//	  break;
//	}
//    }
    
    //-------------------------------------------------------------------------
    // Perform the fit
    // Then fix everything but the SLepton Sector
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      if ( !( !yyFittedVec[i].name.compare("Xtau") || !yyFittedVec[i].name.compare("MSelectronL") 
	      || !yyFittedVec[i].name.compare("MSmuL") || !yyFittedVec[i].name.compare("MStauL")
	      || !yyFittedVec[i].name.compare("MSelectronR") || !yyFittedVec[i].name.compare("MSmuR")
	      || !yyFittedVec[i].name.compare("MStauR") ) ) {
	arguments[0] = i+1;
	cout << "fixing parameter "<< yyFittedVec[i].name << endl; 
	fitter->mnexcm("FIX", arguments, 1,ierr);
      }
    }
    arguments[0] = 2;
    fitter->mnexcm("SET STRATEGY", arguments, 1, ierr);
    arguments[0] = 2000;
    arguments[1] = 0.1;
    fitter->mnexcm("MINIMIZE", arguments, 2,ierr);
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      if ( !( !yyFittedVec[i].name.compare("Xtau") || !yyFittedVec[i].name.compare("MSelectronL") 
	      || !yyFittedVec[i].name.compare("MSmuL") || !yyFittedVec[i].name.compare("MStauL")
	      || !yyFittedVec[i].name.compare("MSelectronR") || !yyFittedVec[i].name.compare("MSmuR")
	      || !yyFittedVec[i].name.compare("MStauR") ) ) {
	arguments[0] = i+1;
	fitter->mnexcm("RELEASE", arguments, 1,ierr);
      }
    }
    
    for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
      yyMeasuredVec[i].temp_nofit = false;
    }  
    
    
    
    for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
      if ( !( !yyMeasuredVec[i].name.compare("massSupL") ||  !yyMeasuredVec[i].name.compare("massSupR") || 
	      !yyMeasuredVec[i].name.compare("massScharmL") ||  !yyMeasuredVec[i].name.compare("massScharmR") || 
	      !yyMeasuredVec[i].name.compare("massStop1") ||  !yyMeasuredVec[i].name.compare("massStop2") ||
	      !yyMeasuredVec[i].name.compare("massSdownL") || !yyMeasuredVec[i].name.compare("massSdownR") || 
	      !yyMeasuredVec[i].name.compare("massSstrangeL") || !yyMeasuredVec[i].name.compare("massSstrangeR") || 
	      !yyMeasuredVec[i].name.compare("massSbottom1") || !yyMeasuredVec[i].name.compare("massSbottom2") || 
	      ((yyMeasuredVec[i].type == xsection) && ((TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SupL"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SupR"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SdownL"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SdownR"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["ScharmL"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["ScharmR"]) || 
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SstrangeL"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["SstrangeR"]) || 
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Stop1"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Stop2"]) || 
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Sbottom1"]) ||
						       (TMath::Abs(yyMeasuredVec[i].products[0])==yyParticleIDs["Sbottom2"]) ) ) ) ) {
	yyMeasuredVec[i].temp_nofit = true;
      }
    } 
    
    
    if (yyScanX) {
      //-------------------------------------------------------------------------
      // scan Xtop
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if (  !yyFittedVec[i].name.compare("Xtop") ) {
	  arguments[0] = i+1;
	  cout << "scanning parameter no " << i+1 << " " << yyFittedVec[i].name << endl;
	  arguments[1] = 100;
	  if (yyBoundsOnX) {
	    arguments[2] = yyFittedVec[i].bound_low;
	    arguments[3] = yyFittedVec[i].bound_up;
	  } else {
	    arguments[2] = yyXscanlow;
	    arguments[3] = yyXscanhigh;	    
	  }
	  fitter->mnexcm("SCAN", arguments, 4,ierr);
	  break;
	}
      }
      
      //-------------------------------------------------------------------------
      // scan Xbottom
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if (  !yyFittedVec[i].name.compare("Xbottom") ) {
	  arguments[0] = i+1;
	  cout << "scanning parameter no " << i+1 << " " << yyFittedVec[i].name << endl;
	  arguments[1] = 100;
	  if (yyBoundsOnX) {
	    arguments[2] = yyFittedVec[i].bound_low;
	    arguments[3] = yyFittedVec[i].bound_up;
	  } else {
	    arguments[2] = yyXscanlow;
	    arguments[3] = yyXscanhigh;	    
	  }
	  fitter->mnexcm("SCAN", arguments, 4,ierr);
	  break;
	}
      }
    }

    //-------------------------------------------------------------------------
    // Perform the fit
    // Then fix everything but the SQuark Sector
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      if ( !( !yyFittedVec[i].name.compare("Xtop") || !yyFittedVec[i].name.compare("Xbottom")   
	      || !yyFittedVec[i].name.compare("MSupL") || !yyFittedVec[i].name.compare("MScharmL")
	      || !yyFittedVec[i].name.compare("MStopL") || !yyFittedVec[i].name.compare("MSupR")
	      || !yyFittedVec[i].name.compare("MScharmR") || !yyFittedVec[i].name.compare("MStopR") 
	      || !yyFittedVec[i].name.compare("MSdownR")
	      || !yyFittedVec[i].name.compare("MSstrangeR") || !yyFittedVec[i].name.compare("MSbottomR") ) ) {
	arguments[0] = i+1;
	cout << "fixing parameter "<< yyFittedVec[i].name << endl; 
	fitter->mnexcm("FIX", arguments, 1,ierr);
      }
    }
    arguments[0] = 2;
    fitter->mnexcm("SET STRATEGY", arguments, 1, ierr);
    arguments[0] = 2000;
    arguments[1] = 0.1;
    fitter->mnexcm("MINIMIZE", arguments, 2,ierr);
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      if ( !( !yyFittedVec[i].name.compare("Xtop") || !yyFittedVec[i].name.compare("Xbottom")  
	      || !yyFittedVec[i].name.compare("MSupL") || !yyFittedVec[i].name.compare("MScharmL")
	      || !yyFittedVec[i].name.compare("MStopL") || !yyFittedVec[i].name.compare("MSupR")
	      || !yyFittedVec[i].name.compare("MScharmR") || !yyFittedVec[i].name.compare("MStopR") 
	      || !yyFittedVec[i].name.compare("MSdownR")
	      || !yyFittedVec[i].name.compare("MSstrangeR") || !yyFittedVec[i].name.compare("MSbottomR") ) ) {
	arguments[0] = i+1;
	fitter->mnexcm("RELEASE", arguments, 1,ierr);
      }
    }
    
    //==================================================================================================
    //==================================================================================================
    // now release all
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      arguments[0] = i+1;
      fitter->mnexcm("RELEASE", arguments, 1,ierr);
    }
    for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
      yyMeasuredVec[i].temp_nofit = false;
    }  
    
    if (yySepFitTanbX) {
      //-------------------------------------------------------------------------
      // Perform the fit
      // Then fix everything but the A's and MStop and TanBeta
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if ( !( !yyFittedVec[i].name.compare("Xtop") || !yyFittedVec[i].name.compare("Xbottom") || !yyFittedVec[i].name.compare("TanBeta") 
		|| !yyFittedVec[i].name.compare("MStopL") || !yyFittedVec[i].name.compare("MStopR") ) ) { 
	  arguments[0] = i+1;
	  cout << "fixing parameter "<< yyFittedVec[i].name << endl; 
	  fitter->mnexcm("FIX", arguments, 1,ierr);
	}
      }
      arguments[0] = 2;
      fitter->mnexcm("SET STRATEGY", arguments, 1, ierr);
      arguments[0] = 20000;
      arguments[1] = 0.1;
      fitter->mnexcm("MINIMIZE", arguments, 2,ierr);
      for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
	if ( !( !yyFittedVec[i].name.compare("Xtop") || !yyFittedVec[i].name.compare("Xbottom") || !yyFittedVec[i].name.compare("TanBeta") 
		|| !yyFittedVec[i].name.compare("MStopL") || !yyFittedVec[i].name.compare("MStopR") ) ) {
	  arguments[0] = i+1;
	  fitter->mnexcm("RELEASE", arguments, 1,ierr);
	}
      }
    }
    
  } // end of if(!yyFitAllDirectly)


  //==================================================================================================
  //==================================================================================================
  //--------------------------------------------------------------------------
  // then release all the fixed parameters
  for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
    arguments[0] = i+1;
    fitter->mnexcm("RELEASE", arguments, 1,ierr);
  }
  for (unsigned int i = 0; i < yyMeasuredVec.size(); i++ ) {
    yyMeasuredVec[i].temp_nofit = false;
  }  


  //--------------------------------------------------------------------------
  // perform the final fit with all parameters free
  arguments[0] = 2;
  fitter->mnexcm("SET STRATEGY", arguments, 1, ierr);
  arguments[0] = 200000;
  arguments[1] = 0.1;
  for (unsigned int i=0; i<yyNumberOfMinimizations; i++) {
    if (i>0) {
      // read the central values and uncertainties
      for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
	fitter->mnpout(j,parname,yyFittedVec[j].value,yyFittedVec[j].error,vlow,vhigh,ierr);
      }      
      //-------------------------------------------------------------------------
      // eventuaslly call simulated annealing
      if (yyUseSimAnnWhile) {
	// open the ntuple file
	TFile *SimAnnNtupFile = new TFile("SimAnnNtupFile.root","UPDATE");
	sprintf ( ntuplename, "ntuple%i", i );
	sprintf ( ntupletext, "path of the simulated annealing No. %i", i );
	sprintf ( ntuplevars, "n:t:f:acc:xopt" );
	for (unsigned int j=0; j < yyFittedVec.size(); j++ ) {
	  sprintf ( ntuplevars, "%s:%s", ntuplevars, yyFittedVec[j].name.c_str() );
	}
	TNtuple *simannntuple = new TNtuple(ntuplename,ntupletext,ntuplevars);
	// check the status of the minimization
	fitter->mnstat(amin, edm, errdef, nvpar, nparx, ierr);	
	if (ierr == 3) {
	  for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
	    yyFittedVec[j].error = 2.*yyFittedVec[j].error;
	  }
	} else {
	  // error matrix is not accurate
	  for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
	    yyFittedVec[j].error = saved_uncertainties[j];
	  }
	}
	simulated_annealing(i,simannntuple);
	// close the ntuple
	SimAnnNtupFile->Write();
	SimAnnNtupFile->Close();
      }
      // reset the uncertainties
      for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
	cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	cout << "resetting parameter errors " << yyFittedVec[j].name << endl;
	cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	fitter->mnparm(j, yyFittedVec[j].name.c_str(),
		       yyFittedVec[j].value, saved_uncertainties[j], 
		       yyFittedVec[j].bound_low, yyFittedVec[j].bound_up,ierr);
      }
    }
    fitter->mnexcm("MINIMIZE", arguments, 2,ierr);
  }
  cout << "returning from MINIMIZE, return value "<< ierr << endl;
  if (yyUseHesse && !yyUseMinos) {
    fitter->mnexcm("HESSE", arguments, 1,ierr);
    cout << "returning from HESSE, return value "<< ierr << endl;
  }
  if (yyUseMinos) {
    fitter->mnexcm("HESSE", arguments, 1,ierr);
    arguments[0] = yyErrDef;
    fitter->mnexcm("SET ERR", arguments, 1,ierr);
    arguments[0] = 200000;
    fitter->mnexcm("MINOS", arguments, 1,ierr);
  }  

  //-------------------------------------------------------------------------
  // Get Results
  fitter->mnstat(amin, edm, errdef, nvpar, nparx, ierr);
  fchisq = amin;
  gchisq = amin;
  gstat  = ierr;
  fFittedCovarianceMatrix = new double[yyFittedVec.size()*yyFittedVec.size()];
  fitter->mnemat(fFittedCovarianceMatrix,yyFittedVec.size());
  //fitter->PrintResults(3, amin);
  
  fSavedFittedCovarianceMatrix = new TMatrixD(yyFittedVec.size(),yyFittedVec.size());
  fSavedFittedCorrelationMatrix = new TMatrixD(yyFittedVec.size(),yyFittedVec.size());
  

  cout << "calculating the covariance and correlation matrix "<<endl;
   //-------------------------------------------------------------------------
  // Calculate Covariance Matrix
  for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
    for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
      (*fSavedFittedCovarianceMatrix)(i,j) = fFittedCovarianceMatrix[i*yyFittedVec.size()+j];
    }
  }
  for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
    for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
      (*fSavedFittedCorrelationMatrix)(i,j) = (*fSavedFittedCovarianceMatrix)(i,j)/
	(TMath::Sqrt((*fSavedFittedCovarianceMatrix)(i,i)*(*fSavedFittedCovarianceMatrix)(j,j)));
    }
  }

  cout << "get the errors" << endl;
  //-------------------------------------------------------------------------
  // Get Errors
  // return values
  for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
    fitter->mnpout(i,parname,yyFittedVec[i].value,yyFittedVec[i].error,vlow,vhigh,ierr);
    if (yyUseMinos) {
      fitter->mnerrs(i, eplus, eminus, eparab, globcc);
      //      yyFittedVec[i].error = TMath::Max(TMath::Abs(eplus), TMath::Abs(eminus));
      yyFittedVec[i].positive_error = eplus;
      yyFittedVec[i].negative_error = eminus;
      // recalculate Error because of yyErrDef
      yyFittedVec[i].error = yyFittedVec[i].error / sqrt(yyErrDef);
      yyFittedVec[i].positive_error = yyFittedVec[i].positive_error / sqrt(yyErrDef);
      yyFittedVec[i].negative_error = yyFittedVec[i].negative_error / sqrt(yyErrDef);
    }
  }

  if (yyUseMinos) {
    // recalculate Covariance Matrix because of yyErrDef
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
	(*fSavedFittedCovarianceMatrix)(i,j) = (*fSavedFittedCovarianceMatrix)(i,j)/yyErrDef;
      }
    }    
  }

  cout << "write the preliminay output file" << endl;
  writeResults("fittino.out.intermediate");

  cout << "eventually get the contours " << endl;
  //-------------------------------------------------------------------------
  // Get the contours with only two parameters free
  if (yyGetContours) {
    TFile* hfile1 = 0;
    hfile1 = new TFile("FitContours_2params_free.root","RECREATE","Fit Contours of the SUSY Paramater Fit with only two params free");

    arguments[0] = 1.0;
    fitter->mnexcm("SET ERR", arguments, 1,ierr);
    //    fitter->mnexcm("MINOS", arguments, 1,ierr);
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
	if (i < j) {
	  // fix all parameters but i and j 
	  for (unsigned int k = 0; k < yyFittedVec.size(); k++ ) {
	    if ( !( (k==i) || (k==j) ) ) {
	      arguments[0] = k+1;
	      cout << "fixing parameter "<< yyFittedVec[k].name << endl; 
	      fitter->mnexcm("FIX", arguments, 1,ierr);
	    }
	  }
	  // now get the contour
	  Int_t ncontpoints = 0;
	  fitter->mncont(i,j,npoints,xarray,yarray,ierr);
	  if (ierr < 4) {
  	      cerr<<"ierr returned from mncont = "<<ierr<<endl;
	      for (unsigned int k = 0; k < yyFittedVec.size(); k++ ) {
		arguments[0] = k+1;
		fitter->mnexcm("RELEASE", arguments, 1,ierr);
	      }
	      continue;
	  }
	  ncontpoints = ierr;
	  for (unsigned int k = 0; k<ncontpoints; k++) {
	    closedxarray[k] = xarray[k];
	    closedyarray[k] = yarray[k];
	    cout<<"x["<<k<<"] = "<<xarray[k]<<"   y["<<k<<"] = "<<yarray[k]<<endl;
	  }
	  closedxarray[ncontpoints] = xarray[0];
	  closedyarray[ncontpoints] = yarray[0];
	  //-------------------------------------------------------------------------
	  // book graphs
	  // create one TGraph object for each pair of variables
	  TGraphTitle = "2Param Fit Contour of ";
	  TGraphTitle.append(yyFittedVec[j].name);
	  TGraphTitle.append(" vs ");
	  TGraphTitle.append(yyFittedVec[i].name);
	  TGraphDummy = new TGraph(ncontpoints+1,closedxarray,closedyarray);
	  TGraphDummy->SetTitle(TGraphTitle.c_str());
	  sprintf(GraphName,"Graph%i",i*yyFittedVec.size()+j);
	  TGraphDummy->SetName(GraphName);
	  TGraphDummy->Write();
	  // free all other parameters
	  for (unsigned int k = 0; k < yyFittedVec.size(); k++ ) {
	    arguments[0] = k+1;
	    fitter->mnexcm("RELEASE", arguments, 1,ierr);
	  }
	}
      }
    }
    cout << "close the output root files " << endl;
    //-------------------------------------------------------------------------
    // Close the Output file
    // close histograms
    hfile1->Write();
    hfile1->Close();
  }

  //-------------------------------------------------------------------------
  // Get Contours
  if (yyGetContours) {
    TFile* hfile = 0;
    hfile = new TFile("FitContours.root","RECREATE","Fit Contours of the SUSY Paramater Fit");

    arguments[0] = 1.0;
    fitter->mnexcm("SET ERR", arguments, 1,ierr);
    //    fitter->mnexcm("MINOS", arguments, 1,ierr);
    for (unsigned int i = 0; i < yyFittedVec.size(); i++ ) {
      for (unsigned int j = 0; j < yyFittedVec.size(); j++ ) {
	if (i < j) {
	  Int_t ncontpoints = 0;
	  fitter->mncont(i,j,npoints,xarray,yarray,ierr);
	  if (ierr < 4) {
  	      cerr<<"ierr returned from mncont = "<<ierr<<endl;
	      continue;
	  }
	  ncontpoints = ierr;
	  for (unsigned int k = 0; k<ncontpoints; k++) {
	    closedxarray[k] = xarray[k];
	    closedyarray[k] = yarray[k];
	    cout<<"x["<<k<<"] = "<<xarray[k]<<"   y["<<k<<"] = "<<yarray[k]<<endl;
	  }
	  closedxarray[ncontpoints] = xarray[0];
	  closedyarray[ncontpoints] = yarray[0];
	  //-------------------------------------------------------------------------
	  // book graphs
	  // create one TGraph object for each pair of variables
	  TGraphTitle = "Fit Contour of ";
	  TGraphTitle.append(yyFittedVec[j].name);
	  TGraphTitle.append(" vs ");
	  TGraphTitle.append(yyFittedVec[i].name);
	  TGraphDummy = new TGraph(ncontpoints+1,closedxarray,closedyarray);
	  TGraphDummy->SetTitle(TGraphTitle.c_str());
	  sprintf(GraphName,"Graph%i",i*yyFittedVec.size()+j);
	  TGraphDummy->SetName(GraphName);
	  TGraphDummy->Write();
	}
      }
    }
    cout << "close the output root files " << endl;
    //-------------------------------------------------------------------------
    // Close the Output file
    // close histograms
    hfile->Write();
    hfile->Close();
  }

  //  arguments[0] = 0;
  //  fitter->mnexcm("RETURN", arguments, 0,ierr);

  cout << "do not delete fitter" << endl;

  // delete fitter;
  // cout << "fitter deleted" << endl;
}


void Fittino::writeResults(const char* filename)
{
    FILE* file = fopen(filename,"w");

    time_t now;
    time(&now);

    tm* tm_now = localtime(&now);
    int sec  = tm_now->tm_sec;
    int min  = tm_now->tm_min;
    int hour = tm_now->tm_hour;
    int mday = tm_now->tm_mday;
    int mon  = tm_now->tm_mon;
    int year = tm_now->tm_year + 1900;
    int wday = tm_now->tm_wday;

    cout << " entering writeResults" << endl;

    string weekday;
    switch (wday) {
	case 0:
	    weekday = "Sunday";
	    break;
	case 1:
	    weekday = "Monday";
	    break;
	case 2:
	    weekday = "Tuesday";
	    break;
	case 3:
	    weekday = "Wednesday";
	    break;
	case 4:
	    weekday = "Thursday";
	    break;
	case 5:
	    weekday = "Friday";
	    break;
	case 6:
	    weekday = "Saterday";
	    break;
	default:
	    cerr<<"Unknown weekday"<<endl;
	    exit(EXIT_FAILURE);
    }

    string month;
    switch (mon) {
	case 0:
	    month = "January";
	    break;
	case 1:
	    month = "February";
	    break;
	case 2:
	    month = "March";
	    break;
	case 3:
	    month = "April";
	    break;
	case 4:
	    month = "May";
	    break;
	case 5:
	    month = "June";
	    break;
	case 6:
	    month = "July";
	    break;
	case 7:
	    month = "August";
	    break;
	case 8:
	    month = "September";
	    break;
	case 9:
	    month = "October";
	    break;
	case 10:
	    month = "November";
	    break;
	case 11:
	    month = "December";
	    break;
	default:
	    cerr<<"Unknown month"<<endl;
	    exit(EXIT_FAILURE);
    }


    fprintf(file,"#################### Fittino Fit Summary ####################\n");
    fprintf(file,"created by Fittino version 0.0.1\n");
    fprintf(file,"on %s, %s %02d, %d at %02d:%02d:%02d\n", weekday.c_str(), month.c_str(),
	      mday, year, hour, min, sec);
    fprintf(file,"\n");
    fprintf(file,"Input values:\n");
    fprintf(file,"=============\n");
    for (int i=0; i<yyMeasuredVec.size(); i++) {
        if (yyMeasuredVec[i].nofit) continue;
        if (yyMeasuredVec[i].name.length() < 20)
 	    fprintf(file,"%20s   %12g +- %12g\n", yyMeasuredVec[i].name.c_str(), yyMeasuredVec[i].value,
		      yyMeasuredVec[i].error);
	else
	  fprintf(file,"%s      %f +- %f\n",yyMeasuredVec[i].name.c_str(),
		  yyMeasuredVec[i].value,yyMeasuredVec[i].error);
    }
    fprintf(file,"\n");
    fprintf(file,"Covariance matrix for input value:\n");
    fprintf(file,"==================================\n");
    fprintf(file,"                     ");
    for (int i=0; i<yyMeasuredVec.size(); i++) {
      if (!yyMeasuredVec[i].nofit) {
	fprintf(file,"%12s", yyMeasuredVec[i].name.c_str());
      }
    }
    fprintf(file,"\n");
    for (int i=0; i<yyMeasuredVec.size(); i++) {
      if (!yyMeasuredVec[i].nofit) {
      fprintf(file,"%20s ", yyMeasuredVec[i].name.c_str());
      for (int j=0; j<yyMeasuredVec.size(); j++) {
	if (!yyMeasuredVec[j].nofit) {
	  fprintf(file,"%12f", fInput->GetMeasuredCorrelationMatrix().GetCovariance(i, j));
	}
      }
      fprintf(file,"\n");
      }
    }
    fprintf(file,"\n");
    fprintf(file,"Fixed values: \n");
    fprintf(file,"=============\n");
    for (int i=0; i<yyFixedPar.size(); i++) {
        if (yyFixedPar[i].name.length() < 20)
 	    fprintf(file,"%20s   %12g\n", yyFixedPar[i].name.c_str(), yyFixedPar[i].value);
	else
	    fprintf(file,"%s      %f\n",yyFixedPar[i].name.c_str(),yyFixedPar[i].value);
    }
    fprintf(file,"\n");
    fprintf(file,"Fitted values:\n");
    fprintf(file,"==============\n");
    for (int i=0; i<yyFittedVec.size(); i++) {
      if (yyFittedVec[i].name.length() < 20) {
	fprintf(file,"%20s   %12g +- %12g", yyFittedVec[i].name.c_str(), yyFittedVec[i].value,
		yyFittedVec[i].error);
	if (yyUseMinos) {
	  fprintf(file," (%12g,%12g)\n",yyFittedVec[i].negative_error,yyFittedVec[i].positive_error);
	} else {
	  fprintf(file,"\n");
	}
      }
      else {
	fprintf(file,"%s      %f +- %f",yyFittedVec[i].name.c_str(),yyFittedVec[i].value,yyFittedVec[i].error);
	if (yyUseMinos) {
	  fprintf(file," (%f,%f)\n",yyFittedVec[i].negative_error,yyFittedVec[i].positive_error);
	} else {
	  fprintf(file,"\n");
	}
      }
    }
    fprintf(file,"\n");
    if (fSavedFittedCovarianceMatrix) {
	fprintf(file,"Covariance matrix for fitted parameters:\n");
	fprintf(file,"========================================\n");
        fprintf(file,"\n");
	fprintf(file,"                   ");
	for (int i=0; i<yyFittedVec.size(); i++) {
	  fprintf(file, "%12s", yyFittedVec[i].name.c_str());
	}
	fprintf(file,"\n");
	for (int i=0; i<yyFittedVec.size(); i++) {
	  fprintf(file,"%20s ", yyFittedVec[i].name.c_str());
	  for (int j=0; j<yyFittedVec.size(); j++) {
	    fprintf(file, "%12g", (*fSavedFittedCovarianceMatrix)(i,j) );
	  }
	  fprintf(file, " \n");
	}
	fprintf(file, " \n");
    }
    fprintf(file,"\n");
    if (fSavedFittedCorrelationMatrix) {
	fprintf(file,"Correlation matrix for fitted parameters:\n");
	fprintf(file,"=========================================\n");
        fprintf(file,"\n");
	fprintf(file,"                   ");
	for (int i=0; i<yyFittedVec.size(); i++) {
	  fprintf(file, "%12s", yyFittedVec[i].name.c_str());
	}
	fprintf(file,"\n");
	for (int i=0; i<yyFittedVec.size(); i++) {
	  fprintf(file,"%20s ", yyFittedVec[i].name.c_str());
	  for (int j=0; j<yyFittedVec.size(); j++) {
	    fprintf(file, "%12g", (*fSavedFittedCorrelationMatrix)(i,j) );
	  }
	  fprintf(file, " \n");
	}
	fprintf(file, " \n");
    }

    cout << "printing chisq result" << endl;

    fprintf(file,"\n");
    fprintf(file,"Chisq of the fit:\n");
    fprintf(file,"=================\n");
    fprintf(file,"\n");
    fprintf(file,"                   chisq = %f\n\n",fchisq);
    fprintf(file,"\n");
    fprintf(file,"Status of the minimization:\n");
    fprintf(file,"===========================\n");
    fprintf(file,"\n");
    if (gstat==0) {
      fprintf(file,"                   Error Matrix not calculated at all\n\n");
    } else if (gstat==1) {
      fprintf(file,"                   Error Matrix available in diagonal approximation only\n\n");
    } else if (gstat==2) {
      fprintf(file,"                   Error Matrix available but forced positive-definite\n\n");
    } else if (gstat==3) {
      fprintf(file,"                   Error Matrix accurate\n\n");
    } else {
      fprintf(file,"                   No information on error matrix accuracy available\n\n");
    }
    fprintf(file, "##################### End of Fit Summary ####################\n");

    cout << "closing output file" << endl;
    fclose(file);
}


// ************************************************************
// fitterFCN
//
// calculates: The chi^2 from the measured values and the SUSPECT prediction
// needs:      all
//
// P. Bechtle, P. Wienemann, 18.09.03
//
// ************************************************************

void fitterFCN(Int_t &, Double_t *, Double_t &f, Double_t *x, Int_t iflag) 
{
  // cout<<"fitterFCN called"<<endl;
  //  niterations++;
  int rc = 0;
  int nobs = 0;
  int ncorr = 0;
  int nbound = 0;

  f = 0.;

  // set values
  yyGeneratorError = false;
  yyParseError     = 0;

  for (unsigned int i = 0; i < yyMeasuredVec.size(); i++) {
    yyMeasuredVec[i].theoset = false;
  }
  // start loop over cross sections and polarisations
  // cout << "having to call Generator " <<  CrossSectionProduction.size() << " times" << endl; 
  // for (unsigned int j = 0; j < CrossSectionProduction.size(); j++) {
    // HERE: WRITE THE LES HOUCHES FILE
  WriteLesHouches(x);
    

//  fitterMassA0.value =  436.2;
//  fitterMassTop.value = x[6];
//  fitterMassSelL.value = 262.5;// x[6];
//  fitterMassSelR.value = 218.6; // x[7];
//  fitterMassSupL.value = 519.1; // x[8];
//  fitterMassSupR.value = 429.6; // x[9];
//  fitterA.value        = -500.0;// x[10];

  if (yyGenerator == SUSPECT) {
    callSuspect();
  }
  else if (yyGenerator == SPHENO) {
    rc = callSPheno();
  }
  if (rc == 2) {
    cout << "SIGINT received in SPheno, exiting" << endl;
    exit (2);
  }
  if (rc > 0) {
    cerr << "Exiting fitterFCN because of problem in SPheno run" << endl;
    f = 111111111111.;
    cout << " f = " << f << endl;
    return;        
  }
  
  // HERE: READ THE LES HOUCHES FILE
  ReadLesHouches();
  
  if (yyGeneratorError) {
    cerr << "Exiting fitterFCN because LesHouches outfile did not exist" << endl;
    f = 111111111111.;
    cout << " f = " << f << endl;
    return;    
  }
  if (yyParseError) {
    cerr << "Exiting fitterFCN because of a parse error in yacc" << endl;
    f = 111111111111.;
    cout << " f = " << f << endl;
    return;
  }
  if (rc > 0) {
    cerr << "Exiting fitterFCN because of return value "<< rc << " from the Generator" << endl;
    f = 111111111111.;
    cout << " f = " << f << endl;
    return;
  }

    // end loop over crosssections and polarisations:
    //  }

  // HERE: FIND OBSERVABLES BELONGING TO THE MEASUREMENTS IN yyMeasuredVec
  // CALCULATE CHISQ
  cout << "calculating chisq" << endl;
  for (unsigned int i = 0; i < yyMeasuredVec.size(); i++) {
    for (unsigned int j = 0; j < yyMeasuredVec.size(); j++) {
      //    cout << "looking at "<< yyMeasuredVec[i].name << endl;
      if ((yyMeasuredVec[i].theoset && yyMeasuredVec[j].theoset) && (!yyMeasuredVec[i].nofit && !yyMeasuredVec[j].nofit )
	  && (!yyMeasuredVec[i].temp_nofit && !yyMeasuredVec[j].temp_nofit ) ) {
	// cout << "theo ist set" << endl;
	if ((yyMeasuredVec[i].error > 0.) && (yyMeasuredVec[j].error > 0.) ) {
	  // find corresponding measurement for prediction... 
	  f += (yyMeasuredVec[i].theovalue-yyMeasuredVec[i].value)
	    * yyMeasuredCorrelationMatrix.GetInverseCovariance(i,j)
	    * (yyMeasuredVec[j].theovalue-yyMeasuredVec[j].value);
	  if ( (i==j) && (yyMeasuredCorrelationMatrix.GetInverseCovariance(i,j) > 1.E-12) ) {
	    nobs++;
	  } else if (yyMeasuredCorrelationMatrix.GetInverseCovariance(i,j) > 1.E-12) {
	    ncorr++;
	  }
	  if ((yyMeasuredCorrelationMatrix.GetInverseCovariance(i,j)>1.E-12) || 
	      (yyMeasuredCorrelationMatrix.GetInverseCovariance(i,j)<-1.E-12)) {
	    cout << i << " " << j << "using obs " << yyMeasuredVec[i].name << " = " << yyMeasuredVec[i].value
		 << "+-" << sqrt(yyMeasuredCorrelationMatrix.GetCovariance(i,j)) 
		 << " (" << (TMath::Abs(yyMeasuredVec[i].value-yyMeasuredVec[i].theovalue))*sqrt(yyMeasuredCorrelationMatrix.GetInverseCovariance(i,j)) << ") " << " at theovalue = " 
		 << yyMeasuredVec[i].theovalue<< endl;
	  }
	} else {
	  if (i == j) {
	    if (yyMeasuredVec[i].theovalue<yyMeasuredVec[i].bound_low ) {
	      f += sqr((yyMeasuredVec[i].theovalue-yyMeasuredVec[i].bound_low)/(0.01*yyMeasuredVec[i].bound_low));
	    } else if (yyMeasuredVec[i].theovalue>yyMeasuredVec[i].bound_up ) {
	      f += sqr((yyMeasuredVec[i].theovalue-yyMeasuredVec[i].bound_up)/(0.01*yyMeasuredVec[i].bound_up));
	    }
	    nbound++;
	  }
	}
      }
    }
  }

  cout << " chisq = " << f << " with " << nobs << " observables (" << ncorr/2 << " correlated), "
       << yyFittedVec.size() << " parameters and " << nbound << " bounds" << endl;

//  f = ( sqr(sRecMassChargino1-fMassChargino1)/sqr(fMassChargino1.error)
//	+ sqr(sRecMassChargino2-fMassChargino2)/sqr(fMassChargino2.error)
//	+ sqr(sRecMassNeutralino1-fMassNeutralino1)/sqr(fMassNeutralino1.error)
//	+ sqr(sRecMassNeutralino2-fMassNeutralino2)/sqr(fMassNeutralino2.error)
//	+ sqr(sRecMassGluino-fMassGluino)/sqr(fMassGluino.error)
//	+ sqr(sRecCos2PhiR-fCos2PhiR)/sqr(fCos2PhiR.error)
//	+ sqr(sRecCos2PhiL-fCos2PhiL)/sqr(fCos2PhiL.error)
//	+ sqr(sRecMassh0-fMassh0)/sqr(fMassh0.error)
//	+ sqr(sRecMassA0-fMassA0)/sqr(fMassA0.error)
//	+ sqr(sRecMassTop-fMassTop)/sqr(fMassTop.error)
//	); 
  

  if ((!(f<0.))&&(!(f>=0.))) {

    f = 10000000.;
    cout << "detected nan!" << endl;

  }

  if (iflag == 10) {
    parameter_t tmpchisq;
    // write individual chisq contributions of each observable to global variable
    for (unsigned int i = 0; i < yyMeasuredVec.size(); i++) {
      if ((yyMeasuredVec[i].theoset ) && (!yyMeasuredVec[i].nofit )
	  && (!yyMeasuredVec[i].temp_nofit ) ) {
	if ((yyMeasuredVec[i].error > 0.) ) {
	  tmpchisq.name = yyMeasuredVec[i].name;
	  tmpchisq.value = sqr(yyMeasuredVec[i].theovalue-yyMeasuredVec[i].value)
	    * yyMeasuredCorrelationMatrix.GetInverseCovariance(i,i);
	  tmpchisq.error = yyMeasuredVec[i].alias;
	  indchisq_vec.push_back(tmpchisq);
	  //	  cout << "Ind Delta Chisq = " << tmpchisq.value << endl;
	}
      }
    }
  }

  if (iflag == 20) {
    // do not use the fitted observables, just test for NANs
    f = sRecMassChargino1 + sRecMassChargino2 + sRecMassNeutralino1 + sRecMassNeutralino2 +
      sRecMassA0 + sRecMassh0 + sRecMassGluino + sRecCos2PhiL + sRecCos2PhiR;
    if ((!(f<0.))&&(!(f>=0.))) {
      f = 10000000.;
      cout << "detected nan in return!" << endl;
    }
  }

  //cout << "f = " << f << endl;

  return;
} 


// ************************************************************
// callSPheno
//
// calculates: The spectrum as predicted by SPheno
// needs:      all
//
// P. Bechtle, P. Wienemann, 18.09.03
//
// ************************************************************

int callSPheno() 
{
  int return_value;

//  FILE * SPheno;
//  SPheno = popen("./SPheno", "r"); 
//  pclose (SPheno);

  if (!yyGeneratorPath.compare(""))
    //return_value = system("./SPheno");
    return_value = checkcall("./SPheno",20);
  else
    //return_value = system(yyGeneratorPath.c_str());
    return_value = checkcall(yyGeneratorPath.c_str(),20);
  if (return_value > 0) {
    return(return_value);
  }
  
  return 0;

}

// ************************************************************
// checkcall
//
// performs:   fork of SPheno, check of runtime of SPheno
// needs:      all
//
// P. Bechtle, P. Wienemann, 21.10.04
//
// ************************************************************
int checkcall(char programpath[1024], unsigned int runtime) 
{
  int rc;

  // preliminary: do nothing
  rc = system (programpath);
  return rc;

}

// ************************************************************
// callSuspect
//
// calculates: The spectrum as predicted by SUSPECT
// needs:      all
//
// P. Bechtle, P. Wienemann, 18.09.03
//
// ************************************************************

void callSuspect() 
{
  
// ==================================================
// This is for Tests with SUSPECT:
//  int imodel = 0; // 0: arbitrary soft terms at low scale
//  int iinput = 0; // 0: ma,mu as input; 1: v_1^2, v_2^2 as input
//  double vintb = fitterTanBeta;
//  double vinmu = fitterMu;
//  double vinma = fitterMassA0;
//  double vinm[3];
//  double vinmsl = fitterMassSelL;
//  double vinmstaur = fitterMassSelR;
//  double vinmsq = fitterMassSupL;
//  double vinmtr = fitterMassSupR;
//  double vina = fitterA;
//  double suspect_mchip[2];
//  double suspect_mchi0[4];
//  double suspect_u[2][2];
//  double suspect_vv[2][2];
//  double suspect_z[4][4];
//  double suspect_xmn[4];
//  double suspect_mgl;
//  double suspect_mh[4];
//  double vinmt = fitterMassTop;
//
//  double PhiR,PhiL;
//
//  vinm[0] = fitterM1;
//  vinm[1] = fitterM2;
//  vinm[2] = fitterAbsM3;
//
//  //cout << "fitterMassA0 " << fitterMassA0 << endl;
//
//  // cout << "calling suspect" << endl;
//  call_suspect_free_(&imodel,&iinput,&vintb,&vinmu,&vinma,
//		    vinm,&vinmsl,&vinmstaur,&vinmsq,&vinmtr,&vina,&vinmt);
//
//  // cout << "getting parameters from suspect" << endl;
//  get_suspect_(suspect_mchip,suspect_mchi0,&suspect_mgl,suspect_mh,suspect_u,suspect_vv,suspect_z,suspect_xmn);
//
//  sRecMassChargino1 = suspect_mchip[0];
//  sRecMassChargino2 = suspect_mchip[1];
//  sRecMassNeutralino1 = suspect_mchi0[0];
//  sRecMassNeutralino2 = suspect_mchi0[1];
//  sRecMassNeutralino3 = suspect_mchi0[2];
//  sRecMassNeutralino4 = suspect_mchi0[3];
//  sRecMassGluino = suspect_mgl;
//  sRecMassh0 = suspect_mh[0];
//  sRecMassH0 = suspect_mh[1];
//  sRecMassA0 = suspect_mh[2];
//  sRecMassHpm = suspect_mh[3];
//  PhiL = TMath::ACos(suspect_u[1][1]);
//  PhiR = TMath::ACos(suspect_vv[1][1]);
//  sRecCos2PhiL = TMath::Cos(2.*PhiL);
//  sRecCos2PhiR = TMath::Cos(2.*PhiR);
//  sRecMassTop = fitterMassTop;
  // cout << "PhiL,PhiR: " << PhiL << " " << PhiR << endl;

  //  cout <<  "suspect_mchip[0]" << suspect_mchip[0] << endl;
  // cout << sRecMassChargino1 << endl;
//  cout << sRecMassChargino2 << endl;
//  cout << sRecMassNeutralino1 << endl;
//  cout << sRecMassNeutralino2 << endl;
//  cout << sRecMassNeutralino3 << endl;
//  cout << sRecMassNeutralino4 << endl;
//  cout << sRecMassGluino << endl;
//  cout << sRecMassh0 << endl;
//  cout << sRecMassH0 << endl;
//  cout << sRecMassA0 << endl;
//  cout << sRecMassHpm << endl;
//    cout << sRecCos2PhiL << endl;
//    cout << sRecCos2PhiR << endl;


  
// ==================================================

}

// ************************************************************
// WriteLesHouches
//
// calculates: nothing
// needs:      x, yyFittedVec
//
// P. Bechtle, P. Wienemann, 06.02.04
//
// ************************************************************
void WriteLesHouches(double* x) 
{
  double local_mu, local_tanb;

  if (FindInFixed("Mu")) {
    local_mu = ReturnFixedValue("Mu")->value;
  }    
  else if (FindInFitted("Mu")) {
    local_mu = x[ReturnFittedPosition("Mu")];
  } 
  else if (FindInUniversality("Mu")) {
    local_mu = x[ReturnFittedPosition(ReturnUniversality("Mu")->universality)];
  }
  if (FindInFixed("TanBeta")) {
    local_tanb = ReturnFixedValue("TanBeta")->value;
  }    
  else if (FindInFitted("TanBeta")) {
    local_tanb = x[ReturnFittedPosition("TanBeta")];
  } 
  else if (FindInUniversality("TanBeta")) {
    local_tanb = x[ReturnFittedPosition(ReturnUniversality("TanBeta")->universality)];
  }
  //  cout << "writeLesHouches: local_tanb = "<< local_tanb << endl;
  //  cout << "writeLesHouches: local_mu = "<< local_mu << endl;

  // open file CrossSections.in
//  fstream CrossSectionsOutfile;
//  CrossSectionsOutfile.open ("CrossSections.in",ofstream::out);
//  //  CrossSectionsOutfile.setf(ios::scientific, ios::floatfield);
//  if (CrossSectionsOutfile.is_open()) {
//    CrossSectionsOutfile << CrossSectionProduction[count][0] << endl;
//    CrossSectionsOutfile << CrossSectionProduction[count][1] << endl;
//    CrossSectionsOutfile << CrossSectionProduction[count][2] << endl;
//    if (yyISR) {
//	CrossSectionsOutfile << ".TRUE." << endl;
//    } else {
//	CrossSectionsOutfile << ".FALSE." << endl;
//    }
//  }
//  CrossSectionsOutfile.close();
  //  system ("cat CrossSections.in");
  // open file LesHouches.in
  fstream LesHouchesOutfile;
  LesHouchesOutfile.open ("LesHouches.in",ofstream::out);
  LesHouchesOutfile.setf(ios::scientific, ios::floatfield);

//  1  -1                  # error level
//  2   0                  # if 1, then SPA conventions are used
// 11   1                  # calculate branching ratios
// 12   1.00000000E-04     # write only branching ratios larger than this value
// 21   1                  # calculate cross section
// 22   5.00000000E+02     # cms energy in GeV
// 23   0.00000000E+00     # polarisation of incoming e- beam
// 24   0.00000000E+00     # polarisation of incoming e+ beam
// 25   0                  # if 0 no ISR is calculated, if 1 ISR is caculated
// 26   1.00000000E-04     # write only cross sections larger than this value [fb]
// 31  -1.00000000E+00     # m_GUT, if < 0 than it determined via g_1=g_2
// 32   0                  # require strict unification g_1=g_2=g_3 if '1' is set
// 33  -1.00000000E+00     # Q_EWSB, if < 0 than  Q_EWSB=sqrt(m_~t1 m_~t2)
// 41   2.49520000E+00     # width of the Z-boson
// 42   2.11800000E+00     # width of the W-boson
// 51   5.10998900E-04     # electron mass
// 52   1.05658357E-01     # muon mass
// 61   2.00000000E+00     # scale where quark masses of first 2 gen. are defined
// 62   3.00000000E-03     # m_u(Q)
// 63   1.20000000E+00     # m_c(Q)
// 64   7.00000000E-03     # m_d(Q)
// 65   1.20000000E-01     # m_s(Q)

  if (LesHouchesOutfile.is_open()) {

    LesHouchesOutfile << "BLOCK MODSEL" << endl;
    LesHouchesOutfile << "    1 0 # general MSSM" << endl;
    LesHouchesOutfile << "BLOCK SPhenoInput" << endl;
    LesHouchesOutfile << "    1  0                  # error level" << endl;
    LesHouchesOutfile << "    2  0                  # if 1, then SPA conventions are used" << endl;
    LesHouchesOutfile << "   11  1                  # calculate branching ratios" << endl;
    LesHouchesOutfile << "   12  1.00000000E-04     # write only branching ratios larger than this value" << endl;
    LesHouchesOutfile << "   21  1                  # calculate cross section" << endl;
    for (unsigned int j = 0; j < CrossSectionProduction.size(); j++) {
      LesHouchesOutfile << "   22  " << CrossSectionProduction[j][0] << "     # cms energy in GeV" << endl;
      LesHouchesOutfile << "   23  " << CrossSectionProduction[j][1] << "     # polarisation of incoming e- beam" << endl;
      LesHouchesOutfile << "   24  " << CrossSectionProduction[j][2] << "     # polarisation of incoming e+ beam" << endl;
      if (!yyISR) {
	LesHouchesOutfile << "   25  0                  # no ISR is calculated" << endl;
      } else {
	LesHouchesOutfile << "   25  1                  # ISR is calculated" << endl;
      }
    }
    LesHouchesOutfile << "   26  1.00000000E-05     # write only cross sections larger than this value [fb]" << endl;
    LesHouchesOutfile << "   31 -1.00000000E+00     # m_GUT, if < 0 than it determined via g_1=g_2" << endl;
    LesHouchesOutfile << "   32  0                  # require strict unification g_1=g_2=g_3 if '1' is set " << endl;
    LesHouchesOutfile << "   33  1000.              #  Q_EWSB, if < 0 than  Q_EWSB=sqrt(m_~t1 m_~t2) " << endl;
    if (FindInFixed("massCharm")) {
      LesHouchesOutfile << "   63 "<<ReturnFixedValue("massCharm")->value<<" # mcharm (fixed)"<<endl;
    }
    else if (FindInFitted("massCharm")) {
      LesHouchesOutfile << "   63 "<<x[ReturnFittedPosition("massCharm")]<<" # mcharm"<<endl;
      cout << "Fitting mCharm " << x[ReturnFittedPosition("massCharm")] << endl;
    } 
    else if (FindInUniversality("massCharm")) {
      LesHouchesOutfile << "   63 "<<x[ReturnFittedPosition(ReturnUniversality("massCharm")->universality)]<<" # massCharm"<<endl;
      cout << "fitting " << ReturnUniversality("massCharm")->universality << " instead of massCharm" << endl;
    }
    else {
      LesHouchesOutfile << "   63 "<<ReturnMeasuredValue("massCharm")->value<<" # mcharm (fixed)"<<endl;
    }
    
    //--------------------------------------------------------------------
    LesHouchesOutfile << "BLOCK SMINPUTS" << endl;
    if (FindInFixed("alphaem")) {
      LesHouchesOutfile << "    1 "<<ReturnFixedValue("alphaem")->value<<" # 1/alpha_em (fixed)"<<endl;
    }
    else if (FindInFitted("alphaem")) {
      LesHouchesOutfile << "    1 "<<x[ReturnFittedPosition("alphaem")]<<" # 1/alpha_em(M_Z)"<<endl;
      cout << "Fitting alphaem " << x[ReturnFittedPosition("alphaem")] << endl;
    } 
    else if (FindInUniversality("alphaem")) {
      LesHouchesOutfile << "    1 "<<x[ReturnFittedPosition(ReturnUniversality("alphaem")->universality)]<<" # 1/alpha_em(M_Z)"<<endl;
      cout << "fiting " << ReturnUniversality("alphaem")->universality << " instead of alphaem" << endl;
    }
    else {
      LesHouchesOutfile << "    1 "<<ReturnMeasuredValue("alphaem")->value<<" # 1/alpha_em(M_Z) (fixed)"<<endl;
    }


    if (FindInFixed("G_F")) {
      LesHouchesOutfile << "    2 "<<ReturnFixedValue("G_F")->value<<" # G_F (fixed)"<<endl;
    }
    else if (FindInFitted("G_F")) {
      LesHouchesOutfile << "    2 "<<x[ReturnFittedPosition("G_F")]<<" # G_F"<<endl;
      cout << "Fitting G_F " << x[ReturnFittedPosition("G_F")] << endl;
    } 
    else if (FindInUniversality("G_F")) {
      LesHouchesOutfile << "    2 "<<x[ReturnFittedPosition(ReturnUniversality("G_F")->universality)]<<" # G_F"<<endl;
      cout << "fitting " << ReturnUniversality("G_F")->universality << " instead of G_F" << endl;
    }
    else {
      LesHouchesOutfile << "    2 "<<ReturnMeasuredValue("G_F")->value<<" # G_F (fixed)"<<endl;
    }

    if (FindInFixed("alphas")) {
      LesHouchesOutfile << "    3 "<<ReturnFixedValue("alphas")->value<<" # alpha_s (fixed)"<<endl;
    }
    else if (FindInFitted("alphas")) {
      LesHouchesOutfile << "    3 "<<x[ReturnFittedPosition("alphas")]<<" # alpha_s"<<endl;
      cout << "Fitting alphas " << x[ReturnFittedPosition("alphas")] << endl;
    } 
    else if (FindInUniversality("alphas")) {
      LesHouchesOutfile << "    3 "<<x[ReturnFittedPosition(ReturnUniversality("alphas")->universality)]<<" # alpha_ss"<<endl;
      cout << "fitting " << ReturnUniversality("alphas")->universality << " instead of alphas" << endl;
    }
    else {
      LesHouchesOutfile << "    3 "<<ReturnMeasuredValue("alphas")->value<<" # alpha_s (fixed)"<<endl;
    }


    if (FindInFixed("massZ")) {
      LesHouchesOutfile << "    4 "<<ReturnFixedValue("massZ")->value<<" # mZ (fixed)"<<endl;
    }
    else if (FindInFitted("massZ")) {
      LesHouchesOutfile << "    4 "<<x[ReturnFittedPosition("massZ")]<<" # mZ"<<endl;
      cout << "Fitting massZ " << x[ReturnFittedPosition("massZ")] << endl;
    } 
    else if (FindInUniversality("massZ")) {
      LesHouchesOutfile << "    4 "<<x[ReturnFittedPosition(ReturnUniversality("massZ")->universality)]<<" # massZ"<<endl;
      cout << "fitting " << ReturnUniversality("massZ")->universality << " instead of massZ" << endl;
    }
    else {
      LesHouchesOutfile << "    4 "<<ReturnMeasuredValue("massZ")->value<<" # mZ (fixed)"<<endl;
    }

    if (FindInFixed("massBottom")) {
      LesHouchesOutfile << "    5 "<<ReturnFixedValue("massBottom")->value<<" # mb(mb) (fixed)"<<endl;
    }
    else if (FindInFitted("massBottom")) {
      LesHouchesOutfile << "    5 "<<x[ReturnFittedPosition("massBottom")]<<" # mb(mb)"<<endl;
      cout << "Fitting mBottom " << x[ReturnFittedPosition("massBottom")] << endl;
    } 
    else if (FindInUniversality("massBottom")) {
      LesHouchesOutfile << "    5 "<<x[ReturnFittedPosition(ReturnUniversality("massBottom")->universality)]<<" # massBottom"<<endl;
      cout << "fitting " << ReturnUniversality("massBottom")->universality << " instead of massBottom" << endl;
    }
    else {
      LesHouchesOutfile << "    5 "<<ReturnMeasuredValue("massBottom")->value<<" # mb(mb) (fixed)"<<endl;
    }


    if (FindInFixed("massTop")) {
      LesHouchesOutfile << "    6 "<<ReturnFixedValue("massTop")->value<<" # mtop (fixed)"<<endl;
    }
    else if (FindInFitted("massTop")) {
      LesHouchesOutfile << "    6 "<<x[ReturnFittedPosition("massTop")]<<" # mtop"<<endl;
      cout << "Fitting mTop " << x[ReturnFittedPosition("massTop")] << endl;
    } 
    else if (FindInUniversality("massTop")) {
      LesHouchesOutfile << "    6 "<<x[ReturnFittedPosition(ReturnUniversality("massTop")->universality)]<<" # massTop"<<endl;
      cout << "fitting " << ReturnUniversality("massTop")->universality << " instead of massTop" << endl;
    }
    else {
      LesHouchesOutfile << "    6 "<<ReturnMeasuredValue("massTop")->value<<" # mtop (fixed)"<<endl;
    }

    if (FindInFixed("massTau")) {
      LesHouchesOutfile << "    7 "<<ReturnFixedValue("massTau")->value<<" # mtau (fixed)"<<endl;
    }
    else if (FindInFitted("massTau")) {
      LesHouchesOutfile << "    7 "<<x[ReturnFittedPosition("massTau")]<<" # mtau"<<endl;
    } 
    else if (FindInUniversality("massTau")) {
      LesHouchesOutfile << "    7 "<<x[ReturnFittedPosition(ReturnUniversality("massTau")->universality)]<<" # massTau"<<endl;
      cout << "fitting " << ReturnUniversality("massTau")->universality << " instead of massTau" << endl;
    }
    else {
      LesHouchesOutfile << "    7 "<<ReturnMeasuredValue("massTau")->value<<" # mtau (fixed)"<<endl;
    }


//    if (FindInFixed("massCharm")) {
//      LesHouchesOutfile << "    8 "<<ReturnFixedValue("massCharm")->value<<" # mcharm (fixed)"<<endl;
//    }
//    else if (FindInFitted("massCharm")) {
//      LesHouchesOutfile << "    8 "<<x[ReturnFittedPosition("massCharm")]<<" # mcharm"<<endl;
//      cout << "Fitting mCharm " << x[ReturnFittedPosition("massCharm")] << endl;
//    } 
//    else if (FindInUniversality("massCharm")) {
//      LesHouchesOutfile << "    8 "<<x[ReturnFittedPosition(ReturnUniversality("massCharm")->universality)]<<" # massCharm"<<endl;
//      cout << "fitting " << ReturnUniversality("massCharm")->universality << " instead of massCharm" << endl;
//    }
//    else {
//      LesHouchesOutfile << "    8 "<<ReturnMeasuredValue("massCharm")->value<<" # mcharm (fixed)"<<endl;
//    }


    // MINPAR
    LesHouchesOutfile << "BLOCK MINPAR" << endl;
    if (FindInFixed("TanBeta")) {
      LesHouchesOutfile << "    3 "<< ReturnFixedValue("TanBeta")->value <<" # tanb (fixed)"<< endl;
    }
    else if (FindInFitted("TanBeta")) {
      cout << "Fitting tanb " << x[ReturnFittedPosition("TanBeta")] << endl;
      LesHouchesOutfile << "    3 "<< x[ReturnFittedPosition("TanBeta")]<<" # tanb"<< endl;
    } 
    else if (FindInUniversality("TanBeta")) {
      LesHouchesOutfile << "    3 "<<x[ReturnFittedPosition(ReturnUniversality("TanBeta")->universality)]<<" # TanBeta"<<endl;
      cout << "fitting " << ReturnUniversality("TanBeta")->universality << " instead of TanBeta" << endl;
    }
    else {
      cerr << "Parameter TanBeta not declared" << endl;
      exit (EXIT_FAILURE);
    }
    // EXTPAR
    LesHouchesOutfile << "BLOCK EXTPAR" << endl;
    //    LesHouchesOutfile << "BLOCK MSOFT" << endl;
    LesHouchesOutfile << "    0 1000. # Input scale for mSUGRA" << endl;
    if (FindInFixed("M1")) {
      LesHouchesOutfile << "    1 "<< ReturnFixedValue("M1")->value <<" # M1 (fixed)"<< endl;
    }    
    else if (FindInFitted("M1")) {
      LesHouchesOutfile << "    1 "<< x[ReturnFittedPosition("M1")]<<" # M1"<< endl;
      cout << "Fitting M1 " << x[ReturnFittedPosition("M1")] << endl;
    } 
    else if (FindInUniversality("M1")) {
      LesHouchesOutfile << "    1 "<<x[ReturnFittedPosition(ReturnUniversality("M1")->universality)]<<" # M1"<<endl;
      cout << "fitting " << ReturnUniversality("M1")->universality << " instead of M1" << endl;
    }
    else {
      cerr << "Parameter M1 not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("M2")) {
      LesHouchesOutfile << "    2 "<< ReturnFixedValue("M2")->value <<" # M2 (fixed)"<< endl;
    }    
    else if (FindInFitted("M2")) {
      LesHouchesOutfile << "    2 "<< x[ReturnFittedPosition("M2")]<<" # M2"<< endl;
      cout << "Fitting M2 " << x[ReturnFittedPosition("M2")] << endl;
    }
    else if (FindInUniversality("M2")) {
      LesHouchesOutfile << "    2 "<<x[ReturnFittedPosition(ReturnUniversality("M2")->universality)]<<" # M2"<<endl;
      cout << "fitting " << ReturnUniversality("M2")->universality << " instead of M2" << endl;
    }
    else {
      cerr << "Parameter M2 not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("M3")) {
      LesHouchesOutfile << "    3 "<< ReturnFixedValue("M3")->value <<" # M3 (fixed)"<< endl;
    }    
    else if (FindInFitted("M3")) {
      LesHouchesOutfile << "    3 "<< x[ReturnFittedPosition("M3")]<<" # M3"<< endl;
      cout << "Fitting M3 " << x[ReturnFittedPosition("M3")] << endl;
    }
    else if (FindInUniversality("M3")) {
      LesHouchesOutfile << "    3 "<<x[ReturnFittedPosition(ReturnUniversality("M3")->universality)]<<" # M3"<<endl;
      cout << "fitting " << ReturnUniversality("M3")->universality << " instead of M3" << endl;
    }
    else {
      cerr << "Parameter M3 not declared" << endl;
      exit (EXIT_FAILURE);
    }


    if (FindInFixed("Xtop")) {
      LesHouchesOutfile << "   11 "<< ReturnFixedValue("Xtop")->value+local_mu/local_tanb <<" # Atop (fixed)"<< endl;
    }    
    else if (FindInFitted("Xtop")) {
      LesHouchesOutfile << "   11 "<< x[ReturnFittedPosition("Xtop")]+local_mu/local_tanb<<" # Atop"<< endl;
      cout << "Fitting Xtop " << x[ReturnFittedPosition("Xtop")] << endl;
    } 
    else if (FindInUniversality("Xtop")) {
      LesHouchesOutfile << "   11 "<<x[ReturnFittedPosition(ReturnUniversality("Xtop")->universality)]+local_mu/local_tanb<<" # Atop"<<endl;
      cout << "fitting " << ReturnUniversality("Xtop")->universality << " instead of Xtop" << endl;
    }
    else {
      cerr << "Parameter Xtop not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("Xbottom")) {
      LesHouchesOutfile << "   12 "<< ReturnFixedValue("Xbottom")->value+local_mu*local_tanb <<" # Abottom (fixed)"<< endl;
    }    
    else if (FindInFitted("Xbottom")) {
      LesHouchesOutfile << "   12 "<< x[ReturnFittedPosition("Xbottom")]+local_mu*local_tanb<<" # Abottom"<< endl;
      cout << "Fitting Xbottom " << x[ReturnFittedPosition("Xbottom")] << endl;
    }  
    else if (FindInUniversality("Xbottom")) {
      LesHouchesOutfile << "   12 "<<x[ReturnFittedPosition(ReturnUniversality("Xbottom")->universality)]+local_mu*local_tanb<<" # Abottom"<<endl;
      cout << "fitting " << ReturnUniversality("Xbottom")->universality << " instead of Xbottom" << endl;
    }
    else {
      cerr << "Parameter Xbottom not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("Xtau")) {
      LesHouchesOutfile << "   13 "<< ReturnFixedValue("Xtau")->value+local_mu*local_tanb <<" # Atau (fixed)"<< endl;
    }    
    else if (FindInFitted("Xtau")) {
      LesHouchesOutfile << "   13 "<< x[ReturnFittedPosition("Xtau")]+local_mu*local_tanb<<" # Atau"<< endl;
      cout << "Fitting Xtau " << x[ReturnFittedPosition("Xtau")] << endl;
    }  
    else if (FindInUniversality("Xtau")) {
      LesHouchesOutfile << "   13 "<<x[ReturnFittedPosition(ReturnUniversality("Xtau")->universality)]+local_mu*local_tanb<<" # Atau"<<endl;
      cout << "fitting " << ReturnUniversality("Xtau")->universality << " instead of Xtau" << endl;
    }
    else {
      cerr << "Parameter Xtau not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("Mu")) {
      LesHouchesOutfile << "   23 "<< ReturnFixedValue("Mu")->value <<" # mu (fixed)"<< endl;
    }    
    else if (FindInFitted("Mu")) {
      LesHouchesOutfile << "   23 "<< x[ReturnFittedPosition("Mu")]<<" # mu"<< endl;
      cout << "Fitting Mu " << x[ReturnFittedPosition("Mu")] << endl;
    } 
    else if (FindInUniversality("Mu")) {
      LesHouchesOutfile << "   23 "<<x[ReturnFittedPosition(ReturnUniversality("Mu")->universality)]<<" # Mu"<<endl;
      cout << "fitting " << ReturnUniversality("Mu")->universality << " instead of Mu" << endl;
    }
    else {
      cerr << "Parameter Mu not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("massA0")) {
      LesHouchesOutfile << "   24 "<< sqr(ReturnFixedValue("massA0")->value) <<" # mA (fixed)"<< endl;
    }    
    else if (FindInFitted("massA0")) {
      LesHouchesOutfile << "   24 "<< sqr(x[ReturnFittedPosition("massA0")])<<" # mA"<< endl;
      cout << "Fitting mA " << x[ReturnFittedPosition("massA0")] << endl;
    } 
    else if (FindInUniversality("massA0")) {
      LesHouchesOutfile << "   24 "<<sqr(x[ReturnFittedPosition(ReturnUniversality("massA0")->universality)])<<" # massA0"<<endl;
      cout << "fitting " << ReturnUniversality("massA0")->universality << " instead of massA0" << endl;
    }
    else {
      cerr << "Parameter massA0 not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("MSelectronL")) {
      LesHouchesOutfile << "   31 "<< ReturnFixedValue("MSelectronL")->value <<" # MSelL (fixed)"<< endl;
    }    
    else if (FindInFitted("MSelectronL")) {
      LesHouchesOutfile << "   31 "<< x[ReturnFittedPosition("MSelectronL")]<<" # MSelL"<< endl;
      cout << "Fitting MSelectronL " << x[ReturnFittedPosition("MSelectronL")] << endl;
    } 
    else if (FindInUniversality("MSelectronL")) {
      LesHouchesOutfile << "   31 "<<x[ReturnFittedPosition(ReturnUniversality("MSelectronL")->universality)]<<" # MSelectronL"<<endl;
      cout << "fitting " << ReturnUniversality("MSelectronL")->universality << " instead of MSelectronL" << endl;
    }
    else {
      cerr << "Parameter MSelectronL not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("MSmuL")) {
      LesHouchesOutfile << "   32 "<< ReturnFixedValue("MSmuL")->value <<" # MSmuL (fixed)"<< endl;
    }    
    else if (FindInFitted("MSmuL")) {
      LesHouchesOutfile << "   32 "<< x[ReturnFittedPosition("MSmuL")]<<" # MSmuL"<< endl;
      cout << "Fitting MSmuL " << x[ReturnFittedPosition("MSmuL")] << endl;
    } 
    else if (FindInUniversality("MSmuL")) {
      LesHouchesOutfile << "   32 "<<x[ReturnFittedPosition(ReturnUniversality("MSmuL")->universality)]<<" # MSmuL"<<endl;
      cout << "fitting " << ReturnUniversality("MSmuL")->universality << " instead of MSmuL" << endl;
    }
    else {
      cerr << "Parameter MSmuL not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("MStauL")) {
      LesHouchesOutfile << "   33 "<< ReturnFixedValue("MStauL")->value <<" # MStauL (fixed)"<< endl;
    }    
    else if (FindInFitted("MStauL")) {
      LesHouchesOutfile << "   33 "<< x[ReturnFittedPosition("MStauL")]<<" # MStauL"<< endl;
      cout << "Fitting MStauL " << x[ReturnFittedPosition("MStauL")] << endl;
    } 
    else if (FindInUniversality("MStauL")) {
      LesHouchesOutfile << "   33 "<<x[ReturnFittedPosition(ReturnUniversality("MStauL")->universality)]<<" # MStauL"<<endl;
      cout << "fitting " << ReturnUniversality("MStauL")->universality << " instead of MStauL" << endl;
    }
    else {
      cerr << "Parameter MStauL not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("MSelectronR")) {
      LesHouchesOutfile << "   34 "<< ReturnFixedValue("MSelectronR")->value <<" # MSelR (fixed)"<< endl;
    }    
    else if (FindInFitted("MSelectronR")) {
      LesHouchesOutfile << "   34 "<< x[ReturnFittedPosition("MSelectronR")]<<" # MSelR"<< endl;
      cout << "Fitting MSelectronR " << x[ReturnFittedPosition("MSelectronR")] << endl;
    } 
    else if (FindInUniversality("MSelectronR")) {
      LesHouchesOutfile << "   34 "<<x[ReturnFittedPosition(ReturnUniversality("MSelectronR")->universality)]<<" # MSelectronR"<<endl;
      cout << "fitting " << ReturnUniversality("MSelectronR")->universality << " instead of MSelectronR" << endl;
    }
    else {
      cerr << "Parameter MSelectronR not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("MSmuR")) {
      LesHouchesOutfile << "   35 "<< ReturnFixedValue("MSmuR")->value <<" # MSmuR (fixed)"<< endl;
    }    
    else if (FindInFitted("MSmuR")) {
      LesHouchesOutfile << "   35 "<< x[ReturnFittedPosition("MSmuR")]<<" # MSmuR"<< endl;
      cout << "Fitting MSmuR " << x[ReturnFittedPosition("MSmuR")] << endl;
    } 
    else if (FindInUniversality("MSmuR")) {
      LesHouchesOutfile << "   35 "<<x[ReturnFittedPosition(ReturnUniversality("MSmuR")->universality)]<<" # MSmuR"<<endl;
      cout << "fitting " << ReturnUniversality("MSmuR")->universality << " instead of MSmuR" << endl;
    }
    else {
      cerr << "Parameter MSmuR not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("MStauR")) {
      LesHouchesOutfile << "   36 "<< ReturnFixedValue("MStauR")->value <<" # MStauR (fixed)"<< endl;
    }    
    else if (FindInFitted("MStauR")) {
      LesHouchesOutfile << "   36 "<< x[ReturnFittedPosition("MStauR")]<<" # MStauR"<< endl;
      cout << "Fitting MStauR " << x[ReturnFittedPosition("MStauR")] << endl;
    } 
    else if (FindInUniversality("MStauR")) {
      LesHouchesOutfile << "   36 "<<x[ReturnFittedPosition(ReturnUniversality("MStauR")->universality)]<<" # MStauR"<<endl;
      cout << "fitting " << ReturnUniversality("MStauR")->universality << " instead of MStauR" << endl;
    }
    else {
      cerr << "Parameter MStauR not declared" << endl;
      exit (EXIT_FAILURE);
    }

    if (FindInFixed("MSupL")) {
      LesHouchesOutfile << "   41 "<< ReturnFixedValue("MSupL")->value <<" # MSupL (fixed)"<< endl;
    }    
    else if (FindInFitted("MSupL")) {
      LesHouchesOutfile << "   41 "<< x[ReturnFittedPosition("MSupL")]<<" # MSupL"<< endl;
      cout << "Fitting MSupL " << x[ReturnFittedPosition("MSupL")] << endl;
    } 
    else if (FindInUniversality("MSupL")) {
      LesHouchesOutfile << "   41 "<<x[ReturnFittedPosition(ReturnUniversality("MSupL")->universality)]<<" # MSupL"<<endl;
      cout << "fitting " << ReturnUniversality("MSupL")->universality << " instead of MSupL" << endl;
    }

    if (FindInFixed("MScharmL")) {
      LesHouchesOutfile << "   42 "<< ReturnFixedValue("MScharmL")->value <<" # MScharmL (fixed)"<< endl;
    }    
    else if (FindInFitted("MScharmL")) {
      LesHouchesOutfile << "   42 "<< x[ReturnFittedPosition("MScharmL")]<<" # MScharmL"<< endl;
      cout << "Fitting MScharmL " << x[ReturnFittedPosition("MScharmL")] << endl;
    } 
    else if (FindInUniversality("MScharmL")) {
      LesHouchesOutfile << "   42 "<<x[ReturnFittedPosition(ReturnUniversality("MScharmL")->universality)]<<" # MScharmL"<<endl;
      cout << "fitting " << ReturnUniversality("MScharmL")->universality << " instead of MScharmL" << endl;
    }

    if (FindInFixed("MStopL")) {
      LesHouchesOutfile << "   43 "<< ReturnFixedValue("MStopL")->value <<" # MStopL (fixed)"<< endl;
    }    
    else if (FindInFitted("MStopL")) {
      LesHouchesOutfile << "   43 "<< x[ReturnFittedPosition("MStopL")]<<" # MStopL"<< endl;
      cout << "Fitting MStopL " << x[ReturnFittedPosition("MStopL")] << endl;
    } 
    else if (FindInUniversality("MStopL")) {
      LesHouchesOutfile << "   43 "<<x[ReturnFittedPosition(ReturnUniversality("MStopL")->universality)]<<" # MStopL"<<endl;
      cout << "fitting " << ReturnUniversality("MStopL")->universality << " instead of MStopL" << endl;
    }


    if (FindInFixed("MSupR")) {
      LesHouchesOutfile << "   44 "<< ReturnFixedValue("MSupR")->value <<" # MSupR (fixed)"<< endl;
    }    
    else if (FindInFitted("MSupR")) {
      LesHouchesOutfile << "   44 "<< x[ReturnFittedPosition("MSupR")]<<" # MSupR"<< endl;
      cout << "Fitting MSupR " << x[ReturnFittedPosition("MSupR")] << endl;
    } 
    else if (FindInUniversality("MSupR")) {
      LesHouchesOutfile << "   44 "<<x[ReturnFittedPosition(ReturnUniversality("MSupR")->universality)]<<" # MSupR"<<endl;
      cout << "fitting " << ReturnUniversality("MSupR")->universality << " instead of MSupR" << endl;
    }

    if (FindInFixed("MScharmR")) {
      LesHouchesOutfile << "   45 "<< ReturnFixedValue("MScharmR")->value <<" # MScharmR (fixed)"<< endl;
    }    
    else if (FindInFitted("MScharmR")) {
      LesHouchesOutfile << "   45 "<< x[ReturnFittedPosition("MScharmR")]<<" # MScharmR"<< endl;
      cout << "Fitting MScharmR " << x[ReturnFittedPosition("MScharmR")] << endl;
    } 
    else if (FindInUniversality("MScharmR")) {
      LesHouchesOutfile << "   45 "<<x[ReturnFittedPosition(ReturnUniversality("MScharmR")->universality)]<<" # MScharmR"<<endl;
      cout << "fitting " << ReturnUniversality("MScharmR")->universality << " instead of MScharmR" << endl;
    }

    if (FindInFixed("MStopR")) {
      LesHouchesOutfile << "   46 "<< ReturnFixedValue("MStopR")->value <<" # MStopR (fixed)"<< endl;
    }    
    else if (FindInFitted("MStopR")) {
      LesHouchesOutfile << "   46 "<< x[ReturnFittedPosition("MStopR")]<<" # MStopR"<< endl;
      cout << "Fitting MStopR " << x[ReturnFittedPosition("MStopR")] << endl;
    } 
    else if (FindInUniversality("MStopR")) {
      LesHouchesOutfile << "   46 "<<x[ReturnFittedPosition(ReturnUniversality("MStopR")->universality)]<<" # MStopR"<<endl;
      cout << "fitting " << ReturnUniversality("MStopR")->universality << " instead of MStopR" << endl;
    }


    if (FindInFixed("MSdownR")) {
      LesHouchesOutfile << "   47 "<< ReturnFixedValue("MSdownR")->value <<" # MSdownR (fixed)"<< endl;
    }    
    else if (FindInFitted("MSdownR")) {
      LesHouchesOutfile << "   47 "<< x[ReturnFittedPosition("MSdownR")]<<" # MSdownR"<< endl;
      cout << "Fitting MSdownR " << x[ReturnFittedPosition("MSdownR")] << endl;
    } 
    else if (FindInUniversality("MSdownR")) {
      LesHouchesOutfile << "   47 "<<x[ReturnFittedPosition(ReturnUniversality("MSdownR")->universality)]<<" # MSdownR"<<endl;
      cout << "fitting " << ReturnUniversality("MSdownR")->universality << " instead of MSdownR" << endl;
    }

    if (FindInFixed("MSstrangeR")) {
      LesHouchesOutfile << "   48 "<< ReturnFixedValue("MSstrangeR")->value <<" # MSstrangeR (fixed)"<< endl;
    }    
    else if (FindInFitted("MSstrangeR")) {
      LesHouchesOutfile << "   48 "<< x[ReturnFittedPosition("MSstrangeR")]<<" # MSstrangeR"<< endl;
      cout << "Fitting MSstrangeR " << x[ReturnFittedPosition("MSstrangeR")] << endl;
    } 
    else if (FindInUniversality("MSstrangeR")) {
      LesHouchesOutfile << "   48 "<<x[ReturnFittedPosition(ReturnUniversality("MSstrangeR")->universality)]<<" # MSstrangeR"<<endl;
      cout << "fitting " << ReturnUniversality("MSstrangeR")->universality << " instead of MSstrangeR" << endl;
    }

    if (FindInFixed("MSbottomR")) {
      LesHouchesOutfile << "   49 "<< ReturnFixedValue("MSbottomR")->value <<" # MSbottomR (fixed)"<< endl;
    }    
    else if (FindInFitted("MSbottomR")) {
      LesHouchesOutfile << "   49 "<< x[ReturnFittedPosition("MSbottomR")]<<" # MSbottomR"<< endl;
      cout << "Fitting MSbottomR " << x[ReturnFittedPosition("MSbottomR")] << endl;
    } 
    else if (FindInUniversality("MSbottomR")) {
      LesHouchesOutfile << "   49 "<<x[ReturnFittedPosition(ReturnUniversality("MSbottomR")->universality)]<<" # MSbottomR"<<endl;
      cout << "fitting " << ReturnUniversality("MSbottomR")->universality << " instead of MSbottomR" << endl;
    }





    // mu L R up charm top


  } else {
    cerr << "error opening LesHouches.in" << endl;
    exit(EXIT_FAILURE);
  }

  //close file
  LesHouchesOutfile.close();
  return;
  
}

MeasuredValue* ReturnMeasuredValue (string name)
{
  for (unsigned int i = 0; i < yyMeasuredVec.size(); i++) {
    if (!yyMeasuredVec[i].name.compare(name)) {
      return &(yyMeasuredVec[i]);
    }
  }
  
  cerr<<"No measured value with name "<<name<<" found"<<endl;
  exit(EXIT_FAILURE);

}

MeasuredValue* ReturnFittedValue (string name)
{
  for (unsigned int i = 0; i < yyFittedVec.size(); i++) {
    if (!yyFittedVec[i].name.compare(name)) {
      return &(yyFittedVec[i]);
    }
  }

  cerr<<"No fitted value with name "<<name<<" found"<<endl;  
  exit(EXIT_FAILURE);

}

MeasuredValue* ReturnUniversality (string name)
{
  for (unsigned int i = 0; i < yyUniversalityVec.size(); i++) {
    if (!yyUniversalityVec[i].name.compare(name)) {
      return &(yyUniversalityVec[i]);
    }
  }

  cerr<<"No UniversalityVec with name "<<name<<" found"<<endl;  
  exit(EXIT_FAILURE);

}

MeasuredValue* ReturnFixedValue (string name)
{
  for (unsigned int i = 0; i < yyFixedVec.size(); i++) {
    if (!yyFixedVec[i].name.compare(name)) {
      return &(yyFixedVec[i]);
    }
  }

  cerr<<"No fixed value with name "<<name<<" found"<<endl;  
  exit(EXIT_FAILURE);

}
int ReturnFittedPosition (string name)
{
  for (unsigned int i = 0; i < yyFittedVec.size(); i++) {
    if (!yyFittedVec[i].name.compare(name)) {
      return i;
    }
  }
  
  cerr<<"No fitted position with name "<<name<<" found"<<endl;
  exit(EXIT_FAILURE);

}

bool FindInFitted (string name)
{
  for (unsigned int i = 0; i < yyFittedVec.size(); i++) {
    if (!yyFittedVec[i].name.compare(name)) {
      //      cout << name << " found in fitted on pos " << i << endl;
      return true;
    }
  }
  //  cout << name << " not found in fitted"<<endl;
  return false;
  
}
int FindInFittedPar (string name)
{
  for (unsigned int i = 0; i < yyFittedPar.size(); i++) {
    if (!yyFittedPar[i].name.compare(name)) {
      //      cout << name << " found in fitted on pos " << i << endl;
      return i;
    }
  }
  //  cout << name << " not found in fitted"<<endl;
  return -1;
  
}
bool FindInFixed (string name)
{
  for (unsigned int i = 0; i < yyFixedVec.size(); i++) {
    if (!yyFixedVec[i].name.compare(name)) {
      //      cout << name << " found in fitted on pos " << i << endl;
      return true;
    }
  }
  //  cout << name << " not found in fitted"<<endl;
  return false;
  
}
bool FindInUniversality (string name)
{
  for (unsigned int i = 0; i < yyUniversalityVec.size(); i++) {
    if (!yyUniversalityVec[i].name.compare(name)) {
      //      cout << name << " found in fitted on pos " << i << endl;
      return true;
    }
  }
  //  cout << name << " not found in fitted"<<endl;
  return false;
  
}

void Fittino::FillFixedParameters()
{
  MeasuredValue tmpValue;

  for (unsigned int i = 0; i < yyFixedPar.size(); i++) {
    for (unsigned int j = 0; j < yyFittedPar.size(); j++) {
      if (!yyFittedPar[j].name.compare(yyFixedPar[i].name)) {
	cerr << "Fixed Parameter " <<  yyFixedPar[i].name << " is also in Fitted Par" << endl;
	exit (EXIT_FAILURE);
      }
    }

    for (unsigned int j = 0; j < yyUniversalityVec.size(); j++) {
      if (!yyUniversalityVec[j].name.compare(yyFixedPar[i].name)) {
	cerr << "Fixed Parameter " <<  yyFixedPar[i].name << " is also in universality vec" << endl;
	exit (EXIT_FAILURE);
      }
    }

    tmpValue.name = yyFixedPar[i].name;
    tmpValue.value = yyFixedPar[i].value;
    tmpValue.error = -1;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    yyFixedVec.push_back(tmpValue);

    cout << "fixed: " << yyFixedPar[i].name << " at " << yyFixedPar[i].value << endl;
  }


  if (!FindInFitted("Xtop") && !FindInFixed("Xtop") && !FindInUniversality("Xtop")) {
    tmpValue.name = "Xtop";
    tmpValue.value = 500.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("Xbottom") && !FindInFixed("Xbottom") && !FindInUniversality("Xbottom")) {  
    tmpValue.name = "Xbottom";
    tmpValue.value = 400.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("Xtau") && !FindInFixed("Xtau") && !FindInUniversality("Xtau")) {  
    tmpValue.name = "Xtau";
    tmpValue.value = 300.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("massA0") && !FindInFixed("massA0") && !FindInUniversality("massA0")) {  
    tmpValue.name = "massA0";
    tmpValue.value = 436.2;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = ID_A;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("Mu") && !FindInFixed("Mu") && !FindInUniversality("Mu")) {  
    tmpValue.name = "Mu";
    tmpValue.value = 250.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MSelectronL") && !FindInFixed("MSelectronL") && !FindInUniversality("MSelectronL")) {  
    tmpValue.name = "MSelectronL";
    tmpValue.value = 550.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MSmuL") && !FindInFixed("MSmuL") && !FindInUniversality("MSmuL")) {  
    tmpValue.name = "MSmuL";
    tmpValue.value = 450.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MStauL") && !FindInFixed("MStauL") && !FindInUniversality("MStauL")) {  
    tmpValue.name = "MStauL";
    tmpValue.value = 400.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MSelectronR") && !FindInFixed("MSelectronR") && !FindInUniversality("MSelectronR")) {  
    tmpValue.name = "MSelectronR";
    tmpValue.value = 560.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MSmuR") && !FindInFixed("MSmuR") && !FindInUniversality("MSmuR")) {  
    tmpValue.name = "MSmuR";
    tmpValue.value = 460.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MStauR") && !FindInFixed("MStauR") && !FindInUniversality("MStauR")) {  
    tmpValue.name = "MStauR";
    tmpValue.value = 410.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }

  if (!FindInFitted("MSupL") && !FindInFixed("MSupL") && !FindInUniversality("MSupL")) {  
    tmpValue.name = "MSupL";
    tmpValue.value = 650.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MScharmL") && !FindInFixed("MScharmL") && !FindInUniversality("MScharmL")) {  
    tmpValue.name = "MScharmL";
    tmpValue.value = 550.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MStopL") && !FindInFixed("MStopL") && !FindInUniversality("MStopL")) {  
    tmpValue.name = "MStopL";
    tmpValue.value = 450.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }

  if (!FindInFitted("MSupR") && !FindInFixed("MSupR") && !FindInUniversality("MSupR")) {  
    tmpValue.name = "MSupR";
    tmpValue.value = 660.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MScharmR") && !FindInFixed("MScharmR") && !FindInUniversality("MScharmR")) {  
    tmpValue.name = "MScharmR";
    tmpValue.value = 560.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MStopR") && !FindInFixed("MStopR") && !FindInUniversality("MStopR")) {  
    tmpValue.name = "MStopR";
    tmpValue.value = 460.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }

  if (!FindInFitted("MSdownR") && !FindInFixed("MSdownR") && !FindInUniversality("MSdownR")) {  
    tmpValue.name = "MSdownR";
    tmpValue.value = 670.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MSstrangeR") && !FindInFixed("MSstrangeR") && !FindInUniversality("MSstrangeR")) {  
    tmpValue.name = "MSstrangeR";
    tmpValue.value = 570.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }
  if (!FindInFitted("MSbottomR") && !FindInFixed("MSbottomR") && !FindInUniversality("MSbottomR")) {  
    tmpValue.name = "MSbottomR";
    tmpValue.value = 470.;
    tmpValue.error = -1.;
    tmpValue.bound_low = -1E+6;
    tmpValue.bound_up = 1E+6;
    tmpValue.id = 0;
    yyFixedVec.push_back(tmpValue);
  }



  return;

}

void   ReadLesHouches()
{
  unsigned int found_prod;
  // first parse the LesHouches outfile and load predictions to maps
  cout << "parsing the les houches file" << endl;
  ParseLesHouches();
  if (yyGeneratorError) {
    cerr << "Error in ParseLesHouches, exiting ReadLesHouches" << endl;
    return;
  }
  vector <unsigned int> used_products;
  bool already_used = false;

  used_products.clear();
  // loop over yyMeasuredVec and fill theovalue
  cout << "interpreting the les houches output " << endl;
  
  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
    if (yyMeasuredVec[i].type == mass) {
      //      cout << "foubd a mass" << endl;
      yyMeasuredVec[i].theovalue = yyMass[yyMeasuredVec[i].id];
      yyMeasuredVec[i].theoset = true;
    }
    else if (yyMeasuredVec[i].type == Pwidth) {
      yyMeasuredVec[i].theovalue = branching_ratios[yyMeasuredVec[i].id].TWidth;
      yyMeasuredVec[i].theoset = true;
    }
    else if (yyMeasuredVec[i].type == xsection) {
      //      cout << " looking for xs " << yyMeasuredVec[i].name << endl; 
      doubleVec_t tmp;
      tmp.clear();
      tmp.push_back(yyMeasuredVec[i].id);
      tmp.push_back(yyMeasuredVec[i].sqrts);
      tmp.push_back(yyMeasuredVec[i].polarisation1);
      tmp.push_back(yyMeasuredVec[i].polarisation2);
      if( cross_sections.find(tmp) != cross_sections.end() ) {
	//	cout << " map element found " << yyMeasuredVec[i].id << " " 
	//	     << yyMeasuredVec[i].polarisation1 << " " << yyMeasuredVec[i].polarisation2 << endl; 
	for (unsigned int j=0; 1; j++) {
	  found_prod = 0;
	  used_products.clear();
	  //          cout << " j = " << j << " products " << (int)give_xs(tmp, j, 2) << " " <<  (int)give_xs(tmp, j, 3) << endl;
	  if (give_xs(tmp, j, 1) == 0) break;
	  //          cout << "looping over cross_sections from 2 to " << give_xs(tmp, j, 1)+1 << endl;
	  for (unsigned int k=2; k<(2+give_xs(tmp, j, 1)); k++) {
	    for (unsigned int m=0;m<yyMeasuredVec[i].products.size();m++) {
	      for (unsigned int n=0;n< used_products.size();n++) {
		if (m==used_products[n]) {
		  already_used = true;
		  break;
		}
	      }
	      if (already_used) {
		already_used = false;
		continue;
	      }
	      //	      cout << " comparing " << yyMeasuredVec[i].products[m] << " with " << (int)give_xs(tmp, j, k) << endl;
	      //	      if (TMath::Abs(yyMeasuredVec[i].products[m]) == TMath::Abs((int)give_xs(tmp, j, k)) ) {
	      if (yyMeasuredVec[i].products[m] == (int)give_xs(tmp, j, k) ) {
		found_prod++;
		used_products.push_back(m);
		break;
	      }
	    }  
	  }
	  //          cout << "found_prod = " << found_prod << endl;
	  if ((found_prod == yyMeasuredVec[i].products.size()) && (found_prod == give_xs(tmp, j, 1))) {
	    //     	    cout << "cross section " << yyMeasuredVec[i].name << " with Ecm " << cross_sections[tmp].ecm << " "  << yyMeasuredVec[i].sqrts << endl;
	    if ((cross_sections[tmp].ecm == yyMeasuredVec[i].sqrts) && 
		(cross_sections[tmp].polarisation[0] == yyMeasuredVec[i].polarisation1) && 
		(cross_sections[tmp].polarisation[1] == yyMeasuredVec[i].polarisation2) ) {
	      yyMeasuredVec[i].theovalue = give_xs(tmp, j, 0);
	      yyMeasuredVec[i].theoset = true;
	      break;
	    }
	  } 
	}
      }
    }
    else if (yyMeasuredVec[i].type == br) {
      for (unsigned int j = 0; 1; j++ ) {
	used_products.clear();
	found_prod = 0;
	if (give_br(yyMeasuredVec[i].id, j, 1) == 0) break;
	for (unsigned int k=2; k<(2+give_br(yyMeasuredVec[i].id, j, 1)); k++) {
	  for (unsigned int m=0;m<yyMeasuredVec[i].daughters.size();m++) {
	    for (unsigned int n=0;n< used_products.size();n++) {
	      if (m==used_products[n]) {
		already_used = true;
		break;
	      }
	    }
	    if (already_used) {
	      already_used = false;
	      continue;
	    }
	    if (TMath::Abs(yyMeasuredVec[i].daughters[m]) == TMath::Abs((int)give_br(yyMeasuredVec[i].id, j, k)) ) {
	      used_products.push_back(m);
	      found_prod++;
	      break;
	    }	    
	  }
	}
	if ((found_prod == yyMeasuredVec[i].daughters.size()) && (found_prod == give_br(yyMeasuredVec[i].id, j, 1))) {
	  yyMeasuredVec[i].theovalue = give_br(yyMeasuredVec[i].id, j, 0);
	  yyMeasuredVec[i].theoset = true;
	  break;
	} 
      }
      if (yyMeasuredVec[i].theoset == false) {
	yyMeasuredVec[i].theovalue = 0.0;
	yyMeasuredVec[i].theoset = true;
      }
    }
    else if (yyMeasuredVec[i].type == Pedge) {
      //      cout << "in theory fitting edge " << yyMeasuredVec[i].name << endl;
      if (yyMeasuredVec[i].id == 1) {
	yyMeasuredVec[i].theovalue = yyMass[yyMeasuredVec[i].daughters[0]] + yyMass[yyMeasuredVec[i].daughters[1]];
	yyMeasuredVec[i].theoset = true;
      }
      else if (yyMeasuredVec[i].id == 2) {
	yyMeasuredVec[i].theovalue = TMath::Abs(yyMass[yyMeasuredVec[i].daughters[0]] - yyMass[yyMeasuredVec[i].daughters[1]]);
	yyMeasuredVec[i].theoset = true;
      }
    }
  }
  // check for xsbr
  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
    if (yyMeasuredVec[i].type == xsbr) {
      yyMeasuredVec[i].theovalue = 1.;
      for (unsigned int j=0; j<yyMeasuredVec[i].daughters.size(); j++) {
	// cout << "multiplying br " << yyMeasuredVec[yyMeasuredVec[i].daughters[j]].name << " = " << 
	//   yyMeasuredVec[yyMeasuredVec[i].daughters[j]].theovalue << endl;
	yyMeasuredVec[i].theovalue *= yyMeasuredVec[yyMeasuredVec[i].daughters[j]].theovalue;
      } 
      for (unsigned int j=0; j<yyMeasuredVec[i].products.size(); j++) {
	// cout << "multiplying xs " << yyMeasuredVec[yyMeasuredVec[i].products[j]].name << " = " << 
	//   yyMeasuredVec[yyMeasuredVec[i].products[j]].theovalue << endl;
	yyMeasuredVec[i].theovalue *= yyMeasuredVec[yyMeasuredVec[i].products[j]].theovalue;
      }        
      yyMeasuredVec[i].theoset = true;
    }
  }  
  // check for brratio
  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
    if (yyMeasuredVec[i].type == brratio) {
      yyMeasuredVec[i].theovalue = yyMeasuredVec[yyMeasuredVec[i].daughters[0]].theovalue / 
	yyMeasuredVec[yyMeasuredVec[i].daughters[1]].theovalue;
      yyMeasuredVec[i].theoset = true;
    }
  }  
  // check whether all has been found
  for (unsigned int i=0; i<yyMeasuredVec.size(); i++) {
    if ((yyMeasuredVec[i].theoset == false) && (yyMeasuredVec[i].nofit == false)) {
      cout << "WARNING: observable " << yyMeasuredVec[i].name << " not found!" << endl;
    }
    if (yyMeasuredVec[i].first == true) {
      yyMeasuredVec[i].first = false;
      if (yyMeasuredVec[i].theoset == true) {
	yyMeasuredVec[i].hasbeenset = true;
      }
    }
    if ((yyMeasuredVec[i].theoset == true) && (yyMeasuredVec[i].hasbeenset == false)) {
      yyMeasuredVec[i].hasbeenset = true;
    }
    if ((yyMeasuredVec[i].hasbeenset == true) && (yyMeasuredVec[i].theoset == false)) {
      yyMeasuredVec[i].theovalue = 11111111111.;
      yyMeasuredVec[i].theoset = true;
    }
  }

  return;

}

void   ParseLesHouches()
{
    static int counter = 0;

    counter++;

    if (yyGenerator == SPHENO) {
      //      system("cp SPheno.spc.test SPheno.spc");
      yyin = fopen("SPheno.spc", "r");
      //	yyin = fopen("leshouches.in", "r");
      if (!yyin) {
	cerr<<"SPheno.spc does not exist"<<endl;
	yyGeneratorError = true;
	return;
      }
      yyparse();
      fclose(yyin);
      system ("rm SPheno.spc");
    }
    else {
      cerr<<"Only SPHENO is implemented"<<endl;
      exit(EXIT_FAILURE);
    }

    cout<<counter<<" ###########################################################"<<endl;


    return;

}


double give_br ( int id, int decay, int element )
{
 
  if( branching_ratios.find( id ) != branching_ratios.end() ) {
    if (branching_ratios[id].decays.size() > (unsigned int)decay ) {
      if (branching_ratios[id].decays[decay].size() > (unsigned int)element) {
        return branching_ratios[id].decays[decay][element];
      } else {
        return 0.;
      }
    } else {
      return 0.;
    }
  } else {
    return 0.;
  }
 
}
 
double give_xs (doubleVec_t initial, int channel, int element )
{
 
  if( cross_sections.find( initial ) != cross_sections.end() ) {
    if (cross_sections[initial].xs.size() > (unsigned int)channel ) {
      if (cross_sections[initial].xs[channel].size() > (unsigned int)element) {
        return cross_sections[initial].xs[channel][element];
      } else {
        return 0.;
      }
    } else {
      return 0.;
    }
  } else {
    return 0.;
  }
 
}


// test simulated annealing

void Fittino::simulated_annealing (int iteration, TNtuple *ntuple)
{

  Double_t xdummy[100];
  int n;
  vector <double> x; 
  vector <double> xvar; 
  bool max = false; 
  double rt; 
  double eps = 0.0001; 
  int ns = 20; 
  int nt; 
  int neps = 4;
  int maxevl; 
  vector <double> lb; 
  vector <double> ub;
  vector <double> c; 
  int iprint = 0;
  int iseed1 = 31327; 
  int iseed2 = 30080; 
  double t = 5.0; // initial temperature
  vector <double> vm;
  vector <double> xopt;
  double fopt = 11111111111.0;
  int nacc = 0; 
  int nfcnev = 0; 
  int nobds = 0;
  int ier = 0; 
  vector <double> fstar; 
  vector <double> xp;
  vector <int> nacp;
  
  /* System generated locals */
  int i1, i2, i3, i4;
  double d1;
  Double_t dummyfloat = 5.;
  Int_t dummyint = 1;
  
  /* Local variables */
  static int nrej;
  static int nnew;
  static bool quit = false;
  bool quit_while_loop = false;
  static double f;
  static int h, i, j, m;
  static double p, ratio;
  static int ndown;
  static int nup;
  static double fp, pp;
  static int lnobds;
  double fvar = 0.;
  double fcubed = 0.;
  double fsum = 0.;
  int nvalid = 0;
  bool firstfalse;
  int accpoint = 0;
  int xoptflag = 0;
  int niter = 0;

  Float_t ntupvars[50];

  // set values
  maxevl = yyMaxCallsSimAnn;
  rt = yyTempRedSimAnn;
  nup = 0;
  nrej = 0;
  nnew = 0;
  ndown = 0;
  lnobds = 0;
  n = yyFittedVec.size();
  if (n<20) {
    nt = 60;
  } else {
    nt = 3*n;
  }
  for (i = 0; i < neps; ++i) {
    fstar.push_back(1e20);
  }

  // fill vector of parameters
  for (unsigned int k = 0; k < yyFittedVec.size(); k++ ) {
    x.push_back(yyFittedVec[k].value);
    xp.push_back(yyFittedVec[k].value);
    xvar.push_back(yyFittedVec[k].value);
    xopt.push_back(yyFittedVec[k].value);
    vm.push_back(yyFittedVec[k].error);
    lb.push_back(yyFittedVec[k].bound_low);
    ub.push_back(yyFittedVec[k].bound_up);
    c.push_back(2.0);
    nacp.push_back(0);
  }

  /*  Evaluate function to be minimized at the starting point */
  for (unsigned int ii = 0; ii < xp.size(); ii++) {
    xdummy[ii] = x[ii];
  }
  fitterFCN(dummyint, &dummyfloat, f, xdummy, 0);
  f = -f;
  ++nfcnev;
  fopt = f;
  fstar[1] = f;

  // ------------------------------------------------------------------
  // if initial temperature is given:
  if (yyInitTempSimAnn>0.) {
    t = yyInitTempSimAnn;
  }
  else {
    // ------------------------------------------------------------------
    // if already very close to the minimum:
    if (TMath::Abs(f)<0.1) {
      t = 0.1;
      if (rt<0.6) {
	rt = 0.6;
      }
    } else if (TMath::Abs(f)<1.) {
      t = 1.;
      if (rt<0.6) {
	rt = 0.6;
      }
    } else if (TMath::Abs(f)<5.) {
      t = 5.;
      if (rt<0.6) {
	rt = 0.6;
      }    
    } else if (TMath::Abs(f)<10.) {
      t = 10.;
      if (rt<0.6) {
	rt = 0.6;
      }    
    } else if (TMath::Abs(f)<100.) {
      t = 100.;
      if (rt<0.6) {
	rt = 0.6;
      }    
    } else if (TMath::Abs(f)<200.) {
      t = 200.;
      if (rt<0.5) {
	rt = 0.5;
      }    
    } else {
      //-------------------------------------------
      // if not very close to the minimum:
      // first adjust the temperature to the variance for variations within vm?
      nvalid = 0;
      fsum = 0.;
      fcubed = 0.;
      for (m = 0; m < 10; m++) {
	for (i = 0; i < n; i++) {
	  xp[i] = xvar[i] + gRandom->Uniform(-1.,1.) * vm[i];
	  //  If variation has led xp out of bounds:
	  while ( (xp[i] < lb[i]) || (xp[i] > ub[i]) ) {
	    xp[i] = gRandom->Uniform(-1.,1.) * vm[i] + yyFittedVec[i].value;
	  }
	  // Then call the function to be minimized at the starting point
	  for (unsigned int ii = 0; ii < xp.size(); ii++) {
	    xdummy[ii] = xp[ii];
	  }
	  fitterFCN(dummyint, &dummyfloat, fp, xdummy, 0);
	  fp = -fp;
	  cout << "fp = " << fp << endl;
	  if (fp > -1e10) {
	    //	fopt = 0.;
	    //      fp = 5.;
	    fsum = fsum + TMath::Abs(fopt - fp);
	    fcubed = fcubed + sqr((fopt - fp));
	    nvalid++;
	  }
	  xp[i] = xvar[i];
	}
	for (i = 0; i < n; i++) {
	  xvar[i] = xp[i];
	}
      }
      // set t to half of the variance
      t = TMath::Sqrt(fcubed/double(nvalid) - sqr(fsum/double(nvalid)))/2.;
    }
  }
  cout << "temperature chosen " << t << endl;

  //-------------------------------------------
  // perform the optimization
  niter = 0;

  // begin new temperature era
  while (!quit_while_loop) {
    nup = 0;
    nrej = 0;
    nnew = 0;
    ndown = 0;
    lnobds = 0;
    niter++;
    // loop over the iterations before temperature reduction:
    for (m = 0; m < nt; m++) {
      // loop over the accepted function evaluations:
      for (j = 0; j < ns; j++) {
	// loop over the variables:
	for (h = 0; h < n; h++) {
	  for (i = 0; i < n; i++) {
	    if (i == h) {
	      xp[i] = x[i] + gRandom->Uniform(-1.,1.) * vm[i];
	    } else {
	      xp[i] = x[i];
	    }
	    //  If variation has led xp out of bounds:
	    firstfalse = true;
	    while ( (xp[i] < lb[i]) || (xp[i] > ub[i]) ) {
	      xp[i] = gRandom->Uniform(-2.,2.) * vm[i] + yyFittedVec[i].value;
	      if (firstfalse) {
		++lnobds;
		++(nobds);
	      }
	      firstfalse = false;
	    }
	  }
	  // Then call the function to be minimized at the starting point
	  for (unsigned int ii = 0; ii < xp.size(); ii++) {
	    xdummy[ii] = xp[ii];
	  }
	  fitterFCN(dummyint, &dummyfloat, fp, xdummy, 0);
	  fp = -fp;
	  ++nfcnev;
	  xoptflag = 0;
	  cout << "running simulated annealing at t = " << t << " m(<"<<nt<<") = " << m << " j(<"<<ns<<") = " << j << " h(<"<<n<<") = " << h << " in niter = " << niter << endl;
	  //  For too many function evaluations, terminate:
	  if (nfcnev >= maxevl) {
	    cout << "terminating simulated annealing prematurely due to number of calls" << endl;
	    fopt = -fopt;
	    ier = 1;
	    for (unsigned int k = 0; k < yyFittedVec.size(); k++ ) {
	      yyFittedVec[k].value = xopt[k];
	    }
	    quit_while_loop = true;
	    return;
	  }
	  //  Accept new point if better chisq
	  if (fp >= f) {
	    for (i = 0; i < n; ++i) {
	      x[i] = xp[i];
	    }
	    f = fp;
	    ++nacc;
	    ++nacp[h];
	    ++nup;
	    accpoint = 1;
	    //  If best chisq so far, record as new best parameter set:
	    if (fp > fopt) {
	      for (i = 0; i < n; ++i) {
		xopt[i] = xp[i];
	      }
	      fopt = fp;
	      ++nnew;
	      xoptflag = 1;
	      cout << "found new optimum at chisq = " << fopt << endl;
	    }
	    //  Metropolis criteria to decide whether point is accepted if chisq is worse than before 
	  } else {
	    d1 = (fp - f) / t;
	    p = TMath::Exp(d1);
	    pp = gRandom->Uniform(1.);
	    if (pp < p) {
	      for (i = 1; i < n; ++i) {
		x[i] = xp[i];
	      }
	      f = fp;
	      ++nacc;
	      ++nacp[h];
	      ++ndown;
	      accpoint = 1;
	    } else {
	      ++nrej;
	      accpoint = 0;
	    }
	  }
	  // write point to ntuple
	  ntupvars[0] = (Float_t)nfcnev;
	  ntupvars[1] = (Float_t)t;
	  ntupvars[2] = -(Float_t)fp;
	  ntupvars[3] = (Float_t)accpoint;
	  ntupvars[4] = (Float_t)xoptflag;
	  for (unsigned int ii = 5; ii < 5+yyFittedVec.size(); ii++) {
	    ntupvars[ii] = xp[ii-5];
	  }
	  ntuple->Fill(ntupvars);
	} // close the outer variable loop
      } // close the loop over the accepted function evaluations
      // go on to adjust vm and t
      // adjust vm to accept between 0.4 and 0.6 of the trial points at the given temperature
      for (i = 0; i < n; ++i) {
	ratio = (double) nacp[i] / (double) (ns);
	if (ratio > 0.6) {
	  vm[i] *= c[i] * (ratio - 0.6) / 0.4 + 1.;
	} else if (ratio < 0.4) {
	  vm[i] /= c[i] * (0.4 - ratio) / 0.4 + 1.;
	}
	if (vm[i] > ub[i] - lb[i]) {
	  vm[i] = ub[i] - lb[i];
	}	
      }
      for (i = 0; i < n; ++i) {
	nacp[i] = 0;
      }
    } // close the loop over the iterations before temperature reduction
    // check whether to terminate
    quit = false;
    fstar[1] = f;
    if (TMath::Abs(fopt - fstar[1]) <= eps) {
      quit = true;
    }
    for (i = 0; i < neps; ++i) {
      if (TMath::Abs(f - fstar[i]) > eps) {
	quit = false;
      }
    }
    //  Terminate simulated annealing if appropriate 
    if (quit) {
      for (i = 0; i < n; ++i) {
	x[i] = xopt[i];
      }
      ier = 0;
      fopt = -fopt;
      quit_while_loop = true;
    } else {
      // If not terminating, prepare t
      cout << "starting new temperature loop at t = " << rt * t <<  endl;
      t = rt * t;
      for (i = neps-1; i >= 1; --i) {
	fstar[i] = fstar[i - 1];
      }
      f = fopt;
      for (i = 0; i < n; ++i) {
	x[i] = xopt[i];
      }
    }
  } // close the while loop 

  //--------------------------------------------
  // write output
    
  for (unsigned int k = 0; k < yyFittedVec.size(); k++ ) {
    yyFittedVec[k].value = xopt[k];
  }
  cout << "optimal function value " << fopt << " after " << nacc << " evaluations" << endl;

}

! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.11.0 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 11:53 on 24.4.2017   
! ----------------------------------------------------------------------  
 
 
Module CouplingsCT_MSSMTriLnV
 
Use Control 
Use Settings 
Use Model_Data_MSSMTriLnV 
Use Mathematics, Only: CompareMatrices, Adjungate 
 
Contains 
 
 Subroutine CalculateCouplingCT(g1,g2,vd,vu,ZH,ZA,ZP,Mu,Yd,Td,ZD,Ye,Te,ZE,             & 
& Yu,Tu,ZU,ZV,L2,T2,L1,T1,TW,g3,UM,UP,ZN,ZDL,ZDR,ZEL,ZER,ZUL,ZUR,pG,dg1,dg2,             & 
& dvd,dvu,dZH,dZA,dZP,dMu,dYd,dTd,dZD,dYe,dTe,dZE,dYu,dTu,dZU,dZV,dL2,dT2,               & 
& dL1,dT1,dSinTW,dCosTW,dTanTW,dg3,dUM,dUP,dZN,dZDL,dZDR,dZEL,dZER,dZUL,dZUR,            & 
& dpG,ctcplAhAhhh,ctcplAhHpmcHpm,ctcplAhSdcSd,ctcplAhSecSe,ctcplAhSucSu,ctcplhhhhhh,     & 
& ctcplhhHpmcHpm,ctcplhhSdcSd,ctcplhhSecSe,ctcplhhSucSu,ctcplhhSvcSv,ctcplHpmSucSd,      & 
& ctcplHpmSvcSe,ctcplSdSvcSd,ctcplSdcHpmcSu,ctcplSdcSdcSv,ctcplSdcSecSu,ctcplSeSucSd,    & 
& ctcplSeSvcSe,ctcplSecHpmcSv,ctcplSecSecSv,ctcplAhhhVZ,ctcplAhHpmcVWm,ctcplAhcHpmVWm,   & 
& ctcplhhHpmcVWm,ctcplhhcHpmVWm,ctcplHpmcHpmVP,ctcplHpmcHpmVZ,ctcplSdcSdVG,              & 
& ctcplSdcSdVP,ctcplSdcSdVZ,ctcplSdcSucVWm,ctcplSecSeVP,ctcplSecSeVZ,ctcplSecSvcVWm,     & 
& ctcplSucSuVG,ctcplSucSuVP,ctcplSucSdVWm,ctcplSucSuVZ,ctcplSvcSeVWm,ctcplSvcSvVZ,       & 
& ctcplhhcVWmVWm,ctcplhhVZVZ,ctcplHpmcVWmVP,ctcplHpmcVWmVZ,ctcplcHpmVPVWm,               & 
& ctcplcHpmVWmVZ,ctcplVGVGVG,ctcplcVWmVPVWm,ctcplcVWmVWmVZ,ctcplcChaChaAhL,              & 
& ctcplcChaChaAhR,ctcplChiChiAhL,ctcplChiChiAhR,ctcplcFdFdAhL,ctcplcFdFdAhR,             & 
& ctcplcFeFeAhL,ctcplcFeFeAhR,ctcplcFuFuAhL,ctcplcFuFuAhR,ctcplChiChacHpmL,              & 
& ctcplChiChacHpmR,ctcplChaFucSdL,ctcplChaFucSdR,ctcplFvChacSeL,ctcplFvChacSeR,          & 
& ctcplcChaChahhL,ctcplcChaChahhR,ctcplcFdChaSuL,ctcplcFdChaSuR,ctcplcFeChaSvL,          & 
& ctcplcFeChaSvR,ctcplChiChihhL,ctcplChiChihhR,ctcplChiFdcSdL,ctcplChiFdcSdR,            & 
& ctcplChiFecSeL,ctcplChiFecSeR,ctcplChiFucSuL,ctcplChiFucSuR,ctcplChiFvSvL,             & 
& ctcplChiFvSvR,ctcplChiFvcSvL,ctcplChiFvcSvR,ctcplcChaChiHpmL,ctcplcChaChiHpmR,         & 
& ctcplcFdChiSdL,ctcplcFdChiSdR,ctcplcFeChiSeL,ctcplcFeChiSeR,ctcplcFuChiSuL,            & 
& ctcplcFuChiSuR,ctcplFvFdcSdL,ctcplFvFdcSdR,ctcplGluFdcSdL,ctcplGluFdcSdR,              & 
& ctcplcFdFdhhL,ctcplcFdFdhhR,ctcplcFdFdSvL,ctcplcFdFdSvR,ctcplcChaFdcSuL,               & 
& ctcplcChaFdcSuR,ctcplcFdFdcSvL,ctcplcFdFdcSvR,ctcplcFeFdcSuL,ctcplcFeFdcSuR,           & 
& ctcplcFuFdcHpmL,ctcplcFuFdcHpmR,ctcplcFuFdcSeL,ctcplcFuFdcSeR,ctcplFeFucSdL,           & 
& ctcplFeFucSdR,ctcplFvFecHpmL,ctcplFvFecHpmR,ctcplFvFecSeL,ctcplFvFecSeR,               & 
& ctcplcFeFehhL,ctcplcFeFehhR,ctcplcFdFeSuL,ctcplcFdFeSuR,ctcplcFeFeSvL,ctcplcFeFeSvR,   & 
& ctcplcChaFecSvL,ctcplcChaFecSvR,ctcplcFeFecSvL,ctcplcFeFecSvR,ctcplGluFucSuL,          & 
& ctcplGluFucSuR,ctcplcFuFuhhL,ctcplcFuFuhhR,ctcplcFdFuHpmL,ctcplcFdFuHpmR,              & 
& ctcplcFdFuSeL,ctcplcFdFuSeR,ctcplcFeFvHpmL,ctcplcFeFvHpmR,ctcplcFdFvSdL,               & 
& ctcplcFdFvSdR,ctcplcChaFvSeL,ctcplcChaFvSeR,ctcplcFeFvSeL,ctcplcFeFvSeR,               & 
& ctcplcFdGluSdL,ctcplcFdGluSdR,ctcplcFuGluSuL,ctcplcFuGluSuR,ctcplcChacFuSdL,           & 
& ctcplcChacFuSdR,ctcplcFecFuSdL,ctcplcFecFuSdR,ctcplChiChacVWmL,ctcplChiChacVWmR,       & 
& ctcplcChaChaVPL,ctcplcChaChaVPR,ctcplcChaChaVZL,ctcplcChaChaVZR,ctcplChiChiVZL,        & 
& ctcplChiChiVZR,ctcplcChaChiVWmL,ctcplcChaChiVWmR,ctcplcFdFdVGL,ctcplcFdFdVGR,          & 
& ctcplcFdFdVPL,ctcplcFdFdVPR,ctcplcFdFdVZL,ctcplcFdFdVZR,ctcplcFuFdcVWmL,               & 
& ctcplcFuFdcVWmR,ctcplFvFecVWmL,ctcplFvFecVWmR,ctcplcFeFeVPL,ctcplcFeFeVPR,             & 
& ctcplcFeFeVZL,ctcplcFeFeVZR,ctcplcFuFuVGL,ctcplcFuFuVGR,ctcplcFuFuVPL,ctcplcFuFuVPR,   & 
& ctcplcFdFuVWmL,ctcplcFdFuVWmR,ctcplcFuFuVZL,ctcplcFuFuVZR,ctcplFvFvVZL,ctcplFvFvVZR,   & 
& ctcplcFeFvVWmL,ctcplcFeFvVWmR,ctcplGluGluVGL,ctcplGluGluVGR)

Implicit None 
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),ZA(2,2),ZP(2,2),TW,g3,dg1,dg2,dvd,dvu,dZH(2,2),dZA(2,2),          & 
& dZP(2,2),dSinTW,dCosTW,dTanTW,dg3

Complex(dp), Intent(in) :: Mu,Yd(3,3),Td(3,3),ZD(6,6),Ye(3,3),Te(3,3),ZE(6,6),Yu(3,3),Tu(3,3),ZU(6,6),           & 
& ZV(3,3),L2(3,3,3),T2(3,3,3),L1(3,3,3),T1(3,3,3),UM(2,2),UP(2,2),ZN(4,4),               & 
& ZDL(3,3),ZDR(3,3),ZEL(3,3),ZER(3,3),ZUL(3,3),ZUR(3,3),pG,dMu,dYd(3,3),dTd(3,3),        & 
& dZD(6,6),dYe(3,3),dTe(3,3),dZE(6,6),dYu(3,3),dTu(3,3),dZU(6,6),dZV(3,3),               & 
& dL2(3,3,3),dT2(3,3,3),dL1(3,3,3),dT1(3,3,3),dUM(2,2),dUP(2,2),dZN(4,4),dZDL(3,3),      & 
& dZDR(3,3),dZEL(3,3),dZER(3,3),dZUL(3,3),dZUR(3,3),dpG

Complex(dp), Intent(out) :: ctcplAhAhhh(2,2,2),ctcplAhHpmcHpm(2,2,2),ctcplAhSdcSd(2,6,6),ctcplAhSecSe(2,6,6),     & 
& ctcplAhSucSu(2,6,6),ctcplhhhhhh(2,2,2),ctcplhhHpmcHpm(2,2,2),ctcplhhSdcSd(2,6,6),      & 
& ctcplhhSecSe(2,6,6),ctcplhhSucSu(2,6,6),ctcplhhSvcSv(2,3,3),ctcplHpmSucSd(2,6,6),      & 
& ctcplHpmSvcSe(2,3,6),ctcplSdSvcSd(6,3,6),ctcplSdcHpmcSu(6,2,6),ctcplSdcSdcSv(6,6,3),   & 
& ctcplSdcSecSu(6,6,6),ctcplSeSucSd(6,6,6),ctcplSeSvcSe(6,3,6),ctcplSecHpmcSv(6,2,3),    & 
& ctcplSecSecSv(6,6,3),ctcplAhhhVZ(2,2),ctcplAhHpmcVWm(2,2),ctcplAhcHpmVWm(2,2),         & 
& ctcplhhHpmcVWm(2,2),ctcplhhcHpmVWm(2,2),ctcplHpmcHpmVP(2,2),ctcplHpmcHpmVZ(2,2),       & 
& ctcplSdcSdVG(6,6),ctcplSdcSdVP(6,6),ctcplSdcSdVZ(6,6),ctcplSdcSucVWm(6,6),             & 
& ctcplSecSeVP(6,6),ctcplSecSeVZ(6,6),ctcplSecSvcVWm(6,3),ctcplSucSuVG(6,6),             & 
& ctcplSucSuVP(6,6),ctcplSucSdVWm(6,6),ctcplSucSuVZ(6,6),ctcplSvcSeVWm(3,6),             & 
& ctcplSvcSvVZ(3,3),ctcplhhcVWmVWm(2),ctcplhhVZVZ(2),ctcplHpmcVWmVP(2),ctcplHpmcVWmVZ(2),& 
& ctcplcHpmVPVWm(2),ctcplcHpmVWmVZ(2),ctcplVGVGVG,ctcplcVWmVPVWm,ctcplcVWmVWmVZ,         & 
& ctcplcChaChaAhL(2,2,2),ctcplcChaChaAhR(2,2,2),ctcplChiChiAhL(4,4,2),ctcplChiChiAhR(4,4,2),& 
& ctcplcFdFdAhL(3,3,2),ctcplcFdFdAhR(3,3,2),ctcplcFeFeAhL(3,3,2),ctcplcFeFeAhR(3,3,2),   & 
& ctcplcFuFuAhL(3,3,2),ctcplcFuFuAhR(3,3,2),ctcplChiChacHpmL(4,2,2),ctcplChiChacHpmR(4,2,2),& 
& ctcplChaFucSdL(2,3,6),ctcplChaFucSdR(2,3,6),ctcplFvChacSeL(3,2,6),ctcplFvChacSeR(3,2,6),& 
& ctcplcChaChahhL(2,2,2),ctcplcChaChahhR(2,2,2),ctcplcFdChaSuL(3,2,6),ctcplcFdChaSuR(3,2,6),& 
& ctcplcFeChaSvL(3,2,3),ctcplcFeChaSvR(3,2,3),ctcplChiChihhL(4,4,2),ctcplChiChihhR(4,4,2),& 
& ctcplChiFdcSdL(4,3,6),ctcplChiFdcSdR(4,3,6),ctcplChiFecSeL(4,3,6),ctcplChiFecSeR(4,3,6),& 
& ctcplChiFucSuL(4,3,6),ctcplChiFucSuR(4,3,6),ctcplChiFvSvL(4,3,3),ctcplChiFvSvR(4,3,3), & 
& ctcplChiFvcSvL(4,3,3),ctcplChiFvcSvR(4,3,3),ctcplcChaChiHpmL(2,4,2),ctcplcChaChiHpmR(2,4,2),& 
& ctcplcFdChiSdL(3,4,6),ctcplcFdChiSdR(3,4,6),ctcplcFeChiSeL(3,4,6),ctcplcFeChiSeR(3,4,6),& 
& ctcplcFuChiSuL(3,4,6),ctcplcFuChiSuR(3,4,6),ctcplFvFdcSdL(3,3,6),ctcplFvFdcSdR(3,3,6), & 
& ctcplGluFdcSdL(3,6),ctcplGluFdcSdR(3,6),ctcplcFdFdhhL(3,3,2),ctcplcFdFdhhR(3,3,2),     & 
& ctcplcFdFdSvL(3,3,3),ctcplcFdFdSvR(3,3,3),ctcplcChaFdcSuL(2,3,6),ctcplcChaFdcSuR(2,3,6),& 
& ctcplcFdFdcSvL(3,3,3),ctcplcFdFdcSvR(3,3,3),ctcplcFeFdcSuL(3,3,6),ctcplcFeFdcSuR(3,3,6),& 
& ctcplcFuFdcHpmL(3,3,2),ctcplcFuFdcHpmR(3,3,2),ctcplcFuFdcSeL(3,3,6),ctcplcFuFdcSeR(3,3,6),& 
& ctcplFeFucSdL(3,3,6),ctcplFeFucSdR(3,3,6),ctcplFvFecHpmL(3,3,2),ctcplFvFecHpmR(3,3,2), & 
& ctcplFvFecSeL(3,3,6),ctcplFvFecSeR(3,3,6),ctcplcFeFehhL(3,3,2),ctcplcFeFehhR(3,3,2),   & 
& ctcplcFdFeSuL(3,3,6),ctcplcFdFeSuR(3,3,6),ctcplcFeFeSvL(3,3,3),ctcplcFeFeSvR(3,3,3),   & 
& ctcplcChaFecSvL(2,3,3),ctcplcChaFecSvR(2,3,3),ctcplcFeFecSvL(3,3,3),ctcplcFeFecSvR(3,3,3),& 
& ctcplGluFucSuL(3,6),ctcplGluFucSuR(3,6),ctcplcFuFuhhL(3,3,2),ctcplcFuFuhhR(3,3,2),     & 
& ctcplcFdFuHpmL(3,3,2),ctcplcFdFuHpmR(3,3,2),ctcplcFdFuSeL(3,3,6),ctcplcFdFuSeR(3,3,6), & 
& ctcplcFeFvHpmL(3,3,2),ctcplcFeFvHpmR(3,3,2),ctcplcFdFvSdL(3,3,6),ctcplcFdFvSdR(3,3,6), & 
& ctcplcChaFvSeL(2,3,6),ctcplcChaFvSeR(2,3,6),ctcplcFeFvSeL(3,3,6),ctcplcFeFvSeR(3,3,6), & 
& ctcplcFdGluSdL(3,6),ctcplcFdGluSdR(3,6),ctcplcFuGluSuL(3,6),ctcplcFuGluSuR(3,6),       & 
& ctcplcChacFuSdL(2,3,6),ctcplcChacFuSdR(2,3,6),ctcplcFecFuSdL(3,3,6),ctcplcFecFuSdR(3,3,6)

Complex(dp), Intent(out) :: ctcplChiChacVWmL(4,2),ctcplChiChacVWmR(4,2),ctcplcChaChaVPL(2,2),ctcplcChaChaVPR(2,2), & 
& ctcplcChaChaVZL(2,2),ctcplcChaChaVZR(2,2),ctcplChiChiVZL(4,4),ctcplChiChiVZR(4,4),     & 
& ctcplcChaChiVWmL(2,4),ctcplcChaChiVWmR(2,4),ctcplcFdFdVGL(3,3),ctcplcFdFdVGR(3,3),     & 
& ctcplcFdFdVPL(3,3),ctcplcFdFdVPR(3,3),ctcplcFdFdVZL(3,3),ctcplcFdFdVZR(3,3),           & 
& ctcplcFuFdcVWmL(3,3),ctcplcFuFdcVWmR(3,3),ctcplFvFecVWmL(3,3),ctcplFvFecVWmR(3,3),     & 
& ctcplcFeFeVPL(3,3),ctcplcFeFeVPR(3,3),ctcplcFeFeVZL(3,3),ctcplcFeFeVZR(3,3),           & 
& ctcplcFuFuVGL(3,3),ctcplcFuFuVGR(3,3),ctcplcFuFuVPL(3,3),ctcplcFuFuVPR(3,3),           & 
& ctcplcFdFuVWmL(3,3),ctcplcFdFuVWmR(3,3),ctcplcFuFuVZL(3,3),ctcplcFuFuVZR(3,3),         & 
& ctcplFvFvVZL(3,3),ctcplFvFvVZR(3,3),ctcplcFeFvVWmL(3,3),ctcplcFeFvVWmR(3,3),           & 
& ctcplGluGluVGL,ctcplGluGluVGR

Integer :: gt1, gt2, gt3, gt4, ct1, ct2, ct3, ct4

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateCouplingCT'
 
ctcplAhAhhh = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
  Do gt3 = 1, 2
Call CT_CouplingAhAhhh(gt1,gt2,gt3,g1,g2,vd,vu,ZH,ZA,dg1,dg2,dvd,dvu,dZH,             & 
& dZA,ctcplAhAhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplAhHpmcHpm = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
  Do gt3 = 1, 2
Call CT_CouplingAhHpmcHpm(gt1,gt2,gt3,g2,vd,vu,ZA,ZP,dg2,dvd,dvu,dZA,dZP,             & 
& ctcplAhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplAhSdcSd = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplingAhSdcSd(gt1,gt2,gt3,Mu,Yd,Td,ZD,ZA,dMu,dYd,dTd,dZD,dZA,               & 
& ctcplAhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplAhSecSe = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplingAhSecSe(gt1,gt2,gt3,Mu,Ye,Te,ZE,ZA,dMu,dYe,dTe,dZE,dZA,               & 
& ctcplAhSecSe(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplAhSucSu = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplingAhSucSu(gt1,gt2,gt3,Mu,Yu,Tu,ZU,ZA,dMu,dYu,dTu,dZU,dZA,               & 
& ctcplAhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplhhhhhh = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
  Do gt3 = 1, 2
Call CT_Couplinghhhhhh(gt1,gt2,gt3,g1,g2,vd,vu,ZH,dg1,dg2,dvd,dvu,dZH,ctcplhhhhhh(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplhhHpmcHpm = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
  Do gt3 = 1, 2
Call CT_CouplinghhHpmcHpm(gt1,gt2,gt3,g1,g2,vd,vu,ZH,ZP,dg1,dg2,dvd,dvu,              & 
& dZH,dZP,ctcplhhHpmcHpm(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplhhSdcSd = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplinghhSdcSd(gt1,gt2,gt3,g1,g2,Mu,Yd,Td,vd,vu,ZD,ZH,dg1,dg2,               & 
& dMu,dYd,dTd,dvd,dvu,dZD,dZH,ctcplhhSdcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplhhSecSe = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplinghhSecSe(gt1,gt2,gt3,g1,g2,Mu,Ye,Te,vd,vu,ZE,ZH,dg1,dg2,               & 
& dMu,dYe,dTe,dvd,dvu,dZE,dZH,ctcplhhSecSe(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplhhSucSu = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplinghhSucSu(gt1,gt2,gt3,g1,g2,Mu,Yu,Tu,vd,vu,ZU,ZH,dg1,dg2,               & 
& dMu,dYu,dTu,dvd,dvu,dZU,dZH,ctcplhhSucSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplhhSvcSv = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplinghhSvcSv(gt1,gt2,gt3,g1,g2,vd,vu,ZH,dg1,dg2,dvd,dvu,dZH,               & 
& ctcplhhSvcSv(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplHpmSucSd = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplingHpmSucSd(gt1,gt2,gt3,g2,Mu,Yd,Td,Yu,Tu,vd,vu,ZD,ZU,ZP,dg2,            & 
& dMu,dYd,dTd,dYu,dTu,dvd,dvu,dZD,dZU,dZP,ctcplHpmSucSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplHpmSvcSe = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingHpmSvcSe(gt1,gt2,gt3,g2,Mu,Ye,Te,vd,vu,ZV,ZE,ZP,dg2,dMu,              & 
& dYe,dTe,dvd,dvu,dZV,dZE,dZP,ctcplHpmSvcSe(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSdSvcSd = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingSdSvcSd(gt1,gt2,gt3,Yd,L2,T2,vd,ZD,ZV,dYd,dL2,dT2,dvd,dZD,            & 
& dZV,ctcplSdSvcSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSdcHpmcSu = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 2
  Do gt3 = 1, 6
Call CT_CouplingSdcHpmcSu(gt1,gt2,gt3,g2,Mu,Yd,Td,Yu,Tu,vd,vu,ZD,ZU,ZP,               & 
& dg2,dMu,dYd,dTd,dYu,dTu,dvd,dvu,dZD,dZU,dZP,ctcplSdcHpmcSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSdcSdcSv = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
  Do gt3 = 1, 3
Call CT_CouplingSdcSdcSv(gt1,gt2,gt3,Yd,L2,T2,vd,ZD,ZV,dYd,dL2,dT2,dvd,               & 
& dZD,dZV,ctcplSdcSdcSv(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSdcSecSu = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplingSdcSecSu(gt1,gt2,gt3,Yd,Ye,L2,T2,Yu,vd,vu,ZD,ZU,ZE,dYd,               & 
& dYe,dL2,dT2,dYu,dvd,dvu,dZD,dZU,dZE,ctcplSdcSecSu(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSeSucSd = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
  Do gt3 = 1, 6
Call CT_CouplingSeSucSd(gt1,gt2,gt3,Yd,Ye,L2,T2,Yu,vd,vu,ZD,ZU,ZE,dYd,dYe,            & 
& dL2,dT2,dYu,dvd,dvu,dZD,dZU,dZE,ctcplSeSucSd(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSeSvcSe = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingSeSvcSe(gt1,gt2,gt3,Ye,L1,T1,vd,ZV,ZE,dYe,dL1,dT1,dvd,dZV,            & 
& dZE,ctcplSeSvcSe(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSecHpmcSv = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 2
  Do gt3 = 1, 3
Call CT_CouplingSecHpmcSv(gt1,gt2,gt3,g2,Mu,Ye,Te,vd,vu,ZV,ZE,ZP,dg2,dMu,             & 
& dYe,dTe,dvd,dvu,dZV,dZE,dZP,ctcplSecHpmcSv(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplSecSecSv = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
  Do gt3 = 1, 3
Call CT_CouplingSecSecSv(gt1,gt2,gt3,Ye,L1,T1,vd,ZV,ZE,dYe,dL1,dT1,dvd,               & 
& dZV,dZE,ctcplSecSecSv(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplAhhhVZ = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplingAhhhVZ(gt1,gt2,g1,g2,ZH,ZA,TW,dg1,dg2,dZH,dZA,dSinTW,dCosTW,          & 
& dTanTW,ctcplAhhhVZ(gt1,gt2))

 End Do 
End Do 


ctcplAhHpmcVWm = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplingAhHpmcVWm(gt1,gt2,g2,ZA,ZP,dg2,dZA,dZP,ctcplAhHpmcVWm(gt1,gt2))

 End Do 
End Do 


ctcplAhcHpmVWm = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplingAhcHpmVWm(gt1,gt2,g2,ZA,ZP,dg2,dZA,dZP,ctcplAhcHpmVWm(gt1,gt2))

 End Do 
End Do 


ctcplhhHpmcVWm = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplinghhHpmcVWm(gt1,gt2,g2,ZH,ZP,dg2,dZH,dZP,ctcplhhHpmcVWm(gt1,gt2))

 End Do 
End Do 


ctcplhhcHpmVWm = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplinghhcHpmVWm(gt1,gt2,g2,ZH,ZP,dg2,dZH,dZP,ctcplhhcHpmVWm(gt1,gt2))

 End Do 
End Do 


ctcplHpmcHpmVP = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplingHpmcHpmVP(gt1,gt2,g1,g2,ZP,TW,dg1,dg2,dZP,dSinTW,dCosTW,              & 
& dTanTW,ctcplHpmcHpmVP(gt1,gt2))

 End Do 
End Do 


ctcplHpmcHpmVZ = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplingHpmcHpmVZ(gt1,gt2,g1,g2,ZP,TW,dg1,dg2,dZP,dSinTW,dCosTW,              & 
& dTanTW,ctcplHpmcHpmVZ(gt1,gt2))

 End Do 
End Do 


ctcplSdcSdVG = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSdcSdVG(gt1,gt2,g3,dg3,ctcplSdcSdVG(gt1,gt2))

 End Do 
End Do 


ctcplSdcSdVP = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSdcSdVP(gt1,gt2,g1,g2,ZD,TW,dg1,dg2,dZD,dSinTW,dCosTW,dTanTW,         & 
& ctcplSdcSdVP(gt1,gt2))

 End Do 
End Do 


ctcplSdcSdVZ = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSdcSdVZ(gt1,gt2,g1,g2,ZD,TW,dg1,dg2,dZD,dSinTW,dCosTW,dTanTW,         & 
& ctcplSdcSdVZ(gt1,gt2))

 End Do 
End Do 


ctcplSdcSucVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSdcSucVWm(gt1,gt2,g2,ZD,ZU,dg2,dZD,dZU,ctcplSdcSucVWm(gt1,gt2))

 End Do 
End Do 


ctcplSecSeVP = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSecSeVP(gt1,gt2,g1,g2,ZE,TW,dg1,dg2,dZE,dSinTW,dCosTW,dTanTW,         & 
& ctcplSecSeVP(gt1,gt2))

 End Do 
End Do 


ctcplSecSeVZ = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSecSeVZ(gt1,gt2,g1,g2,ZE,TW,dg1,dg2,dZE,dSinTW,dCosTW,dTanTW,         & 
& ctcplSecSeVZ(gt1,gt2))

 End Do 
End Do 


ctcplSecSvcVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 3
Call CT_CouplingSecSvcVWm(gt1,gt2,g2,ZV,ZE,dg2,dZV,dZE,ctcplSecSvcVWm(gt1,gt2))

 End Do 
End Do 


ctcplSucSuVG = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSucSuVG(gt1,gt2,g3,dg3,ctcplSucSuVG(gt1,gt2))

 End Do 
End Do 


ctcplSucSuVP = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSucSuVP(gt1,gt2,g1,g2,ZU,TW,dg1,dg2,dZU,dSinTW,dCosTW,dTanTW,         & 
& ctcplSucSuVP(gt1,gt2))

 End Do 
End Do 


ctcplSucSdVWm = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSucSdVWm(gt1,gt2,g2,ZD,ZU,dg2,dZD,dZU,ctcplSucSdVWm(gt1,gt2))

 End Do 
End Do 


ctcplSucSuVZ = 0._dp 
Do gt1 = 1, 6
 Do gt2 = 1, 6
Call CT_CouplingSucSuVZ(gt1,gt2,g1,g2,ZU,TW,dg1,dg2,dZU,dSinTW,dCosTW,dTanTW,         & 
& ctcplSucSuVZ(gt1,gt2))

 End Do 
End Do 


ctcplSvcSeVWm = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 6
Call CT_CouplingSvcSeVWm(gt1,gt2,g2,ZV,ZE,dg2,dZV,dZE,ctcplSvcSeVWm(gt1,gt2))

 End Do 
End Do 


ctcplSvcSvVZ = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingSvcSvVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplSvcSvVZ(gt1,gt2))

 End Do 
End Do 


ctcplhhcVWmVWm = 0._dp 
Do gt1 = 1, 2
Call CT_CouplinghhcVWmVWm(gt1,g2,vd,vu,ZH,dg2,dvd,dvu,dZH,ctcplhhcVWmVWm(gt1))

End Do 


ctcplhhVZVZ = 0._dp 
Do gt1 = 1, 2
Call CT_CouplinghhVZVZ(gt1,g1,g2,vd,vu,ZH,TW,dg1,dg2,dvd,dvu,dZH,dSinTW,              & 
& dCosTW,dTanTW,ctcplhhVZVZ(gt1))

End Do 


ctcplHpmcVWmVP = 0._dp 
Do gt1 = 1, 2
Call CT_CouplingHpmcVWmVP(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,dSinTW,           & 
& dCosTW,dTanTW,ctcplHpmcVWmVP(gt1))

End Do 


ctcplHpmcVWmVZ = 0._dp 
Do gt1 = 1, 2
Call CT_CouplingHpmcVWmVZ(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,dSinTW,           & 
& dCosTW,dTanTW,ctcplHpmcVWmVZ(gt1))

End Do 


ctcplcHpmVPVWm = 0._dp 
Do gt1 = 1, 2
Call CT_CouplingcHpmVPVWm(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,dSinTW,           & 
& dCosTW,dTanTW,ctcplcHpmVPVWm(gt1))

End Do 


ctcplcHpmVWmVZ = 0._dp 
Do gt1 = 1, 2
Call CT_CouplingcHpmVWmVZ(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,dSinTW,           & 
& dCosTW,dTanTW,ctcplcHpmVWmVZ(gt1))

End Do 


ctcplVGVGVG = 0._dp 
Call CT_CouplingVGVGVG(g3,dg3,ctcplVGVGVG)



ctcplcVWmVPVWm = 0._dp 
Call CT_CouplingcVWmVPVWm(g2,TW,dg2,dSinTW,dCosTW,dTanTW,ctcplcVWmVPVWm)



ctcplcVWmVWmVZ = 0._dp 
Call CT_CouplingcVWmVWmVZ(g2,TW,dg2,dSinTW,dCosTW,dTanTW,ctcplcVWmVWmVZ)



ctcplcChaChaAhL = 0._dp 
ctcplcChaChaAhR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
  Do gt3 = 1, 2
Call CT_CouplingcChaChaAh(gt1,gt2,gt3,g2,ZA,UM,UP,dg2,dZA,dUM,dUP,ctcplcChaChaAhL(gt1,gt2,gt3)& 
& ,ctcplcChaChaAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiChiAhL = 0._dp 
ctcplChiChiAhR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 4
  Do gt3 = 1, 2
Call CT_CouplingChiChiAh(gt1,gt2,gt3,g1,g2,ZA,ZN,dg1,dg2,dZA,dZN,ctcplChiChiAhL(gt1,gt2,gt3)& 
& ,ctcplChiChiAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdFdAhL = 0._dp 
ctcplcFdFdAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFdFdAh(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,dYd,dZA,dZDL,dZDR,ctcplcFdFdAhL(gt1,gt2,gt3)& 
& ,ctcplcFdFdAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeFeAhL = 0._dp 
ctcplcFeFeAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFeFeAh(gt1,gt2,gt3,Ye,ZA,ZEL,ZER,dYe,dZA,dZEL,dZER,ctcplcFeFeAhL(gt1,gt2,gt3)& 
& ,ctcplcFeFeAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFuFuAhL = 0._dp 
ctcplcFuFuAhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFuFuAh(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,dYu,dZA,dZUL,dZUR,ctcplcFuFuAhL(gt1,gt2,gt3)& 
& ,ctcplcFuFuAhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiChacHpmL = 0._dp 
ctcplChiChacHpmR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 2
  Do gt3 = 1, 2
Call CT_CouplingChiChacHpm(gt1,gt2,gt3,g1,g2,ZP,ZN,UM,UP,dg1,dg2,dZP,dZN,             & 
& dUM,dUP,ctcplChiChacHpmL(gt1,gt2,gt3),ctcplChiChacHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChaFucSdL = 0._dp 
ctcplChaFucSdR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingChaFucSd(gt1,gt2,gt3,g2,Yd,Yu,ZD,UM,UP,ZUL,ZUR,dg2,dYd,               & 
& dYu,dZD,dUM,dUP,dZUL,dZUR,ctcplChaFucSdL(gt1,gt2,gt3),ctcplChaFucSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplFvChacSeL = 0._dp 
ctcplFvChacSeR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 2
  Do gt3 = 1, 6
Call CT_CouplingFvChacSe(gt1,gt2,gt3,g2,Ye,ZE,UM,dg2,dYe,dZE,dUM,ctcplFvChacSeL(gt1,gt2,gt3)& 
& ,ctcplFvChacSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcChaChahhL = 0._dp 
ctcplcChaChahhR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
  Do gt3 = 1, 2
Call CT_CouplingcChaChahh(gt1,gt2,gt3,g2,ZH,UM,UP,dg2,dZH,dUM,dUP,ctcplcChaChahhL(gt1,gt2,gt3)& 
& ,ctcplcChaChahhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdChaSuL = 0._dp 
ctcplcFdChaSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 2
  Do gt3 = 1, 6
Call CT_CouplingcFdChaSu(gt1,gt2,gt3,g2,Yd,Yu,ZU,UM,UP,ZDL,ZDR,dg2,dYd,               & 
& dYu,dZU,dUM,dUP,dZDL,dZDR,ctcplcFdChaSuL(gt1,gt2,gt3),ctcplcFdChaSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeChaSvL = 0._dp 
ctcplcFeChaSvR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 2
  Do gt3 = 1, 3
Call CT_CouplingcFeChaSv(gt1,gt2,gt3,g2,Ye,ZV,UM,UP,ZEL,ZER,dg2,dYe,dZV,              & 
& dUM,dUP,dZEL,dZER,ctcplcFeChaSvL(gt1,gt2,gt3),ctcplcFeChaSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiChihhL = 0._dp 
ctcplChiChihhR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 4
  Do gt3 = 1, 2
Call CT_CouplingChiChihh(gt1,gt2,gt3,g1,g2,ZH,ZN,dg1,dg2,dZH,dZN,ctcplChiChihhL(gt1,gt2,gt3)& 
& ,ctcplChiChihhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiFdcSdL = 0._dp 
ctcplChiFdcSdR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingChiFdcSd(gt1,gt2,gt3,g1,g2,Yd,ZD,ZN,ZDL,ZDR,dg1,dg2,dYd,              & 
& dZD,dZN,dZDL,dZDR,ctcplChiFdcSdL(gt1,gt2,gt3),ctcplChiFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiFecSeL = 0._dp 
ctcplChiFecSeR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingChiFecSe(gt1,gt2,gt3,g1,g2,Ye,ZE,ZN,ZEL,ZER,dg1,dg2,dYe,              & 
& dZE,dZN,dZEL,dZER,ctcplChiFecSeL(gt1,gt2,gt3),ctcplChiFecSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiFucSuL = 0._dp 
ctcplChiFucSuR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingChiFucSu(gt1,gt2,gt3,g1,g2,Yu,ZU,ZN,ZUL,ZUR,dg1,dg2,dYu,              & 
& dZU,dZN,dZUL,dZUR,ctcplChiFucSuL(gt1,gt2,gt3),ctcplChiFucSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiFvSvL = 0._dp 
ctcplChiFvSvR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplingChiFvSv(gt1,gt2,gt3,g1,g2,ZV,ZN,dg1,dg2,dZV,dZN,ctcplChiFvSvL(gt1,gt2,gt3)& 
& ,ctcplChiFvSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiFvcSvL = 0._dp 
ctcplChiFvcSvR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplingChiFvcSv(gt1,gt2,gt3,g1,g2,ZV,ZN,dg1,dg2,dZV,dZN,ctcplChiFvcSvL(gt1,gt2,gt3)& 
& ,ctcplChiFvcSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcChaChiHpmL = 0._dp 
ctcplcChaChiHpmR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 4
  Do gt3 = 1, 2
Call CT_CouplingcChaChiHpm(gt1,gt2,gt3,g1,g2,ZP,ZN,UM,UP,dg1,dg2,dZP,dZN,             & 
& dUM,dUP,ctcplcChaChiHpmL(gt1,gt2,gt3),ctcplcChaChiHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdChiSdL = 0._dp 
ctcplcFdChiSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 4
  Do gt3 = 1, 6
Call CT_CouplingcFdChiSd(gt1,gt2,gt3,g1,g2,Yd,ZD,ZN,ZDL,ZDR,dg1,dg2,dYd,              & 
& dZD,dZN,dZDL,dZDR,ctcplcFdChiSdL(gt1,gt2,gt3),ctcplcFdChiSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeChiSeL = 0._dp 
ctcplcFeChiSeR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 4
  Do gt3 = 1, 6
Call CT_CouplingcFeChiSe(gt1,gt2,gt3,g1,g2,Ye,ZE,ZN,ZEL,ZER,dg1,dg2,dYe,              & 
& dZE,dZN,dZEL,dZER,ctcplcFeChiSeL(gt1,gt2,gt3),ctcplcFeChiSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFuChiSuL = 0._dp 
ctcplcFuChiSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 4
  Do gt3 = 1, 6
Call CT_CouplingcFuChiSu(gt1,gt2,gt3,g1,g2,Yu,ZU,ZN,ZUL,ZUR,dg1,dg2,dYu,              & 
& dZU,dZN,dZUL,dZUR,ctcplcFuChiSuL(gt1,gt2,gt3),ctcplcFuChiSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplFvFdcSdL = 0._dp 
ctcplFvFdcSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingFvFdcSd(gt1,gt2,gt3,L2,ZD,ZDL,ZDR,dL2,dZD,dZDL,dZDR,ctcplFvFdcSdL(gt1,gt2,gt3)& 
& ,ctcplFvFdcSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplGluFdcSdL = 0._dp 
ctcplGluFdcSdR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CT_CouplingGluFdcSd(gt2,gt3,g3,pG,ZD,ZDL,ZDR,dg3,dpG,dZD,dZDL,dZDR,              & 
& ctcplGluFdcSdL(gt2,gt3),ctcplGluFdcSdR(gt2,gt3))

 End Do 
End Do 


ctcplcFdFdhhL = 0._dp 
ctcplcFdFdhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFdFdhh(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,dYd,dZH,dZDL,dZDR,ctcplcFdFdhhL(gt1,gt2,gt3)& 
& ,ctcplcFdFdhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdFdSvL = 0._dp 
ctcplcFdFdSvR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplingcFdFdSv(gt1,gt2,gt3,L2,ZV,ZDL,ZDR,dL2,dZV,dZDL,dZDR,ctcplcFdFdSvL(gt1,gt2,gt3)& 
& ,ctcplcFdFdSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcChaFdcSuL = 0._dp 
ctcplcChaFdcSuR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcChaFdcSu(gt1,gt2,gt3,g2,Yd,Yu,ZU,UM,UP,ZDL,ZDR,dg2,dYd,              & 
& dYu,dZU,dUM,dUP,dZDL,dZDR,ctcplcChaFdcSuL(gt1,gt2,gt3),ctcplcChaFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdFdcSvL = 0._dp 
ctcplcFdFdcSvR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplingcFdFdcSv(gt1,gt2,gt3,L2,ZV,ZDL,ZDR,dL2,dZV,dZDL,dZDR,ctcplcFdFdcSvL(gt1,gt2,gt3)& 
& ,ctcplcFdFdcSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeFdcSuL = 0._dp 
ctcplcFeFdcSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcFeFdcSu(gt1,gt2,gt3,L2,ZU,ZEL,ZDR,dL2,dZU,dZEL,dZDR,ctcplcFeFdcSuL(gt1,gt2,gt3)& 
& ,ctcplcFeFdcSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFuFdcHpmL = 0._dp 
ctcplcFuFdcHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFuFdcHpm(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,dYd,dYu,               & 
& dZP,dZDL,dZDR,dZUL,dZUR,ctcplcFuFdcHpmL(gt1,gt2,gt3),ctcplcFuFdcHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFuFdcSeL = 0._dp 
ctcplcFuFdcSeR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcFuFdcSe(gt1,gt2,gt3,L2,ZE,ZDR,ZUL,dL2,dZE,dZDR,dZUL,ctcplcFuFdcSeL(gt1,gt2,gt3)& 
& ,ctcplcFuFdcSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplFeFucSdL = 0._dp 
ctcplFeFucSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingFeFucSd(gt1,gt2,gt3,L2,ZD,ZEL,ZUL,dL2,dZD,dZEL,dZUL,ctcplFeFucSdL(gt1,gt2,gt3)& 
& ,ctcplFeFucSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplFvFecHpmL = 0._dp 
ctcplFvFecHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingFvFecHpm(gt1,gt2,gt3,Ye,ZP,ZER,dYe,dZP,dZER,ctcplFvFecHpmL(gt1,gt2,gt3)& 
& ,ctcplFvFecHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplFvFecSeL = 0._dp 
ctcplFvFecSeR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingFvFecSe(gt1,gt2,gt3,L1,ZE,ZEL,ZER,dL1,dZE,dZEL,dZER,ctcplFvFecSeL(gt1,gt2,gt3)& 
& ,ctcplFvFecSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeFehhL = 0._dp 
ctcplcFeFehhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFeFehh(gt1,gt2,gt3,Ye,ZH,ZEL,ZER,dYe,dZH,dZEL,dZER,ctcplcFeFehhL(gt1,gt2,gt3)& 
& ,ctcplcFeFehhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdFeSuL = 0._dp 
ctcplcFdFeSuR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcFdFeSu(gt1,gt2,gt3,L2,ZU,ZEL,ZDR,dL2,dZU,dZEL,dZDR,ctcplcFdFeSuL(gt1,gt2,gt3)& 
& ,ctcplcFdFeSuR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeFeSvL = 0._dp 
ctcplcFeFeSvR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplingcFeFeSv(gt1,gt2,gt3,L1,ZV,ZEL,ZER,dL1,dZV,dZEL,dZER,ctcplcFeFeSvL(gt1,gt2,gt3)& 
& ,ctcplcFeFeSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcChaFecSvL = 0._dp 
ctcplcChaFecSvR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplingcChaFecSv(gt1,gt2,gt3,g2,Ye,ZV,UM,UP,ZEL,ZER,dg2,dYe,dZV,             & 
& dUM,dUP,dZEL,dZER,ctcplcChaFecSvL(gt1,gt2,gt3),ctcplcChaFecSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeFecSvL = 0._dp 
ctcplcFeFecSvR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 3
Call CT_CouplingcFeFecSv(gt1,gt2,gt3,L1,ZV,ZEL,ZER,dL1,dZV,dZEL,dZER,ctcplcFeFecSvL(gt1,gt2,gt3)& 
& ,ctcplcFeFecSvR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplGluFucSuL = 0._dp 
ctcplGluFucSuR = 0._dp 
Do gt2 = 1, 3
 Do gt3 = 1, 6
Call CT_CouplingGluFucSu(gt2,gt3,g3,pG,ZU,ZUL,ZUR,dg3,dpG,dZU,dZUL,dZUR,              & 
& ctcplGluFucSuL(gt2,gt3),ctcplGluFucSuR(gt2,gt3))

 End Do 
End Do 


ctcplcFuFuhhL = 0._dp 
ctcplcFuFuhhR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFuFuhh(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,dYu,dZH,dZUL,dZUR,ctcplcFuFuhhL(gt1,gt2,gt3)& 
& ,ctcplcFuFuhhR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdFuHpmL = 0._dp 
ctcplcFdFuHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFdFuHpm(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,dYd,dYu,dZP,            & 
& dZDL,dZDR,dZUL,dZUR,ctcplcFdFuHpmL(gt1,gt2,gt3),ctcplcFdFuHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdFuSeL = 0._dp 
ctcplcFdFuSeR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcFdFuSe(gt1,gt2,gt3,L2,ZE,ZDR,ZUL,dL2,dZE,dZDR,dZUL,ctcplcFdFuSeL(gt1,gt2,gt3)& 
& ,ctcplcFdFuSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeFvHpmL = 0._dp 
ctcplcFeFvHpmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 2
Call CT_CouplingcFeFvHpm(gt1,gt2,gt3,Ye,ZP,ZER,dYe,dZP,dZER,ctcplcFeFvHpmL(gt1,gt2,gt3)& 
& ,ctcplcFeFvHpmR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdFvSdL = 0._dp 
ctcplcFdFvSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcFdFvSd(gt1,gt2,gt3,L2,ZD,ZDL,ZDR,dL2,dZD,dZDL,dZDR,ctcplcFdFvSdL(gt1,gt2,gt3)& 
& ,ctcplcFdFvSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcChaFvSeL = 0._dp 
ctcplcChaFvSeR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcChaFvSe(gt1,gt2,gt3,g2,Ye,ZE,UM,dg2,dYe,dZE,dUM,ctcplcChaFvSeL(gt1,gt2,gt3)& 
& ,ctcplcChaFvSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFeFvSeL = 0._dp 
ctcplcFeFvSeR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcFeFvSe(gt1,gt2,gt3,L1,ZE,ZEL,ZER,dL1,dZE,dZEL,dZER,ctcplcFeFvSeL(gt1,gt2,gt3)& 
& ,ctcplcFeFvSeR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFdGluSdL = 0._dp 
ctcplcFdGluSdR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CT_CouplingcFdGluSd(gt1,gt3,g3,pG,ZD,ZDL,ZDR,dg3,dpG,dZD,dZDL,dZDR,              & 
& ctcplcFdGluSdL(gt1,gt3),ctcplcFdGluSdR(gt1,gt3))

 End Do 
End Do 


ctcplcFuGluSuL = 0._dp 
ctcplcFuGluSuR = 0._dp 
Do gt1 = 1, 3
 Do gt3 = 1, 6
Call CT_CouplingcFuGluSu(gt1,gt3,g3,pG,ZU,ZUL,ZUR,dg3,dpG,dZU,dZUL,dZUR,              & 
& ctcplcFuGluSuL(gt1,gt3),ctcplcFuGluSuR(gt1,gt3))

 End Do 
End Do 


ctcplcChacFuSdL = 0._dp 
ctcplcChacFuSdR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcChacFuSd(gt1,gt2,gt3,g2,Yd,Yu,ZD,UM,UP,ZUL,ZUR,dg2,dYd,              & 
& dYu,dZD,dUM,dUP,dZUL,dZUR,ctcplcChacFuSdL(gt1,gt2,gt3),ctcplcChacFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplcFecFuSdL = 0._dp 
ctcplcFecFuSdR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
  Do gt3 = 1, 6
Call CT_CouplingcFecFuSd(gt1,gt2,gt3,L2,ZD,ZEL,ZUL,dL2,dZD,dZEL,dZUL,ctcplcFecFuSdL(gt1,gt2,gt3)& 
& ,ctcplcFecFuSdR(gt1,gt2,gt3))

  End Do 
 End Do 
End Do 


ctcplChiChacVWmL = 0._dp 
ctcplChiChacVWmR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 2
Call CT_CouplingChiChacVWm(gt1,gt2,g2,ZN,UM,UP,dg2,dZN,dUM,dUP,ctcplChiChacVWmL(gt1,gt2)& 
& ,ctcplChiChacVWmR(gt1,gt2))

 End Do 
End Do 


ctcplcChaChaVPL = 0._dp 
ctcplcChaChaVPR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplingcChaChaVP(gt1,gt2,g1,g2,UM,UP,TW,dg1,dg2,dUM,dUP,dSinTW,              & 
& dCosTW,dTanTW,ctcplcChaChaVPL(gt1,gt2),ctcplcChaChaVPR(gt1,gt2))

 End Do 
End Do 


ctcplcChaChaVZL = 0._dp 
ctcplcChaChaVZR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 2
Call CT_CouplingcChaChaVZ(gt1,gt2,g1,g2,UM,UP,TW,dg1,dg2,dUM,dUP,dSinTW,              & 
& dCosTW,dTanTW,ctcplcChaChaVZL(gt1,gt2),ctcplcChaChaVZR(gt1,gt2))

 End Do 
End Do 


ctcplChiChiVZL = 0._dp 
ctcplChiChiVZR = 0._dp 
Do gt1 = 1, 4
 Do gt2 = 1, 4
Call CT_CouplingChiChiVZ(gt1,gt2,g1,g2,ZN,TW,dg1,dg2,dZN,dSinTW,dCosTW,               & 
& dTanTW,ctcplChiChiVZL(gt1,gt2),ctcplChiChiVZR(gt1,gt2))

 End Do 
End Do 


ctcplcChaChiVWmL = 0._dp 
ctcplcChaChiVWmR = 0._dp 
Do gt1 = 1, 2
 Do gt2 = 1, 4
Call CT_CouplingcChaChiVWm(gt1,gt2,g2,ZN,UM,UP,dg2,dZN,dUM,dUP,ctcplcChaChiVWmL(gt1,gt2)& 
& ,ctcplcChaChiVWmR(gt1,gt2))

 End Do 
End Do 


ctcplcFdFdVGL = 0._dp 
ctcplcFdFdVGR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFdFdVG(gt1,gt2,g3,dg3,ctcplcFdFdVGL(gt1,gt2),ctcplcFdFdVGR(gt1,gt2))

 End Do 
End Do 


ctcplcFdFdVPL = 0._dp 
ctcplcFdFdVPR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFdFdVP(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplcFdFdVPL(gt1,gt2)& 
& ,ctcplcFdFdVPR(gt1,gt2))

 End Do 
End Do 


ctcplcFdFdVZL = 0._dp 
ctcplcFdFdVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFdFdVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplcFdFdVZL(gt1,gt2)& 
& ,ctcplcFdFdVZR(gt1,gt2))

 End Do 
End Do 


ctcplcFuFdcVWmL = 0._dp 
ctcplcFuFdcVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFuFdcVWm(gt1,gt2,g2,ZDL,ZUL,dg2,dZDL,dZUL,ctcplcFuFdcVWmL(gt1,gt2)   & 
& ,ctcplcFuFdcVWmR(gt1,gt2))

 End Do 
End Do 


ctcplFvFecVWmL = 0._dp 
ctcplFvFecVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingFvFecVWm(gt1,gt2,g2,ZEL,dg2,dZEL,ctcplFvFecVWmL(gt1,gt2)              & 
& ,ctcplFvFecVWmR(gt1,gt2))

 End Do 
End Do 


ctcplcFeFeVPL = 0._dp 
ctcplcFeFeVPR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFeFeVP(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplcFeFeVPL(gt1,gt2)& 
& ,ctcplcFeFeVPR(gt1,gt2))

 End Do 
End Do 


ctcplcFeFeVZL = 0._dp 
ctcplcFeFeVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFeFeVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplcFeFeVZL(gt1,gt2)& 
& ,ctcplcFeFeVZR(gt1,gt2))

 End Do 
End Do 


ctcplcFuFuVGL = 0._dp 
ctcplcFuFuVGR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFuFuVG(gt1,gt2,g3,dg3,ctcplcFuFuVGL(gt1,gt2),ctcplcFuFuVGR(gt1,gt2))

 End Do 
End Do 


ctcplcFuFuVPL = 0._dp 
ctcplcFuFuVPR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFuFuVP(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplcFuFuVPL(gt1,gt2)& 
& ,ctcplcFuFuVPR(gt1,gt2))

 End Do 
End Do 


ctcplcFdFuVWmL = 0._dp 
ctcplcFdFuVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFdFuVWm(gt1,gt2,g2,ZDL,ZUL,dg2,dZDL,dZUL,ctcplcFdFuVWmL(gt1,gt2)     & 
& ,ctcplcFdFuVWmR(gt1,gt2))

 End Do 
End Do 


ctcplcFuFuVZL = 0._dp 
ctcplcFuFuVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFuFuVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplcFuFuVZL(gt1,gt2)& 
& ,ctcplcFuFuVZR(gt1,gt2))

 End Do 
End Do 


ctcplFvFvVZL = 0._dp 
ctcplFvFvVZR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingFvFvVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,ctcplFvFvVZL(gt1,gt2)& 
& ,ctcplFvFvVZR(gt1,gt2))

 End Do 
End Do 


ctcplcFeFvVWmL = 0._dp 
ctcplcFeFvVWmR = 0._dp 
Do gt1 = 1, 3
 Do gt2 = 1, 3
Call CT_CouplingcFeFvVWm(gt1,gt2,g2,ZEL,dg2,dZEL,ctcplcFeFvVWmL(gt1,gt2)              & 
& ,ctcplcFeFvVWmR(gt1,gt2))

 End Do 
End Do 


ctcplGluGluVGL = 0._dp 
ctcplGluGluVGR = 0._dp 
Call CT_CouplingGluGluVG(g3,pG,dg3,dpG,ctcplGluGluVGL,ctcplGluGluVGR)



Iname = Iname - 1 
End Subroutine CalculateCouplingCT

Subroutine CT_CouplingAhAhhh(gt1,gt2,gt3,g1,g2,vd,vu,ZH,ZA,dg1,dg2,dvd,               & 
& dvu,dZH,dZA,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),ZA(2,2),dg1,dg2,dvd,dvu,dZH(2,2),dZA(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhAhhh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g1**2*vd*dZH(gt3,1)*ZA(gt1,1)*ZA(gt2,1))/4._dp
res = res-(g2**2*vd*dZH(gt3,1)*ZA(gt1,1)*ZA(gt2,1))/4._dp
res = res+(g1**2*vu*dZH(gt3,2)*ZA(gt1,1)*ZA(gt2,1))/4._dp
res = res+(g2**2*vu*dZH(gt3,2)*ZA(gt1,1)*ZA(gt2,1))/4._dp
res = res+(g1**2*vd*dZH(gt3,1)*ZA(gt1,2)*ZA(gt2,2))/4._dp
res = res+(g2**2*vd*dZH(gt3,1)*ZA(gt1,2)*ZA(gt2,2))/4._dp
res = res-(g1**2*vu*dZH(gt3,2)*ZA(gt1,2)*ZA(gt2,2))/4._dp
res = res-(g2**2*vu*dZH(gt3,2)*ZA(gt1,2)*ZA(gt2,2))/4._dp
res = res-(g1**2*vd*dZA(gt2,1)*ZA(gt1,1)*ZH(gt3,1))/4._dp
res = res-(g2**2*vd*dZA(gt2,1)*ZA(gt1,1)*ZH(gt3,1))/4._dp
res = res+(g1**2*vd*dZA(gt2,2)*ZA(gt1,2)*ZH(gt3,1))/4._dp
res = res+(g2**2*vd*dZA(gt2,2)*ZA(gt1,2)*ZH(gt3,1))/4._dp
res = res-(g1**2*vd*dZA(gt1,1)*ZA(gt2,1)*ZH(gt3,1))/4._dp
res = res-(g2**2*vd*dZA(gt1,1)*ZA(gt2,1)*ZH(gt3,1))/4._dp
res = res-(dvd*g1**2*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,1))/4._dp
res = res-(dvd*g2**2*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,1))/4._dp
res = res-(dg1*g1*vd*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,1))/2._dp
res = res-(dg2*g2*vd*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,1))/2._dp
res = res+(g1**2*vd*dZA(gt1,2)*ZA(gt2,2)*ZH(gt3,1))/4._dp
res = res+(g2**2*vd*dZA(gt1,2)*ZA(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dvd*g1**2*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dvd*g2**2*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dg1*g1*vd*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,1))/2._dp
res = res+(dg2*g2*vd*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,1))/2._dp
res = res+(g1**2*vu*dZA(gt2,1)*ZA(gt1,1)*ZH(gt3,2))/4._dp
res = res+(g2**2*vu*dZA(gt2,1)*ZA(gt1,1)*ZH(gt3,2))/4._dp
res = res-(g1**2*vu*dZA(gt2,2)*ZA(gt1,2)*ZH(gt3,2))/4._dp
res = res-(g2**2*vu*dZA(gt2,2)*ZA(gt1,2)*ZH(gt3,2))/4._dp
res = res+(g1**2*vu*dZA(gt1,1)*ZA(gt2,1)*ZH(gt3,2))/4._dp
res = res+(g2**2*vu*dZA(gt1,1)*ZA(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dvu*g1**2*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dvu*g2**2*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dg1*g1*vu*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,2))/2._dp
res = res+(dg2*g2*vu*ZA(gt1,1)*ZA(gt2,1)*ZH(gt3,2))/2._dp
res = res-(g1**2*vu*dZA(gt1,2)*ZA(gt2,2)*ZH(gt3,2))/4._dp
res = res-(g2**2*vu*dZA(gt1,2)*ZA(gt2,2)*ZH(gt3,2))/4._dp
res = res-(dvu*g1**2*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,2))/4._dp
res = res-(dvu*g2**2*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,2))/4._dp
res = res-(dg1*g1*vu*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,2))/2._dp
res = res-(dg2*g2*vu*ZA(gt1,2)*ZA(gt2,2)*ZH(gt3,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhAhhh  
 
 
Subroutine CT_CouplingAhHpmcHpm(gt1,gt2,gt3,g2,vd,vu,ZA,ZP,dg2,dvd,dvu,               & 
& dZA,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,vd,vu,ZA(2,2),ZP(2,2),dg2,dvd,dvu,dZA(2,2),dZP(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhHpmcHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g2**2*vu*dZP(gt3,2)*ZA(gt1,1)*ZP(gt2,1))/4._dp
res = res-(g2**2*vd*dZP(gt3,2)*ZA(gt1,2)*ZP(gt2,1))/4._dp
res = res+(g2**2*vu*dZP(gt3,1)*ZA(gt1,1)*ZP(gt2,2))/4._dp
res = res+(g2**2*vd*dZP(gt3,1)*ZA(gt1,2)*ZP(gt2,2))/4._dp
res = res+(g2**2*vu*dZP(gt2,2)*ZA(gt1,1)*ZP(gt3,1))/4._dp
res = res+(g2**2*vd*dZP(gt2,2)*ZA(gt1,2)*ZP(gt3,1))/4._dp
res = res+(g2**2*vu*dZA(gt1,1)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res+(g2**2*vd*dZA(gt1,2)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res+(dvu*g2**2*ZA(gt1,1)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res+(dg2*g2*vu*ZA(gt1,1)*ZP(gt2,2)*ZP(gt3,1))/2._dp
res = res+(dvd*g2**2*ZA(gt1,2)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res+(dg2*g2*vd*ZA(gt1,2)*ZP(gt2,2)*ZP(gt3,1))/2._dp
res = res-(g2**2*vu*dZP(gt2,1)*ZA(gt1,1)*ZP(gt3,2))/4._dp
res = res-(g2**2*vd*dZP(gt2,1)*ZA(gt1,2)*ZP(gt3,2))/4._dp
res = res-(g2**2*vu*dZA(gt1,1)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(g2**2*vd*dZA(gt1,2)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(dvu*g2**2*ZA(gt1,1)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(dg2*g2*vu*ZA(gt1,1)*ZP(gt2,1)*ZP(gt3,2))/2._dp
res = res-(dvd*g2**2*ZA(gt1,2)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(dg2*g2*vd*ZA(gt1,2)*ZP(gt2,1)*ZP(gt3,2))/2._dp
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhHpmcHpm  
 
 
Subroutine CT_CouplingAhSdcSd(gt1,gt2,gt3,Mu,Yd,Td,ZD,ZA,dMu,dYd,dTd,dZD,dZA,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZA(2,2),dZA(2,2)

Complex(dp), Intent(in) :: Mu,Yd(3,3),Td(3,3),ZD(6,6),dMu,dYd(3,3),dTd(3,3),dZD(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhSdcSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZD(gt2,3 + j1))*Conjg(Td(j1,j2))*dZD(gt3,j2)*ZA(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(Yd(j1,j2))*Conjg(ZD(gt2,3 + j1))*dZD(gt3,j2)*ZA(gt1,2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZD(gt2,j2))*dZD(gt3,3 + j1)*Yd(j1,j2)*ZA(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZD(gt2,j2))*dZA(gt1,2)*Yd(j1,j2)*ZD(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZD(gt2,j2))*dTd(j1,j2)*ZA(gt1,1)*ZD(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZD(gt2,j2))*dYd(j1,j2)*ZA(gt1,2)*ZD(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(dZD(gt2,j2))*Yd(j1,j2)*ZA(gt1,2)*ZD(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dMu)*Conjg(ZD(gt2,j2))*Yd(j1,j2)*ZA(gt1,2)*ZD(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZD(gt2,3 + j1))*Conjg(Td(j1,j2))*dZA(gt1,1)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(Yd(j1,j2))*Conjg(ZD(gt2,3 + j1))*dZA(gt1,2)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dTd(j1,j2))*Conjg(ZD(gt2,3 + j1))*ZA(gt1,1)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZD(gt2,3 + j1))*Conjg(Td(j1,j2))*ZA(gt1,1)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(dZD(gt2,3 + j1))*Conjg(Yd(j1,j2))*ZA(gt1,2)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(dYd(j1,j2))*Conjg(ZD(gt2,3 + j1))*ZA(gt1,2)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((dMu*Conjg(Yd(j1,j2))*Conjg(ZD(gt2,3 + j1))*ZA(gt1,2)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZD(gt2,j2))*dZD(gt3,3 + j1)*ZA(gt1,1)*Td(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZD(gt2,j2))*dZA(gt1,1)*ZD(gt3,3 + j1)*Td(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dZD(gt2,j2))*ZA(gt1,1)*ZD(gt3,3 + j1)*Td(j1,j2))/sqrt(2._dp)
End Do 
End Do 
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhSdcSd  
 
 
Subroutine CT_CouplingAhSecSe(gt1,gt2,gt3,Mu,Ye,Te,ZE,ZA,dMu,dYe,dTe,dZE,dZA,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZA(2,2),dZA(2,2)

Complex(dp), Intent(in) :: Mu,Ye(3,3),Te(3,3),ZE(6,6),dMu,dYe(3,3),dTe(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhSecSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZE(gt2,3 + j1))*Conjg(Te(j1,j2))*dZE(gt3,j2)*ZA(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(Ye(j1,j2))*Conjg(ZE(gt2,3 + j1))*dZE(gt3,j2)*ZA(gt1,2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZE(gt2,j2))*dZE(gt3,3 + j1)*Ye(j1,j2)*ZA(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZE(gt2,j2))*dZA(gt1,2)*Ye(j1,j2)*ZE(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZE(gt2,j2))*dTe(j1,j2)*ZA(gt1,1)*ZE(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZE(gt2,j2))*dYe(j1,j2)*ZA(gt1,2)*ZE(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(dZE(gt2,j2))*Ye(j1,j2)*ZA(gt1,2)*ZE(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dMu)*Conjg(ZE(gt2,j2))*Ye(j1,j2)*ZA(gt1,2)*ZE(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZE(gt2,3 + j1))*Conjg(Te(j1,j2))*dZA(gt1,1)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(Ye(j1,j2))*Conjg(ZE(gt2,3 + j1))*dZA(gt1,2)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dTe(j1,j2))*Conjg(ZE(gt2,3 + j1))*ZA(gt1,1)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZE(gt2,3 + j1))*Conjg(Te(j1,j2))*ZA(gt1,1)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(dZE(gt2,3 + j1))*Conjg(Ye(j1,j2))*ZA(gt1,2)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(dYe(j1,j2))*Conjg(ZE(gt2,3 + j1))*ZA(gt1,2)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((dMu*Conjg(Ye(j1,j2))*Conjg(ZE(gt2,3 + j1))*ZA(gt1,2)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZE(gt2,j2))*dZE(gt3,3 + j1)*ZA(gt1,1)*Te(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZE(gt2,j2))*dZA(gt1,1)*ZE(gt3,3 + j1)*Te(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dZE(gt2,j2))*ZA(gt1,1)*ZE(gt3,3 + j1)*Te(j1,j2))/sqrt(2._dp)
End Do 
End Do 
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhSecSe  
 
 
Subroutine CT_CouplingAhSucSu(gt1,gt2,gt3,Mu,Yu,Tu,ZU,ZA,dMu,dYu,dTu,dZU,dZA,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZA(2,2),dZA(2,2)

Complex(dp), Intent(in) :: Mu,Yu(3,3),Tu(3,3),ZU(6,6),dMu,dYu(3,3),dTu(3,3),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhSucSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*dZU(gt3,j2)*ZA(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZU(gt2,j2))*dZU(gt3,3 + j1)*Yu(j1,j2)*ZA(gt1,1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*dZU(gt3,j2)*ZA(gt1,2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZU(gt2,j2))*dZA(gt1,1)*Yu(j1,j2)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZU(gt2,j2))*dYu(j1,j2)*ZA(gt1,1)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(dZU(gt2,j2))*Yu(j1,j2)*ZA(gt1,1)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dMu)*Conjg(ZU(gt2,j2))*Yu(j1,j2)*ZA(gt1,1)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZU(gt2,j2))*dTu(j1,j2)*ZA(gt1,2)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*dZA(gt1,1)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*dZA(gt1,2)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(dZU(gt2,3 + j1))*Conjg(Yu(j1,j2))*ZA(gt1,1)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Mu*Conjg(dYu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZA(gt1,1)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((dMu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZA(gt1,1)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dTu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZA(gt1,2)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*ZA(gt1,2)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZU(gt2,j2))*dZU(gt3,3 + j1)*ZA(gt1,2)*Tu(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZU(gt2,j2))*dZA(gt1,2)*ZU(gt3,3 + j1)*Tu(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dZU(gt2,j2))*ZA(gt1,2)*ZU(gt3,3 + j1)*Tu(j1,j2))/sqrt(2._dp)
End Do 
End Do 
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhSucSu  
 
 
Subroutine CT_Couplinghhhhhh(gt1,gt2,gt3,g1,g2,vd,vu,ZH,dg1,dg2,dvd,dvu,dZH,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),dg1,dg2,dvd,dvu,dZH(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_Couplinghhhhhh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res+(-3*g1**2*vd*dZH(gt3,1)*ZH(gt1,1)*ZH(gt2,1))/4._dp
res = res+(-3*g2**2*vd*dZH(gt3,1)*ZH(gt1,1)*ZH(gt2,1))/4._dp
res = res+(g1**2*vu*dZH(gt3,2)*ZH(gt1,1)*ZH(gt2,1))/4._dp
res = res+(g2**2*vu*dZH(gt3,2)*ZH(gt1,1)*ZH(gt2,1))/4._dp
res = res+(g1**2*vu*dZH(gt3,1)*ZH(gt1,2)*ZH(gt2,1))/4._dp
res = res+(g2**2*vu*dZH(gt3,1)*ZH(gt1,2)*ZH(gt2,1))/4._dp
res = res+(g1**2*vd*dZH(gt3,2)*ZH(gt1,2)*ZH(gt2,1))/4._dp
res = res+(g2**2*vd*dZH(gt3,2)*ZH(gt1,2)*ZH(gt2,1))/4._dp
res = res+(g1**2*vu*dZH(gt3,1)*ZH(gt1,1)*ZH(gt2,2))/4._dp
res = res+(g2**2*vu*dZH(gt3,1)*ZH(gt1,1)*ZH(gt2,2))/4._dp
res = res+(g1**2*vd*dZH(gt3,2)*ZH(gt1,1)*ZH(gt2,2))/4._dp
res = res+(g2**2*vd*dZH(gt3,2)*ZH(gt1,1)*ZH(gt2,2))/4._dp
res = res+(g1**2*vd*dZH(gt3,1)*ZH(gt1,2)*ZH(gt2,2))/4._dp
res = res+(g2**2*vd*dZH(gt3,1)*ZH(gt1,2)*ZH(gt2,2))/4._dp
res = res+(-3*g1**2*vu*dZH(gt3,2)*ZH(gt1,2)*ZH(gt2,2))/4._dp
res = res+(-3*g2**2*vu*dZH(gt3,2)*ZH(gt1,2)*ZH(gt2,2))/4._dp
res = res+(-3*g1**2*vd*dZH(gt2,1)*ZH(gt1,1)*ZH(gt3,1))/4._dp
res = res+(-3*g2**2*vd*dZH(gt2,1)*ZH(gt1,1)*ZH(gt3,1))/4._dp
res = res+(g1**2*vu*dZH(gt2,2)*ZH(gt1,1)*ZH(gt3,1))/4._dp
res = res+(g2**2*vu*dZH(gt2,2)*ZH(gt1,1)*ZH(gt3,1))/4._dp
res = res+(g1**2*vu*dZH(gt2,1)*ZH(gt1,2)*ZH(gt3,1))/4._dp
res = res+(g2**2*vu*dZH(gt2,1)*ZH(gt1,2)*ZH(gt3,1))/4._dp
res = res+(g1**2*vd*dZH(gt2,2)*ZH(gt1,2)*ZH(gt3,1))/4._dp
res = res+(g2**2*vd*dZH(gt2,2)*ZH(gt1,2)*ZH(gt3,1))/4._dp
res = res+(-3*g1**2*vd*dZH(gt1,1)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(-3*g2**2*vd*dZH(gt1,1)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(g1**2*vu*dZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(g2**2*vu*dZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(-3*dvd*g1**2*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(-3*dvd*g2**2*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(-3*dg1*g1*vd*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,1))/2._dp
res = res+(-3*dg2*g2*vd*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,1))/2._dp
res = res+(dvu*g1**2*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(dvu*g2**2*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1))/4._dp
res = res+(dg1*g1*vu*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1))/2._dp
res = res+(dg2*g2*vu*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,1))/2._dp
res = res+(g1**2*vu*dZH(gt1,1)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(g2**2*vu*dZH(gt1,1)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(g1**2*vd*dZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(g2**2*vd*dZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dvu*g1**2*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dvu*g2**2*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dg1*g1*vu*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,1))/2._dp
res = res+(dg2*g2*vu*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,1))/2._dp
res = res+(dvd*g1**2*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dvd*g2**2*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1))/4._dp
res = res+(dg1*g1*vd*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1))/2._dp
res = res+(dg2*g2*vd*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,1))/2._dp
res = res+(g1**2*vu*dZH(gt2,1)*ZH(gt1,1)*ZH(gt3,2))/4._dp
res = res+(g2**2*vu*dZH(gt2,1)*ZH(gt1,1)*ZH(gt3,2))/4._dp
res = res+(g1**2*vd*dZH(gt2,2)*ZH(gt1,1)*ZH(gt3,2))/4._dp
res = res+(g2**2*vd*dZH(gt2,2)*ZH(gt1,1)*ZH(gt3,2))/4._dp
res = res+(g1**2*vd*dZH(gt2,1)*ZH(gt1,2)*ZH(gt3,2))/4._dp
res = res+(g2**2*vd*dZH(gt2,1)*ZH(gt1,2)*ZH(gt3,2))/4._dp
res = res+(-3*g1**2*vu*dZH(gt2,2)*ZH(gt1,2)*ZH(gt3,2))/4._dp
res = res+(-3*g2**2*vu*dZH(gt2,2)*ZH(gt1,2)*ZH(gt3,2))/4._dp
res = res+(g1**2*vu*dZH(gt1,1)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(g2**2*vu*dZH(gt1,1)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(g1**2*vd*dZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(g2**2*vd*dZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dvu*g1**2*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dvu*g2**2*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dg1*g1*vu*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,2))/2._dp
res = res+(dg2*g2*vu*ZH(gt1,1)*ZH(gt2,1)*ZH(gt3,2))/2._dp
res = res+(dvd*g1**2*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dvd*g2**2*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2))/4._dp
res = res+(dg1*g1*vd*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2))/2._dp
res = res+(dg2*g2*vd*ZH(gt1,2)*ZH(gt2,1)*ZH(gt3,2))/2._dp
res = res+(g1**2*vd*dZH(gt1,1)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(g2**2*vd*dZH(gt1,1)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(-3*g1**2*vu*dZH(gt1,2)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(-3*g2**2*vu*dZH(gt1,2)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(dvd*g1**2*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(dvd*g2**2*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(dg1*g1*vd*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,2))/2._dp
res = res+(dg2*g2*vd*ZH(gt1,1)*ZH(gt2,2)*ZH(gt3,2))/2._dp
res = res+(-3*dvu*g1**2*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(-3*dvu*g2**2*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,2))/4._dp
res = res+(-3*dg1*g1*vu*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,2))/2._dp
res = res+(-3*dg2*g2*vu*ZH(gt1,2)*ZH(gt2,2)*ZH(gt3,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_Couplinghhhhhh  
 
 
Subroutine CT_CouplinghhHpmcHpm(gt1,gt2,gt3,g1,g2,vd,vu,ZH,ZP,dg1,dg2,dvd,            & 
& dvu,dZH,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),ZP(2,2),dg1,dg2,dvd,dvu,dZH(2,2),dZP(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhHpmcHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g1**2*vd*dZP(gt3,1)*ZH(gt1,1)*ZP(gt2,1))/4._dp
res = res-(g2**2*vd*dZP(gt3,1)*ZH(gt1,1)*ZP(gt2,1))/4._dp
res = res-(g2**2*vu*dZP(gt3,2)*ZH(gt1,1)*ZP(gt2,1))/4._dp
res = res+(g1**2*vu*dZP(gt3,1)*ZH(gt1,2)*ZP(gt2,1))/4._dp
res = res-(g2**2*vu*dZP(gt3,1)*ZH(gt1,2)*ZP(gt2,1))/4._dp
res = res-(g2**2*vd*dZP(gt3,2)*ZH(gt1,2)*ZP(gt2,1))/4._dp
res = res-(g2**2*vu*dZP(gt3,1)*ZH(gt1,1)*ZP(gt2,2))/4._dp
res = res+(g1**2*vd*dZP(gt3,2)*ZH(gt1,1)*ZP(gt2,2))/4._dp
res = res-(g2**2*vd*dZP(gt3,2)*ZH(gt1,1)*ZP(gt2,2))/4._dp
res = res-(g2**2*vd*dZP(gt3,1)*ZH(gt1,2)*ZP(gt2,2))/4._dp
res = res-(g1**2*vu*dZP(gt3,2)*ZH(gt1,2)*ZP(gt2,2))/4._dp
res = res-(g2**2*vu*dZP(gt3,2)*ZH(gt1,2)*ZP(gt2,2))/4._dp
res = res-(g1**2*vd*dZP(gt2,1)*ZH(gt1,1)*ZP(gt3,1))/4._dp
res = res-(g2**2*vd*dZP(gt2,1)*ZH(gt1,1)*ZP(gt3,1))/4._dp
res = res-(g2**2*vu*dZP(gt2,2)*ZH(gt1,1)*ZP(gt3,1))/4._dp
res = res+(g1**2*vu*dZP(gt2,1)*ZH(gt1,2)*ZP(gt3,1))/4._dp
res = res-(g2**2*vu*dZP(gt2,1)*ZH(gt1,2)*ZP(gt3,1))/4._dp
res = res-(g2**2*vd*dZP(gt2,2)*ZH(gt1,2)*ZP(gt3,1))/4._dp
res = res-(g1**2*vd*dZH(gt1,1)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res-(g2**2*vd*dZH(gt1,1)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res+(g1**2*vu*dZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res-(g2**2*vu*dZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res-(dvd*g1**2*ZH(gt1,1)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res-(dvd*g2**2*ZH(gt1,1)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res-(dg1*g1*vd*ZH(gt1,1)*ZP(gt2,1)*ZP(gt3,1))/2._dp
res = res-(dg2*g2*vd*ZH(gt1,1)*ZP(gt2,1)*ZP(gt3,1))/2._dp
res = res+(dvu*g1**2*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res-(dvu*g2**2*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1))/4._dp
res = res+(dg1*g1*vu*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1))/2._dp
res = res-(dg2*g2*vu*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,1))/2._dp
res = res-(g2**2*vu*dZH(gt1,1)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res-(g2**2*vd*dZH(gt1,2)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res-(dvu*g2**2*ZH(gt1,1)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res-(dg2*g2*vu*ZH(gt1,1)*ZP(gt2,2)*ZP(gt3,1))/2._dp
res = res-(dvd*g2**2*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,1))/4._dp
res = res-(dg2*g2*vd*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,1))/2._dp
res = res-(g2**2*vu*dZP(gt2,1)*ZH(gt1,1)*ZP(gt3,2))/4._dp
res = res+(g1**2*vd*dZP(gt2,2)*ZH(gt1,1)*ZP(gt3,2))/4._dp
res = res-(g2**2*vd*dZP(gt2,2)*ZH(gt1,1)*ZP(gt3,2))/4._dp
res = res-(g2**2*vd*dZP(gt2,1)*ZH(gt1,2)*ZP(gt3,2))/4._dp
res = res-(g1**2*vu*dZP(gt2,2)*ZH(gt1,2)*ZP(gt3,2))/4._dp
res = res-(g2**2*vu*dZP(gt2,2)*ZH(gt1,2)*ZP(gt3,2))/4._dp
res = res-(g2**2*vu*dZH(gt1,1)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(g2**2*vd*dZH(gt1,2)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(dvu*g2**2*ZH(gt1,1)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(dg2*g2*vu*ZH(gt1,1)*ZP(gt2,1)*ZP(gt3,2))/2._dp
res = res-(dvd*g2**2*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,2))/4._dp
res = res-(dg2*g2*vd*ZH(gt1,2)*ZP(gt2,1)*ZP(gt3,2))/2._dp
res = res+(g1**2*vd*dZH(gt1,1)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res-(g2**2*vd*dZH(gt1,1)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res-(g1**2*vu*dZH(gt1,2)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res-(g2**2*vu*dZH(gt1,2)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res+(dvd*g1**2*ZH(gt1,1)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res-(dvd*g2**2*ZH(gt1,1)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res+(dg1*g1*vd*ZH(gt1,1)*ZP(gt2,2)*ZP(gt3,2))/2._dp
res = res-(dg2*g2*vd*ZH(gt1,1)*ZP(gt2,2)*ZP(gt3,2))/2._dp
res = res-(dvu*g1**2*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res-(dvu*g2**2*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,2))/4._dp
res = res-(dg1*g1*vu*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,2))/2._dp
res = res-(dg2*g2*vu*ZH(gt1,2)*ZP(gt2,2)*ZP(gt3,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhHpmcHpm  
 
 
Subroutine CT_CouplinghhSdcSd(gt1,gt2,gt3,g1,g2,Mu,Yd,Td,vd,vu,ZD,ZH,dg1,             & 
& dg2,dMu,dYd,dTd,dvd,dvu,dZD,dZH,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),dg1,dg2,dvd,dvu,dZH(2,2)

Complex(dp), Intent(in) :: Mu,Yd(3,3),Td(3,3),ZD(6,6),dMu,dYd(3,3),dTd(3,3),dZD(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhSdcSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZD(gt2,j1))*dZH(gt1,1)*ZD(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vd*Conjg(ZD(gt2,j1))*dZH(gt1,1)*ZD(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZD(gt2,j1))*dZH(gt1,2)*ZD(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZD(gt2,j1))*dZH(gt1,2)*ZD(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZD(gt2,3 + j1))*dZH(gt1,1)*ZD(gt3,3 + j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZD(gt2,3 + j1))*dZH(gt1,2)*ZD(gt3,3 + j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZD(gt2,j1))*dZD(gt3,j1)*ZH(gt1,1))/12._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vd*Conjg(ZD(gt2,j1))*dZD(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZD(gt2,3 + j1))*dZD(gt3,3 + j1)*ZH(gt1,1))/6._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(dZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,1))/12._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vd*Conjg(dZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res+(dvd*g1**2*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,1))/12._dp
End Do 
Do j1 = 1,3
res = res+(dvd*g2**2*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res+(dg1*g1*vd*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dg2*g2*vd*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,1))/2._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(dZD(gt2,3 + j1))*ZD(gt3,3 + j1)*ZH(gt1,1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dvd*g1**2*Conjg(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)*ZH(gt1,1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dg1*g1*vd*Conjg(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)*ZH(gt1,1))/3._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZD(gt2,j1))*dZD(gt3,j1)*ZH(gt1,2))/12._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZD(gt2,j1))*dZD(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZD(gt2,3 + j1))*dZD(gt3,3 + j1)*ZH(gt1,2))/6._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(dZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,2))/12._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(dZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res-(dvu*g1**2*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,2))/12._dp
End Do 
Do j1 = 1,3
res = res-(dvu*g2**2*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res-(dg1*g1*vu*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,2))/6._dp
End Do 
Do j1 = 1,3
res = res-(dg2*g2*vu*Conjg(ZD(gt2,j1))*ZD(gt3,j1)*ZH(gt1,2))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(dZD(gt2,3 + j1))*ZD(gt3,3 + j1)*ZH(gt1,2))/6._dp
End Do 
Do j1 = 1,3
res = res-(dvu*g1**2*Conjg(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)*ZH(gt1,2))/6._dp
End Do 
Do j1 = 1,3
res = res-(dg1*g1*vu*Conjg(ZD(gt2,3 + j1))*ZD(gt3,3 + j1)*ZH(gt1,2))/3._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZD(gt2,j2))*dZH(gt1,2)*Yd(j1,j2)*ZD(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZD(gt2,3 + j1))*Conjg(Td(j1,j2))*dZH(gt1,1)*ZD(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(Yd(j1,j2))*Conjg(ZD(gt2,3 + j1))*dZH(gt1,2)*ZD(gt3,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZD(gt2,3 + j1))*Conjg(Td(j1,j2))*dZD(gt3,j2)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZD(gt2,j2))*dTd(j1,j2)*ZD(gt3,3 + j1)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dTd(j1,j2))*Conjg(ZD(gt2,3 + j1))*ZD(gt3,j2)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZD(gt2,3 + j1))*Conjg(Td(j1,j2))*ZD(gt3,j2)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(Yd(j1,j2))*Conjg(ZD(gt2,3 + j1))*dZD(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZD(gt2,j2))*dZD(gt3,3 + j1)*Yd(j1,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZD(gt2,j2))*dYd(j1,j2)*ZD(gt3,3 + j1)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(dZD(gt2,j2))*Yd(j1,j2)*ZD(gt3,3 + j1)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dMu)*Conjg(ZD(gt2,j2))*Yd(j1,j2)*ZD(gt3,3 + j1)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(dZD(gt2,3 + j1))*Conjg(Yd(j1,j2))*ZD(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(dYd(j1,j2))*Conjg(ZD(gt2,3 + j1))*ZD(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(dMu*Conjg(Yd(j1,j2))*Conjg(ZD(gt2,3 + j1))*ZD(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZD(gt2,j2))*dZH(gt1,1)*ZD(gt3,3 + j1)*Td(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZD(gt2,j2))*dZD(gt3,3 + j1)*ZH(gt1,1)*Td(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZD(gt2,j2))*ZD(gt3,3 + j1)*ZH(gt1,1)*Td(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Yd(j3,j1))*Conjg(ZD(gt2,3 + j3))*dZH(gt1,1)*Yd(j2,j1)*ZD(gt3,3 + j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Yd(j1,j3))*Conjg(ZD(gt2,j2))*dZH(gt1,1)*Yd(j1,j2)*ZD(gt3,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Yd(j1,j3))*Conjg(ZD(gt2,j2))*dZD(gt3,j3)*Yd(j1,j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Yd(j3,j1))*Conjg(ZD(gt2,3 + j3))*dZD(gt3,3 + j2)*Yd(j2,j1)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Yd(j3,j1))*Conjg(ZD(gt2,3 + j3))*dYd(j2,j1)*ZD(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZD(gt2,3 + j3))*Conjg(Yd(j3,j1))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dYd(j3,j1))*Conjg(ZD(gt2,3 + j3))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(Yd(j3,j1))*Conjg(ZD(gt2,3 + j3))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Yd(j1,j3))*Conjg(ZD(gt2,j2))*dYd(j1,j2)*ZD(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZD(gt2,j2))*Conjg(Yd(j1,j3))*Yd(j1,j2)*ZD(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dYd(j1,j3))*Conjg(ZD(gt2,j2))*Yd(j1,j2)*ZD(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(Yd(j1,j3))*Conjg(ZD(gt2,j2))*Yd(j1,j2)*ZD(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhSdcSd  
 
 
Subroutine CT_CouplinghhSecSe(gt1,gt2,gt3,g1,g2,Mu,Ye,Te,vd,vu,ZE,ZH,dg1,             & 
& dg2,dMu,dYe,dTe,dvd,dvu,dZE,dZH,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),dg1,dg2,dvd,dvu,dZH(2,2)

Complex(dp), Intent(in) :: Mu,Ye(3,3),Te(3,3),ZE(6,6),dMu,dYe(3,3),dTe(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhSecSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g1**2*vd*Conjg(ZE(gt2,j1))*dZH(gt1,1)*ZE(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vd*Conjg(ZE(gt2,j1))*dZH(gt1,1)*ZE(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vu*Conjg(ZE(gt2,j1))*dZH(gt1,2)*ZE(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZE(gt2,j1))*dZH(gt1,2)*ZE(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZE(gt2,3 + j1))*dZH(gt1,1)*ZE(gt3,3 + j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZE(gt2,3 + j1))*dZH(gt1,2)*ZE(gt3,3 + j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vd*Conjg(ZE(gt2,j1))*dZE(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vd*Conjg(ZE(gt2,j1))*dZE(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZE(gt2,3 + j1))*dZE(gt3,3 + j1)*ZH(gt1,1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vd*Conjg(dZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vd*Conjg(dZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res-(dvd*g1**2*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res+(dvd*g2**2*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res-(dg1*g1*vd*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,1))/2._dp
End Do 
Do j1 = 1,3
res = res+(dg2*g2*vd*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,1))/2._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(dZE(gt2,3 + j1))*ZE(gt3,3 + j1)*ZH(gt1,1))/2._dp
End Do 
Do j1 = 1,3
res = res+(dvd*g1**2*Conjg(ZE(gt2,3 + j1))*ZE(gt3,3 + j1)*ZH(gt1,1))/2._dp
End Do 
Do j1 = 1,3
res = res+dg1*g1*vd*Conjg(ZE(gt2,3 + j1))*ZE(gt3,3 + j1)*ZH(gt1,1)
End Do 
Do j1 = 1,3
res = res+(g1**2*vu*Conjg(ZE(gt2,j1))*dZE(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZE(gt2,j1))*dZE(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZE(gt2,3 + j1))*dZE(gt3,3 + j1)*ZH(gt1,2))/2._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vu*Conjg(dZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(dZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res+(dvu*g1**2*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res-(dvu*g2**2*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res+(dg1*g1*vu*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,2))/2._dp
End Do 
Do j1 = 1,3
res = res-(dg2*g2*vu*Conjg(ZE(gt2,j1))*ZE(gt3,j1)*ZH(gt1,2))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(dZE(gt2,3 + j1))*ZE(gt3,3 + j1)*ZH(gt1,2))/2._dp
End Do 
Do j1 = 1,3
res = res-(dvu*g1**2*Conjg(ZE(gt2,3 + j1))*ZE(gt3,3 + j1)*ZH(gt1,2))/2._dp
End Do 
Do j1 = 1,3
res = res-(dg1*g1*vu*Conjg(ZE(gt2,3 + j1))*ZE(gt3,3 + j1)*ZH(gt1,2))
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZE(gt2,j2))*dZH(gt1,2)*Ye(j1,j2)*ZE(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZE(gt2,3 + j1))*Conjg(Te(j1,j2))*dZH(gt1,1)*ZE(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(Ye(j1,j2))*Conjg(ZE(gt2,3 + j1))*dZH(gt1,2)*ZE(gt3,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZE(gt2,3 + j1))*Conjg(Te(j1,j2))*dZE(gt3,j2)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZE(gt2,j2))*dTe(j1,j2)*ZE(gt3,3 + j1)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dTe(j1,j2))*Conjg(ZE(gt2,3 + j1))*ZE(gt3,j2)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZE(gt2,3 + j1))*Conjg(Te(j1,j2))*ZE(gt3,j2)*ZH(gt1,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(Ye(j1,j2))*Conjg(ZE(gt2,3 + j1))*dZE(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZE(gt2,j2))*dZE(gt3,3 + j1)*Ye(j1,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZE(gt2,j2))*dYe(j1,j2)*ZE(gt3,3 + j1)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(dZE(gt2,j2))*Ye(j1,j2)*ZE(gt3,3 + j1)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dMu)*Conjg(ZE(gt2,j2))*Ye(j1,j2)*ZE(gt3,3 + j1)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(dZE(gt2,3 + j1))*Conjg(Ye(j1,j2))*ZE(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(dYe(j1,j2))*Conjg(ZE(gt2,3 + j1))*ZE(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(dMu*Conjg(Ye(j1,j2))*Conjg(ZE(gt2,3 + j1))*ZE(gt3,j2)*ZH(gt1,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZE(gt2,j2))*dZH(gt1,1)*ZE(gt3,3 + j1)*Te(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZE(gt2,j2))*dZE(gt3,3 + j1)*ZH(gt1,1)*Te(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZE(gt2,j2))*ZE(gt3,3 + j1)*ZH(gt1,1)*Te(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j3,j1))*Conjg(ZE(gt2,3 + j3))*dZH(gt1,1)*Ye(j2,j1)*ZE(gt3,3 + j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j1,j3))*Conjg(ZE(gt2,j2))*dZH(gt1,1)*Ye(j1,j2)*ZE(gt3,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j1,j3))*Conjg(ZE(gt2,j2))*dZE(gt3,j3)*Ye(j1,j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j3,j1))*Conjg(ZE(gt2,3 + j3))*dZE(gt3,3 + j2)*Ye(j2,j1)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j3,j1))*Conjg(ZE(gt2,3 + j3))*dYe(j2,j1)*ZE(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZE(gt2,3 + j3))*Conjg(Ye(j3,j1))*Ye(j2,j1)*ZE(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dYe(j3,j1))*Conjg(ZE(gt2,3 + j3))*Ye(j2,j1)*ZE(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(Ye(j3,j1))*Conjg(ZE(gt2,3 + j3))*Ye(j2,j1)*ZE(gt3,3 + j2)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j1,j3))*Conjg(ZE(gt2,j2))*dYe(j1,j2)*ZE(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZE(gt2,j2))*Conjg(Ye(j1,j3))*Ye(j1,j2)*ZE(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dYe(j1,j3))*Conjg(ZE(gt2,j2))*Ye(j1,j2)*ZE(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(Ye(j1,j3))*Conjg(ZE(gt2,j2))*Ye(j1,j2)*ZE(gt3,j3)*ZH(gt1,1))
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhSecSe  
 
 
Subroutine CT_CouplinghhSucSu(gt1,gt2,gt3,g1,g2,Mu,Yu,Tu,vd,vu,ZU,ZH,dg1,             & 
& dg2,dMu,dYu,dTu,dvd,dvu,dZU,dZH,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),dg1,dg2,dvd,dvu,dZH(2,2)

Complex(dp), Intent(in) :: Mu,Yu(3,3),Tu(3,3),ZU(6,6),dMu,dYu(3,3),dTu(3,3),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhSucSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZU(gt2,j1))*dZU(gt3,j1)*ZH(gt1,1))/12._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZU(gt2,j1))*dZU(gt3,j1)*ZH(gt1,1))/4._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vd*Conjg(ZU(gt2,3 + j1))*dZU(gt3,3 + j1)*ZH(gt1,1))/3._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZU(gt2,j1))*dZU(gt3,j1)*ZH(gt1,2))/12._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vu*Conjg(ZU(gt2,j1))*dZU(gt3,j1)*ZH(gt1,2))/4._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vu*Conjg(ZU(gt2,3 + j1))*dZU(gt3,3 + j1)*ZH(gt1,2))/3._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(ZU(gt2,j1))*dZH(gt1,1)*ZU(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZU(gt2,j1))*dZH(gt1,1)*ZU(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(ZU(gt2,j1))*dZH(gt1,2)*ZU(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vu*Conjg(ZU(gt2,j1))*dZH(gt1,2)*ZU(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vd*Conjg(dZU(gt2,j1))*ZH(gt1,1)*ZU(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(dZU(gt2,j1))*ZH(gt1,1)*ZU(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res+(dvd*g1**2*Conjg(ZU(gt2,j1))*ZH(gt1,1)*ZU(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res-(dvd*g2**2*Conjg(ZU(gt2,j1))*ZH(gt1,1)*ZU(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res+(dg1*g1*vd*Conjg(ZU(gt2,j1))*ZH(gt1,1)*ZU(gt3,j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(dg2*g2*vd*Conjg(ZU(gt2,j1))*ZH(gt1,1)*ZU(gt3,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vu*Conjg(dZU(gt2,j1))*ZH(gt1,2)*ZU(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res+(g2**2*vu*Conjg(dZU(gt2,j1))*ZH(gt1,2)*ZU(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res-(dvu*g1**2*Conjg(ZU(gt2,j1))*ZH(gt1,2)*ZU(gt3,j1))/12._dp
End Do 
Do j1 = 1,3
res = res+(dvu*g2**2*Conjg(ZU(gt2,j1))*ZH(gt1,2)*ZU(gt3,j1))/4._dp
End Do 
Do j1 = 1,3
res = res-(dg1*g1*vu*Conjg(ZU(gt2,j1))*ZH(gt1,2)*ZU(gt3,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dg2*g2*vu*Conjg(ZU(gt2,j1))*ZH(gt1,2)*ZU(gt3,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vd*Conjg(ZU(gt2,3 + j1))*dZH(gt1,1)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vu*Conjg(ZU(gt2,3 + j1))*dZH(gt1,2)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res-(g1**2*vd*Conjg(dZU(gt2,3 + j1))*ZH(gt1,1)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res-(dvd*g1**2*Conjg(ZU(gt2,3 + j1))*ZH(gt1,1)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(-2*dg1*g1*vd*Conjg(ZU(gt2,3 + j1))*ZH(gt1,1)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(g1**2*vu*Conjg(dZU(gt2,3 + j1))*ZH(gt1,2)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(dvu*g1**2*Conjg(ZU(gt2,3 + j1))*ZH(gt1,2)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(2*dg1*g1*vu*Conjg(ZU(gt2,3 + j1))*ZH(gt1,2)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*dZU(gt3,j2)*ZH(gt1,1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZU(gt2,j2))*dZU(gt3,3 + j1)*Yu(j1,j2)*ZH(gt1,1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*dZU(gt3,j2)*ZH(gt1,2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZU(gt2,j2))*dZH(gt1,1)*Yu(j1,j2)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(ZU(gt2,j2))*dYu(j1,j2)*ZH(gt1,1)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(Mu)*Conjg(dZU(gt2,j2))*Yu(j1,j2)*ZH(gt1,1)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dMu)*Conjg(ZU(gt2,j2))*Yu(j1,j2)*ZH(gt1,1)*ZU(gt3,3 + j1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZU(gt2,j2))*dTu(j1,j2)*ZH(gt1,2)*ZU(gt3,3 + j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*dZH(gt1,1)*ZU(gt3,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*dZH(gt1,2)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(dZU(gt2,3 + j1))*Conjg(Yu(j1,j2))*ZH(gt1,1)*ZU(gt3,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(Mu*Conjg(dYu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZH(gt1,1)*ZU(gt3,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+(dMu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZH(gt1,1)*ZU(gt3,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dTu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZH(gt1,2)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*ZH(gt1,2)*ZU(gt3,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZU(gt2,j2))*dZU(gt3,3 + j1)*ZH(gt1,2)*Tu(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(ZU(gt2,j2))*dZH(gt1,2)*ZU(gt3,3 + j1)*Tu(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res-((Conjg(dZU(gt2,j2))*ZH(gt1,2)*ZU(gt3,3 + j1)*Tu(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*dZU(gt3,j3)*Yu(j1,j2)*ZH(gt1,2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dZU(gt3,3 + j2)*Yu(j2,j1)*ZH(gt1,2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dZH(gt1,2)*Yu(j2,j1)*ZU(gt3,3 + j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dYu(j2,j1)*ZH(gt1,2)*ZU(gt3,3 + j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(dZU(gt2,3 + j3))*Conjg(Yu(j3,j1))*Yu(j2,j1)*ZH(gt1,2)*ZU(gt3,3 + j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(dYu(j3,j1))*Conjg(ZU(gt2,3 + j3))*Yu(j2,j1)*ZH(gt1,2)*ZU(gt3,3 + j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*Yu(j2,j1)*ZH(gt1,2)*ZU(gt3,3 + j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*dZH(gt1,2)*Yu(j1,j2)*ZU(gt3,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*dYu(j1,j2)*ZH(gt1,2)*ZU(gt3,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(dZU(gt2,j2))*Conjg(Yu(j1,j3))*Yu(j1,j2)*ZH(gt1,2)*ZU(gt3,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vu*Conjg(dYu(j1,j3))*Conjg(ZU(gt2,j2))*Yu(j1,j2)*ZH(gt1,2)*ZU(gt3,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*Yu(j1,j2)*ZH(gt1,2)*ZU(gt3,j3))
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhSucSu  
 
 
Subroutine CT_CouplinghhSvcSv(gt1,gt2,gt3,g1,g2,vd,vu,ZH,dg1,dg2,dvd,dvu,dZH,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),dg1,dg2,dvd,dvu,dZH(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhSvcSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
If ((gt2.eq.gt3)) Then 
res = res-(g1**2*vd*dZH(gt1,1))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res-(g2**2*vd*dZH(gt1,1))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res+(g1**2*vu*dZH(gt1,2))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res+(g2**2*vu*dZH(gt1,2))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res-(dvd*g1**2*ZH(gt1,1))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res-(dvd*g2**2*ZH(gt1,1))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res-(dg1*g1*vd*ZH(gt1,1))/2._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res-(dg2*g2*vd*ZH(gt1,1))/2._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res+(dvu*g1**2*ZH(gt1,2))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res+(dvu*g2**2*ZH(gt1,2))/4._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res+(dg1*g1*vu*ZH(gt1,2))/2._dp
End If 
If ((gt2.eq.gt3)) Then 
res = res+(dg2*g2*vu*ZH(gt1,2))/2._dp
End If 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhSvcSv  
 
 
Subroutine CT_CouplingHpmSucSd(gt1,gt2,gt3,g2,Mu,Yd,Td,Yu,Tu,vd,vu,ZD,ZU,             & 
& ZP,dg2,dMu,dYd,dTd,dYu,dTu,dvd,dvu,dZD,dZU,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,vd,vu,ZP(2,2),dg2,dvd,dvu,dZP(2,2)

Complex(dp), Intent(in) :: Mu,Yd(3,3),Td(3,3),Yu(3,3),Tu(3,3),ZD(6,6),ZU(6,6),dMu,dYd(3,3),dTd(3,3),             & 
& dYu(3,3),dTu(3,3),dZD(6,6),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingHpmSucSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZU(gt2,j1))*dZP(gt1,1)*ZD(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZU(gt2,j1))*dZP(gt1,2)*ZD(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZU(gt2,j1))*dZD(gt3,j1)*ZP(gt1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(dZU(gt2,j1))*ZD(gt3,j1)*ZP(gt1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvd*g2**2*Conjg(ZU(gt2,j1))*ZD(gt3,j1)*ZP(gt1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vd*Conjg(ZU(gt2,j1))*ZD(gt3,j1)*ZP(gt1,1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZU(gt2,j1))*dZD(gt3,j1)*ZP(gt1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(dZU(gt2,j1))*ZD(gt3,j1)*ZP(gt1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvu*g2**2*Conjg(ZU(gt2,j1))*ZD(gt3,j1)*ZP(gt1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vu*Conjg(ZU(gt2,j1))*ZD(gt3,j1)*ZP(gt1,2))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZU(gt2,j2))*dZP(gt1,2)*Yd(j1,j2)*ZD(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*dZP(gt1,1)*ZD(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*dZP(gt1,2)*ZD(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*dZD(gt3,j2)*ZP(gt1,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZU(gt2,j2))*dTd(j1,j2)*ZD(gt3,3 + j1)*ZP(gt1,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(dZU(gt2,3 + j1))*Conjg(Yu(j1,j2))*ZD(gt3,j2)*ZP(gt1,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(dYu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZD(gt3,j2)*ZP(gt1,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+dMu*Conjg(Yu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZD(gt3,j2)*ZP(gt1,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*dZD(gt3,j2)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZU(gt2,j2))*dZD(gt3,3 + j1)*Yd(j1,j2)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZU(gt2,j2))*dYd(j1,j2)*ZD(gt3,3 + j1)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(dZU(gt2,j2))*Yd(j1,j2)*ZD(gt3,3 + j1)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dMu)*Conjg(ZU(gt2,j2))*Yd(j1,j2)*ZD(gt3,3 + j1)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dTu(j1,j2))*Conjg(ZU(gt2,3 + j1))*ZD(gt3,j2)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZU(gt2,3 + j1))*Conjg(Tu(j1,j2))*ZD(gt3,j2)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZU(gt2,j2))*dZP(gt1,1)*ZD(gt3,3 + j1)*Td(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZU(gt2,j2))*dZD(gt3,3 + j1)*ZP(gt1,1)*Td(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZU(gt2,j2))*ZD(gt3,3 + j1)*ZP(gt1,1)*Td(j1,j2)
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dZP(gt1,1)*Yd(j2,j1)*ZD(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dZP(gt1,2)*Yd(j2,j1)*ZD(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j3))*Conjg(ZU(gt2,j2))*dZP(gt1,1)*Yd(j1,j2)*ZD(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*dZP(gt1,2)*Yu(j1,j2)*ZD(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j3))*Conjg(ZU(gt2,j2))*dZD(gt3,j3)*Yd(j1,j2)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dZD(gt3,3 + j2)*Yd(j2,j1)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dYd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dZU(gt2,3 + j3))*Conjg(Yu(j3,j1))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dYu(j3,j1))*Conjg(ZU(gt2,3 + j3))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvu*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j3))*Conjg(ZU(gt2,j2))*dYd(j1,j2)*ZD(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZU(gt2,j2))*Conjg(Yd(j1,j3))*Yd(j1,j2)*ZD(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYd(j1,j3))*Conjg(ZU(gt2,j2))*Yd(j1,j2)*ZD(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Yd(j1,j3))*Conjg(ZU(gt2,j2))*Yd(j1,j2)*ZD(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dZD(gt3,3 + j2)*Yd(j2,j1)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*dZD(gt3,j3)*Yu(j1,j2)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*dYd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZU(gt2,3 + j3))*Conjg(Yu(j3,j1))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYu(j3,j1))*Conjg(ZU(gt2,3 + j3))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Yu(j3,j1))*Conjg(ZU(gt2,3 + j3))*Yd(j2,j1)*ZD(gt3,3 + j2)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*dYu(j1,j2)*ZD(gt3,j3)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dZU(gt2,j2))*Conjg(Yu(j1,j3))*Yu(j1,j2)*ZD(gt3,j3)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dYu(j1,j3))*Conjg(ZU(gt2,j2))*Yu(j1,j2)*ZD(gt3,j3)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvu*Conjg(Yu(j1,j3))*Conjg(ZU(gt2,j2))*Yu(j1,j2)*ZD(gt3,j3)*ZP(gt1,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingHpmSucSd  
 
 
Subroutine CT_CouplingHpmSvcSe(gt1,gt2,gt3,g2,Mu,Ye,Te,vd,vu,ZV,ZE,ZP,dg2,            & 
& dMu,dYe,dTe,dvd,dvu,dZV,dZE,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,vd,vu,ZP(2,2),dg2,dvd,dvu,dZP(2,2)

Complex(dp), Intent(in) :: Mu,Ye(3,3),Te(3,3),ZV(3,3),ZE(6,6),dMu,dYe(3,3),dTe(3,3),dZV(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingHpmSvcSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZV(gt2,j1))*dZP(gt1,1)*ZE(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZV(gt2,j1))*dZP(gt1,2)*ZE(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZV(gt2,j1))*dZE(gt3,j1)*ZP(gt1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(dZV(gt2,j1))*ZE(gt3,j1)*ZP(gt1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvd*g2**2*Conjg(ZV(gt2,j1))*ZE(gt3,j1)*ZP(gt1,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vd*Conjg(ZV(gt2,j1))*ZE(gt3,j1)*ZP(gt1,1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZV(gt2,j1))*dZE(gt3,j1)*ZP(gt1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(dZV(gt2,j1))*ZE(gt3,j1)*ZP(gt1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvu*g2**2*Conjg(ZV(gt2,j1))*ZE(gt3,j1)*ZP(gt1,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vu*Conjg(ZV(gt2,j1))*ZE(gt3,j1)*ZP(gt1,2))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZV(gt2,j2))*dZP(gt1,2)*Ye(j1,j2)*ZE(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZV(gt2,j2))*dTe(j1,j2)*ZE(gt3,3 + j1)*ZP(gt1,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZV(gt2,j2))*dZE(gt3,3 + j1)*Ye(j1,j2)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZV(gt2,j2))*dYe(j1,j2)*ZE(gt3,3 + j1)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(dZV(gt2,j2))*Ye(j1,j2)*ZE(gt3,3 + j1)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dMu)*Conjg(ZV(gt2,j2))*Ye(j1,j2)*ZE(gt3,3 + j1)*ZP(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZV(gt2,j2))*dZP(gt1,1)*ZE(gt3,3 + j1)*Te(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZV(gt2,j2))*dZE(gt3,3 + j1)*ZP(gt1,1)*Te(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZV(gt2,j2))*ZE(gt3,3 + j1)*ZP(gt1,1)*Te(j1,j2)
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j3))*Conjg(ZV(gt2,j2))*dZP(gt1,1)*Ye(j1,j2)*ZE(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j3))*Conjg(ZV(gt2,j2))*dZE(gt3,j3)*Ye(j1,j2)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j3))*Conjg(ZV(gt2,j2))*dYe(j1,j2)*ZE(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZV(gt2,j2))*Conjg(Ye(j1,j3))*Ye(j1,j2)*ZE(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYe(j1,j3))*Conjg(ZV(gt2,j2))*Ye(j1,j2)*ZE(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Ye(j1,j3))*Conjg(ZV(gt2,j2))*Ye(j1,j2)*ZE(gt3,j3)*ZP(gt1,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingHpmSvcSe  
 
 
Subroutine CT_CouplingSdSvcSd(gt1,gt2,gt3,Yd,L2,T2,vd,ZD,ZV,dYd,dL2,dT2,              & 
& dvd,dZD,dZV,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: vd,dvd

Complex(dp), Intent(in) :: Yd(3,3),L2(3,3,3),T2(3,3,3),ZD(6,6),ZV(3,3),dYd(3,3),dL2(3,3,3),dT2(3,3,3),           & 
& dZD(6,6),dZV(3,3)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdSvcSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZD(gt1,j2))*Conjg(ZV(gt2,j1))*dT2(j1,j2,j3)*ZD(gt3,3 + j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZD(gt1,j2))*Conjg(ZV(gt2,j1))*dZD(gt3,3 + j3)*T2(j1,j2,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dZV(gt2,j1))*Conjg(ZD(gt1,j2))*ZD(gt3,3 + j3)*T2(j1,j2,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dZD(gt1,j2))*Conjg(ZV(gt2,j1))*ZD(gt3,3 + j3)*T2(j1,j2,j3))
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(Yd(j4,j1))*Conjg(ZD(gt1,3 + j4))*Conjg(ZV(gt2,j2))*dZD(gt3,3 + j3)*L2(j2,j1,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(Yd(j1,j4))*Conjg(ZD(gt1,j3))*Conjg(ZV(gt2,j2))*dZD(gt3,j4)*L2(j2,j3,j1))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(Yd(j4,j1))*Conjg(ZD(gt1,3 + j4))*Conjg(ZV(gt2,j2))*dL2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dZV(gt2,j2))*Conjg(Yd(j4,j1))*Conjg(ZD(gt1,3 + j4))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dZD(gt1,3 + j4))*Conjg(Yd(j4,j1))*Conjg(ZV(gt2,j2))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dYd(j4,j1))*Conjg(ZD(gt1,3 + j4))*Conjg(ZV(gt2,j2))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((dvd*Conjg(Yd(j4,j1))*Conjg(ZD(gt1,3 + j4))*Conjg(ZV(gt2,j2))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(Yd(j1,j4))*Conjg(ZD(gt1,j3))*Conjg(ZV(gt2,j2))*dL2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dZV(gt2,j2))*Conjg(Yd(j1,j4))*Conjg(ZD(gt1,j3))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dZD(gt1,j3))*Conjg(Yd(j1,j4))*Conjg(ZV(gt2,j2))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dYd(j1,j4))*Conjg(ZD(gt1,j3))*Conjg(ZV(gt2,j2))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((dvd*Conjg(Yd(j1,j4))*Conjg(ZD(gt1,j3))*Conjg(ZV(gt2,j2))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdSvcSd  
 
 
Subroutine CT_CouplingSdcHpmcSu(gt1,gt2,gt3,g2,Mu,Yd,Td,Yu,Tu,vd,vu,ZD,               & 
& ZU,ZP,dg2,dMu,dYd,dTd,dYu,dTu,dvd,dvu,dZD,dZU,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,vd,vu,ZP(2,2),dg2,dvd,dvu,dZP(2,2)

Complex(dp), Intent(in) :: Mu,Yd(3,3),Td(3,3),Yu(3,3),Tu(3,3),ZD(6,6),ZU(6,6),dMu,dYd(3,3),dTd(3,3),             & 
& dYu(3,3),dTu(3,3),dZD(6,6),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdcHpmcSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZD(gt1,j1))*dZU(gt3,j1)*ZP(gt2,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZD(gt1,j1))*dZU(gt3,j1)*ZP(gt2,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZD(gt1,j1))*dZP(gt2,1)*ZU(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZD(gt1,j1))*dZP(gt2,2)*ZU(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(dZD(gt1,j1))*ZP(gt2,1)*ZU(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvd*g2**2*Conjg(ZD(gt1,j1))*ZP(gt2,1)*ZU(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vd*Conjg(ZD(gt1,j1))*ZP(gt2,1)*ZU(gt3,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(dZD(gt1,j1))*ZP(gt2,2)*ZU(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvu*g2**2*Conjg(ZD(gt1,j1))*ZP(gt2,2)*ZU(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vu*Conjg(ZD(gt1,j1))*ZP(gt2,2)*ZU(gt3,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZD(gt1,3 + j1))*Conjg(Td(j1,j2))*dZU(gt3,j2)*ZP(gt2,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZD(gt1,j2))*dZU(gt3,3 + j1)*Yu(j1,j2)*ZP(gt2,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(Yd(j1,j2))*Conjg(ZD(gt1,3 + j1))*dZU(gt3,j2)*ZP(gt2,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZD(gt1,j2))*dZP(gt2,1)*Yu(j1,j2)*ZU(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(ZD(gt1,j2))*dYu(j1,j2)*ZP(gt2,1)*ZU(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(Mu)*Conjg(dZD(gt1,j2))*Yu(j1,j2)*ZP(gt2,1)*ZU(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dMu)*Conjg(ZD(gt1,j2))*Yu(j1,j2)*ZP(gt2,1)*ZU(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZD(gt1,j2))*dTu(j1,j2)*ZP(gt2,2)*ZU(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZD(gt1,3 + j1))*Conjg(Td(j1,j2))*dZP(gt2,1)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(Yd(j1,j2))*Conjg(ZD(gt1,3 + j1))*dZP(gt2,2)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dTd(j1,j2))*Conjg(ZD(gt1,3 + j1))*ZP(gt2,1)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZD(gt1,3 + j1))*Conjg(Td(j1,j2))*ZP(gt2,1)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(dZD(gt1,3 + j1))*Conjg(Yd(j1,j2))*ZP(gt2,2)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(dYd(j1,j2))*Conjg(ZD(gt1,3 + j1))*ZP(gt2,2)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+dMu*Conjg(Yd(j1,j2))*Conjg(ZD(gt1,3 + j1))*ZP(gt2,2)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZD(gt1,j2))*dZU(gt3,3 + j1)*ZP(gt2,2)*Tu(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZD(gt1,j2))*dZP(gt2,2)*ZU(gt3,3 + j1)*Tu(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZD(gt1,j2))*ZP(gt2,2)*ZU(gt3,3 + j1)*Tu(j1,j2)
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j3))*Conjg(ZD(gt1,j2))*dZU(gt3,j3)*Yd(j1,j2)*ZP(gt2,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*dZU(gt3,3 + j2)*Yu(j2,j1)*ZP(gt2,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j1,j3))*Conjg(ZD(gt1,j2))*dZU(gt3,j3)*Yu(j1,j2)*ZP(gt2,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*dZU(gt3,3 + j2)*Yu(j2,j1)*ZP(gt2,2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*dZP(gt2,1)*Yu(j2,j1)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*dZP(gt2,2)*Yu(j2,j1)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*dYu(j2,j1)*ZP(gt2,1)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dZD(gt1,3 + j3))*Conjg(Yd(j3,j1))*Yu(j2,j1)*ZP(gt2,1)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dYd(j3,j1))*Conjg(ZD(gt1,3 + j3))*Yu(j2,j1)*ZP(gt2,1)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvu*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*Yu(j2,j1)*ZP(gt2,1)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*dYu(j2,j1)*ZP(gt2,2)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZD(gt1,3 + j3))*Conjg(Yd(j3,j1))*Yu(j2,j1)*ZP(gt2,2)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYd(j3,j1))*Conjg(ZD(gt1,3 + j3))*Yu(j2,j1)*ZP(gt2,2)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Yd(j3,j1))*Conjg(ZD(gt1,3 + j3))*Yu(j2,j1)*ZP(gt2,2)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j3))*Conjg(ZD(gt1,j2))*dZP(gt2,1)*Yd(j1,j2)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j1,j3))*Conjg(ZD(gt1,j2))*dZP(gt2,2)*Yu(j1,j2)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j3))*Conjg(ZD(gt1,j2))*dYd(j1,j2)*ZP(gt2,1)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZD(gt1,j2))*Conjg(Yd(j1,j3))*Yd(j1,j2)*ZP(gt2,1)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYd(j1,j3))*Conjg(ZD(gt1,j2))*Yd(j1,j2)*ZP(gt2,1)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Yd(j1,j3))*Conjg(ZD(gt1,j2))*Yd(j1,j2)*ZP(gt2,1)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j1,j3))*Conjg(ZD(gt1,j2))*dYu(j1,j2)*ZP(gt2,2)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dZD(gt1,j2))*Conjg(Yu(j1,j3))*Yu(j1,j2)*ZP(gt2,2)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dYu(j1,j3))*Conjg(ZD(gt1,j2))*Yu(j1,j2)*ZP(gt2,2)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvu*Conjg(Yu(j1,j3))*Conjg(ZD(gt1,j2))*Yu(j1,j2)*ZP(gt2,2)*ZU(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdcHpmcSu  
 
 
Subroutine CT_CouplingSdcSdcSv(gt1,gt2,gt3,Yd,L2,T2,vd,ZD,ZV,dYd,dL2,dT2,             & 
& dvd,dZD,dZV,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: vd,dvd

Complex(dp), Intent(in) :: Yd(3,3),L2(3,3,3),T2(3,3,3),ZD(6,6),ZV(3,3),dYd(3,3),dL2(3,3,3),dT2(3,3,3),           & 
& dZD(6,6),dZV(3,3)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdcSdcSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZD(gt1,3 + j3))*Conjg(T2(j1,j2,j3))*dZV(gt3,j1)*ZD(gt2,j2))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZD(gt1,3 + j3))*Conjg(T2(j1,j2,j3))*dZD(gt2,j2)*ZV(gt3,j1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dT2(j1,j2,j3))*Conjg(ZD(gt1,3 + j3))*ZD(gt2,j2)*ZV(gt3,j1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dZD(gt1,3 + j3))*Conjg(T2(j1,j2,j3))*ZD(gt2,j2)*ZV(gt3,j1))
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*dZV(gt3,j3)*Yd(j2,j1)*ZD(gt2,3 + j2))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*dZV(gt3,j3)*Yd(j1,j2)*ZD(gt2,j4))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*dZD(gt2,j4)*Yd(j1,j2)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*dZD(gt2,3 + j2)*Yd(j2,j1)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*dYd(j2,j1)*ZD(gt2,3 + j2)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dZD(gt1,3 + j4))*Conjg(L2(j3,j1,j4))*Yd(j2,j1)*ZD(gt2,3 + j2)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dL2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*Yd(j2,j1)*ZD(gt2,3 + j2)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((dvd*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*Yd(j2,j1)*ZD(gt2,3 + j2)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*dYd(j1,j2)*ZD(gt2,j4)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dZD(gt1,j2))*Conjg(L2(j3,j4,j1))*Yd(j1,j2)*ZD(gt2,j4)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((vd*Conjg(dL2(j3,j4,j1))*Conjg(ZD(gt1,j2))*Yd(j1,j2)*ZD(gt2,j4)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-((dvd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*Yd(j1,j2)*ZD(gt2,j4)*ZV(gt3,j3))/sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdcSdcSv  
 
 
Subroutine CT_CouplingSdcSecSu(gt1,gt2,gt3,Yd,Ye,L2,T2,Yu,vd,vu,ZD,ZU,ZE,             & 
& dYd,dYe,dL2,dT2,dYu,dvd,dvu,dZD,dZU,dZE,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: vd,vu,dvd,dvu

Complex(dp), Intent(in) :: Yd(3,3),Ye(3,3),L2(3,3,3),T2(3,3,3),Yu(3,3),ZD(6,6),ZU(6,6),ZE(6,6),dYd(3,3),         & 
& dYe(3,3),dL2(3,3,3),dT2(3,3,3),dYu(3,3),dZD(6,6),dZU(6,6),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdcSecSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZD(gt1,3 + j3))*Conjg(T2(j1,j2,j3))*dZU(gt3,j2)*ZE(gt2,j1)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZD(gt1,3 + j3))*Conjg(T2(j1,j2,j3))*dZE(gt2,j1)*ZU(gt3,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dT2(j1,j2,j3))*Conjg(ZD(gt1,3 + j3))*ZE(gt2,j1)*ZU(gt3,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZD(gt1,3 + j3))*Conjg(T2(j1,j2,j3))*ZE(gt2,j1)*ZU(gt3,j2)
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L2(j1,j4,j3))*Conjg(ZD(gt1,3 + j3))*dZU(gt3,j4)*Ye(j2,j1)*ZE(gt2,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*dZU(gt3,j4)*Yd(j1,j2)*ZE(gt2,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*dZU(gt3,3 + j2)*Yu(j2,j1)*ZE(gt2,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*dZE(gt2,j3)*Yu(j2,j1)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*dYu(j2,j1)*ZE(gt2,j3)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dZD(gt1,3 + j4))*Conjg(L2(j3,j1,j4))*Yu(j2,j1)*ZE(gt2,j3)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dL2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*Yu(j2,j1)*ZE(gt2,j3)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvu*Conjg(L2(j3,j1,j4))*Conjg(ZD(gt1,3 + j4))*Yu(j2,j1)*ZE(gt2,j3)*ZU(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*dZE(gt2,j3)*Yd(j1,j2)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L2(j1,j4,j3))*Conjg(ZD(gt1,3 + j3))*dZE(gt2,3 + j2)*Ye(j2,j1)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L2(j1,j4,j3))*Conjg(ZD(gt1,3 + j3))*dYe(j2,j1)*ZE(gt2,3 + j2)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZD(gt1,3 + j3))*Conjg(L2(j1,j4,j3))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dL2(j1,j4,j3))*Conjg(ZD(gt1,3 + j3))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(L2(j1,j4,j3))*Conjg(ZD(gt1,3 + j3))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*dYd(j1,j2)*ZE(gt2,j3)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZD(gt1,j2))*Conjg(L2(j3,j4,j1))*Yd(j1,j2)*ZE(gt2,j3)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dL2(j3,j4,j1))*Conjg(ZD(gt1,j2))*Yd(j1,j2)*ZE(gt2,j3)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(L2(j3,j4,j1))*Conjg(ZD(gt1,j2))*Yd(j1,j2)*ZE(gt2,j3)*ZU(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdcSecSu  
 
 
Subroutine CT_CouplingSeSucSd(gt1,gt2,gt3,Yd,Ye,L2,T2,Yu,vd,vu,ZD,ZU,ZE,              & 
& dYd,dYe,dL2,dT2,dYu,dvd,dvu,dZD,dZU,dZE,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: vd,vu,dvd,dvu

Complex(dp), Intent(in) :: Yd(3,3),Ye(3,3),L2(3,3,3),T2(3,3,3),Yu(3,3),ZD(6,6),ZU(6,6),ZE(6,6),dYd(3,3),         & 
& dYe(3,3),dL2(3,3,3),dT2(3,3,3),dYu(3,3),dZD(6,6),dZU(6,6),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSeSucSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZE(gt1,j1))*Conjg(ZU(gt2,j2))*dT2(j1,j2,j3)*ZD(gt3,3 + j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZE(gt1,j1))*Conjg(ZU(gt2,j2))*dZD(gt3,3 + j3)*T2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZU(gt2,j2))*Conjg(ZE(gt1,j1))*ZD(gt3,3 + j3)*T2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZE(gt1,j1))*Conjg(ZU(gt2,j2))*ZD(gt3,3 + j3)*T2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZU(gt2,j3))*dZD(gt3,3 + j2)*L2(j1,j3,j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j4,j1))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,3 + j4))*dZD(gt3,3 + j3)*L2(j2,j1,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,j3))*dZD(gt3,j4)*L2(j2,j3,j1))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZU(gt2,j3))*dL2(j1,j3,j2)*ZD(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZU(gt2,j3))*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*L2(j1,j3,j2)*ZD(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZE(gt1,3 + j4))*Conjg(Ye(j4,j1))*Conjg(ZU(gt2,j3))*L2(j1,j3,j2)*ZD(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYe(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZU(gt2,j3))*L2(j1,j3,j2)*ZD(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZU(gt2,j3))*L2(j1,j3,j2)*ZD(gt3,3 + j2))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(Yu(j4,j1))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,3 + j4))*dL2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dZU(gt2,3 + j4))*Conjg(Yu(j4,j1))*Conjg(ZE(gt1,j2))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dZE(gt1,j2))*Conjg(Yu(j4,j1))*Conjg(ZU(gt2,3 + j4))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vu*Conjg(dYu(j4,j1))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,3 + j4))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvu*Conjg(Yu(j4,j1))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,3 + j4))*L2(j2,j1,j3)*ZD(gt3,3 + j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Yd(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,j3))*dL2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZU(gt2,j3))*Conjg(Yd(j1,j4))*Conjg(ZE(gt1,j2))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZE(gt1,j2))*Conjg(Yd(j1,j4))*Conjg(ZU(gt2,j3))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYd(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,j3))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Yd(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZU(gt2,j3))*L2(j2,j3,j1)*ZD(gt3,j4))/sqrt(2._dp)
End Do 
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSeSucSd  
 
 
Subroutine CT_CouplingSeSvcSe(gt1,gt2,gt3,Ye,L1,T1,vd,ZV,ZE,dYe,dL1,dT1,              & 
& dvd,dZV,dZE,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: vd,dvd

Complex(dp), Intent(in) :: Ye(3,3),L1(3,3,3),T1(3,3,3),ZV(3,3),ZE(6,6),dYe(3,3),dL1(3,3,3),dT1(3,3,3),           & 
& dZV(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSeSvcSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZE(gt1,j2))*Conjg(ZV(gt2,j1))*dT1(j1,j2,j3)*ZE(gt3,3 + j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZE(gt1,j1))*Conjg(ZV(gt2,j2))*dT1(j1,j2,j3)*ZE(gt3,3 + j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZE(gt1,j2))*Conjg(ZV(gt2,j1))*dZE(gt3,3 + j3)*T1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZE(gt1,j1))*Conjg(ZV(gt2,j2))*dZE(gt3,3 + j3)*T1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dZV(gt2,j2))*Conjg(ZE(gt1,j1))*ZE(gt3,3 + j3)*T1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dZV(gt2,j1))*Conjg(ZE(gt1,j2))*ZE(gt3,3 + j3)*T1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dZE(gt1,j2))*Conjg(ZV(gt2,j1))*ZE(gt3,3 + j3)*T1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dZE(gt1,j1))*Conjg(ZV(gt2,j2))*ZE(gt3,3 + j3)*T1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j3))*dZE(gt3,3 + j2)*L1(j1,j3,j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j2))*dZE(gt3,3 + j3)*L1(j2,j1,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j3))*Conjg(ZV(gt2,j2))*dZE(gt3,j4)*L1(j2,j3,j1))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZV(gt2,j3))*dZE(gt3,j4)*L1(j2,j3,j1))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j3))*dL1(j1,j3,j2)*ZE(gt3,3 + j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZV(gt2,j3))*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*L1(j1,j3,j2)*ZE(gt3,3 + j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZE(gt1,3 + j4))*Conjg(Ye(j4,j1))*Conjg(ZV(gt2,j3))*L1(j1,j3,j2)*ZE(gt3,3 + j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYe(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j3))*L1(j1,j3,j2)*ZE(gt3,3 + j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j3))*L1(j1,j3,j2)*ZE(gt3,3 + j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j2))*dL1(j2,j1,j3)*ZE(gt3,3 + j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZV(gt2,j2))*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*L1(j2,j1,j3)*ZE(gt3,3 + j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZE(gt1,3 + j4))*Conjg(Ye(j4,j1))*Conjg(ZV(gt2,j2))*L1(j2,j1,j3)*ZE(gt3,3 + j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dYe(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j2))*L1(j2,j1,j3)*ZE(gt3,3 + j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(Ye(j4,j1))*Conjg(ZE(gt1,3 + j4))*Conjg(ZV(gt2,j2))*L1(j2,j1,j3)*ZE(gt3,3 + j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j3))*Conjg(ZV(gt2,j2))*dL1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZV(gt2,j3))*dL1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZV(gt2,j3))*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j2))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZV(gt2,j2))*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j3))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZE(gt1,j3))*Conjg(Ye(j1,j4))*Conjg(ZV(gt2,j2))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dYe(j1,j4))*Conjg(ZE(gt1,j3))*Conjg(ZV(gt2,j2))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j3))*Conjg(ZV(gt2,j2))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZE(gt1,j2))*Conjg(Ye(j1,j4))*Conjg(ZV(gt2,j3))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYe(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZV(gt2,j3))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Ye(j1,j4))*Conjg(ZE(gt1,j2))*Conjg(ZV(gt2,j3))*L1(j2,j3,j1)*ZE(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSeSvcSe  
 
 
Subroutine CT_CouplingSecHpmcSv(gt1,gt2,gt3,g2,Mu,Ye,Te,vd,vu,ZV,ZE,ZP,               & 
& dg2,dMu,dYe,dTe,dvd,dvu,dZV,dZE,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,vd,vu,ZP(2,2),dg2,dvd,dvu,dZP(2,2)

Complex(dp), Intent(in) :: Mu,Ye(3,3),Te(3,3),ZV(3,3),ZE(6,6),dMu,dYe(3,3),dTe(3,3),dZV(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSecHpmcSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZE(gt1,j1))*dZV(gt3,j1)*ZP(gt2,1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZE(gt1,j1))*dZV(gt3,j1)*ZP(gt2,2))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(ZE(gt1,j1))*dZP(gt2,1)*ZV(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(ZE(gt1,j1))*dZP(gt2,2)*ZV(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vd*Conjg(dZE(gt1,j1))*ZP(gt2,1)*ZV(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvd*g2**2*Conjg(ZE(gt1,j1))*ZP(gt2,1)*ZV(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vd*Conjg(ZE(gt1,j1))*ZP(gt2,1)*ZV(gt3,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(g2**2*vu*Conjg(dZE(gt1,j1))*ZP(gt2,2)*ZV(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-(dvu*g2**2*Conjg(ZE(gt1,j1))*ZP(gt2,2)*ZV(gt3,j1))/(2._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
res = res-((dg2*g2*vu*Conjg(ZE(gt1,j1))*ZP(gt2,2)*ZV(gt3,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZE(gt1,3 + j1))*Conjg(Te(j1,j2))*dZV(gt3,j2)*ZP(gt2,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(Ye(j1,j2))*Conjg(ZE(gt1,3 + j1))*dZV(gt3,j2)*ZP(gt2,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(ZE(gt1,3 + j1))*Conjg(Te(j1,j2))*dZP(gt2,1)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(Ye(j1,j2))*Conjg(ZE(gt1,3 + j1))*dZP(gt2,2)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dTe(j1,j2))*Conjg(ZE(gt1,3 + j1))*ZP(gt2,1)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Conjg(dZE(gt1,3 + j1))*Conjg(Te(j1,j2))*ZP(gt2,1)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(dZE(gt1,3 + j1))*Conjg(Ye(j1,j2))*ZP(gt2,2)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+Mu*Conjg(dYe(j1,j2))*Conjg(ZE(gt1,3 + j1))*ZP(gt2,2)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
res = res+dMu*Conjg(Ye(j1,j2))*Conjg(ZE(gt1,3 + j1))*ZP(gt2,2)*ZV(gt3,j2)
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j3))*Conjg(ZE(gt1,j2))*dZV(gt3,j3)*Ye(j1,j2)*ZP(gt2,1))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j3))*Conjg(ZE(gt1,j2))*dZP(gt2,1)*Ye(j1,j2)*ZV(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(Ye(j1,j3))*Conjg(ZE(gt1,j2))*dYe(j1,j2)*ZP(gt2,1)*ZV(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZE(gt1,j2))*Conjg(Ye(j1,j3))*Ye(j1,j2)*ZP(gt2,1)*ZV(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dYe(j1,j3))*Conjg(ZE(gt1,j2))*Ye(j1,j2)*ZP(gt2,1)*ZV(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(Ye(j1,j3))*Conjg(ZE(gt1,j2))*Ye(j1,j2)*ZP(gt2,1)*ZV(gt3,j3))/sqrt(2._dp)
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSecHpmcSv  
 
 
Subroutine CT_CouplingSecSecSv(gt1,gt2,gt3,Ye,L1,T1,vd,ZV,ZE,dYe,dL1,dT1,             & 
& dvd,dZV,dZE,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: vd,dvd

Complex(dp), Intent(in) :: Ye(3,3),L1(3,3,3),T1(3,3,3),ZV(3,3),ZE(6,6),dYe(3,3),dL1(3,3,3),dT1(3,3,3),           & 
& dZV(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSecSecSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZE(gt1,3 + j3))*Conjg(T1(j1,j2,j3))*dZV(gt3,j2)*ZE(gt2,j1))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZE(gt1,3 + j3))*Conjg(T1(j1,j2,j3))*dZV(gt3,j1)*ZE(gt2,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(ZE(gt1,3 + j3))*Conjg(T1(j1,j2,j3))*dZE(gt2,j2)*ZV(gt3,j1))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dT1(j1,j2,j3))*Conjg(ZE(gt1,3 + j3))*ZE(gt2,j2)*ZV(gt3,j1))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(Conjg(dZE(gt1,3 + j3))*Conjg(T1(j1,j2,j3))*ZE(gt2,j2)*ZV(gt3,j1))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(ZE(gt1,3 + j3))*Conjg(T1(j1,j2,j3))*dZE(gt2,j1)*ZV(gt3,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dT1(j1,j2,j3))*Conjg(ZE(gt1,3 + j3))*ZE(gt2,j1)*ZV(gt3,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(Conjg(dZE(gt1,3 + j3))*Conjg(T1(j1,j2,j3))*ZE(gt2,j1)*ZV(gt3,j2))/2._dp
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(L1(j3,j1,j4))*Conjg(ZE(gt1,3 + j4))*dZV(gt3,j3)*Ye(j2,j1)*ZE(gt2,3 + j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L1(j1,j4,j3))*Conjg(ZE(gt1,3 + j3))*dZV(gt3,j4)*Ye(j2,j1)*ZE(gt2,3 + j2))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*dZV(gt3,j4)*Ye(j1,j2)*ZE(gt2,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*dZV(gt3,j3)*Ye(j1,j2)*ZE(gt2,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*dZE(gt2,j4)*Ye(j1,j2)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(L1(j3,j1,j4))*Conjg(ZE(gt1,3 + j4))*dZE(gt2,3 + j2)*Ye(j2,j1)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(L1(j3,j1,j4))*Conjg(ZE(gt1,3 + j4))*dYe(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZE(gt1,3 + j4))*Conjg(L1(j3,j1,j4))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dL1(j3,j1,j4))*Conjg(ZE(gt1,3 + j4))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(L1(j3,j1,j4))*Conjg(ZE(gt1,3 + j4))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*dYe(j1,j2)*ZE(gt2,j4)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dZE(gt1,j2))*Conjg(L1(j3,j4,j1))*Ye(j1,j2)*ZE(gt2,j4)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(vd*Conjg(dL1(j3,j4,j1))*Conjg(ZE(gt1,j2))*Ye(j1,j2)*ZE(gt2,j4)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res-(dvd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*Ye(j1,j2)*ZE(gt2,j4)*ZV(gt3,j3))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*dZE(gt2,j3)*Ye(j1,j2)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L1(j1,j4,j3))*Conjg(ZE(gt1,3 + j3))*dZE(gt2,3 + j2)*Ye(j2,j1)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L1(j1,j4,j3))*Conjg(ZE(gt1,3 + j3))*dYe(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZE(gt1,3 + j3))*Conjg(L1(j1,j4,j3))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dL1(j1,j4,j3))*Conjg(ZE(gt1,3 + j3))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(L1(j1,j4,j3))*Conjg(ZE(gt1,3 + j3))*Ye(j2,j1)*ZE(gt2,3 + j2)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*dYe(j1,j2)*ZE(gt2,j3)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dZE(gt1,j2))*Conjg(L1(j3,j4,j1))*Ye(j1,j2)*ZE(gt2,j3)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(vd*Conjg(dL1(j3,j4,j1))*Conjg(ZE(gt1,j2))*Ye(j1,j2)*ZE(gt2,j3)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
Do j4 = 1,3
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
res = res+(dvd*Conjg(L1(j3,j4,j1))*Conjg(ZE(gt1,j2))*Ye(j1,j2)*ZE(gt2,j3)*ZV(gt3,j4))/(2._dp*sqrt(2._dp))
End Do 
End Do 
End Do 
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSecSecSv  
 
 
Subroutine CT_CouplingAhhhVZ(gt1,gt2,g1,g2,ZH,ZA,TW,dg1,dg2,dZH,dZA,dSinTW,           & 
& dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,ZH(2,2),ZA(2,2),TW,dg1,dg2,dZH(2,2),dZA(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhhhVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g2*Cos(TW)*dZH(gt2,1)*ZA(gt1,1))/2._dp
res = res-(g1*dZH(gt2,1)*Sin(TW)*ZA(gt1,1))/2._dp
res = res+(g2*Cos(TW)*dZH(gt2,2)*ZA(gt1,2))/2._dp
res = res+(g1*dZH(gt2,2)*Sin(TW)*ZA(gt1,2))/2._dp
res = res-(g2*Cos(TW)*dZA(gt1,1)*ZH(gt2,1))/2._dp
res = res-(g1*dZA(gt1,1)*Sin(TW)*ZH(gt2,1))/2._dp
res = res-(dSinTW*g1*ZA(gt1,1)*ZH(gt2,1))/2._dp
res = res-(dCosTW*g2*ZA(gt1,1)*ZH(gt2,1))/2._dp
res = res-(dg2*Cos(TW)*ZA(gt1,1)*ZH(gt2,1))/2._dp
res = res-(dg1*Sin(TW)*ZA(gt1,1)*ZH(gt2,1))/2._dp
res = res+(g2*Cos(TW)*dZA(gt1,2)*ZH(gt2,2))/2._dp
res = res+(g1*dZA(gt1,2)*Sin(TW)*ZH(gt2,2))/2._dp
res = res+(dSinTW*g1*ZA(gt1,2)*ZH(gt2,2))/2._dp
res = res+(dCosTW*g2*ZA(gt1,2)*ZH(gt2,2))/2._dp
res = res+(dg2*Cos(TW)*ZA(gt1,2)*ZH(gt2,2))/2._dp
res = res+(dg1*Sin(TW)*ZA(gt1,2)*ZH(gt2,2))/2._dp
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhhhVZ  
 
 
Subroutine CT_CouplingAhHpmcVWm(gt1,gt2,g2,ZA,ZP,dg2,dZA,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,ZA(2,2),ZP(2,2),dg2,dZA(2,2),dZP(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhHpmcVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g2*dZP(gt2,1)*ZA(gt1,1))/2._dp
res = res-(g2*dZP(gt2,2)*ZA(gt1,2))/2._dp
res = res-(g2*dZA(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dg2*ZA(gt1,1)*ZP(gt2,1))/2._dp
res = res-(g2*dZA(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dg2*ZA(gt1,2)*ZP(gt2,2))/2._dp
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhHpmcVWm  
 
 
Subroutine CT_CouplingAhcHpmVWm(gt1,gt2,g2,ZA,ZP,dg2,dZA,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,ZA(2,2),ZP(2,2),dg2,dZA(2,2),dZP(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingAhcHpmVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g2*dZP(gt2,1)*ZA(gt1,1))/2._dp
res = res-(g2*dZP(gt2,2)*ZA(gt1,2))/2._dp
res = res-(g2*dZA(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dg2*ZA(gt1,1)*ZP(gt2,1))/2._dp
res = res-(g2*dZA(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dg2*ZA(gt1,2)*ZP(gt2,2))/2._dp
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingAhcHpmVWm  
 
 
Subroutine CT_CouplinghhHpmcVWm(gt1,gt2,g2,ZH,ZP,dg2,dZH,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,ZH(2,2),ZP(2,2),dg2,dZH(2,2),dZP(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhHpmcVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g2*dZP(gt2,1)*ZH(gt1,1))/2._dp
res = res+(g2*dZP(gt2,2)*ZH(gt1,2))/2._dp
res = res-(g2*dZH(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dg2*ZH(gt1,1)*ZP(gt2,1))/2._dp
res = res+(g2*dZH(gt1,2)*ZP(gt2,2))/2._dp
res = res+(dg2*ZH(gt1,2)*ZP(gt2,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhHpmcVWm  
 
 
Subroutine CT_CouplinghhcHpmVWm(gt1,gt2,g2,ZH,ZP,dg2,dZH,dZP,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,ZH(2,2),ZP(2,2),dg2,dZH(2,2),dZP(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhcHpmVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res+(g2*dZP(gt2,1)*ZH(gt1,1))/2._dp
res = res-(g2*dZP(gt2,2)*ZH(gt1,2))/2._dp
res = res+(g2*dZH(gt1,1)*ZP(gt2,1))/2._dp
res = res+(dg2*ZH(gt1,1)*ZP(gt2,1))/2._dp
res = res-(g2*dZH(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dg2*ZH(gt1,2)*ZP(gt2,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhcHpmVWm  
 
 
Subroutine CT_CouplingHpmcHpmVP(gt1,gt2,g1,g2,ZP,TW,dg1,dg2,dZP,dSinTW,               & 
& dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,ZP(2,2),TW,dg1,dg2,dZP(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingHpmcHpmVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g1*Cos(TW)*dZP(gt2,1)*ZP(gt1,1))/2._dp
res = res-(g2*dZP(gt2,1)*Sin(TW)*ZP(gt1,1))/2._dp
res = res-(g1*Cos(TW)*dZP(gt2,2)*ZP(gt1,2))/2._dp
res = res-(g2*dZP(gt2,2)*Sin(TW)*ZP(gt1,2))/2._dp
res = res-(g1*Cos(TW)*dZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(g2*dZP(gt1,1)*Sin(TW)*ZP(gt2,1))/2._dp
res = res-(dCosTW*g1*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dSinTW*g2*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dg1*Cos(TW)*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dg2*Sin(TW)*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(g1*Cos(TW)*dZP(gt1,2)*ZP(gt2,2))/2._dp
res = res-(g2*dZP(gt1,2)*Sin(TW)*ZP(gt2,2))/2._dp
res = res-(dCosTW*g1*ZP(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dSinTW*g2*ZP(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dg1*Cos(TW)*ZP(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dg2*Sin(TW)*ZP(gt1,2)*ZP(gt2,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingHpmcHpmVP  
 
 
Subroutine CT_CouplingHpmcHpmVZ(gt1,gt2,g1,g2,ZP,TW,dg1,dg2,dZP,dSinTW,               & 
& dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,ZP(2,2),TW,dg1,dg2,dZP(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingHpmcHpmVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g2*Cos(TW)*dZP(gt2,1)*ZP(gt1,1))/2._dp
res = res+(g1*dZP(gt2,1)*Sin(TW)*ZP(gt1,1))/2._dp
res = res-(g2*Cos(TW)*dZP(gt2,2)*ZP(gt1,2))/2._dp
res = res+(g1*dZP(gt2,2)*Sin(TW)*ZP(gt1,2))/2._dp
res = res-(g2*Cos(TW)*dZP(gt1,1)*ZP(gt2,1))/2._dp
res = res+(g1*dZP(gt1,1)*Sin(TW)*ZP(gt2,1))/2._dp
res = res+(dSinTW*g1*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dCosTW*g2*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(dg2*Cos(TW)*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res+(dg1*Sin(TW)*ZP(gt1,1)*ZP(gt2,1))/2._dp
res = res-(g2*Cos(TW)*dZP(gt1,2)*ZP(gt2,2))/2._dp
res = res+(g1*dZP(gt1,2)*Sin(TW)*ZP(gt2,2))/2._dp
res = res+(dSinTW*g1*ZP(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dCosTW*g2*ZP(gt1,2)*ZP(gt2,2))/2._dp
res = res-(dg2*Cos(TW)*ZP(gt1,2)*ZP(gt2,2))/2._dp
res = res+(dg1*Sin(TW)*ZP(gt1,2)*ZP(gt2,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingHpmcHpmVZ  
 
 
Subroutine CT_CouplingSdcSdVG(gt1,gt2,g3,dg3,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdcSdVG' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
If ((gt1.eq.gt2)) Then 
res = res+dg3
End If 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdcSdVG  
 
 
Subroutine CT_CouplingSdcSdVP(gt1,gt2,g1,g2,ZD,TW,dg1,dg2,dZD,dSinTW,dCosTW,          & 
& dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: ZD(6,6),dZD(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdcSdVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g1*Conjg(ZD(gt1,j1))*Cos(TW)*dZD(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(ZD(gt1,3 + j1))*Cos(TW)*dZD(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res-(g2*Conjg(ZD(gt1,j1))*dZD(gt2,j1)*Sin(TW))/2._dp
End Do 
Do j1 = 1,3
res = res+(dCosTW*g1*Conjg(ZD(gt1,j1))*ZD(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(dSinTW*g2*Conjg(ZD(gt1,j1))*ZD(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(g1*Conjg(dZD(gt1,j1))*Cos(TW)*ZD(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dg1*Conjg(ZD(gt1,j1))*Cos(TW)*ZD(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(g2*Conjg(dZD(gt1,j1))*Sin(TW)*ZD(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dg2*Conjg(ZD(gt1,j1))*Sin(TW)*ZD(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dCosTW*g1*Conjg(ZD(gt1,3 + j1))*ZD(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(dZD(gt1,3 + j1))*Cos(TW)*ZD(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res-(dg1*Conjg(ZD(gt1,3 + j1))*Cos(TW)*ZD(gt2,3 + j1))/3._dp
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdcSdVP  
 
 
Subroutine CT_CouplingSdcSdVZ(gt1,gt2,g1,g2,ZD,TW,dg1,dg2,dZD,dSinTW,dCosTW,          & 
& dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: ZD(6,6),dZD(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdcSdVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g2*Conjg(ZD(gt1,j1))*Cos(TW)*dZD(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(ZD(gt1,j1))*dZD(gt2,j1)*Sin(TW))/6._dp
End Do 
Do j1 = 1,3
res = res+(g1*Conjg(ZD(gt1,3 + j1))*dZD(gt2,3 + j1)*Sin(TW))/3._dp
End Do 
Do j1 = 1,3
res = res-(dSinTW*g1*Conjg(ZD(gt1,j1))*ZD(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(dCosTW*g2*Conjg(ZD(gt1,j1))*ZD(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g2*Conjg(dZD(gt1,j1))*Cos(TW)*ZD(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dg2*Conjg(ZD(gt1,j1))*Cos(TW)*ZD(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(dZD(gt1,j1))*Sin(TW)*ZD(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(dg1*Conjg(ZD(gt1,j1))*Sin(TW)*ZD(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dSinTW*g1*Conjg(ZD(gt1,3 + j1))*ZD(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(g1*Conjg(dZD(gt1,3 + j1))*Sin(TW)*ZD(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(dg1*Conjg(ZD(gt1,3 + j1))*Sin(TW)*ZD(gt2,3 + j1))/3._dp
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdcSdVZ  
 
 
Subroutine CT_CouplingSdcSucVWm(gt1,gt2,g2,ZD,ZU,dg2,dZD,dZU,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZD(6,6),ZU(6,6),dZD(6,6),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSdcSucVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g2*Conjg(ZD(gt1,j1))*dZU(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(g2*Conjg(dZD(gt1,j1))*ZU(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(dg2*Conjg(ZD(gt1,j1))*ZU(gt2,j1))/sqrt(2._dp)
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSdcSucVWm  
 
 
Subroutine CT_CouplingSecSeVP(gt1,gt2,g1,g2,ZE,TW,dg1,dg2,dZE,dSinTW,dCosTW,          & 
& dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: ZE(6,6),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSecSeVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g1*Conjg(ZE(gt1,j1))*Cos(TW)*dZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(ZE(gt1,3 + j1))*Cos(TW)*dZE(gt2,3 + j1))
End Do 
Do j1 = 1,3
res = res-(g2*Conjg(ZE(gt1,j1))*dZE(gt2,j1)*Sin(TW))/2._dp
End Do 
Do j1 = 1,3
res = res-(dCosTW*g1*Conjg(ZE(gt1,j1))*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dSinTW*g2*Conjg(ZE(gt1,j1))*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(dZE(gt1,j1))*Cos(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dg1*Conjg(ZE(gt1,j1))*Cos(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g2*Conjg(dZE(gt1,j1))*Sin(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dg2*Conjg(ZE(gt1,j1))*Sin(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dCosTW*g1*Conjg(ZE(gt1,3 + j1))*ZE(gt2,3 + j1))
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(dZE(gt1,3 + j1))*Cos(TW)*ZE(gt2,3 + j1))
End Do 
Do j1 = 1,3
res = res-(dg1*Conjg(ZE(gt1,3 + j1))*Cos(TW)*ZE(gt2,3 + j1))
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSecSeVP  
 
 
Subroutine CT_CouplingSecSeVZ(gt1,gt2,g1,g2,ZE,TW,dg1,dg2,dZE,dSinTW,dCosTW,          & 
& dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: ZE(6,6),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSecSeVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res-(g2*Conjg(ZE(gt1,j1))*Cos(TW)*dZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(g1*Conjg(ZE(gt1,j1))*dZE(gt2,j1)*Sin(TW))/2._dp
End Do 
Do j1 = 1,3
res = res+g1*Conjg(ZE(gt1,3 + j1))*dZE(gt2,3 + j1)*Sin(TW)
End Do 
Do j1 = 1,3
res = res+(dSinTW*g1*Conjg(ZE(gt1,j1))*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dCosTW*g2*Conjg(ZE(gt1,j1))*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g2*Conjg(dZE(gt1,j1))*Cos(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(dg2*Conjg(ZE(gt1,j1))*Cos(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(g1*Conjg(dZE(gt1,j1))*Sin(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(dg1*Conjg(ZE(gt1,j1))*Sin(TW)*ZE(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+dSinTW*g1*Conjg(ZE(gt1,3 + j1))*ZE(gt2,3 + j1)
End Do 
Do j1 = 1,3
res = res+g1*Conjg(dZE(gt1,3 + j1))*Sin(TW)*ZE(gt2,3 + j1)
End Do 
Do j1 = 1,3
res = res+dg1*Conjg(ZE(gt1,3 + j1))*Sin(TW)*ZE(gt2,3 + j1)
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSecSeVZ  
 
 
Subroutine CT_CouplingSecSvcVWm(gt1,gt2,g2,ZV,ZE,dg2,dZV,dZE,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZV(3,3),ZE(6,6),dZV(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSecSvcVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g2*Conjg(ZE(gt1,j1))*dZV(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(g2*Conjg(dZE(gt1,j1))*ZV(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(dg2*Conjg(ZE(gt1,j1))*ZV(gt2,j1))/sqrt(2._dp)
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSecSvcVWm  
 
 
Subroutine CT_CouplingSucSuVG(gt1,gt2,g3,dg3,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSucSuVG' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
If ((gt1.eq.gt2)) Then 
res = res+dg3
End If 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSucSuVG  
 
 
Subroutine CT_CouplingSucSuVP(gt1,gt2,g1,g2,ZU,TW,dg1,dg2,dZU,dSinTW,dCosTW,          & 
& dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: ZU(6,6),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSucSuVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g1*Conjg(ZU(gt1,j1))*Cos(TW)*dZU(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(2*g1*Conjg(ZU(gt1,3 + j1))*Cos(TW)*dZU(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(g2*Conjg(ZU(gt1,j1))*dZU(gt2,j1)*Sin(TW))/2._dp
End Do 
Do j1 = 1,3
res = res+(dCosTW*g1*Conjg(ZU(gt1,j1))*ZU(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dSinTW*g2*Conjg(ZU(gt1,j1))*ZU(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(g1*Conjg(dZU(gt1,j1))*Cos(TW)*ZU(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dg1*Conjg(ZU(gt1,j1))*Cos(TW)*ZU(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(g2*Conjg(dZU(gt1,j1))*Sin(TW)*ZU(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(dg2*Conjg(ZU(gt1,j1))*Sin(TW)*ZU(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(2*dCosTW*g1*Conjg(ZU(gt1,3 + j1))*ZU(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(2*g1*Conjg(dZU(gt1,3 + j1))*Cos(TW)*ZU(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(2*dg1*Conjg(ZU(gt1,3 + j1))*Cos(TW)*ZU(gt2,3 + j1))/3._dp
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSucSuVP  
 
 
Subroutine CT_CouplingSucSdVWm(gt1,gt2,g2,ZD,ZU,dg2,dZD,dZU,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZD(6,6),ZU(6,6),dZD(6,6),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSucSdVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g2*Conjg(ZU(gt1,j1))*dZD(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(g2*Conjg(dZU(gt1,j1))*ZD(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(dg2*Conjg(ZU(gt1,j1))*ZD(gt2,j1))/sqrt(2._dp)
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSucSdVWm  
 
 
Subroutine CT_CouplingSucSuVZ(gt1,gt2,g1,g2,ZU,TW,dg1,dg2,dZU,dSinTW,dCosTW,          & 
& dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: ZU(6,6),dZU(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSucSuVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g2*Conjg(ZU(gt1,j1))*Cos(TW)*dZU(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(ZU(gt1,j1))*dZU(gt2,j1)*Sin(TW))/6._dp
End Do 
Do j1 = 1,3
res = res+(-2*g1*Conjg(ZU(gt1,3 + j1))*dZU(gt2,3 + j1)*Sin(TW))/3._dp
End Do 
Do j1 = 1,3
res = res-(dSinTW*g1*Conjg(ZU(gt1,j1))*ZU(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(dCosTW*g2*Conjg(ZU(gt1,j1))*ZU(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(g2*Conjg(dZU(gt1,j1))*Cos(TW)*ZU(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res+(dg2*Conjg(ZU(gt1,j1))*Cos(TW)*ZU(gt2,j1))/2._dp
End Do 
Do j1 = 1,3
res = res-(g1*Conjg(dZU(gt1,j1))*Sin(TW)*ZU(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res-(dg1*Conjg(ZU(gt1,j1))*Sin(TW)*ZU(gt2,j1))/6._dp
End Do 
Do j1 = 1,3
res = res+(-2*dSinTW*g1*Conjg(ZU(gt1,3 + j1))*ZU(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(-2*g1*Conjg(dZU(gt1,3 + j1))*Sin(TW)*ZU(gt2,3 + j1))/3._dp
End Do 
Do j1 = 1,3
res = res+(-2*dg1*Conjg(ZU(gt1,3 + j1))*Sin(TW)*ZU(gt2,3 + j1))/3._dp
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSucSuVZ  
 
 
Subroutine CT_CouplingSvcSeVWm(gt1,gt2,g2,ZV,ZE,dg2,dZV,dZE,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZV(3,3),ZE(6,6),dZV(3,3),dZE(6,6)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSvcSeVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
Do j1 = 1,3
res = res+(g2*Conjg(ZV(gt1,j1))*dZE(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(g2*Conjg(dZV(gt1,j1))*ZE(gt2,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
res = res+(dg2*Conjg(ZV(gt1,j1))*ZE(gt2,j1))/sqrt(2._dp)
End Do 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSvcSeVWm  
 
 
Subroutine CT_CouplingSvcSvVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingSvcSvVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

res = 0._dp 
If ((gt1.eq.gt2)) Then 
res = res+(dSinTW*g1)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
res = res+(dCosTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
res = res+(dg2*Cos(TW))/2._dp
End If 
If ((gt1.eq.gt2)) Then 
res = res+(dg1*Sin(TW))/2._dp
End If 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingSvcSvVZ  
 
 
Subroutine CT_CouplinghhcVWmVWm(gt1,g2,vd,vu,ZH,dg2,dvd,dvu,dZH,res)

Implicit None 

Integer, Intent(in) :: gt1
Real(dp), Intent(in) :: g2,vd,vu,ZH(2,2),dg2,dvd,dvu,dZH(2,2)

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhcVWmVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res+(g2**2*vd*dZH(gt1,1))/2._dp
res = res+(g2**2*vu*dZH(gt1,2))/2._dp
res = res+(dvd*g2**2*ZH(gt1,1))/2._dp
res = res+dg2*g2*vd*ZH(gt1,1)
res = res+(dvu*g2**2*ZH(gt1,2))/2._dp
res = res+dg2*g2*vu*ZH(gt1,2)
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhcVWmVWm  
 
 
Subroutine CT_CouplinghhVZVZ(gt1,g1,g2,vd,vu,ZH,TW,dg1,dg2,dvd,dvu,dZH,               & 
& dSinTW,dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1
Real(dp), Intent(in) :: g1,g2,vd,vu,ZH(2,2),TW,dg1,dg2,dvd,dvu,dZH(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplinghhVZVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res+(g1**2*vd*dZH(gt1,1))/4._dp
res = res+(g2**2*vd*dZH(gt1,1))/4._dp
res = res-(g1**2*vd*Cos(TW)**2*dZH(gt1,1))/4._dp
res = res+(g2**2*vd*Cos(TW)**2*dZH(gt1,1))/4._dp
res = res+(g1**2*vu*dZH(gt1,2))/4._dp
res = res+(g2**2*vu*dZH(gt1,2))/4._dp
res = res-(g1**2*vu*Cos(TW)**2*dZH(gt1,2))/4._dp
res = res+(g2**2*vu*Cos(TW)**2*dZH(gt1,2))/4._dp
res = res+g1*g2*vd*Cos(TW)*dZH(gt1,1)*Sin(TW)
res = res+g1*g2*vu*Cos(TW)*dZH(gt1,2)*Sin(TW)
res = res+(g1**2*vd*dZH(gt1,1)*Sin(TW)**2)/4._dp
res = res-(g2**2*vd*dZH(gt1,1)*Sin(TW)**2)/4._dp
res = res+(g1**2*vu*dZH(gt1,2)*Sin(TW)**2)/4._dp
res = res-(g2**2*vu*dZH(gt1,2)*Sin(TW)**2)/4._dp
res = res+(dvd*g1**2*ZH(gt1,1))/4._dp
res = res+(dvd*g2**2*ZH(gt1,1))/4._dp
res = res+(dg1*g1*vd*ZH(gt1,1))/2._dp
res = res+(dg2*g2*vd*ZH(gt1,1))/2._dp
res = res-(dCosTW*g1**2*vd*Cos(TW)*ZH(gt1,1))/2._dp
res = res+dSinTW*g1*g2*vd*Cos(TW)*ZH(gt1,1)
res = res+(dCosTW*g2**2*vd*Cos(TW)*ZH(gt1,1))/2._dp
res = res-(dvd*g1**2*Cos(TW)**2*ZH(gt1,1))/4._dp
res = res+(dvd*g2**2*Cos(TW)**2*ZH(gt1,1))/4._dp
res = res-(dg1*g1*vd*Cos(TW)**2*ZH(gt1,1))/2._dp
res = res+(dg2*g2*vd*Cos(TW)**2*ZH(gt1,1))/2._dp
res = res+(dSinTW*g1**2*vd*Sin(TW)*ZH(gt1,1))/2._dp
res = res+dCosTW*g1*g2*vd*Sin(TW)*ZH(gt1,1)
res = res-(dSinTW*g2**2*vd*Sin(TW)*ZH(gt1,1))/2._dp
res = res+dvd*g1*g2*Cos(TW)*Sin(TW)*ZH(gt1,1)
res = res+dg2*g1*vd*Cos(TW)*Sin(TW)*ZH(gt1,1)
res = res+dg1*g2*vd*Cos(TW)*Sin(TW)*ZH(gt1,1)
res = res+(dvd*g1**2*Sin(TW)**2*ZH(gt1,1))/4._dp
res = res-(dvd*g2**2*Sin(TW)**2*ZH(gt1,1))/4._dp
res = res+(dg1*g1*vd*Sin(TW)**2*ZH(gt1,1))/2._dp
res = res-(dg2*g2*vd*Sin(TW)**2*ZH(gt1,1))/2._dp
res = res+(dvu*g1**2*ZH(gt1,2))/4._dp
res = res+(dvu*g2**2*ZH(gt1,2))/4._dp
res = res+(dg1*g1*vu*ZH(gt1,2))/2._dp
res = res+(dg2*g2*vu*ZH(gt1,2))/2._dp
res = res-(dCosTW*g1**2*vu*Cos(TW)*ZH(gt1,2))/2._dp
res = res+dSinTW*g1*g2*vu*Cos(TW)*ZH(gt1,2)
res = res+(dCosTW*g2**2*vu*Cos(TW)*ZH(gt1,2))/2._dp
res = res-(dvu*g1**2*Cos(TW)**2*ZH(gt1,2))/4._dp
res = res+(dvu*g2**2*Cos(TW)**2*ZH(gt1,2))/4._dp
res = res-(dg1*g1*vu*Cos(TW)**2*ZH(gt1,2))/2._dp
res = res+(dg2*g2*vu*Cos(TW)**2*ZH(gt1,2))/2._dp
res = res+(dSinTW*g1**2*vu*Sin(TW)*ZH(gt1,2))/2._dp
res = res+dCosTW*g1*g2*vu*Sin(TW)*ZH(gt1,2)
res = res-(dSinTW*g2**2*vu*Sin(TW)*ZH(gt1,2))/2._dp
res = res+dvu*g1*g2*Cos(TW)*Sin(TW)*ZH(gt1,2)
res = res+dg2*g1*vu*Cos(TW)*Sin(TW)*ZH(gt1,2)
res = res+dg1*g2*vu*Cos(TW)*Sin(TW)*ZH(gt1,2)
res = res+(dvu*g1**2*Sin(TW)**2*ZH(gt1,2))/4._dp
res = res-(dvu*g2**2*Sin(TW)**2*ZH(gt1,2))/4._dp
res = res+(dg1*g1*vu*Sin(TW)**2*ZH(gt1,2))/2._dp
res = res-(dg2*g2*vu*Sin(TW)**2*ZH(gt1,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplinghhVZVZ  
 
 
Subroutine CT_CouplingHpmcVWmVP(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,            & 
& dSinTW,dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1
Real(dp), Intent(in) :: g1,g2,vd,vu,ZP(2,2),TW,dg1,dg2,dvd,dvu,dZP(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingHpmcVWmVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g1*g2*vd*Cos(TW)*dZP(gt1,1))/2._dp
res = res+(g1*g2*vu*Cos(TW)*dZP(gt1,2))/2._dp
res = res-(dCosTW*g1*g2*vd*ZP(gt1,1))/2._dp
res = res-(dvd*g1*g2*Cos(TW)*ZP(gt1,1))/2._dp
res = res-(dg2*g1*vd*Cos(TW)*ZP(gt1,1))/2._dp
res = res-(dg1*g2*vd*Cos(TW)*ZP(gt1,1))/2._dp
res = res+(dCosTW*g1*g2*vu*ZP(gt1,2))/2._dp
res = res+(dvu*g1*g2*Cos(TW)*ZP(gt1,2))/2._dp
res = res+(dg2*g1*vu*Cos(TW)*ZP(gt1,2))/2._dp
res = res+(dg1*g2*vu*Cos(TW)*ZP(gt1,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingHpmcVWmVP  
 
 
Subroutine CT_CouplingHpmcVWmVZ(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,            & 
& dSinTW,dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1
Real(dp), Intent(in) :: g1,g2,vd,vu,ZP(2,2),TW,dg1,dg2,dvd,dvu,dZP(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingHpmcVWmVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res+(g1*g2*vd*dZP(gt1,1)*Sin(TW))/2._dp
res = res-(g1*g2*vu*dZP(gt1,2)*Sin(TW))/2._dp
res = res+(dSinTW*g1*g2*vd*ZP(gt1,1))/2._dp
res = res+(dvd*g1*g2*Sin(TW)*ZP(gt1,1))/2._dp
res = res+(dg2*g1*vd*Sin(TW)*ZP(gt1,1))/2._dp
res = res+(dg1*g2*vd*Sin(TW)*ZP(gt1,1))/2._dp
res = res-(dSinTW*g1*g2*vu*ZP(gt1,2))/2._dp
res = res-(dvu*g1*g2*Sin(TW)*ZP(gt1,2))/2._dp
res = res-(dg2*g1*vu*Sin(TW)*ZP(gt1,2))/2._dp
res = res-(dg1*g2*vu*Sin(TW)*ZP(gt1,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingHpmcVWmVZ  
 
 
Subroutine CT_CouplingcHpmVPVWm(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,            & 
& dSinTW,dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1
Real(dp), Intent(in) :: g1,g2,vd,vu,ZP(2,2),TW,dg1,dg2,dvd,dvu,dZP(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcHpmVPVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res-(g1*g2*vd*Cos(TW)*dZP(gt1,1))/2._dp
res = res+(g1*g2*vu*Cos(TW)*dZP(gt1,2))/2._dp
res = res-(dCosTW*g1*g2*vd*ZP(gt1,1))/2._dp
res = res-(dvd*g1*g2*Cos(TW)*ZP(gt1,1))/2._dp
res = res-(dg2*g1*vd*Cos(TW)*ZP(gt1,1))/2._dp
res = res-(dg1*g2*vd*Cos(TW)*ZP(gt1,1))/2._dp
res = res+(dCosTW*g1*g2*vu*ZP(gt1,2))/2._dp
res = res+(dvu*g1*g2*Cos(TW)*ZP(gt1,2))/2._dp
res = res+(dg2*g1*vu*Cos(TW)*ZP(gt1,2))/2._dp
res = res+(dg1*g2*vu*Cos(TW)*ZP(gt1,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcHpmVPVWm  
 
 
Subroutine CT_CouplingcHpmVWmVZ(gt1,g1,g2,vd,vu,ZP,TW,dg1,dg2,dvd,dvu,dZP,            & 
& dSinTW,dCosTW,dTanTW,res)

Implicit None 

Integer, Intent(in) :: gt1
Real(dp), Intent(in) :: g1,g2,vd,vu,ZP(2,2),TW,dg1,dg2,dvd,dvu,dZP(2,2),dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcHpmVWmVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

res = 0._dp 
res = res+(g1*g2*vd*dZP(gt1,1)*Sin(TW))/2._dp
res = res-(g1*g2*vu*dZP(gt1,2)*Sin(TW))/2._dp
res = res+(dSinTW*g1*g2*vd*ZP(gt1,1))/2._dp
res = res+(dvd*g1*g2*Sin(TW)*ZP(gt1,1))/2._dp
res = res+(dg2*g1*vd*Sin(TW)*ZP(gt1,1))/2._dp
res = res+(dg1*g2*vd*Sin(TW)*ZP(gt1,1))/2._dp
res = res-(dSinTW*g1*g2*vu*ZP(gt1,2))/2._dp
res = res-(dvu*g1*g2*Sin(TW)*ZP(gt1,2))/2._dp
res = res-(dg2*g1*vu*Sin(TW)*ZP(gt1,2))/2._dp
res = res-(dg1*g2*vu*Sin(TW)*ZP(gt1,2))/2._dp
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcHpmVWmVZ  
 
 
Subroutine CT_CouplingVGVGVG(g3,dg3,res)

Implicit None 

Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingVGVGVG' 
 
res = 0._dp 
res = res+dg3
res = -(0.,1.)*res 
 
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingVGVGVG  
 
 
Subroutine CT_CouplingcVWmVPVWm(g2,TW,dg2,dSinTW,dCosTW,dTanTW,res)

Implicit None 

Real(dp), Intent(in) :: g2,TW,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcVWmVPVWm' 
 
res = 0._dp 
res = res+dSinTW*g2
res = res+dg2*Sin(TW)
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcVWmVPVWm  
 
 
Subroutine CT_CouplingcVWmVWmVZ(g2,TW,dg2,dSinTW,dCosTW,dTanTW,res)

Implicit None 

Real(dp), Intent(in) :: g2,TW,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: res 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcVWmVWmVZ' 
 
res = 0._dp 
res = res-(dCosTW*g2)
res = res-(dg2*Cos(TW))
If (Real(res,dp).ne.Real(res,dp)) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcVWmVWmVZ  
 
 
Subroutine CT_CouplingcChaChaAh(gt1,gt2,gt3,g2,ZA,UM,UP,dg2,dZA,dUM,dUP,              & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,ZA(2,2),dg2,dZA(2,2)

Complex(dp), Intent(in) :: UM(2,2),UP(2,2),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaChaAh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL-((g2*Conjg(UM(gt2,2))*Conjg(UP(gt1,1))*dZA(gt3,1))/sqrt(2._dp))
resL = resL-((g2*Conjg(UM(gt2,1))*Conjg(UP(gt1,2))*dZA(gt3,2))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUP(gt1,1))*Conjg(UM(gt2,2))*ZA(gt3,1))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUM(gt2,2))*Conjg(UP(gt1,1))*ZA(gt3,1))/sqrt(2._dp))
resL = resL-((dg2*Conjg(UM(gt2,2))*Conjg(UP(gt1,1))*ZA(gt3,1))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUP(gt1,2))*Conjg(UM(gt2,1))*ZA(gt3,2))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUM(gt2,1))*Conjg(UP(gt1,2))*ZA(gt3,2))/sqrt(2._dp))
resL = resL-((dg2*Conjg(UM(gt2,1))*Conjg(UP(gt1,2))*ZA(gt3,2))/sqrt(2._dp))
resR = 0._dp 
resR = resR+(g2*dZA(gt3,1)*UM(gt1,2)*UP(gt2,1))/sqrt(2._dp)
resR = resR+(g2*dZA(gt3,2)*UM(gt1,1)*UP(gt2,2))/sqrt(2._dp)
resR = resR+(g2*dUP(gt2,1)*UM(gt1,2)*ZA(gt3,1))/sqrt(2._dp)
resR = resR+(g2*dUM(gt1,2)*UP(gt2,1)*ZA(gt3,1))/sqrt(2._dp)
resR = resR+(dg2*UM(gt1,2)*UP(gt2,1)*ZA(gt3,1))/sqrt(2._dp)
resR = resR+(g2*dUP(gt2,2)*UM(gt1,1)*ZA(gt3,2))/sqrt(2._dp)
resR = resR+(g2*dUM(gt1,1)*UP(gt2,2)*ZA(gt3,2))/sqrt(2._dp)
resR = resR+(dg2*UM(gt1,1)*UP(gt2,2)*ZA(gt3,2))/sqrt(2._dp)
resL = -(0.,1.)*resL 
 
resR = -(0.,1.)*resR 
 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaChaAh  
 
 
Subroutine CT_CouplingChiChiAh(gt1,gt2,gt3,g1,g2,ZA,ZN,dg1,dg2,dZA,dZN,               & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,ZA(2,2),dg1,dg2,dZA(2,2)

Complex(dp), Intent(in) :: ZN(4,4),dZN(4,4)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiChiAh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL+(g1*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,1))*dZA(gt3,1))/2._dp
resL = resL-(g2*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,2))*dZA(gt3,1))/2._dp
resL = resL+(g1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,3))*dZA(gt3,1))/2._dp
resL = resL-(g2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,3))*dZA(gt3,1))/2._dp
resL = resL-(g1*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,1))*dZA(gt3,2))/2._dp
resL = resL+(g2*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,2))*dZA(gt3,2))/2._dp
resL = resL-(g1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,4))*dZA(gt3,2))/2._dp
resL = resL+(g2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,4))*dZA(gt3,2))/2._dp
resL = resL+(g1*Conjg(dZN(gt2,3))*Conjg(ZN(gt1,1))*ZA(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt2,3))*Conjg(ZN(gt1,2))*ZA(gt3,1))/2._dp
resL = resL+(g1*Conjg(dZN(gt2,1))*Conjg(ZN(gt1,3))*ZA(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt2,2))*Conjg(ZN(gt1,3))*ZA(gt3,1))/2._dp
resL = resL+(g1*Conjg(dZN(gt1,3))*Conjg(ZN(gt2,1))*ZA(gt3,1))/2._dp
resL = resL+(dg1*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,1))*ZA(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt1,3))*Conjg(ZN(gt2,2))*ZA(gt3,1))/2._dp
resL = resL-(dg2*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,2))*ZA(gt3,1))/2._dp
resL = resL+(g1*Conjg(dZN(gt1,1))*Conjg(ZN(gt2,3))*ZA(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt1,2))*Conjg(ZN(gt2,3))*ZA(gt3,1))/2._dp
resL = resL+(dg1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,3))*ZA(gt3,1))/2._dp
resL = resL-(dg2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,3))*ZA(gt3,1))/2._dp
resL = resL-(g1*Conjg(dZN(gt2,4))*Conjg(ZN(gt1,1))*ZA(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt2,4))*Conjg(ZN(gt1,2))*ZA(gt3,2))/2._dp
resL = resL-(g1*Conjg(dZN(gt2,1))*Conjg(ZN(gt1,4))*ZA(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt2,2))*Conjg(ZN(gt1,4))*ZA(gt3,2))/2._dp
resL = resL-(g1*Conjg(dZN(gt1,4))*Conjg(ZN(gt2,1))*ZA(gt3,2))/2._dp
resL = resL-(dg1*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,1))*ZA(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt1,4))*Conjg(ZN(gt2,2))*ZA(gt3,2))/2._dp
resL = resL+(dg2*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,2))*ZA(gt3,2))/2._dp
resL = resL-(g1*Conjg(dZN(gt1,1))*Conjg(ZN(gt2,4))*ZA(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt1,2))*Conjg(ZN(gt2,4))*ZA(gt3,2))/2._dp
resL = resL-(dg1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,4))*ZA(gt3,2))/2._dp
resL = resL+(dg2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,4))*ZA(gt3,2))/2._dp
resR = 0._dp 
resR = resR-(g1*dZN(gt2,3)*ZA(gt3,1)*ZN(gt1,1))/2._dp
resR = resR+(g1*dZN(gt2,4)*ZA(gt3,2)*ZN(gt1,1))/2._dp
resR = resR+(g2*dZN(gt2,3)*ZA(gt3,1)*ZN(gt1,2))/2._dp
resR = resR-(g2*dZN(gt2,4)*ZA(gt3,2)*ZN(gt1,2))/2._dp
resR = resR-(g1*dZN(gt2,1)*ZA(gt3,1)*ZN(gt1,3))/2._dp
resR = resR+(g2*dZN(gt2,2)*ZA(gt3,1)*ZN(gt1,3))/2._dp
resR = resR+(g1*dZN(gt2,1)*ZA(gt3,2)*ZN(gt1,4))/2._dp
resR = resR-(g2*dZN(gt2,2)*ZA(gt3,2)*ZN(gt1,4))/2._dp
resR = resR-(g1*dZN(gt1,3)*ZA(gt3,1)*ZN(gt2,1))/2._dp
resR = resR+(g1*dZN(gt1,4)*ZA(gt3,2)*ZN(gt2,1))/2._dp
resR = resR-(g1*dZA(gt3,1)*ZN(gt1,3)*ZN(gt2,1))/2._dp
resR = resR-(dg1*ZA(gt3,1)*ZN(gt1,3)*ZN(gt2,1))/2._dp
resR = resR+(g1*dZA(gt3,2)*ZN(gt1,4)*ZN(gt2,1))/2._dp
resR = resR+(dg1*ZA(gt3,2)*ZN(gt1,4)*ZN(gt2,1))/2._dp
resR = resR+(g2*dZN(gt1,3)*ZA(gt3,1)*ZN(gt2,2))/2._dp
resR = resR-(g2*dZN(gt1,4)*ZA(gt3,2)*ZN(gt2,2))/2._dp
resR = resR+(g2*dZA(gt3,1)*ZN(gt1,3)*ZN(gt2,2))/2._dp
resR = resR+(dg2*ZA(gt3,1)*ZN(gt1,3)*ZN(gt2,2))/2._dp
resR = resR-(g2*dZA(gt3,2)*ZN(gt1,4)*ZN(gt2,2))/2._dp
resR = resR-(dg2*ZA(gt3,2)*ZN(gt1,4)*ZN(gt2,2))/2._dp
resR = resR-(g1*dZN(gt1,1)*ZA(gt3,1)*ZN(gt2,3))/2._dp
resR = resR+(g2*dZN(gt1,2)*ZA(gt3,1)*ZN(gt2,3))/2._dp
resR = resR-(g1*dZA(gt3,1)*ZN(gt1,1)*ZN(gt2,3))/2._dp
resR = resR-(dg1*ZA(gt3,1)*ZN(gt1,1)*ZN(gt2,3))/2._dp
resR = resR+(g2*dZA(gt3,1)*ZN(gt1,2)*ZN(gt2,3))/2._dp
resR = resR+(dg2*ZA(gt3,1)*ZN(gt1,2)*ZN(gt2,3))/2._dp
resR = resR+(g1*dZN(gt1,1)*ZA(gt3,2)*ZN(gt2,4))/2._dp
resR = resR-(g2*dZN(gt1,2)*ZA(gt3,2)*ZN(gt2,4))/2._dp
resR = resR+(g1*dZA(gt3,2)*ZN(gt1,1)*ZN(gt2,4))/2._dp
resR = resR+(dg1*ZA(gt3,2)*ZN(gt1,1)*ZN(gt2,4))/2._dp
resR = resR-(g2*dZA(gt3,2)*ZN(gt1,2)*ZN(gt2,4))/2._dp
resR = resR-(dg2*ZA(gt3,2)*ZN(gt1,2)*ZN(gt2,4))/2._dp
resL = -(0.,1.)*resL 
 
resR = -(0.,1.)*resR 
 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiChiAh  
 
 
Subroutine CT_CouplingcFdFdAh(gt1,gt2,gt3,Yd,ZA,ZDL,ZDR,dYd,dZA,dZDL,dZDR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZA(2,2),dZA(2,2)

Complex(dp), Intent(in) :: Yd(3,3),ZDL(3,3),ZDR(3,3),dYd(3,3),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFdAh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZDL(gt2,j2))*Conjg(ZDR(gt1,j1))*dZA(gt3,1)*Yd(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZDL(gt2,j2))*Conjg(ZDR(gt1,j1))*dYd(j1,j2)*ZA(gt3,1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZDR(gt1,j1))*Conjg(ZDL(gt2,j2))*Yd(j1,j2)*ZA(gt3,1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZDL(gt2,j2))*Conjg(ZDR(gt1,j1))*Yd(j1,j2)*ZA(gt3,1))/sqrt(2._dp)
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yd(j1,j2))*dZDR(gt2,j1)*ZA(gt3,1)*ZDL(gt1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yd(j1,j2))*dZDL(gt1,j2)*ZA(gt3,1)*ZDR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yd(j1,j2))*dZA(gt3,1)*ZDL(gt1,j2)*ZDR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(dYd(j1,j2))*ZA(gt3,1)*ZDL(gt1,j2)*ZDR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
resL = -(0.,1.)*resL 
 
resR = -(0.,1.)*resR 
 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFdAh  
 
 
Subroutine CT_CouplingcFeFeAh(gt1,gt2,gt3,Ye,ZA,ZEL,ZER,dYe,dZA,dZEL,dZER,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZA(2,2),dZA(2,2)

Complex(dp), Intent(in) :: Ye(3,3),ZEL(3,3),ZER(3,3),dYe(3,3),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFeAh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZEL(gt2,j2))*Conjg(ZER(gt1,j1))*dZA(gt3,1)*Ye(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZEL(gt2,j2))*Conjg(ZER(gt1,j1))*dYe(j1,j2)*ZA(gt3,1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZER(gt1,j1))*Conjg(ZEL(gt2,j2))*Ye(j1,j2)*ZA(gt3,1))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZEL(gt2,j2))*Conjg(ZER(gt1,j1))*Ye(j1,j2)*ZA(gt3,1))/sqrt(2._dp)
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Ye(j1,j2))*dZER(gt2,j1)*ZA(gt3,1)*ZEL(gt1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Ye(j1,j2))*dZEL(gt1,j2)*ZA(gt3,1)*ZER(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Ye(j1,j2))*dZA(gt3,1)*ZEL(gt1,j2)*ZER(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(dYe(j1,j2))*ZA(gt3,1)*ZEL(gt1,j2)*ZER(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
resL = -(0.,1.)*resL 
 
resR = -(0.,1.)*resR 
 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFeAh  
 
 
Subroutine CT_CouplingcFuFuAh(gt1,gt2,gt3,Yu,ZA,ZUL,ZUR,dYu,dZA,dZUL,dZUR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZA(2,2),dZA(2,2)

Complex(dp), Intent(in) :: Yu(3,3),ZUL(3,3),ZUR(3,3),dYu(3,3),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFuAh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZUL(gt2,j2))*Conjg(ZUR(gt1,j1))*dZA(gt3,2)*Yu(j1,j2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZUL(gt2,j2))*Conjg(ZUR(gt1,j1))*dYu(j1,j2)*ZA(gt3,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZUR(gt1,j1))*Conjg(ZUL(gt2,j2))*Yu(j1,j2)*ZA(gt3,2))/sqrt(2._dp)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZUL(gt2,j2))*Conjg(ZUR(gt1,j1))*Yu(j1,j2)*ZA(gt3,2))/sqrt(2._dp)
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yu(j1,j2))*dZUR(gt2,j1)*ZA(gt3,2)*ZUL(gt1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yu(j1,j2))*dZUL(gt1,j2)*ZA(gt3,2)*ZUR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yu(j1,j2))*dZA(gt3,2)*ZUL(gt1,j2)*ZUR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(dYu(j1,j2))*ZA(gt3,2)*ZUL(gt1,j2)*ZUR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
resL = -(0.,1.)*resL 
 
resR = -(0.,1.)*resR 
 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFuAh  
 
 
Subroutine CT_CouplingChiChacHpm(gt1,gt2,gt3,g1,g2,ZP,ZN,UM,UP,dg1,dg2,               & 
& dZP,dZN,dUM,dUP,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,ZP(2,2),dg1,dg2,dZP(2,2)

Complex(dp), Intent(in) :: ZN(4,4),UM(2,2),UP(2,2),dZN(4,4),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiChacHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL+(g1*Conjg(UM(gt2,2))*Conjg(ZN(gt1,1))*dZP(gt3,1))/sqrt(2._dp)
resL = resL+(g2*Conjg(UM(gt2,2))*Conjg(ZN(gt1,2))*dZP(gt3,1))/sqrt(2._dp)
resL = resL-(g2*Conjg(UM(gt2,1))*Conjg(ZN(gt1,3))*dZP(gt3,1))
resL = resL-(g2*Conjg(dZN(gt1,3))*Conjg(UM(gt2,1))*ZP(gt3,1))
resL = resL+(g1*Conjg(dZN(gt1,1))*Conjg(UM(gt2,2))*ZP(gt3,1))/sqrt(2._dp)
resL = resL+(g2*Conjg(dZN(gt1,2))*Conjg(UM(gt2,2))*ZP(gt3,1))/sqrt(2._dp)
resL = resL+(g1*Conjg(dUM(gt2,2))*Conjg(ZN(gt1,1))*ZP(gt3,1))/sqrt(2._dp)
resL = resL+(dg1*Conjg(UM(gt2,2))*Conjg(ZN(gt1,1))*ZP(gt3,1))/sqrt(2._dp)
resL = resL+(g2*Conjg(dUM(gt2,2))*Conjg(ZN(gt1,2))*ZP(gt3,1))/sqrt(2._dp)
resL = resL+(dg2*Conjg(UM(gt2,2))*Conjg(ZN(gt1,2))*ZP(gt3,1))/sqrt(2._dp)
resL = resL-(g2*Conjg(dUM(gt2,1))*Conjg(ZN(gt1,3))*ZP(gt3,1))
resL = resL-(dg2*Conjg(UM(gt2,1))*Conjg(ZN(gt1,3))*ZP(gt3,1))
resR = 0._dp 
resR = resR-((g1*dZP(gt3,2)*UP(gt2,2)*ZN(gt1,1))/sqrt(2._dp))
resR = resR-((g2*dZP(gt3,2)*UP(gt2,2)*ZN(gt1,2))/sqrt(2._dp))
resR = resR-(g2*dZP(gt3,2)*UP(gt2,1)*ZN(gt1,4))
resR = resR-(g2*dZN(gt1,4)*UP(gt2,1)*ZP(gt3,2))
resR = resR-((g1*dZN(gt1,1)*UP(gt2,2)*ZP(gt3,2))/sqrt(2._dp))
resR = resR-((g2*dZN(gt1,2)*UP(gt2,2)*ZP(gt3,2))/sqrt(2._dp))
resR = resR-((g1*dUP(gt2,2)*ZN(gt1,1)*ZP(gt3,2))/sqrt(2._dp))
resR = resR-((dg1*UP(gt2,2)*ZN(gt1,1)*ZP(gt3,2))/sqrt(2._dp))
resR = resR-((g2*dUP(gt2,2)*ZN(gt1,2)*ZP(gt3,2))/sqrt(2._dp))
resR = resR-((dg2*UP(gt2,2)*ZN(gt1,2)*ZP(gt3,2))/sqrt(2._dp))
resR = resR-(g2*dUP(gt2,1)*ZN(gt1,4)*ZP(gt3,2))
resR = resR-(dg2*UP(gt2,1)*ZN(gt1,4)*ZP(gt3,2))
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiChacHpm  
 
 
Subroutine CT_CouplingChaFucSd(gt1,gt2,gt3,g2,Yd,Yu,ZD,UM,UP,ZUL,ZUR,dg2,             & 
& dYd,dYu,dZD,dUM,dUP,dZUL,dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Yd(3,3),Yu(3,3),ZD(6,6),UM(2,2),UP(2,2),ZUL(3,3),ZUR(3,3),dYd(3,3),dYu(3,3),          & 
& dZD(6,6),dUM(2,2),dUP(2,2),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChaFucSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(g2*Conjg(UM(gt1,1))*Conjg(ZUL(gt2,j1))*dZD(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(g2*Conjg(dZUL(gt2,j1))*Conjg(UM(gt1,1))*ZD(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(g2*Conjg(dUM(gt1,1))*Conjg(ZUL(gt2,j1))*ZD(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(dg2*Conjg(UM(gt1,1))*Conjg(ZUL(gt2,j1))*ZD(gt3,j1))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(UM(gt1,2))*Conjg(ZUL(gt2,j2))*dZD(gt3,3 + j1)*Yd(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(UM(gt1,2))*Conjg(ZUL(gt2,j2))*dYd(j1,j2)*ZD(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZUL(gt2,j2))*Conjg(UM(gt1,2))*Yd(j1,j2)*ZD(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dUM(gt1,2))*Conjg(ZUL(gt2,j2))*Yd(j1,j2)*ZD(gt3,3 + j1)
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*dZUR(gt2,j1)*UP(gt1,2)*ZD(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*dZD(gt3,j2)*UP(gt1,2)*ZUR(gt2,j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*dUP(gt1,2)*ZD(gt3,j2)*ZUR(gt2,j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dYu(j1,j2))*UP(gt1,2)*ZD(gt3,j2)*ZUR(gt2,j1)
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChaFucSd  
 
 
Subroutine CT_CouplingFvChacSe(gt1,gt2,gt3,g2,Ye,ZE,UM,dg2,dYe,dZE,dUM,               & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Ye(3,3),ZE(6,6),UM(2,2),dYe(3,3),dZE(6,6),dUM(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingFvChacSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL+Conjg(UM(gt2,2))*dZE(gt3,3 + j1)*Ye(j1,gt1)
End Do 
Do j1 = 1,3
resL = resL+Conjg(UM(gt2,2))*dYe(j1,gt1)*ZE(gt3,3 + j1)
End Do 
Do j1 = 1,3
resL = resL+Conjg(dUM(gt2,2))*Ye(j1,gt1)*ZE(gt3,3 + j1)
End Do 
If ((gt1.le.3).And.(gt1.ge.1)) Then 
resL = resL-(g2*Conjg(UM(gt2,1))*dZE(gt3,gt1))
End If 
If ((gt1.le.3).And.(gt1.ge.1)) Then 
resL = resL-(g2*Conjg(dUM(gt2,1))*ZE(gt3,gt1))
End If 
If ((gt1.le.3).And.(gt1.ge.1)) Then 
resL = resL-(dg2*Conjg(UM(gt2,1))*ZE(gt3,gt1))
End If 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingFvChacSe  
 
 
Subroutine CT_CouplingcChaChahh(gt1,gt2,gt3,g2,ZH,UM,UP,dg2,dZH,dUM,dUP,              & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,ZH(2,2),dg2,dZH(2,2)

Complex(dp), Intent(in) :: UM(2,2),UP(2,2),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaChahh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL-((g2*Conjg(UM(gt2,2))*Conjg(UP(gt1,1))*dZH(gt3,1))/sqrt(2._dp))
resL = resL-((g2*Conjg(UM(gt2,1))*Conjg(UP(gt1,2))*dZH(gt3,2))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUP(gt1,1))*Conjg(UM(gt2,2))*ZH(gt3,1))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUM(gt2,2))*Conjg(UP(gt1,1))*ZH(gt3,1))/sqrt(2._dp))
resL = resL-((dg2*Conjg(UM(gt2,2))*Conjg(UP(gt1,1))*ZH(gt3,1))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUP(gt1,2))*Conjg(UM(gt2,1))*ZH(gt3,2))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUM(gt2,1))*Conjg(UP(gt1,2))*ZH(gt3,2))/sqrt(2._dp))
resL = resL-((dg2*Conjg(UM(gt2,1))*Conjg(UP(gt1,2))*ZH(gt3,2))/sqrt(2._dp))
resR = 0._dp 
resR = resR-((g2*dZH(gt3,1)*UM(gt1,2)*UP(gt2,1))/sqrt(2._dp))
resR = resR-((g2*dZH(gt3,2)*UM(gt1,1)*UP(gt2,2))/sqrt(2._dp))
resR = resR-((g2*dUP(gt2,1)*UM(gt1,2)*ZH(gt3,1))/sqrt(2._dp))
resR = resR-((g2*dUM(gt1,2)*UP(gt2,1)*ZH(gt3,1))/sqrt(2._dp))
resR = resR-((dg2*UM(gt1,2)*UP(gt2,1)*ZH(gt3,1))/sqrt(2._dp))
resR = resR-((g2*dUP(gt2,2)*UM(gt1,1)*ZH(gt3,2))/sqrt(2._dp))
resR = resR-((g2*dUM(gt1,1)*UP(gt2,2)*ZH(gt3,2))/sqrt(2._dp))
resR = resR-((dg2*UM(gt1,1)*UP(gt2,2)*ZH(gt3,2))/sqrt(2._dp))
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaChahh  
 
 
Subroutine CT_CouplingcFdChaSu(gt1,gt2,gt3,g2,Yd,Yu,ZU,UM,UP,ZDL,ZDR,dg2,             & 
& dYd,dYu,dZU,dUM,dUP,dZDL,dZDR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Yd(3,3),Yu(3,3),ZU(6,6),UM(2,2),UP(2,2),ZDL(3,3),ZDR(3,3),dYd(3,3),dYu(3,3),          & 
& dZU(6,6),dUM(2,2),dUP(2,2),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdChaSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(UM(gt2,2))*Conjg(ZDR(gt1,j1))*Conjg(ZU(gt3,j2))*dYd(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZU(gt3,j2))*Conjg(UM(gt2,2))*Conjg(ZDR(gt1,j1))*Yd(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZDR(gt1,j1))*Conjg(UM(gt2,2))*Conjg(ZU(gt3,j2))*Yd(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dUM(gt2,2))*Conjg(ZDR(gt1,j1))*Conjg(ZU(gt3,j2))*Yd(j1,j2)
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(g2*Conjg(ZU(gt3,j1))*dZDL(gt1,j1)*UP(gt2,1))
End Do 
Do j1 = 1,3
resR = resR-(g2*Conjg(ZU(gt3,j1))*dUP(gt2,1)*ZDL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(g2*Conjg(dZU(gt3,j1))*UP(gt2,1)*ZDL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(dg2*Conjg(ZU(gt3,j1))*UP(gt2,1)*ZDL(gt1,j1))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*Conjg(ZU(gt3,3 + j1))*dZDL(gt1,j2)*UP(gt2,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*Conjg(ZU(gt3,3 + j1))*dUP(gt2,2)*ZDL(gt1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dZU(gt3,3 + j1))*Conjg(Yu(j1,j2))*UP(gt2,2)*ZDL(gt1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dYu(j1,j2))*Conjg(ZU(gt3,3 + j1))*UP(gt2,2)*ZDL(gt1,j2)
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdChaSu  
 
 
Subroutine CT_CouplingcFeChaSv(gt1,gt2,gt3,g2,Ye,ZV,UM,UP,ZEL,ZER,dg2,dYe,            & 
& dZV,dUM,dUP,dZEL,dZER,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Ye(3,3),ZV(3,3),UM(2,2),UP(2,2),ZEL(3,3),ZER(3,3),dYe(3,3),dZV(3,3),dUM(2,2),         & 
& dUP(2,2),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeChaSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(UM(gt2,2))*Conjg(ZER(gt1,j1))*Conjg(ZV(gt3,j2))*dYe(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZV(gt3,j2))*Conjg(UM(gt2,2))*Conjg(ZER(gt1,j1))*Ye(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZER(gt1,j1))*Conjg(UM(gt2,2))*Conjg(ZV(gt3,j2))*Ye(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dUM(gt2,2))*Conjg(ZER(gt1,j1))*Conjg(ZV(gt3,j2))*Ye(j1,j2)
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(g2*Conjg(ZV(gt3,j1))*dZEL(gt1,j1)*UP(gt2,1))
End Do 
Do j1 = 1,3
resR = resR-(g2*Conjg(ZV(gt3,j1))*dUP(gt2,1)*ZEL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(g2*Conjg(dZV(gt3,j1))*UP(gt2,1)*ZEL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(dg2*Conjg(ZV(gt3,j1))*UP(gt2,1)*ZEL(gt1,j1))
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeChaSv  
 
 
Subroutine CT_CouplingChiChihh(gt1,gt2,gt3,g1,g2,ZH,ZN,dg1,dg2,dZH,dZN,               & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,ZH(2,2),dg1,dg2,dZH(2,2)

Complex(dp), Intent(in) :: ZN(4,4),dZN(4,4)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiChihh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL+(g1*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,1))*dZH(gt3,1))/2._dp
resL = resL-(g2*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,2))*dZH(gt3,1))/2._dp
resL = resL+(g1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,3))*dZH(gt3,1))/2._dp
resL = resL-(g2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,3))*dZH(gt3,1))/2._dp
resL = resL-(g1*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,1))*dZH(gt3,2))/2._dp
resL = resL+(g2*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,2))*dZH(gt3,2))/2._dp
resL = resL-(g1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,4))*dZH(gt3,2))/2._dp
resL = resL+(g2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,4))*dZH(gt3,2))/2._dp
resL = resL+(g1*Conjg(dZN(gt2,3))*Conjg(ZN(gt1,1))*ZH(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt2,3))*Conjg(ZN(gt1,2))*ZH(gt3,1))/2._dp
resL = resL+(g1*Conjg(dZN(gt2,1))*Conjg(ZN(gt1,3))*ZH(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt2,2))*Conjg(ZN(gt1,3))*ZH(gt3,1))/2._dp
resL = resL+(g1*Conjg(dZN(gt1,3))*Conjg(ZN(gt2,1))*ZH(gt3,1))/2._dp
resL = resL+(dg1*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,1))*ZH(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt1,3))*Conjg(ZN(gt2,2))*ZH(gt3,1))/2._dp
resL = resL-(dg2*Conjg(ZN(gt1,3))*Conjg(ZN(gt2,2))*ZH(gt3,1))/2._dp
resL = resL+(g1*Conjg(dZN(gt1,1))*Conjg(ZN(gt2,3))*ZH(gt3,1))/2._dp
resL = resL-(g2*Conjg(dZN(gt1,2))*Conjg(ZN(gt2,3))*ZH(gt3,1))/2._dp
resL = resL+(dg1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,3))*ZH(gt3,1))/2._dp
resL = resL-(dg2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,3))*ZH(gt3,1))/2._dp
resL = resL-(g1*Conjg(dZN(gt2,4))*Conjg(ZN(gt1,1))*ZH(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt2,4))*Conjg(ZN(gt1,2))*ZH(gt3,2))/2._dp
resL = resL-(g1*Conjg(dZN(gt2,1))*Conjg(ZN(gt1,4))*ZH(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt2,2))*Conjg(ZN(gt1,4))*ZH(gt3,2))/2._dp
resL = resL-(g1*Conjg(dZN(gt1,4))*Conjg(ZN(gt2,1))*ZH(gt3,2))/2._dp
resL = resL-(dg1*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,1))*ZH(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt1,4))*Conjg(ZN(gt2,2))*ZH(gt3,2))/2._dp
resL = resL+(dg2*Conjg(ZN(gt1,4))*Conjg(ZN(gt2,2))*ZH(gt3,2))/2._dp
resL = resL-(g1*Conjg(dZN(gt1,1))*Conjg(ZN(gt2,4))*ZH(gt3,2))/2._dp
resL = resL+(g2*Conjg(dZN(gt1,2))*Conjg(ZN(gt2,4))*ZH(gt3,2))/2._dp
resL = resL-(dg1*Conjg(ZN(gt1,1))*Conjg(ZN(gt2,4))*ZH(gt3,2))/2._dp
resL = resL+(dg2*Conjg(ZN(gt1,2))*Conjg(ZN(gt2,4))*ZH(gt3,2))/2._dp
resR = 0._dp 
resR = resR+(g1*dZN(gt2,3)*ZH(gt3,1)*ZN(gt1,1))/2._dp
resR = resR-(g1*dZN(gt2,4)*ZH(gt3,2)*ZN(gt1,1))/2._dp
resR = resR-(g2*dZN(gt2,3)*ZH(gt3,1)*ZN(gt1,2))/2._dp
resR = resR+(g2*dZN(gt2,4)*ZH(gt3,2)*ZN(gt1,2))/2._dp
resR = resR+(g1*dZN(gt2,1)*ZH(gt3,1)*ZN(gt1,3))/2._dp
resR = resR-(g2*dZN(gt2,2)*ZH(gt3,1)*ZN(gt1,3))/2._dp
resR = resR-(g1*dZN(gt2,1)*ZH(gt3,2)*ZN(gt1,4))/2._dp
resR = resR+(g2*dZN(gt2,2)*ZH(gt3,2)*ZN(gt1,4))/2._dp
resR = resR+(g1*dZN(gt1,3)*ZH(gt3,1)*ZN(gt2,1))/2._dp
resR = resR-(g1*dZN(gt1,4)*ZH(gt3,2)*ZN(gt2,1))/2._dp
resR = resR+(g1*dZH(gt3,1)*ZN(gt1,3)*ZN(gt2,1))/2._dp
resR = resR+(dg1*ZH(gt3,1)*ZN(gt1,3)*ZN(gt2,1))/2._dp
resR = resR-(g1*dZH(gt3,2)*ZN(gt1,4)*ZN(gt2,1))/2._dp
resR = resR-(dg1*ZH(gt3,2)*ZN(gt1,4)*ZN(gt2,1))/2._dp
resR = resR-(g2*dZN(gt1,3)*ZH(gt3,1)*ZN(gt2,2))/2._dp
resR = resR+(g2*dZN(gt1,4)*ZH(gt3,2)*ZN(gt2,2))/2._dp
resR = resR-(g2*dZH(gt3,1)*ZN(gt1,3)*ZN(gt2,2))/2._dp
resR = resR-(dg2*ZH(gt3,1)*ZN(gt1,3)*ZN(gt2,2))/2._dp
resR = resR+(g2*dZH(gt3,2)*ZN(gt1,4)*ZN(gt2,2))/2._dp
resR = resR+(dg2*ZH(gt3,2)*ZN(gt1,4)*ZN(gt2,2))/2._dp
resR = resR+(g1*dZN(gt1,1)*ZH(gt3,1)*ZN(gt2,3))/2._dp
resR = resR-(g2*dZN(gt1,2)*ZH(gt3,1)*ZN(gt2,3))/2._dp
resR = resR+(g1*dZH(gt3,1)*ZN(gt1,1)*ZN(gt2,3))/2._dp
resR = resR+(dg1*ZH(gt3,1)*ZN(gt1,1)*ZN(gt2,3))/2._dp
resR = resR-(g2*dZH(gt3,1)*ZN(gt1,2)*ZN(gt2,3))/2._dp
resR = resR-(dg2*ZH(gt3,1)*ZN(gt1,2)*ZN(gt2,3))/2._dp
resR = resR-(g1*dZN(gt1,1)*ZH(gt3,2)*ZN(gt2,4))/2._dp
resR = resR+(g2*dZN(gt1,2)*ZH(gt3,2)*ZN(gt2,4))/2._dp
resR = resR-(g1*dZH(gt3,2)*ZN(gt1,1)*ZN(gt2,4))/2._dp
resR = resR-(dg1*ZH(gt3,2)*ZN(gt1,1)*ZN(gt2,4))/2._dp
resR = resR+(g2*dZH(gt3,2)*ZN(gt1,2)*ZN(gt2,4))/2._dp
resR = resR+(dg2*ZH(gt3,2)*ZN(gt1,2)*ZN(gt2,4))/2._dp
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiChihh  
 
 
Subroutine CT_CouplingChiFdcSd(gt1,gt2,gt3,g1,g2,Yd,ZD,ZN,ZDL,ZDR,dg1,dg2,            & 
& dYd,dZD,dZN,dZDL,dZDR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: Yd(3,3),ZD(6,6),ZN(4,4),ZDL(3,3),ZDR(3,3),dYd(3,3),dZD(6,6),dZN(4,4),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiFdcSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(g1*Conjg(ZDL(gt2,j1))*Conjg(ZN(gt1,1))*dZD(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL+(g2*Conjg(ZDL(gt2,j1))*Conjg(ZN(gt1,2))*dZD(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL-(g1*Conjg(dZN(gt1,1))*Conjg(ZDL(gt2,j1))*ZD(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL+(g2*Conjg(dZN(gt1,2))*Conjg(ZDL(gt2,j1))*ZD(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL-(g1*Conjg(dZDL(gt2,j1))*Conjg(ZN(gt1,1))*ZD(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-(dg1*Conjg(ZDL(gt2,j1))*Conjg(ZN(gt1,1))*ZD(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL+(g2*Conjg(dZDL(gt2,j1))*Conjg(ZN(gt1,2))*ZD(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(dg2*Conjg(ZDL(gt2,j1))*Conjg(ZN(gt1,2))*ZD(gt3,j1))/sqrt(2._dp)
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZDL(gt2,j2))*Conjg(ZN(gt1,3))*dZD(gt3,3 + j1)*Yd(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZDL(gt2,j2))*Conjg(ZN(gt1,3))*dYd(j1,j2)*ZD(gt3,3 + j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZN(gt1,3))*Conjg(ZDL(gt2,j2))*Yd(j1,j2)*ZD(gt3,3 + j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZDL(gt2,j2))*Conjg(ZN(gt1,3))*Yd(j1,j2)*ZD(gt3,3 + j1))
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g1*dZN(gt1,1)*ZD(gt3,3 + j1)*ZDR(gt2,j1))/3._dp
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g1*dZDR(gt2,j1)*ZD(gt3,3 + j1)*ZN(gt1,1))/3._dp
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g1*dZD(gt3,3 + j1)*ZDR(gt2,j1)*ZN(gt1,1))/3._dp
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*dg1*ZD(gt3,3 + j1)*ZDR(gt2,j1)*ZN(gt1,1))/3._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yd(j1,j2))*dZN(gt1,3)*ZD(gt3,j2)*ZDR(gt2,j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yd(j1,j2))*dZDR(gt2,j1)*ZD(gt3,j2)*ZN(gt1,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yd(j1,j2))*dZD(gt3,j2)*ZDR(gt2,j1)*ZN(gt1,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dYd(j1,j2))*ZD(gt3,j2)*ZDR(gt2,j1)*ZN(gt1,3))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiFdcSd  
 
 
Subroutine CT_CouplingChiFecSe(gt1,gt2,gt3,g1,g2,Ye,ZE,ZN,ZEL,ZER,dg1,dg2,            & 
& dYe,dZE,dZN,dZEL,dZER,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: Ye(3,3),ZE(6,6),ZN(4,4),ZEL(3,3),ZER(3,3),dYe(3,3),dZE(6,6),dZN(4,4),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiFecSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL+(g1*Conjg(ZEL(gt2,j1))*Conjg(ZN(gt1,1))*dZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(g2*Conjg(ZEL(gt2,j1))*Conjg(ZN(gt1,2))*dZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(g1*Conjg(dZN(gt1,1))*Conjg(ZEL(gt2,j1))*ZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(g2*Conjg(dZN(gt1,2))*Conjg(ZEL(gt2,j1))*ZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(g1*Conjg(dZEL(gt2,j1))*Conjg(ZN(gt1,1))*ZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(dg1*Conjg(ZEL(gt2,j1))*Conjg(ZN(gt1,1))*ZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(g2*Conjg(dZEL(gt2,j1))*Conjg(ZN(gt1,2))*ZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resL = resL+(dg2*Conjg(ZEL(gt2,j1))*Conjg(ZN(gt1,2))*ZE(gt3,j1))/sqrt(2._dp)
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZEL(gt2,j2))*Conjg(ZN(gt1,3))*dZE(gt3,3 + j1)*Ye(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZEL(gt2,j2))*Conjg(ZN(gt1,3))*dYe(j1,j2)*ZE(gt3,3 + j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZN(gt1,3))*Conjg(ZEL(gt2,j2))*Ye(j1,j2)*ZE(gt3,3 + j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZEL(gt2,j2))*Conjg(ZN(gt1,3))*Ye(j1,j2)*ZE(gt3,3 + j1))
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g1*dZN(gt1,1)*ZE(gt3,3 + j1)*ZER(gt2,j1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g1*dZER(gt2,j1)*ZE(gt3,3 + j1)*ZN(gt1,1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g1*dZE(gt3,3 + j1)*ZER(gt2,j1)*ZN(gt1,1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*dg1*ZE(gt3,3 + j1)*ZER(gt2,j1)*ZN(gt1,1))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Ye(j1,j2))*dZN(gt1,3)*ZE(gt3,j2)*ZER(gt2,j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Ye(j1,j2))*dZER(gt2,j1)*ZE(gt3,j2)*ZN(gt1,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Ye(j1,j2))*dZE(gt3,j2)*ZER(gt2,j1)*ZN(gt1,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dYe(j1,j2))*ZE(gt3,j2)*ZER(gt2,j1)*ZN(gt1,3))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiFecSe  
 
 
Subroutine CT_CouplingChiFucSu(gt1,gt2,gt3,g1,g2,Yu,ZU,ZN,ZUL,ZUR,dg1,dg2,            & 
& dYu,dZU,dZN,dZUL,dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: Yu(3,3),ZU(6,6),ZN(4,4),ZUL(3,3),ZUR(3,3),dYu(3,3),dZU(6,6),dZN(4,4),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiFucSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(g1*Conjg(ZN(gt1,1))*Conjg(ZUL(gt2,j1))*dZU(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((g2*Conjg(ZN(gt1,2))*Conjg(ZUL(gt2,j1))*dZU(gt3,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-(g1*Conjg(dZUL(gt2,j1))*Conjg(ZN(gt1,1))*ZU(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((g2*Conjg(dZUL(gt2,j1))*Conjg(ZN(gt1,2))*ZU(gt3,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-(g1*Conjg(dZN(gt1,1))*Conjg(ZUL(gt2,j1))*ZU(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((g2*Conjg(dZN(gt1,2))*Conjg(ZUL(gt2,j1))*ZU(gt3,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-(dg1*Conjg(ZN(gt1,1))*Conjg(ZUL(gt2,j1))*ZU(gt3,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((dg2*Conjg(ZN(gt1,2))*Conjg(ZUL(gt2,j1))*ZU(gt3,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZN(gt1,4))*Conjg(ZUL(gt2,j2))*dZU(gt3,3 + j1)*Yu(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZN(gt1,4))*Conjg(ZUL(gt2,j2))*dYu(j1,j2)*ZU(gt3,3 + j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZUL(gt2,j2))*Conjg(ZN(gt1,4))*Yu(j1,j2)*ZU(gt3,3 + j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZN(gt1,4))*Conjg(ZUL(gt2,j2))*Yu(j1,j2)*ZU(gt3,3 + j1))
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR+(2*sqrt(2._dp)*g1*dZUR(gt2,j1)*ZN(gt1,1)*ZU(gt3,3 + j1))/3._dp
End Do 
Do j1 = 1,3
resR = resR+(2*sqrt(2._dp)*g1*dZU(gt3,3 + j1)*ZN(gt1,1)*ZUR(gt2,j1))/3._dp
End Do 
Do j1 = 1,3
resR = resR+(2*sqrt(2._dp)*g1*dZN(gt1,1)*ZU(gt3,3 + j1)*ZUR(gt2,j1))/3._dp
End Do 
Do j1 = 1,3
resR = resR+(2*sqrt(2._dp)*dg1*ZN(gt1,1)*ZU(gt3,3 + j1)*ZUR(gt2,j1))/3._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yu(j1,j2))*dZUR(gt2,j1)*ZN(gt1,4)*ZU(gt3,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yu(j1,j2))*dZU(gt3,j2)*ZN(gt1,4)*ZUR(gt2,j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yu(j1,j2))*dZN(gt1,4)*ZU(gt3,j2)*ZUR(gt2,j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dYu(j1,j2))*ZN(gt1,4)*ZU(gt3,j2)*ZUR(gt2,j1))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiFucSu  
 
 
Subroutine CT_CouplingChiFvSv(gt1,gt2,gt3,g1,g2,ZV,ZN,dg1,dg2,dZV,dZN,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: ZV(3,3),ZN(4,4),dZV(3,3),dZN(4,4)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiFvSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR+(g1*Conjg(ZV(gt3,gt2))*dZN(gt1,1))/sqrt(2._dp)
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR-((g2*Conjg(ZV(gt3,gt2))*dZN(gt1,2))/sqrt(2._dp))
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR+(g1*Conjg(dZV(gt3,gt2))*ZN(gt1,1))/sqrt(2._dp)
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR+(dg1*Conjg(ZV(gt3,gt2))*ZN(gt1,1))/sqrt(2._dp)
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR-((g2*Conjg(dZV(gt3,gt2))*ZN(gt1,2))/sqrt(2._dp))
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR-((dg2*Conjg(ZV(gt3,gt2))*ZN(gt1,2))/sqrt(2._dp))
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiFvSv  
 
 
Subroutine CT_CouplingChiFvcSv(gt1,gt2,gt3,g1,g2,ZV,ZN,dg1,dg2,dZV,dZN,               & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: ZV(3,3),ZN(4,4),dZV(3,3),dZN(4,4)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiFvcSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL+(g1*Conjg(ZN(gt1,1))*dZV(gt3,gt2))/sqrt(2._dp)
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL-((g2*Conjg(ZN(gt1,2))*dZV(gt3,gt2))/sqrt(2._dp))
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL+(g1*Conjg(dZN(gt1,1))*ZV(gt3,gt2))/sqrt(2._dp)
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL-((g2*Conjg(dZN(gt1,2))*ZV(gt3,gt2))/sqrt(2._dp))
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL+(dg1*Conjg(ZN(gt1,1))*ZV(gt3,gt2))/sqrt(2._dp)
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL-((dg2*Conjg(ZN(gt1,2))*ZV(gt3,gt2))/sqrt(2._dp))
End If 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiFvcSv  
 
 
Subroutine CT_CouplingcChaChiHpm(gt1,gt2,gt3,g1,g2,ZP,ZN,UM,UP,dg1,dg2,               & 
& dZP,dZN,dUM,dUP,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,ZP(2,2),dg1,dg2,dZP(2,2)

Complex(dp), Intent(in) :: ZN(4,4),UM(2,2),UP(2,2),dZN(4,4),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaChiHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL-((g1*Conjg(UP(gt1,2))*Conjg(ZN(gt2,1))*dZP(gt3,2))/sqrt(2._dp))
resL = resL-((g2*Conjg(UP(gt1,2))*Conjg(ZN(gt2,2))*dZP(gt3,2))/sqrt(2._dp))
resL = resL-(g2*Conjg(UP(gt1,1))*Conjg(ZN(gt2,4))*dZP(gt3,2))
resL = resL-(g2*Conjg(dZN(gt2,4))*Conjg(UP(gt1,1))*ZP(gt3,2))
resL = resL-((g1*Conjg(dZN(gt2,1))*Conjg(UP(gt1,2))*ZP(gt3,2))/sqrt(2._dp))
resL = resL-((g2*Conjg(dZN(gt2,2))*Conjg(UP(gt1,2))*ZP(gt3,2))/sqrt(2._dp))
resL = resL-((g1*Conjg(dUP(gt1,2))*Conjg(ZN(gt2,1))*ZP(gt3,2))/sqrt(2._dp))
resL = resL-((dg1*Conjg(UP(gt1,2))*Conjg(ZN(gt2,1))*ZP(gt3,2))/sqrt(2._dp))
resL = resL-((g2*Conjg(dUP(gt1,2))*Conjg(ZN(gt2,2))*ZP(gt3,2))/sqrt(2._dp))
resL = resL-((dg2*Conjg(UP(gt1,2))*Conjg(ZN(gt2,2))*ZP(gt3,2))/sqrt(2._dp))
resL = resL-(g2*Conjg(dUP(gt1,1))*Conjg(ZN(gt2,4))*ZP(gt3,2))
resL = resL-(dg2*Conjg(UP(gt1,1))*Conjg(ZN(gt2,4))*ZP(gt3,2))
resR = 0._dp 
resR = resR+(g1*dZP(gt3,1)*UM(gt1,2)*ZN(gt2,1))/sqrt(2._dp)
resR = resR+(g2*dZP(gt3,1)*UM(gt1,2)*ZN(gt2,2))/sqrt(2._dp)
resR = resR-(g2*dZP(gt3,1)*UM(gt1,1)*ZN(gt2,3))
resR = resR-(g2*dZN(gt2,3)*UM(gt1,1)*ZP(gt3,1))
resR = resR+(g1*dZN(gt2,1)*UM(gt1,2)*ZP(gt3,1))/sqrt(2._dp)
resR = resR+(g2*dZN(gt2,2)*UM(gt1,2)*ZP(gt3,1))/sqrt(2._dp)
resR = resR+(g1*dUM(gt1,2)*ZN(gt2,1)*ZP(gt3,1))/sqrt(2._dp)
resR = resR+(dg1*UM(gt1,2)*ZN(gt2,1)*ZP(gt3,1))/sqrt(2._dp)
resR = resR+(g2*dUM(gt1,2)*ZN(gt2,2)*ZP(gt3,1))/sqrt(2._dp)
resR = resR+(dg2*UM(gt1,2)*ZN(gt2,2)*ZP(gt3,1))/sqrt(2._dp)
resR = resR-(g2*dUM(gt1,1)*ZN(gt2,3)*ZP(gt3,1))
resR = resR-(dg2*UM(gt1,1)*ZN(gt2,3)*ZP(gt3,1))
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaChiHpm  
 
 
Subroutine CT_CouplingcFdChiSd(gt1,gt2,gt3,g1,g2,Yd,ZD,ZN,ZDL,ZDR,dg1,dg2,            & 
& dYd,dZD,dZN,dZDL,dZDR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: Yd(3,3),ZD(6,6),ZN(4,4),ZDL(3,3),ZDR(3,3),dYd(3,3),dZD(6,6),dZN(4,4),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdChiSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g1*Conjg(dZN(gt2,1))*Conjg(ZD(gt3,3 + j1))*Conjg(ZDR(gt1,j1)))/3._dp
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g1*Conjg(dZDR(gt1,j1))*Conjg(ZD(gt3,3 + j1))*Conjg(ZN(gt2,1)))/3._dp
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g1*Conjg(dZD(gt3,3 + j1))*Conjg(ZDR(gt1,j1))*Conjg(ZN(gt2,1)))/3._dp
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*dg1*Conjg(ZD(gt3,3 + j1))*Conjg(ZDR(gt1,j1))*Conjg(ZN(gt2,1)))/3._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZD(gt3,j2))*Conjg(ZDR(gt1,j1))*Conjg(ZN(gt2,3))*dYd(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZN(gt2,3))*Conjg(ZD(gt3,j2))*Conjg(ZDR(gt1,j1))*Yd(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZDR(gt1,j1))*Conjg(ZD(gt3,j2))*Conjg(ZN(gt2,3))*Yd(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZD(gt3,j2))*Conjg(ZDR(gt1,j1))*Conjg(ZN(gt2,3))*Yd(j1,j2))
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(g1*Conjg(ZD(gt3,j1))*dZN(gt2,1)*ZDL(gt1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR+(g2*Conjg(ZD(gt3,j1))*dZN(gt2,2)*ZDL(gt1,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR-(g1*Conjg(ZD(gt3,j1))*dZDL(gt1,j1)*ZN(gt2,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-(g1*Conjg(dZD(gt3,j1))*ZDL(gt1,j1)*ZN(gt2,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-(dg1*Conjg(ZD(gt3,j1))*ZDL(gt1,j1)*ZN(gt2,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR+(g2*Conjg(ZD(gt3,j1))*dZDL(gt1,j1)*ZN(gt2,2))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(g2*Conjg(dZD(gt3,j1))*ZDL(gt1,j1)*ZN(gt2,2))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(dg2*Conjg(ZD(gt3,j1))*ZDL(gt1,j1)*ZN(gt2,2))/sqrt(2._dp)
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yd(j1,j2))*Conjg(ZD(gt3,3 + j1))*dZN(gt2,3)*ZDL(gt1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yd(j1,j2))*Conjg(ZD(gt3,3 + j1))*dZDL(gt1,j2)*ZN(gt2,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dZD(gt3,3 + j1))*Conjg(Yd(j1,j2))*ZDL(gt1,j2)*ZN(gt2,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dYd(j1,j2))*Conjg(ZD(gt3,3 + j1))*ZDL(gt1,j2)*ZN(gt2,3))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdChiSd  
 
 
Subroutine CT_CouplingcFeChiSe(gt1,gt2,gt3,g1,g2,Ye,ZE,ZN,ZEL,ZER,dg1,dg2,            & 
& dYe,dZE,dZN,dZEL,dZER,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: Ye(3,3),ZE(6,6),ZN(4,4),ZEL(3,3),ZER(3,3),dYe(3,3),dZE(6,6),dZN(4,4),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeChiSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g1*Conjg(dZN(gt2,1))*Conjg(ZE(gt3,3 + j1))*Conjg(ZER(gt1,j1)))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g1*Conjg(dZER(gt1,j1))*Conjg(ZE(gt3,3 + j1))*Conjg(ZN(gt2,1)))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g1*Conjg(dZE(gt3,3 + j1))*Conjg(ZER(gt1,j1))*Conjg(ZN(gt2,1)))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*dg1*Conjg(ZE(gt3,3 + j1))*Conjg(ZER(gt1,j1))*Conjg(ZN(gt2,1)))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZE(gt3,j2))*Conjg(ZER(gt1,j1))*Conjg(ZN(gt2,3))*dYe(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZN(gt2,3))*Conjg(ZE(gt3,j2))*Conjg(ZER(gt1,j1))*Ye(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZER(gt1,j1))*Conjg(ZE(gt3,j2))*Conjg(ZN(gt2,3))*Ye(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZE(gt3,j2))*Conjg(ZER(gt1,j1))*Conjg(ZN(gt2,3))*Ye(j1,j2))
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR+(g1*Conjg(ZE(gt3,j1))*dZN(gt2,1)*ZEL(gt1,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(g2*Conjg(ZE(gt3,j1))*dZN(gt2,2)*ZEL(gt1,j1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(g1*Conjg(ZE(gt3,j1))*dZEL(gt1,j1)*ZN(gt2,1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(g1*Conjg(dZE(gt3,j1))*ZEL(gt1,j1)*ZN(gt2,1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(dg1*Conjg(ZE(gt3,j1))*ZEL(gt1,j1)*ZN(gt2,1))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(g2*Conjg(ZE(gt3,j1))*dZEL(gt1,j1)*ZN(gt2,2))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(g2*Conjg(dZE(gt3,j1))*ZEL(gt1,j1)*ZN(gt2,2))/sqrt(2._dp)
End Do 
Do j1 = 1,3
resR = resR+(dg2*Conjg(ZE(gt3,j1))*ZEL(gt1,j1)*ZN(gt2,2))/sqrt(2._dp)
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Ye(j1,j2))*Conjg(ZE(gt3,3 + j1))*dZN(gt2,3)*ZEL(gt1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Ye(j1,j2))*Conjg(ZE(gt3,3 + j1))*dZEL(gt1,j2)*ZN(gt2,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dZE(gt3,3 + j1))*Conjg(Ye(j1,j2))*ZEL(gt1,j2)*ZN(gt2,3))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dYe(j1,j2))*Conjg(ZE(gt3,3 + j1))*ZEL(gt1,j2)*ZN(gt2,3))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeChiSe  
 
 
Subroutine CT_CouplingcFuChiSu(gt1,gt2,gt3,g1,g2,Yu,ZU,ZN,ZUL,ZUR,dg1,dg2,            & 
& dYu,dZU,dZN,dZUL,dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g1,g2,dg1,dg2

Complex(dp), Intent(in) :: Yu(3,3),ZU(6,6),ZN(4,4),ZUL(3,3),ZUR(3,3),dYu(3,3),dZU(6,6),dZN(4,4),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuChiSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL+(2*sqrt(2._dp)*g1*Conjg(dZUR(gt1,j1))*Conjg(ZN(gt2,1))*Conjg(ZU(gt3,3 + j1)))/3._dp
End Do 
Do j1 = 1,3
resL = resL+(2*sqrt(2._dp)*g1*Conjg(dZU(gt3,3 + j1))*Conjg(ZN(gt2,1))*Conjg(ZUR(gt1,j1)))/3._dp
End Do 
Do j1 = 1,3
resL = resL+(2*sqrt(2._dp)*g1*Conjg(dZN(gt2,1))*Conjg(ZU(gt3,3 + j1))*Conjg(ZUR(gt1,j1)))/3._dp
End Do 
Do j1 = 1,3
resL = resL+(2*sqrt(2._dp)*dg1*Conjg(ZN(gt2,1))*Conjg(ZU(gt3,3 + j1))*Conjg(ZUR(gt1,j1)))/3._dp
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZN(gt2,4))*Conjg(ZU(gt3,j2))*Conjg(ZUR(gt1,j1))*dYu(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZUR(gt1,j1))*Conjg(ZN(gt2,4))*Conjg(ZU(gt3,j2))*Yu(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZU(gt3,j2))*Conjg(ZN(gt2,4))*Conjg(ZUR(gt1,j1))*Yu(j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZN(gt2,4))*Conjg(ZU(gt3,j2))*Conjg(ZUR(gt1,j1))*Yu(j1,j2))
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(g1*Conjg(ZU(gt3,j1))*dZUL(gt1,j1)*ZN(gt2,1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-((g2*Conjg(ZU(gt3,j1))*dZUL(gt1,j1)*ZN(gt2,2))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-(g1*Conjg(ZU(gt3,j1))*dZN(gt2,1)*ZUL(gt1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-((g2*Conjg(ZU(gt3,j1))*dZN(gt2,2)*ZUL(gt1,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-(g1*Conjg(dZU(gt3,j1))*ZN(gt2,1)*ZUL(gt1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-(dg1*Conjg(ZU(gt3,j1))*ZN(gt2,1)*ZUL(gt1,j1))/(3._dp*sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-((g2*Conjg(dZU(gt3,j1))*ZN(gt2,2)*ZUL(gt1,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resR = resR-((dg2*Conjg(ZU(gt3,j1))*ZN(gt2,2)*ZUL(gt1,j1))/sqrt(2._dp))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yu(j1,j2))*Conjg(ZU(gt3,3 + j1))*dZUL(gt1,j2)*ZN(gt2,4))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(Yu(j1,j2))*Conjg(ZU(gt3,3 + j1))*dZN(gt2,4)*ZUL(gt1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dZU(gt3,3 + j1))*Conjg(Yu(j1,j2))*ZN(gt2,4)*ZUL(gt1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dYu(j1,j2))*Conjg(ZU(gt3,3 + j1))*ZN(gt2,4)*ZUL(gt1,j2))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuChiSu  
 
 
Subroutine CT_CouplingFvFdcSd(gt1,gt2,gt3,L2,ZD,ZDL,ZDR,dL2,dZD,dZDL,dZDR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZD(6,6),ZDL(3,3),ZDR(3,3),dL2(3,3,3),dZD(6,6),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingFvFdcSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZDL(gt2,j1))*dZD(gt3,3 + j2)*L2(gt1,j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZDL(gt2,j1))*dL2(gt1,j1,j2)*ZD(gt3,3 + j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZDL(gt2,j1))*L2(gt1,j1,j2)*ZD(gt3,3 + j2))
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L2(gt1,j1,j2))*dZDR(gt2,j2)*ZD(gt3,j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L2(gt1,j1,j2))*dZD(gt3,j1)*ZDR(gt2,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dL2(gt1,j1,j2))*ZD(gt3,j1)*ZDR(gt2,j2))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingFvFdcSd  
 
 
Subroutine CT_CouplingGluFdcSd(gt2,gt3,g3,pG,ZD,ZDL,ZDR,dg3,dpG,dZD,dZDL,             & 
& dZDR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt2,gt3
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(in) :: pG,ZD(6,6),ZDL(3,3),ZDR(3,3),dpG,dZD(6,6),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingGluFdcSd' 
 
If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g3*pG*Conjg(ZDL(gt2,j1))*dZD(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g3*pG*Conjg(dZDL(gt2,j1))*ZD(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*dpG*g3*Conjg(ZDL(gt2,j1))*ZD(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*dg3*pG*Conjg(ZDL(gt2,j1))*ZD(gt3,j1))
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*g3*Conjg(pG)*dZDR(gt2,j1)*ZD(gt3,3 + j1)
End Do 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*g3*Conjg(pG)*dZD(gt3,3 + j1)*ZDR(gt2,j1)
End Do 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*g3*Conjg(dpG)*ZD(gt3,3 + j1)*ZDR(gt2,j1)
End Do 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*dg3*Conjg(pG)*ZD(gt3,3 + j1)*ZDR(gt2,j1)
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingGluFdcSd  
 
 
Subroutine CT_CouplingcFdFdhh(gt1,gt2,gt3,Yd,ZH,ZDL,ZDR,dYd,dZH,dZDL,dZDR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZH(2,2),dZH(2,2)

Complex(dp), Intent(in) :: Yd(3,3),ZDL(3,3),ZDR(3,3),dYd(3,3),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFdhh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(ZDL(gt2,j2))*Conjg(ZDR(gt1,j1))*dZH(gt3,1)*Yd(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(ZDL(gt2,j2))*Conjg(ZDR(gt1,j1))*dYd(j1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(dZDR(gt1,j1))*Conjg(ZDL(gt2,j2))*Yd(j1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(dZDL(gt2,j2))*Conjg(ZDR(gt1,j1))*Yd(j1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yd(j1,j2))*dZH(gt3,1)*ZDL(gt1,j2)*ZDR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yd(j1,j2))*dZDR(gt2,j1)*ZDL(gt1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yd(j1,j2))*dZDL(gt1,j2)*ZDR(gt2,j1)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(dYd(j1,j2))*ZDL(gt1,j2)*ZDR(gt2,j1)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFdhh  
 
 
Subroutine CT_CouplingcFdFdSv(gt1,gt2,gt3,L2,ZV,ZDL,ZDR,dL2,dZV,dZDL,dZDR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZV(3,3),ZDL(3,3),ZDR(3,3),dL2(3,3,3),dZV(3,3),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFdSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZDL(gt2,j2))*Conjg(ZDR(gt1,j3))*Conjg(ZV(gt3,j1))*dL2(j1,j2,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZV(gt3,j1))*Conjg(ZDL(gt2,j2))*Conjg(ZDR(gt1,j3))*L2(j1,j2,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZDR(gt1,j3))*Conjg(ZDL(gt2,j2))*Conjg(ZV(gt3,j1))*L2(j1,j2,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZDL(gt2,j2))*Conjg(ZDR(gt1,j3))*Conjg(ZV(gt3,j1))*L2(j1,j2,j3))
End Do 
End Do 
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFdSv  
 
 
Subroutine CT_CouplingcChaFdcSu(gt1,gt2,gt3,g2,Yd,Yu,ZU,UM,UP,ZDL,ZDR,dg2,            & 
& dYd,dYu,dZU,dUM,dUP,dZDL,dZDR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Yd(3,3),Yu(3,3),ZU(6,6),UM(2,2),UP(2,2),ZDL(3,3),ZDR(3,3),dYd(3,3),dYu(3,3),          & 
& dZU(6,6),dUM(2,2),dUP(2,2),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaFdcSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(g2*Conjg(UP(gt1,1))*Conjg(ZDL(gt2,j1))*dZU(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(g2*Conjg(dZDL(gt2,j1))*Conjg(UP(gt1,1))*ZU(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(g2*Conjg(dUP(gt1,1))*Conjg(ZDL(gt2,j1))*ZU(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(dg2*Conjg(UP(gt1,1))*Conjg(ZDL(gt2,j1))*ZU(gt3,j1))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(UP(gt1,2))*Conjg(ZDL(gt2,j2))*dZU(gt3,3 + j1)*Yu(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(UP(gt1,2))*Conjg(ZDL(gt2,j2))*dYu(j1,j2)*ZU(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZDL(gt2,j2))*Conjg(UP(gt1,2))*Yu(j1,j2)*ZU(gt3,3 + j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dUP(gt1,2))*Conjg(ZDL(gt2,j2))*Yu(j1,j2)*ZU(gt3,3 + j1)
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*dZU(gt3,j2)*UM(gt1,2)*ZDR(gt2,j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*dZDR(gt2,j1)*UM(gt1,2)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*dUM(gt1,2)*ZDR(gt2,j1)*ZU(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dYd(j1,j2))*UM(gt1,2)*ZDR(gt2,j1)*ZU(gt3,j2)
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaFdcSu  
 
 
Subroutine CT_CouplingcFdFdcSv(gt1,gt2,gt3,L2,ZV,ZDL,ZDR,dL2,dZV,dZDL,dZDR,           & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZV(3,3),ZDL(3,3),ZDR(3,3),dL2(3,3,3),dZV(3,3),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFdcSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L2(j1,j2,j3))*dZV(gt3,j1)*ZDL(gt1,j2)*ZDR(gt2,j3))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L2(j1,j2,j3))*dZDR(gt2,j3)*ZDL(gt1,j2)*ZV(gt3,j1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L2(j1,j2,j3))*dZDL(gt1,j2)*ZDR(gt2,j3)*ZV(gt3,j1))
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dL2(j1,j2,j3))*ZDL(gt1,j2)*ZDR(gt2,j3)*ZV(gt3,j1))
End Do 
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFdcSv  
 
 
Subroutine CT_CouplingcFeFdcSu(gt1,gt2,gt3,L2,ZU,ZEL,ZDR,dL2,dZU,dZEL,dZDR,           & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZU(6,6),ZEL(3,3),ZDR(3,3),dL2(3,3,3),dZU(6,6),dZEL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFdcSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*dZU(gt3,j2)*ZDR(gt2,j3)*ZEL(gt1,j1)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*dZEL(gt1,j1)*ZDR(gt2,j3)*ZU(gt3,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*dZDR(gt2,j3)*ZEL(gt1,j1)*ZU(gt3,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dL2(j1,j2,j3))*ZDR(gt2,j3)*ZEL(gt1,j1)*ZU(gt3,j2)
End Do 
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFdcSu  
 
 
Subroutine CT_CouplingcFuFdcHpm(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,dYd,             & 
& dYu,dZP,dZDL,dZDR,dZUL,dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZP(2,2),dZP(2,2)

Complex(dp), Intent(in) :: Yd(3,3),Yu(3,3),ZDL(3,3),ZDR(3,3),ZUL(3,3),ZUR(3,3),dYd(3,3),dYu(3,3),dZDL(3,3),      & 
& dZDR(3,3),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFdcHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZDL(gt2,j2))*Conjg(ZUR(gt1,j1))*dZP(gt3,2)*Yu(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZDL(gt2,j2))*Conjg(ZUR(gt1,j1))*dYu(j1,j2)*ZP(gt3,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZUR(gt1,j1))*Conjg(ZDL(gt2,j2))*Yu(j1,j2)*ZP(gt3,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZDL(gt2,j2))*Conjg(ZUR(gt1,j1))*Yu(j1,j2)*ZP(gt3,2)
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*dZUL(gt1,j2)*ZDR(gt2,j1)*ZP(gt3,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*dZP(gt3,1)*ZDR(gt2,j1)*ZUL(gt1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*dZDR(gt2,j1)*ZP(gt3,1)*ZUL(gt1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dYd(j1,j2))*ZDR(gt2,j1)*ZP(gt3,1)*ZUL(gt1,j2)
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFdcHpm  
 
 
Subroutine CT_CouplingcFuFdcSe(gt1,gt2,gt3,L2,ZE,ZDR,ZUL,dL2,dZE,dZDR,dZUL,           & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZE(6,6),ZDR(3,3),ZUL(3,3),dL2(3,3,3),dZE(6,6),dZDR(3,3),dZUL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFdcSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*dZUL(gt1,j2)*ZDR(gt2,j3)*ZE(gt3,j1)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*dZE(gt3,j1)*ZDR(gt2,j3)*ZUL(gt1,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*dZDR(gt2,j3)*ZE(gt3,j1)*ZUL(gt1,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dL2(j1,j2,j3))*ZDR(gt2,j3)*ZE(gt3,j1)*ZUL(gt1,j2)
End Do 
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFdcSe  
 
 
Subroutine CT_CouplingFeFucSd(gt1,gt2,gt3,L2,ZD,ZEL,ZUL,dL2,dZD,dZEL,dZUL,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZD(6,6),ZEL(3,3),ZUL(3,3),dL2(3,3,3),dZD(6,6),dZEL(3,3),dZUL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingFeFucSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZEL(gt1,j1))*Conjg(ZUL(gt2,j2))*dZD(gt3,3 + j3)*L2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZEL(gt1,j1))*Conjg(ZUL(gt2,j2))*dL2(j1,j2,j3)*ZD(gt3,3 + j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZUL(gt2,j2))*Conjg(ZEL(gt1,j1))*L2(j1,j2,j3)*ZD(gt3,3 + j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZEL(gt1,j1))*Conjg(ZUL(gt2,j2))*L2(j1,j2,j3)*ZD(gt3,3 + j3)
End Do 
End Do 
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingFeFucSd  
 
 
Subroutine CT_CouplingFvFecHpm(gt1,gt2,gt3,Ye,ZP,ZER,dYe,dZP,dZER,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZP(2,2),dZP(2,2)

Complex(dp), Intent(in) :: Ye(3,3),ZER(3,3),dYe(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingFvFecHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
Do j1 = 1,3
resR = resR+Conjg(Ye(j1,gt1))*dZP(gt3,1)*ZER(gt2,j1)
End Do 
Do j1 = 1,3
resR = resR+Conjg(Ye(j1,gt1))*dZER(gt2,j1)*ZP(gt3,1)
End Do 
Do j1 = 1,3
resR = resR+Conjg(dYe(j1,gt1))*ZER(gt2,j1)*ZP(gt3,1)
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingFvFecHpm  
 
 
Subroutine CT_CouplingFvFecSe(gt1,gt2,gt3,L1,ZE,ZEL,ZER,dL1,dZE,dZEL,dZER,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L1(3,3,3),ZE(6,6),ZEL(3,3),ZER(3,3),dL1(3,3,3),dZE(6,6),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingFvFecSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZEL(gt2,j1))*dZE(gt3,3 + j2)*L1(gt1,j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZEL(gt2,j1))*dZE(gt3,3 + j2)*L1(j1,gt1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZEL(gt2,j1))*dL1(gt1,j1,j2)*ZE(gt3,3 + j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZEL(gt2,j1))*dL1(j1,gt1,j2)*ZE(gt3,3 + j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZEL(gt2,j1))*L1(gt1,j1,j2)*ZE(gt3,3 + j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZEL(gt2,j1))*L1(j1,gt1,j2)*ZE(gt3,3 + j2))/2._dp
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L1(gt1,j1,j2))*dZER(gt2,j2)*ZE(gt3,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(L1(j1,gt1,j2))*dZER(gt2,j2)*ZE(gt3,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L1(gt1,j1,j2))*dZE(gt3,j1)*ZER(gt2,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(L1(j1,gt1,j2))*dZE(gt3,j1)*ZER(gt2,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dL1(gt1,j1,j2))*ZE(gt3,j1)*ZER(gt2,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(dL1(j1,gt1,j2))*ZE(gt3,j1)*ZER(gt2,j2))/2._dp
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingFvFecSe  
 
 
Subroutine CT_CouplingcFeFehh(gt1,gt2,gt3,Ye,ZH,ZEL,ZER,dYe,dZH,dZEL,dZER,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZH(2,2),dZH(2,2)

Complex(dp), Intent(in) :: Ye(3,3),ZEL(3,3),ZER(3,3),dYe(3,3),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFehh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(ZEL(gt2,j2))*Conjg(ZER(gt1,j1))*dZH(gt3,1)*Ye(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(ZEL(gt2,j2))*Conjg(ZER(gt1,j1))*dYe(j1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(dZER(gt1,j1))*Conjg(ZEL(gt2,j2))*Ye(j1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(dZEL(gt2,j2))*Conjg(ZER(gt1,j1))*Ye(j1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Ye(j1,j2))*dZH(gt3,1)*ZEL(gt1,j2)*ZER(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Ye(j1,j2))*dZER(gt2,j1)*ZEL(gt1,j2)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Ye(j1,j2))*dZEL(gt1,j2)*ZER(gt2,j1)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(dYe(j1,j2))*ZEL(gt1,j2)*ZER(gt2,j1)*ZH(gt3,1))/sqrt(2._dp))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFehh  
 
 
Subroutine CT_CouplingcFdFeSu(gt1,gt2,gt3,L2,ZU,ZEL,ZDR,dL2,dZU,dZEL,dZDR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZU(6,6),ZEL(3,3),ZDR(3,3),dL2(3,3,3),dZU(6,6),dZEL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFeSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZDR(gt1,j3))*Conjg(ZEL(gt2,j1))*Conjg(ZU(gt3,j2))*dL2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZU(gt3,j2))*Conjg(ZDR(gt1,j3))*Conjg(ZEL(gt2,j1))*L2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZEL(gt2,j1))*Conjg(ZDR(gt1,j3))*Conjg(ZU(gt3,j2))*L2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZDR(gt1,j3))*Conjg(ZEL(gt2,j1))*Conjg(ZU(gt3,j2))*L2(j1,j2,j3)
End Do 
End Do 
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFeSu  
 
 
Subroutine CT_CouplingcFeFeSv(gt1,gt2,gt3,L1,ZV,ZEL,ZER,dL1,dZV,dZEL,dZER,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L1(3,3,3),ZV(3,3),ZEL(3,3),ZER(3,3),dL1(3,3,3),dZV(3,3),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFeSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZEL(gt2,j2))*Conjg(ZER(gt1,j3))*Conjg(ZV(gt3,j1))*dL1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZEL(gt2,j1))*Conjg(ZER(gt1,j3))*Conjg(ZV(gt3,j2))*dL1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZV(gt3,j2))*Conjg(ZEL(gt2,j1))*Conjg(ZER(gt1,j3))*L1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZV(gt3,j1))*Conjg(ZEL(gt2,j2))*Conjg(ZER(gt1,j3))*L1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZER(gt1,j3))*Conjg(ZEL(gt2,j2))*Conjg(ZV(gt3,j1))*L1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZEL(gt2,j2))*Conjg(ZER(gt1,j3))*Conjg(ZV(gt3,j1))*L1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZER(gt1,j3))*Conjg(ZEL(gt2,j1))*Conjg(ZV(gt3,j2))*L1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZEL(gt2,j1))*Conjg(ZER(gt1,j3))*Conjg(ZV(gt3,j2))*L1(j1,j2,j3))/2._dp
End Do 
End Do 
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFeSv  
 
 
Subroutine CT_CouplingcChaFecSv(gt1,gt2,gt3,g2,Ye,ZV,UM,UP,ZEL,ZER,dg2,               & 
& dYe,dZV,dUM,dUP,dZEL,dZER,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Ye(3,3),ZV(3,3),UM(2,2),UP(2,2),ZEL(3,3),ZER(3,3),dYe(3,3),dZV(3,3),dUM(2,2),         & 
& dUP(2,2),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaFecSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(g2*Conjg(UP(gt1,1))*Conjg(ZEL(gt2,j1))*dZV(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(g2*Conjg(dZEL(gt2,j1))*Conjg(UP(gt1,1))*ZV(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(g2*Conjg(dUP(gt1,1))*Conjg(ZEL(gt2,j1))*ZV(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(dg2*Conjg(UP(gt1,1))*Conjg(ZEL(gt2,j1))*ZV(gt3,j1))
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Ye(j1,j2))*dZV(gt3,j2)*UM(gt1,2)*ZER(gt2,j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Ye(j1,j2))*dZER(gt2,j1)*UM(gt1,2)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Ye(j1,j2))*dUM(gt1,2)*ZER(gt2,j1)*ZV(gt3,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dYe(j1,j2))*UM(gt1,2)*ZER(gt2,j1)*ZV(gt3,j2)
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaFecSv  
 
 
Subroutine CT_CouplingcFeFecSv(gt1,gt2,gt3,L1,ZV,ZEL,ZER,dL1,dZV,dZEL,dZER,           & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L1(3,3,3),ZV(3,3),ZEL(3,3),ZER(3,3),dL1(3,3,3),dZV(3,3),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFecSv' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(L1(j1,j2,j3))*dZV(gt3,j2)*ZEL(gt1,j1)*ZER(gt2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L1(j1,j2,j3))*dZV(gt3,j1)*ZEL(gt1,j2)*ZER(gt2,j3))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L1(j1,j2,j3))*dZER(gt2,j3)*ZEL(gt1,j2)*ZV(gt3,j1))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L1(j1,j2,j3))*dZEL(gt1,j2)*ZER(gt2,j3)*ZV(gt3,j1))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dL1(j1,j2,j3))*ZEL(gt1,j2)*ZER(gt2,j3)*ZV(gt3,j1))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(L1(j1,j2,j3))*dZER(gt2,j3)*ZEL(gt1,j1)*ZV(gt3,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(L1(j1,j2,j3))*dZEL(gt1,j1)*ZER(gt2,j3)*ZV(gt3,j2))/2._dp
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(dL1(j1,j2,j3))*ZEL(gt1,j1)*ZER(gt2,j3)*ZV(gt3,j2))/2._dp
End Do 
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFecSv  
 
 
Subroutine CT_CouplingGluFucSu(gt2,gt3,g3,pG,ZU,ZUL,ZUR,dg3,dpG,dZU,dZUL,             & 
& dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt2,gt3
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(in) :: pG,ZU(6,6),ZUL(3,3),ZUR(3,3),dpG,dZU(6,6),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingGluFucSu' 
 
If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g3*pG*Conjg(ZUL(gt2,j1))*dZU(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*g3*pG*Conjg(dZUL(gt2,j1))*ZU(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*dpG*g3*Conjg(ZUL(gt2,j1))*ZU(gt3,j1))
End Do 
Do j1 = 1,3
resL = resL-(sqrt(2._dp)*dg3*pG*Conjg(ZUL(gt2,j1))*ZU(gt3,j1))
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*g3*Conjg(pG)*dZUR(gt2,j1)*ZU(gt3,3 + j1)
End Do 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*g3*Conjg(pG)*dZU(gt3,3 + j1)*ZUR(gt2,j1)
End Do 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*g3*Conjg(dpG)*ZU(gt3,3 + j1)*ZUR(gt2,j1)
End Do 
Do j1 = 1,3
resR = resR+sqrt(2._dp)*dg3*Conjg(pG)*ZU(gt3,3 + j1)*ZUR(gt2,j1)
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingGluFucSu  
 
 
Subroutine CT_CouplingcFuFuhh(gt1,gt2,gt3,Yu,ZH,ZUL,ZUR,dYu,dZH,dZUL,dZUR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZH(2,2),dZH(2,2)

Complex(dp), Intent(in) :: Yu(3,3),ZUL(3,3),ZUR(3,3),dYu(3,3),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFuhh' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(ZUL(gt2,j2))*Conjg(ZUR(gt1,j1))*dZH(gt3,2)*Yu(j1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(ZUL(gt2,j2))*Conjg(ZUR(gt1,j1))*dYu(j1,j2)*ZH(gt3,2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(dZUR(gt1,j1))*Conjg(ZUL(gt2,j2))*Yu(j1,j2)*ZH(gt3,2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-((Conjg(dZUL(gt2,j2))*Conjg(ZUR(gt1,j1))*Yu(j1,j2)*ZH(gt3,2))/sqrt(2._dp))
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yu(j1,j2))*dZUR(gt2,j1)*ZH(gt3,2)*ZUL(gt1,j2))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yu(j1,j2))*dZUL(gt1,j2)*ZH(gt3,2)*ZUR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(Yu(j1,j2))*dZH(gt3,2)*ZUL(gt1,j2)*ZUR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-((Conjg(dYu(j1,j2))*ZH(gt3,2)*ZUL(gt1,j2)*ZUR(gt2,j1))/sqrt(2._dp))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFuhh  
 
 
Subroutine CT_CouplingcFdFuHpm(gt1,gt2,gt3,Yd,Yu,ZP,ZDL,ZDR,ZUL,ZUR,dYd,              & 
& dYu,dZP,dZDL,dZDR,dZUL,dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZP(2,2),dZP(2,2)

Complex(dp), Intent(in) :: Yd(3,3),Yu(3,3),ZDL(3,3),ZDR(3,3),ZUL(3,3),ZUR(3,3),dYd(3,3),dYu(3,3),dZDL(3,3),      & 
& dZDR(3,3),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFuHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZDR(gt1,j1))*Conjg(ZUL(gt2,j2))*dZP(gt3,1)*Yd(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZDR(gt1,j1))*Conjg(ZUL(gt2,j2))*dYd(j1,j2)*ZP(gt3,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZUL(gt2,j2))*Conjg(ZDR(gt1,j1))*Yd(j1,j2)*ZP(gt3,1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZDR(gt1,j1))*Conjg(ZUL(gt2,j2))*Yd(j1,j2)*ZP(gt3,1)
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*dZUR(gt2,j1)*ZDL(gt1,j2)*ZP(gt3,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*dZP(gt3,2)*ZDL(gt1,j2)*ZUR(gt2,j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yu(j1,j2))*dZDL(gt1,j2)*ZP(gt3,2)*ZUR(gt2,j1)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dYu(j1,j2))*ZDL(gt1,j2)*ZP(gt3,2)*ZUR(gt2,j1)
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFuHpm  
 
 
Subroutine CT_CouplingcFdFuSe(gt1,gt2,gt3,L2,ZE,ZDR,ZUL,dL2,dZE,dZDR,dZUL,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZE(6,6),ZDR(3,3),ZUL(3,3),dL2(3,3,3),dZE(6,6),dZDR(3,3),dZUL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFuSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(ZDR(gt1,j3))*Conjg(ZE(gt3,j1))*Conjg(ZUL(gt2,j2))*dL2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZUL(gt2,j2))*Conjg(ZDR(gt1,j3))*Conjg(ZE(gt3,j1))*L2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZE(gt3,j1))*Conjg(ZDR(gt1,j3))*Conjg(ZUL(gt2,j2))*L2(j1,j2,j3)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZDR(gt1,j3))*Conjg(ZE(gt3,j1))*Conjg(ZUL(gt2,j2))*L2(j1,j2,j3)
End Do 
End Do 
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFuSe  
 
 
Subroutine CT_CouplingcFeFvHpm(gt1,gt2,gt3,Ye,ZP,ZER,dYe,dZP,dZER,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: ZP(2,2),dZP(2,2)

Complex(dp), Intent(in) :: Ye(3,3),ZER(3,3),dYe(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFvHpm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL+Conjg(ZER(gt1,j1))*dZP(gt3,1)*Ye(j1,gt2)
End Do 
Do j1 = 1,3
resL = resL+Conjg(ZER(gt1,j1))*dYe(j1,gt2)*ZP(gt3,1)
End Do 
Do j1 = 1,3
resL = resL+Conjg(dZER(gt1,j1))*Ye(j1,gt2)*ZP(gt3,1)
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFvHpm  
 
 
Subroutine CT_CouplingcFdFvSd(gt1,gt2,gt3,L2,ZD,ZDL,ZDR,dL2,dZD,dZDL,dZDR,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZD(6,6),ZDL(3,3),ZDR(3,3),dL2(3,3,3),dZD(6,6),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFvSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZD(gt3,j1))*Conjg(ZDR(gt1,j2))*dL2(gt2,j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZDR(gt1,j2))*Conjg(ZD(gt3,j1))*L2(gt2,j1,j2))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZD(gt3,j1))*Conjg(ZDR(gt1,j2))*L2(gt2,j1,j2))
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L2(gt2,j1,j2))*Conjg(ZD(gt3,3 + j2))*dZDL(gt1,j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dZD(gt3,3 + j2))*Conjg(L2(gt2,j1,j2))*ZDL(gt1,j1))
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dL2(gt2,j1,j2))*Conjg(ZD(gt3,3 + j2))*ZDL(gt1,j1))
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFvSd  
 
 
Subroutine CT_CouplingcChaFvSe(gt1,gt2,gt3,g2,Ye,ZE,UM,dg2,dYe,dZE,dUM,               & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Ye(3,3),ZE(6,6),UM(2,2),dYe(3,3),dZE(6,6),dUM(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaFvSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
Do j1 = 1,3
resR = resR+Conjg(Ye(j1,gt2))*Conjg(ZE(gt3,3 + j1))*dUM(gt1,2)
End Do 
Do j1 = 1,3
resR = resR+Conjg(dZE(gt3,3 + j1))*Conjg(Ye(j1,gt2))*UM(gt1,2)
End Do 
Do j1 = 1,3
resR = resR+Conjg(dYe(j1,gt2))*Conjg(ZE(gt3,3 + j1))*UM(gt1,2)
End Do 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR-(g2*Conjg(ZE(gt3,gt2))*dUM(gt1,1))
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR-(g2*Conjg(dZE(gt3,gt2))*UM(gt1,1))
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resR = resR-(dg2*Conjg(ZE(gt3,gt2))*UM(gt1,1))
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaFvSe  
 
 
Subroutine CT_CouplingcFeFvSe(gt1,gt2,gt3,L1,ZE,ZEL,ZER,dL1,dZE,dZEL,dZER,            & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L1(3,3,3),ZE(6,6),ZEL(3,3),ZER(3,3),dL1(3,3,3),dZE(6,6),dZEL(3,3),dZER(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFvSe' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(ZE(gt3,j1))*Conjg(ZER(gt1,j2))*dL1(gt2,j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(ZE(gt3,j1))*Conjg(ZER(gt1,j2))*dL1(j1,gt2,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZER(gt1,j2))*Conjg(ZE(gt3,j1))*L1(gt2,j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL-(Conjg(dZE(gt3,j1))*Conjg(ZER(gt1,j2))*L1(gt2,j1,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZER(gt1,j2))*Conjg(ZE(gt3,j1))*L1(j1,gt2,j2))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+(Conjg(dZE(gt3,j1))*Conjg(ZER(gt1,j2))*L1(j1,gt2,j2))/2._dp
End Do 
End Do 
resR = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(L1(gt2,j1,j2))*Conjg(ZE(gt3,3 + j2))*dZEL(gt1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(L1(j1,gt2,j2))*Conjg(ZE(gt3,3 + j2))*dZEL(gt1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dZE(gt3,3 + j2))*Conjg(L1(gt2,j1,j2))*ZEL(gt1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(dZE(gt3,3 + j2))*Conjg(L1(j1,gt2,j2))*ZEL(gt1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR-(Conjg(dL1(gt2,j1,j2))*Conjg(ZE(gt3,3 + j2))*ZEL(gt1,j1))/2._dp
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+(Conjg(dL1(j1,gt2,j2))*Conjg(ZE(gt3,3 + j2))*ZEL(gt1,j1))/2._dp
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFvSe  
 
 
Subroutine CT_CouplingcFdGluSd(gt1,gt3,g3,pG,ZD,ZDL,ZDR,dg3,dpG,dZD,dZDL,             & 
& dZDR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt3
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(in) :: pG,ZD(6,6),ZDL(3,3),ZDR(3,3),dpG,dZD(6,6),dZDL(3,3),dZDR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdGluSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*g3*pG*Conjg(dZDR(gt1,j1))*Conjg(ZD(gt3,3 + j1))
End Do 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*g3*pG*Conjg(dZD(gt3,3 + j1))*Conjg(ZDR(gt1,j1))
End Do 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*dpG*g3*Conjg(ZD(gt3,3 + j1))*Conjg(ZDR(gt1,j1))
End Do 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*dg3*pG*Conjg(ZD(gt3,3 + j1))*Conjg(ZDR(gt1,j1))
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g3*Conjg(pG)*Conjg(ZD(gt3,j1))*dZDL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g3*Conjg(pG)*Conjg(dZD(gt3,j1))*ZDL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g3*Conjg(dpG)*Conjg(ZD(gt3,j1))*ZDL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*dg3*Conjg(pG)*Conjg(ZD(gt3,j1))*ZDL(gt1,j1))
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdGluSd  
 
 
Subroutine CT_CouplingcFuGluSu(gt1,gt3,g3,pG,ZU,ZUL,ZUR,dg3,dpG,dZU,dZUL,             & 
& dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt3
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(in) :: pG,ZU(6,6),ZUL(3,3),ZUR(3,3),dpG,dZU(6,6),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuGluSu' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*g3*pG*Conjg(dZUR(gt1,j1))*Conjg(ZU(gt3,3 + j1))
End Do 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*g3*pG*Conjg(dZU(gt3,3 + j1))*Conjg(ZUR(gt1,j1))
End Do 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*dpG*g3*Conjg(ZU(gt3,3 + j1))*Conjg(ZUR(gt1,j1))
End Do 
Do j1 = 1,3
resL = resL+sqrt(2._dp)*dg3*pG*Conjg(ZU(gt3,3 + j1))*Conjg(ZUR(gt1,j1))
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g3*Conjg(pG)*Conjg(ZU(gt3,j1))*dZUL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g3*Conjg(pG)*Conjg(dZU(gt3,j1))*ZUL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*g3*Conjg(dpG)*Conjg(ZU(gt3,j1))*ZUL(gt1,j1))
End Do 
Do j1 = 1,3
resR = resR-(sqrt(2._dp)*dg3*Conjg(pG)*Conjg(ZU(gt3,j1))*ZUL(gt1,j1))
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuGluSu  
 
 
Subroutine CT_CouplingcChacFuSd(gt1,gt2,gt3,g2,Yd,Yu,ZD,UM,UP,ZUL,ZUR,dg2,            & 
& dYd,dYu,dZD,dUM,dUP,dZUL,dZUR,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: Yd(3,3),Yu(3,3),ZD(6,6),UM(2,2),UP(2,2),ZUL(3,3),ZUR(3,3),dYd(3,3),dYu(3,3),          & 
& dZD(6,6),dUM(2,2),dUP(2,2),dZUL(3,3),dZUR(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChacFuSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(UP(gt1,2))*Conjg(ZD(gt3,j2))*Conjg(ZUR(gt2,j1))*dYu(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZUR(gt2,j1))*Conjg(UP(gt1,2))*Conjg(ZD(gt3,j2))*Yu(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dZD(gt3,j2))*Conjg(UP(gt1,2))*Conjg(ZUR(gt2,j1))*Yu(j1,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resL = resL+Conjg(dUP(gt1,2))*Conjg(ZD(gt3,j2))*Conjg(ZUR(gt2,j1))*Yu(j1,j2)
End Do 
End Do 
resR = 0._dp 
Do j1 = 1,3
resR = resR-(g2*Conjg(ZD(gt3,j1))*dZUL(gt2,j1)*UM(gt1,1))
End Do 
Do j1 = 1,3
resR = resR-(g2*Conjg(ZD(gt3,j1))*dUM(gt1,1)*ZUL(gt2,j1))
End Do 
Do j1 = 1,3
resR = resR-(g2*Conjg(dZD(gt3,j1))*UM(gt1,1)*ZUL(gt2,j1))
End Do 
Do j1 = 1,3
resR = resR-(dg2*Conjg(ZD(gt3,j1))*UM(gt1,1)*ZUL(gt2,j1))
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*Conjg(ZD(gt3,3 + j1))*dZUL(gt2,j2)*UM(gt1,2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(Yd(j1,j2))*Conjg(ZD(gt3,3 + j1))*dUM(gt1,2)*ZUL(gt2,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dZD(gt3,3 + j1))*Conjg(Yd(j1,j2))*UM(gt1,2)*ZUL(gt2,j2)
End Do 
End Do 
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dYd(j1,j2))*Conjg(ZD(gt3,3 + j1))*UM(gt1,2)*ZUL(gt2,j2)
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChacFuSd  
 
 
Subroutine CT_CouplingcFecFuSd(gt1,gt2,gt3,L2,ZD,ZEL,ZUL,dL2,dZD,dZEL,dZUL,           & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2,gt3
Complex(dp), Intent(in) :: L2(3,3,3),ZD(6,6),ZEL(3,3),ZUL(3,3),dL2(3,3,3),dZD(6,6),dZEL(3,3),dZUL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFecFuSd' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

If ((gt3.Lt.1).Or.(gt3.Gt.6)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt3 out of range', gt3 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt3 out of range', gt3 
  Call TerminateProgram 
End If 

resL = 0._dp 
resR = 0._dp 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*Conjg(ZD(gt3,3 + j3))*dZUL(gt2,j2)*ZEL(gt1,j1)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(L2(j1,j2,j3))*Conjg(ZD(gt3,3 + j3))*dZEL(gt1,j1)*ZUL(gt2,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dZD(gt3,3 + j3))*Conjg(L2(j1,j2,j3))*ZEL(gt1,j1)*ZUL(gt2,j2)
End Do 
End Do 
End Do 
Do j3 = 1,3
Do j2 = 1,3
Do j1 = 1,3
resR = resR+Conjg(dL2(j1,j2,j3))*Conjg(ZD(gt3,3 + j3))*ZEL(gt1,j1)*ZUL(gt2,j2)
End Do 
End Do 
End Do 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFecFuSd  
 
 
Subroutine CT_CouplingChiChacVWm(gt1,gt2,g2,ZN,UM,UP,dg2,dZN,dUM,dUP,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZN(4,4),UM(2,2),UP(2,2),dZN(4,4),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiChacVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL-(g2*Conjg(UM(gt2,1))*dZN(gt1,2))
resL = resL-((g2*Conjg(UM(gt2,2))*dZN(gt1,3))/sqrt(2._dp))
resL = resL-(g2*Conjg(dUM(gt2,1))*ZN(gt1,2))
resL = resL-(dg2*Conjg(UM(gt2,1))*ZN(gt1,2))
resL = resL-((g2*Conjg(dUM(gt2,2))*ZN(gt1,3))/sqrt(2._dp))
resL = resL-((dg2*Conjg(UM(gt2,2))*ZN(gt1,3))/sqrt(2._dp))
resR = 0._dp 
resR = resR-(g2*Conjg(ZN(gt1,2))*dUP(gt2,1))
resR = resR+(g2*Conjg(ZN(gt1,4))*dUP(gt2,2))/sqrt(2._dp)
resR = resR-(g2*Conjg(dZN(gt1,2))*UP(gt2,1))
resR = resR-(dg2*Conjg(ZN(gt1,2))*UP(gt2,1))
resR = resR+(g2*Conjg(dZN(gt1,4))*UP(gt2,2))/sqrt(2._dp)
resR = resR+(dg2*Conjg(ZN(gt1,4))*UP(gt2,2))/sqrt(2._dp)
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiChacVWm  
 
 
Subroutine CT_CouplingcChaChaVP(gt1,gt2,g1,g2,UM,UP,TW,dg1,dg2,dUM,dUP,               & 
& dSinTW,dCosTW,dTanTW,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: UM(2,2),UP(2,2),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaChaVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL+(g1*Conjg(UM(gt2,2))*Cos(TW)*dUM(gt1,2))/2._dp
resL = resL+g2*Conjg(UM(gt2,1))*dUM(gt1,1)*Sin(TW)
resL = resL+(g2*Conjg(UM(gt2,2))*dUM(gt1,2)*Sin(TW))/2._dp
resL = resL+dSinTW*g2*Conjg(UM(gt2,1))*UM(gt1,1)
resL = resL+g2*Conjg(dUM(gt2,1))*Sin(TW)*UM(gt1,1)
resL = resL+dg2*Conjg(UM(gt2,1))*Sin(TW)*UM(gt1,1)
resL = resL+(dCosTW*g1*Conjg(UM(gt2,2))*UM(gt1,2))/2._dp
resL = resL+(dSinTW*g2*Conjg(UM(gt2,2))*UM(gt1,2))/2._dp
resL = resL+(g1*Conjg(dUM(gt2,2))*Cos(TW)*UM(gt1,2))/2._dp
resL = resL+(dg1*Conjg(UM(gt2,2))*Cos(TW)*UM(gt1,2))/2._dp
resL = resL+(g2*Conjg(dUM(gt2,2))*Sin(TW)*UM(gt1,2))/2._dp
resL = resL+(dg2*Conjg(UM(gt2,2))*Sin(TW)*UM(gt1,2))/2._dp
resR = 0._dp 
resR = resR+(g1*Conjg(UP(gt1,2))*Cos(TW)*dUP(gt2,2))/2._dp
resR = resR+g2*Conjg(UP(gt1,1))*dUP(gt2,1)*Sin(TW)
resR = resR+(g2*Conjg(UP(gt1,2))*dUP(gt2,2)*Sin(TW))/2._dp
resR = resR+dSinTW*g2*Conjg(UP(gt1,1))*UP(gt2,1)
resR = resR+g2*Conjg(dUP(gt1,1))*Sin(TW)*UP(gt2,1)
resR = resR+dg2*Conjg(UP(gt1,1))*Sin(TW)*UP(gt2,1)
resR = resR+(dCosTW*g1*Conjg(UP(gt1,2))*UP(gt2,2))/2._dp
resR = resR+(dSinTW*g2*Conjg(UP(gt1,2))*UP(gt2,2))/2._dp
resR = resR+(g1*Conjg(dUP(gt1,2))*Cos(TW)*UP(gt2,2))/2._dp
resR = resR+(dg1*Conjg(UP(gt1,2))*Cos(TW)*UP(gt2,2))/2._dp
resR = resR+(g2*Conjg(dUP(gt1,2))*Sin(TW)*UP(gt2,2))/2._dp
resR = resR+(dg2*Conjg(UP(gt1,2))*Sin(TW)*UP(gt2,2))/2._dp
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaChaVP  
 
 
Subroutine CT_CouplingcChaChaVZ(gt1,gt2,g1,g2,UM,UP,TW,dg1,dg2,dUM,dUP,               & 
& dSinTW,dCosTW,dTanTW,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: UM(2,2),UP(2,2),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaChaVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL+g2*Conjg(UM(gt2,1))*Cos(TW)*dUM(gt1,1)
resL = resL+(g2*Conjg(UM(gt2,2))*Cos(TW)*dUM(gt1,2))/2._dp
resL = resL-(g1*Conjg(UM(gt2,2))*dUM(gt1,2)*Sin(TW))/2._dp
resL = resL+dCosTW*g2*Conjg(UM(gt2,1))*UM(gt1,1)
resL = resL+g2*Conjg(dUM(gt2,1))*Cos(TW)*UM(gt1,1)
resL = resL+dg2*Conjg(UM(gt2,1))*Cos(TW)*UM(gt1,1)
resL = resL-(dSinTW*g1*Conjg(UM(gt2,2))*UM(gt1,2))/2._dp
resL = resL+(dCosTW*g2*Conjg(UM(gt2,2))*UM(gt1,2))/2._dp
resL = resL+(g2*Conjg(dUM(gt2,2))*Cos(TW)*UM(gt1,2))/2._dp
resL = resL+(dg2*Conjg(UM(gt2,2))*Cos(TW)*UM(gt1,2))/2._dp
resL = resL-(g1*Conjg(dUM(gt2,2))*Sin(TW)*UM(gt1,2))/2._dp
resL = resL-(dg1*Conjg(UM(gt2,2))*Sin(TW)*UM(gt1,2))/2._dp
resR = 0._dp 
resR = resR+g2*Conjg(UP(gt1,1))*Cos(TW)*dUP(gt2,1)
resR = resR+(g2*Conjg(UP(gt1,2))*Cos(TW)*dUP(gt2,2))/2._dp
resR = resR-(g1*Conjg(UP(gt1,2))*dUP(gt2,2)*Sin(TW))/2._dp
resR = resR+dCosTW*g2*Conjg(UP(gt1,1))*UP(gt2,1)
resR = resR+g2*Conjg(dUP(gt1,1))*Cos(TW)*UP(gt2,1)
resR = resR+dg2*Conjg(UP(gt1,1))*Cos(TW)*UP(gt2,1)
resR = resR-(dSinTW*g1*Conjg(UP(gt1,2))*UP(gt2,2))/2._dp
resR = resR+(dCosTW*g2*Conjg(UP(gt1,2))*UP(gt2,2))/2._dp
resR = resR+(g2*Conjg(dUP(gt1,2))*Cos(TW)*UP(gt2,2))/2._dp
resR = resR+(dg2*Conjg(UP(gt1,2))*Cos(TW)*UP(gt2,2))/2._dp
resR = resR-(g1*Conjg(dUP(gt1,2))*Sin(TW)*UP(gt2,2))/2._dp
resR = resR-(dg1*Conjg(UP(gt1,2))*Sin(TW)*UP(gt2,2))/2._dp
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaChaVZ  
 
 
Subroutine CT_CouplingChiChiVZ(gt1,gt2,g1,g2,ZN,TW,dg1,dg2,dZN,dSinTW,dCosTW,         & 
& dTanTW,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(in) :: ZN(4,4),dZN(4,4)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingChiChiVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL-(g2*Conjg(ZN(gt2,3))*Cos(TW)*dZN(gt1,3))/2._dp
resL = resL+(g2*Conjg(ZN(gt2,4))*Cos(TW)*dZN(gt1,4))/2._dp
resL = resL-(g1*Conjg(ZN(gt2,3))*dZN(gt1,3)*Sin(TW))/2._dp
resL = resL+(g1*Conjg(ZN(gt2,4))*dZN(gt1,4)*Sin(TW))/2._dp
resL = resL-(dSinTW*g1*Conjg(ZN(gt2,3))*ZN(gt1,3))/2._dp
resL = resL-(dCosTW*g2*Conjg(ZN(gt2,3))*ZN(gt1,3))/2._dp
resL = resL-(g2*Conjg(dZN(gt2,3))*Cos(TW)*ZN(gt1,3))/2._dp
resL = resL-(dg2*Conjg(ZN(gt2,3))*Cos(TW)*ZN(gt1,3))/2._dp
resL = resL-(g1*Conjg(dZN(gt2,3))*Sin(TW)*ZN(gt1,3))/2._dp
resL = resL-(dg1*Conjg(ZN(gt2,3))*Sin(TW)*ZN(gt1,3))/2._dp
resL = resL+(dSinTW*g1*Conjg(ZN(gt2,4))*ZN(gt1,4))/2._dp
resL = resL+(dCosTW*g2*Conjg(ZN(gt2,4))*ZN(gt1,4))/2._dp
resL = resL+(g2*Conjg(dZN(gt2,4))*Cos(TW)*ZN(gt1,4))/2._dp
resL = resL+(dg2*Conjg(ZN(gt2,4))*Cos(TW)*ZN(gt1,4))/2._dp
resL = resL+(g1*Conjg(dZN(gt2,4))*Sin(TW)*ZN(gt1,4))/2._dp
resL = resL+(dg1*Conjg(ZN(gt2,4))*Sin(TW)*ZN(gt1,4))/2._dp
resR = 0._dp 
resR = resR+(g2*Conjg(ZN(gt1,3))*Cos(TW)*dZN(gt2,3))/2._dp
resR = resR-(g2*Conjg(ZN(gt1,4))*Cos(TW)*dZN(gt2,4))/2._dp
resR = resR+(g1*Conjg(ZN(gt1,3))*dZN(gt2,3)*Sin(TW))/2._dp
resR = resR-(g1*Conjg(ZN(gt1,4))*dZN(gt2,4)*Sin(TW))/2._dp
resR = resR+(dSinTW*g1*Conjg(ZN(gt1,3))*ZN(gt2,3))/2._dp
resR = resR+(dCosTW*g2*Conjg(ZN(gt1,3))*ZN(gt2,3))/2._dp
resR = resR+(g2*Conjg(dZN(gt1,3))*Cos(TW)*ZN(gt2,3))/2._dp
resR = resR+(dg2*Conjg(ZN(gt1,3))*Cos(TW)*ZN(gt2,3))/2._dp
resR = resR+(g1*Conjg(dZN(gt1,3))*Sin(TW)*ZN(gt2,3))/2._dp
resR = resR+(dg1*Conjg(ZN(gt1,3))*Sin(TW)*ZN(gt2,3))/2._dp
resR = resR-(dSinTW*g1*Conjg(ZN(gt1,4))*ZN(gt2,4))/2._dp
resR = resR-(dCosTW*g2*Conjg(ZN(gt1,4))*ZN(gt2,4))/2._dp
resR = resR-(g2*Conjg(dZN(gt1,4))*Cos(TW)*ZN(gt2,4))/2._dp
resR = resR-(dg2*Conjg(ZN(gt1,4))*Cos(TW)*ZN(gt2,4))/2._dp
resR = resR-(g1*Conjg(dZN(gt1,4))*Sin(TW)*ZN(gt2,4))/2._dp
resR = resR-(dg1*Conjg(ZN(gt1,4))*Sin(TW)*ZN(gt2,4))/2._dp
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingChiChiVZ  
 
 
Subroutine CT_CouplingcChaChiVWm(gt1,gt2,g2,ZN,UM,UP,dg2,dZN,dUM,dUP,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZN(4,4),UM(2,2),UP(2,2),dZN(4,4),dUM(2,2),dUP(2,2)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcChaChiVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.2)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.4)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
resL = resL-(g2*Conjg(ZN(gt2,2))*dUM(gt1,1))
resL = resL-((g2*Conjg(ZN(gt2,3))*dUM(gt1,2))/sqrt(2._dp))
resL = resL-(g2*Conjg(dZN(gt2,2))*UM(gt1,1))
resL = resL-(dg2*Conjg(ZN(gt2,2))*UM(gt1,1))
resL = resL-((g2*Conjg(dZN(gt2,3))*UM(gt1,2))/sqrt(2._dp))
resL = resL-((dg2*Conjg(ZN(gt2,3))*UM(gt1,2))/sqrt(2._dp))
resR = 0._dp 
resR = resR-(g2*Conjg(UP(gt1,1))*dZN(gt2,2))
resR = resR+(g2*Conjg(UP(gt1,2))*dZN(gt2,4))/sqrt(2._dp)
resR = resR-(g2*Conjg(dUP(gt1,1))*ZN(gt2,2))
resR = resR-(dg2*Conjg(UP(gt1,1))*ZN(gt2,2))
resR = resR+(g2*Conjg(dUP(gt1,2))*ZN(gt2,4))/sqrt(2._dp)
resR = resR+(dg2*Conjg(UP(gt1,2))*ZN(gt2,4))/sqrt(2._dp)
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcChaChiVWm  
 
 
Subroutine CT_CouplingcFdFdVG(gt1,gt2,g3,dg3,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFdVG' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL-1._dp*(dg3)
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR-1._dp*(dg3)
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFdVG  
 
 
Subroutine CT_CouplingcFdFdVP(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,          & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFdVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL-(dCosTW*g1)/6._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dSinTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dg1*Cos(TW))/6._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dg2*Sin(TW))/2._dp
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR+(dCosTW*g1)/3._dp
End If 
If ((gt1.eq.gt2)) Then 
resR = resR+(dg1*Cos(TW))/3._dp
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFdVP  
 
 
Subroutine CT_CouplingcFdFdVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,          & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFdVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL+(dSinTW*g1)/6._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dCosTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dg2*Cos(TW))/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dg1*Sin(TW))/6._dp
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR-(dSinTW*g1)/3._dp
End If 
If ((gt1.eq.gt2)) Then 
resR = resR-(dg1*Sin(TW))/3._dp
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFdVZ  
 
 
Subroutine CT_CouplingcFuFdcVWm(gt1,gt2,g2,ZDL,ZUL,dg2,dZDL,dZUL,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZDL(3,3),ZUL(3,3),dZDL(3,3),dZUL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFdcVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-((g2*Conjg(ZDL(gt2,j1))*dZUL(gt1,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((g2*Conjg(dZDL(gt2,j1))*ZUL(gt1,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((dg2*Conjg(ZDL(gt2,j1))*ZUL(gt1,j1))/sqrt(2._dp))
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFdcVWm  
 
 
Subroutine CT_CouplingFvFecVWm(gt1,gt2,g2,ZEL,dg2,dZEL,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZEL(3,3),dZEL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingFvFecVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.le.3).And.(gt1.ge.1)) Then 
resL = resL-((g2*Conjg(dZEL(gt2,gt1)))/sqrt(2._dp))
End If 
If ((gt1.le.3).And.(gt1.ge.1)) Then 
resL = resL-((dg2*Conjg(ZEL(gt2,gt1)))/sqrt(2._dp))
End If 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingFvFecVWm  
 
 
Subroutine CT_CouplingcFeFeVP(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,          & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFeVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL+(dCosTW*g1)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dSinTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dg1*Cos(TW))/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dg2*Sin(TW))/2._dp
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR+dCosTW*g1
End If 
If ((gt1.eq.gt2)) Then 
resR = resR+dg1*Cos(TW)
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFeVP  
 
 
Subroutine CT_CouplingcFeFeVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,          & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFeVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL-(dSinTW*g1)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dCosTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dg2*Cos(TW))/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dg1*Sin(TW))/2._dp
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR-(dSinTW*g1)
End If 
If ((gt1.eq.gt2)) Then 
resR = resR-(dg1*Sin(TW))
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFeVZ  
 
 
Subroutine CT_CouplingcFuFuVG(gt1,gt2,g3,dg3,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFuVG' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL-1._dp*(dg3)
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR-1._dp*(dg3)
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFuVG  
 
 
Subroutine CT_CouplingcFuFuVP(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,          & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFuVP' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL-(dCosTW*g1)/6._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dSinTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dg1*Cos(TW))/6._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dg2*Sin(TW))/2._dp
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR+(-2*dCosTW*g1)/3._dp
End If 
If ((gt1.eq.gt2)) Then 
resR = resR+(-2*dg1*Cos(TW))/3._dp
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFuVP  
 
 
Subroutine CT_CouplingcFdFuVWm(gt1,gt2,g2,ZDL,ZUL,dg2,dZDL,dZUL,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZDL(3,3),ZUL(3,3),dZDL(3,3),dZUL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFdFuVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
Do j1 = 1,3
resL = resL-((g2*Conjg(ZUL(gt2,j1))*dZDL(gt1,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((g2*Conjg(dZUL(gt2,j1))*ZDL(gt1,j1))/sqrt(2._dp))
End Do 
Do j1 = 1,3
resL = resL-((dg2*Conjg(ZUL(gt2,j1))*ZDL(gt1,j1))/sqrt(2._dp))
End Do 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFdFuVWm  
 
 
Subroutine CT_CouplingcFuFuVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,          & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFuFuVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL+(dSinTW*g1)/6._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dCosTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dg2*Cos(TW))/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL+(dg1*Sin(TW))/6._dp
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR+(2*dSinTW*g1)/3._dp
End If 
If ((gt1.eq.gt2)) Then 
resR = resR+(2*dg1*Sin(TW))/3._dp
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFuFuVZ  
 
 
Subroutine CT_CouplingFvFvVZ(gt1,gt2,g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW,           & 
& resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g1,g2,TW,dg1,dg2,dSinTW,dCosTW,dTanTW

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingFvFvVZ' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt1.eq.gt2)) Then 
resL = resL-(dSinTW*g1)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dCosTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dg2*Cos(TW))/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resL = resL-(dg1*Sin(TW))/2._dp
End If 
resR = 0._dp 
If ((gt1.eq.gt2)) Then 
resR = resR+(dSinTW*g1)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resR = resR+(dCosTW*g2)/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resR = resR+(dg2*Cos(TW))/2._dp
End If 
If ((gt1.eq.gt2)) Then 
resR = resR+(dg1*Sin(TW))/2._dp
End If 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingFvFvVZ  
 
 
Subroutine CT_CouplingcFeFvVWm(gt1,gt2,g2,ZEL,dg2,dZEL,resL,resR)

Implicit None 

Integer, Intent(in) :: gt1,gt2
Real(dp), Intent(in) :: g2,dg2

Complex(dp), Intent(in) :: ZEL(3,3),dZEL(3,3)

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingcFeFvVWm' 
 
If ((gt1.Lt.1).Or.(gt1.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt1 out of range', gt1 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt1 out of range', gt1 
  Call TerminateProgram 
End If 

If ((gt2.Lt.1).Or.(gt2.Gt.3)) Then 
  Write (ErrCan,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (ErrCan,*) 'index gt2 out of range', gt2 
  Write (*,*) 'Problem in Subroutine ',NameOfUnit(Iname) 
  Write (*,*) 'index gt2 out of range', gt2 
  Call TerminateProgram 
End If 

resL = 0._dp 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL-((g2*dZEL(gt1,gt2))/sqrt(2._dp))
End If 
If ((gt2.le.3).And.(gt2.ge.1)) Then 
resL = resL-((dg2*ZEL(gt1,gt2))/sqrt(2._dp))
End If 
resR = 0._dp 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingcFeFvVWm  
 
 
Subroutine CT_CouplingGluGluVG(g3,pG,dg3,dpG,resL,resR)

Implicit None 

Real(dp), Intent(in) :: g3,dg3

Complex(dp), Intent(in) :: pG,dpG

Complex(dp), Intent(out) :: resL, resR 
 
Integer :: j1,j2,j3,j4,j5,j6, j7, j8, j9, j10, j11, j12 
Iname = Iname +1 
NameOfUnit(Iname) = 'CT_CouplingGluGluVG' 
 
resL = 0._dp 
resL = resL-(g3*pG*Conjg(dpG))
resL = resL-(dpG*g3*Conjg(pG))
resL = resL-(dg3*pG*Conjg(pG))
resR = 0._dp 
resR = resR-(g3*pG*Conjg(dpG))
resR = resR-(dpG*g3*Conjg(pG))
resR = resR-(dg3*pG*Conjg(pG))
resL = -(0.,1.)*resL 
 
resR = -(0.,1.)*resR 
 
If ((Real(resL,dp).ne.Real(resL,dp)).or.(Real(resR,dp).ne.Real(resR,dp))) Then 
 Write(*,*) "NaN appearing in ",NameOfUnit(Iname) 
 Call TerminateProgram 
End If 


Iname = Iname - 1 
 
End Subroutine CT_CouplingGluGluVG  
 
 
End Module CouplingsCT_MSSMTriLnV 

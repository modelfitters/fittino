! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.12.3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 13:03 on 7.2.2018   
! ----------------------------------------------------------------------  
 
 
Module LoopMasses_MSSMTriLnV 
 
Use Control 
Use Settings 
Use Couplings_MSSMTriLnV 
Use LoopFunctions 
Use AddLoopFunctions 
Use Mathematics 
Use MathematicsQP 
Use Model_Data_MSSMTriLnV 
Use StandardModel 
Use Tadpoles_MSSMTriLnV 
 Use EffectivePotential_MSSMTriLnV 
 Use Pole2L_MSSMTriLnV 
 Use TreeLevelMasses_MSSMTriLnV 
 
Use TwoLoopHiggsMass_SARAH 
 
Real(dp), Private :: MSd_1L(6), MSd2_1L(6)  
Complex(dp), Private :: ZD_1L(6,6)  
Real(dp), Private :: MSv_1L(3), MSv2_1L(3)  
Complex(dp), Private :: ZV_1L(3,3)  
Real(dp), Private :: MSu_1L(6), MSu2_1L(6)  
Complex(dp), Private :: ZU_1L(6,6)  
Real(dp), Private :: MSe_1L(6), MSe2_1L(6)  
Complex(dp), Private :: ZE_1L(6,6)  
Real(dp), Private :: Mhh_1L(2), Mhh2_1L(2)  
Complex(dp), Private :: ZH_1L(2,2)  
Real(dp), Private :: MAh_1L(2), MAh2_1L(2)  
Complex(dp), Private :: ZA_1L(2,2)  
Real(dp), Private :: MHpm_1L(2), MHpm2_1L(2)  
Complex(dp), Private :: ZP_1L(2,2)  
Real(dp), Private :: MChi_1L(4), MChi2_1L(4)  
Complex(dp), Private :: ZN_1L(4,4)  
Real(dp), Private :: MCha_1L(2), MCha2_1L(2)  
Complex(dp), Private :: UM_1L(2,2),UP_1L(2,2)
Real(dp), Private :: MGlu_1L, MGlu2_1L  
Real(dp), Private :: MVZ_1L, MVZ2_1L  
Real(dp), Private :: MVWm_1L, MVWm2_1L  
Real(dp) :: pi2A0  
Real(dp) :: ti_ep2L(2)  
Real(dp) :: pi_ep2L(2,2)
Real(dp) :: Pi2S_EffPot(2,2)
Real(dp) :: PiP2S_EffPot(2,2)
Contains 
 
Subroutine OneLoopMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,            & 
& MFu,MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,           & 
& MVWm,MVWm2,MVZ,MVZ2,pG,TW,UM,UP,v,ZA,ZD,ZDL,ZDR,ZE,ZEL,ZER,ZH,ZN,ZP,ZU,ZUL,            & 
& ZUR,ZV,ZW,ZZ,alphaH,betaH,vd,vu,g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,T2,Tu,             & 
& Bmu,mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,kont)

Implicit None 
Real(dp),Intent(inout) :: g1,g2,g3,mHd2,mHu2

Complex(dp),Intent(inout) :: Yd(3,3),Ye(3,3),L1(3,3,3),L2(3,3,3),Yu(3,3),Mu,Td(3,3),Te(3,3),T1(3,3,3),             & 
& T2(3,3,3),Tu(3,3),Bmu,mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3),M1,M2,M3

Real(dp),Intent(inout) :: MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),               & 
& MFe2(3),MFu(3),MFu2(3),MGlu,MGlu2,Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MSd(6),              & 
& MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3),MVWm,MVWm2,MVZ,MVZ2,              & 
& TW,v,ZA(2,2),ZH(2,2),ZP(2,2),ZZ(2,2),alphaH,betaH

Complex(dp),Intent(inout) :: pG,UM(2,2),UP(2,2),ZD(6,6),ZDL(3,3),ZDR(3,3),ZE(6,6),ZEL(3,3),ZER(3,3),               & 
& ZN(4,4),ZU(6,6),ZUL(3,3),ZUR(3,3),ZV(3,3),ZW(2,2)

Real(dp),Intent(inout) :: vd,vu

Complex(dp) :: cplAhAhcVWmVWm(2,2),cplAhAhUhh(2,2,2),cplAhAhUhhUhh(2,2,2,2),cplAhAhUHpmcUHpm(2,2,2,2),& 
& cplAhAhUSdcUSd(2,2,6,6),cplAhAhUSecUSe(2,2,6,6),cplAhAhUSucUSu(2,2,6,6),               & 
& cplAhAhUSvcUSv(2,2,3,3),cplAhAhVZVZ(2,2),cplAhcUHpmVWm(2,2),cplAhhhVZ(2,2),            & 
& cplAhHpmcUHpm(2,2,2),cplAhHpmcVWm(2,2),cplAhSdcUSd(2,6,6),cplAhSecUSe(2,6,6),          & 
& cplAhSucUSu(2,6,6),cplAhUhhVZ(2,2),cplcChaChaUAhL(2,2,2),cplcChaChaUAhR(2,2,2),        & 
& cplcChaChaUhhL(2,2,2),cplcChaChaUhhR(2,2,2),cplcChaChaVPL(2,2),cplcChaChaVPR(2,2),     & 
& cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),cplcChacUFuSdL(2,3,6),cplcChacUFuSdR(2,3,6),     & 
& cplcChaFdcUSuL(2,3,6),cplcChaFdcUSuR(2,3,6),cplcChaFecUSvL(2,3,3),cplcChaFecUSvR(2,3,3),& 
& cplcChaUChiHpmL(2,4,2),cplcChaUChiHpmR(2,4,2),cplcChaUChiVWmL(2,4),cplcChaUChiVWmR(2,4),& 
& cplcChaUFvSeL(2,3,6),cplcChaUFvSeR(2,3,6),cplcFdFdcUSvL(3,3,3),cplcFdFdcUSvR(3,3,3),   & 
& cplcFdFdUAhL(3,3,2),cplcFdFdUAhR(3,3,2),cplcFdFdUhhL(3,3,2),cplcFdFdUhhR(3,3,2),       & 
& cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),cplcFdFdVPL(3,3),cplcFdFdVPR(3,3),cplcFdFdVZL(3,3),  & 
& cplcFdFdVZR(3,3),cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFdUChiSdL(3,4,6),             & 
& cplcFdUChiSdR(3,4,6),cplcFdUFvSdL(3,3,6),cplcFdUFvSdR(3,3,6),cplcFecUFuSdL(3,3,6),     & 
& cplcFecUFuSdR(3,3,6),cplcFeFdcUSuL(3,3,6),cplcFeFdcUSuR(3,3,6),cplcFeFecUSvL(3,3,3),   & 
& cplcFeFecUSvR(3,3,3),cplcFeFeUAhL(3,3,2),cplcFeFeUAhR(3,3,2),cplcFeFeUhhL(3,3,2),      & 
& cplcFeFeUhhR(3,3,2),cplcFeFeVPL(3,3),cplcFeFeVPR(3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),& 
& cplcFeUChiSeL(3,4,6),cplcFeUChiSeR(3,4,6),cplcFeUFvHpmL(3,3,2),cplcFeUFvHpmR(3,3,2),   & 
& cplcFeUFvSeL(3,3,6),cplcFeUFvSeR(3,3,6),cplcFeUFvVWmL(3,3),cplcFeUFvVWmR(3,3),         & 
& cplcFuFdcUHpmL(3,3,2),cplcFuFdcUHpmR(3,3,2),cplcFuFdcUSeL(3,3,6),cplcFuFdcUSeR(3,3,6), & 
& cplcFuFdcVWmL(3,3),cplcFuFdcVWmR(3,3),cplcFuFuUAhL(3,3,2),cplcFuFuUAhR(3,3,2),         & 
& cplcFuFuUhhL(3,3,2),cplcFuFuUhhR(3,3,2),cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),             & 
& cplcFuFuVPL(3,3),cplcFuFuVPR(3,3),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcFuGluSuL(3,6), & 
& cplcFuGluSuR(3,6),cplcFuUChiSuL(3,4,6),cplcFuUChiSuR(3,4,6),cplcgAgWmcVWm,             & 
& cplcgGgGVG,cplcgWmgWmUAh(2),cplcgWmgWmUhh(2),cplcgWmgWmVP,cplcgWmgWmVZ,cplcgWmgZUHpm(2),& 
& cplcgWpCgAcVWm,cplcgWpCgWpCUAh(2),cplcgWpCgWpCUhh(2),cplcgWpCgWpCVP,cplcgWpCgWpCVZ,    & 
& cplcgWpCgZcUHpm(2),cplcgWpCgZcVWm,cplcgZgWmcUHpm(2),cplcgZgWmcVWm,cplcgZgWpCUHpm(2),   & 
& cplcgZgZUhh(2),cplChaFucUSdL(2,3,6),cplChaFucUSdR(2,3,6),cplChiChacUHpmL(4,2,2),       & 
& cplChiChacUHpmR(4,2,2),cplChiChacVWmL(4,2),cplChiChacVWmR(4,2),cplChiChiUAhL(4,4,2),   & 
& cplChiChiUAhR(4,4,2),cplChiChiUhhL(4,4,2),cplChiChiUhhR(4,4,2),cplChiChiVZL(4,4),      & 
& cplChiChiVZR(4,4),cplChiFdcUSdL(4,3,6),cplChiFdcUSdR(4,3,6),cplChiFecUSeL(4,3,6),      & 
& cplChiFecUSeR(4,3,6),cplChiFucUSuL(4,3,6),cplChiFucUSuR(4,3,6),cplChiFvcUSvL(4,3,3),   & 
& cplChiFvcUSvR(4,3,3),cplChiUFvcSvL(4,3,3),cplChiUFvcSvR(4,3,3),cplChiUFvSvL(4,3,3),    & 
& cplChiUFvSvR(4,3,3),cplcHpmVPVWm(2),cplcHpmVWmVZ(2),cplcUChacFuSdL(2,3,6),             & 
& cplcUChacFuSdR(2,3,6),cplcUChaChaAhL(2,2,2),cplcUChaChaAhR(2,2,2),cplcUChaChahhL(2,2,2),& 
& cplcUChaChahhR(2,2,2),cplcUChaChaVPL(2,2),cplcUChaChaVPR(2,2),cplcUChaChaVZL(2,2),     & 
& cplcUChaChaVZR(2,2),cplcUChaChiHpmL(2,4,2),cplcUChaChiHpmR(2,4,2),cplcUChaChiVWmL(2,4)

Complex(dp) :: cplcUChaChiVWmR(2,4),cplcUChaFdcSuL(2,3,6),cplcUChaFdcSuR(2,3,6),cplcUChaFecSvL(2,3,3),& 
& cplcUChaFecSvR(2,3,3),cplcUChaFvSeL(2,3,6),cplcUChaFvSeR(2,3,6),cplcUFdChaSuL(3,2,6),  & 
& cplcUFdChaSuR(3,2,6),cplcUFdChiSdL(3,4,6),cplcUFdChiSdR(3,4,6),cplcUFdFdAhL(3,3,2),    & 
& cplcUFdFdAhR(3,3,2),cplcUFdFdcSvL(3,3,3),cplcUFdFdcSvR(3,3,3),cplcUFdFdhhL(3,3,2),     & 
& cplcUFdFdhhR(3,3,2),cplcUFdFdVGL(3,3),cplcUFdFdVGR(3,3),cplcUFdFdVPL(3,3),             & 
& cplcUFdFdVPR(3,3),cplcUFdFdVZL(3,3),cplcUFdFdVZR(3,3),cplcUFdFeSuL(3,3,6),             & 
& cplcUFdFeSuR(3,3,6),cplcUFdFuHpmL(3,3,2),cplcUFdFuHpmR(3,3,2),cplcUFdFuSeL(3,3,6),     & 
& cplcUFdFuSeR(3,3,6),cplcUFdFuVWmL(3,3),cplcUFdFuVWmR(3,3),cplcUFdFvSdL(3,3,6),         & 
& cplcUFdFvSdR(3,3,6),cplcUFdGluSdL(3,6),cplcUFdGluSdR(3,6),cplcUFecFuSdL(3,3,6),        & 
& cplcUFecFuSdR(3,3,6),cplcUFeChaSvL(3,2,3),cplcUFeChaSvR(3,2,3),cplcUFeChiSeL(3,4,6),   & 
& cplcUFeChiSeR(3,4,6),cplcUFeFdcSuL(3,3,6),cplcUFeFdcSuR(3,3,6),cplcUFeFeAhL(3,3,2),    & 
& cplcUFeFeAhR(3,3,2),cplcUFeFecSvL(3,3,3),cplcUFeFecSvR(3,3,3),cplcUFeFehhL(3,3,2),     & 
& cplcUFeFehhR(3,3,2),cplcUFeFeVPL(3,3),cplcUFeFeVPR(3,3),cplcUFeFeVZL(3,3),             & 
& cplcUFeFeVZR(3,3),cplcUFeFvHpmL(3,3,2),cplcUFeFvHpmR(3,3,2),cplcUFeFvSeL(3,3,6),       & 
& cplcUFeFvSeR(3,3,6),cplcUFeFvVWmL(3,3),cplcUFeFvVWmR(3,3),cplcUFuChiSuL(3,4,6),        & 
& cplcUFuChiSuR(3,4,6),cplcUFuFdcHpmL(3,3,2),cplcUFuFdcHpmR(3,3,2),cplcUFuFdcSeL(3,3,6), & 
& cplcUFuFdcSeR(3,3,6),cplcUFuFdcVWmL(3,3),cplcUFuFdcVWmR(3,3),cplcUFuFuAhL(3,3,2),      & 
& cplcUFuFuAhR(3,3,2),cplcUFuFuhhL(3,3,2),cplcUFuFuhhR(3,3,2),cplcUFuFuVGL(3,3),         & 
& cplcUFuFuVGR(3,3),cplcUFuFuVPL(3,3),cplcUFuFuVPR(3,3),cplcUFuFuVZL(3,3),               & 
& cplcUFuFuVZR(3,3),cplcUFuGluSuL(3,6),cplcUFuGluSuR(3,6),cplcUHpmVPVWm(2),              & 
& cplcUHpmVWmVZ(2),cplcVWmcVWmVWmVWm1,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,             & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmVPVPVWm3,cplcVWmVPVWm,cplcVWmVPVWmVZ1,          & 
& cplcVWmVPVWmVZ2,cplcVWmVPVWmVZ3,cplcVWmVWmVZ,cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,          & 
& cplcVWmVWmVZVZ3,cplFeFucUSdL(3,3,6),cplFeFucUSdR(3,3,6),cplFvChacUSeL(3,2,6),          & 
& cplFvChacUSeR(3,2,6),cplFvFdcUSdL(3,3,6),cplFvFdcUSdR(3,3,6),cplFvFecUHpmL(3,3,2),     & 
& cplFvFecUHpmR(3,3,2),cplFvFecUSeL(3,3,6),cplFvFecUSeR(3,3,6),cplFvFecVWmL(3,3),        & 
& cplFvFecVWmR(3,3),cplFvFvVZL(3,3),cplFvFvVZR(3,3),cplGluFdcSdL(3,6),cplGluFdcSdR(3,6), & 
& cplGluFdcUSdL(3,6),cplGluFdcUSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6),             & 
& cplGluFucUSuL(3,6),cplGluFucUSuR(3,6),cplGluGluVGL,cplGluGluVGR,cplhhcUHpmVWm(2,2),    & 
& cplhhcVWmVWm(2),cplhhhhcVWmVWm(2,2),cplhhhhUHpmcUHpm(2,2,2,2),cplhhhhUSdcUSd(2,2,6,6), & 
& cplhhhhUSecUSe(2,2,6,6),cplhhhhUSucUSu(2,2,6,6),cplhhhhUSvcUSv(2,2,3,3),               & 
& cplhhhhVZVZ(2,2),cplhhHpmcUHpm(2,2,2),cplhhHpmcVWm(2,2),cplhhSdcUSd(2,6,6),            & 
& cplhhSecUSe(2,6,6),cplhhSucUSu(2,6,6),cplhhSvcUSv(2,3,3),cplhhVZVZ(2),cplHpmcHpmcVWmVWm(2,2),& 
& cplHpmcHpmVP(2,2),cplHpmcHpmVPVP(2,2),cplHpmcHpmVPVZ(2,2),cplHpmcHpmVZ(2,2),           & 
& cplHpmcHpmVZVZ(2,2),cplHpmcUHpmVP(2,2),cplHpmcUHpmVZ(2,2),cplHpmcVWmVP(2),             & 
& cplHpmcVWmVZ(2),cplHpmSucUSd(2,6,6),cplHpmSvcUSe(2,3,6),cplHpmUSdcHpmcUSd(2,6,2,6),    & 
& cplHpmUSecHpmcUSe(2,6,2,6),cplHpmUSucHpmcUSu(2,6,2,6),cplHpmUSvcHpmcUSv(2,3,2,3)

Complex(dp) :: cplSdcHpmcUSu(6,2,6),cplSdcSdcUSv(6,6,3),cplSdcSdcVWmVWm(6,6),cplSdcSdVG(6,6),         & 
& cplSdcSdVGVG(6,6),cplSdcSdVP(6,6),cplSdcSdVPVP(6,6),cplSdcSdVPVZ(6,6),cplSdcSdVZ(6,6), & 
& cplSdcSdVZVZ(6,6),cplSdcSecUSu(6,6,6),cplSdcSucVWm(6,6),cplSdcUHpmcSu(6,2,6),          & 
& cplSdcUSdcSv(6,6,3),cplSdcUSdVG(6,6),cplSdcUSdVP(6,6),cplSdcUSdVZ(6,6),cplSdcUSecSu(6,6,6),& 
& cplSdcUSucVWm(6,6),cplSdUSecSdcUSe(6,6,6,6),cplSdUSucSdcUSu(6,6,6,6),cplSdUSvcSdcUSv(6,3,6,3),& 
& cplSecHpmcUSv(6,2,3),cplSecSecUSv(6,6,3),cplSecSecVWmVWm(6,6),cplSecSeVP(6,6),         & 
& cplSecSeVPVP(6,6),cplSecSeVPVZ(6,6),cplSecSeVZ(6,6),cplSecSeVZVZ(6,6),cplSecSvcVWm(6,3),& 
& cplSecUHpmcSv(6,2,3),cplSecUSecSv(6,6,3),cplSecUSeVP(6,6),cplSecUSeVZ(6,6),            & 
& cplSecUSvcVWm(6,3),cplSeSucUSd(6,6,6),cplSeUSucSecUSu(6,6,6,6),cplSeUSvcSecUSv(6,3,6,3),& 
& cplSucSucVWmVWm(6,6),cplSucSuVG(6,6),cplSucSuVGVG(6,6),cplSucSuVP(6,6),cplSucSuVPVP(6,6),& 
& cplSucSuVPVZ(6,6),cplSucSuVZ(6,6),cplSucSuVZVZ(6,6),cplSucUSdVWm(6,6),cplSucUSuVG(6,6),& 
& cplSucUSuVP(6,6),cplSucUSuVZ(6,6),cplSuUSvcSucUSv(6,3,6,3),cplSvcSvcVWmVWm(3,3),       & 
& cplSvcSvVZ(3,3),cplSvcSvVZVZ(3,3),cplSvcUSeVWm(3,6),cplSvcUSvVZ(3,3),cplUAhAhhh(2,2,2),& 
& cplUAhhhVZ(2,2),cplUAhHpmcHpm(2,2,2),cplUAhHpmcVWm(2,2),cplUAhSdcSd(2,6,6),            & 
& cplUAhSecSe(2,6,6),cplUAhSucSu(2,6,6),cplUAhUAhAhAh(2,2,2,2),cplUAhUAhcVWmVWm(2,2),    & 
& cplUAhUAhhhhh(2,2,2,2),cplUAhUAhHpmcHpm(2,2,2,2),cplUAhUAhSdcSd(2,2,6,6),              & 
& cplUAhUAhSecSe(2,2,6,6),cplUAhUAhSucSu(2,2,6,6),cplUAhUAhSvcSv(2,2,3,3),               & 
& cplUAhUAhVZVZ(2,2),cplUChiChacHpmL(4,2,2),cplUChiChacHpmR(4,2,2),cplUChiChacVWmL(4,2), & 
& cplUChiChacVWmR(4,2),cplUChiChiAhL(4,4,2),cplUChiChiAhR(4,4,2),cplUChiChihhL(4,4,2),   & 
& cplUChiChihhR(4,4,2),cplUChiChiVZL(4,4),cplUChiChiVZR(4,4),cplUChiFdcSdL(4,3,6),       & 
& cplUChiFdcSdR(4,3,6),cplUChiFecSeL(4,3,6),cplUChiFecSeR(4,3,6),cplUChiFucSuL(4,3,6),   & 
& cplUChiFucSuR(4,3,6),cplUChiFvcSvL(4,3,3),cplUChiFvcSvR(4,3,3),cplUChiFvSvL(4,3,3),    & 
& cplUChiFvSvR(4,3,3),cplUFvChacSeL(3,2,6),cplUFvChacSeR(3,2,6),cplUFvFdcSdL(3,3,6),     & 
& cplUFvFdcSdR(3,3,6),cplUFvFecHpmL(3,3,2),cplUFvFecHpmR(3,3,2),cplUFvFecSeL(3,3,6),     & 
& cplUFvFecSeR(3,3,6),cplUFvFecVWmL(3,3),cplUFvFecVWmR(3,3),cplUFvFvVZL(3,3),            & 
& cplUFvFvVZR(3,3),cplUhhcVWmVWm(2),cplUhhhhhh(2,2,2),cplUhhHpmcHpm(2,2,2),              & 
& cplUhhHpmcVWm(2,2),cplUhhSdcSd(2,6,6),cplUhhSecSe(2,6,6),cplUhhSucSu(2,6,6),           & 
& cplUhhSvcSv(2,3,3),cplUhhUhhcVWmVWm(2,2),cplUhhUhhhhhh(2,2,2,2),cplUhhUhhHpmcHpm(2,2,2,2),& 
& cplUhhUhhSdcSd(2,2,6,6),cplUhhUhhSecSe(2,2,6,6),cplUhhUhhSucSu(2,2,6,6),               & 
& cplUhhUhhSvcSv(2,2,3,3),cplUhhUhhVZVZ(2,2),cplUhhVZVZ(2),cplUHpmcUHpmcVWmVWm(2,2),     & 
& cplUHpmcUHpmVPVP(2,2),cplUHpmcUHpmVZVZ(2,2),cplUHpmHpmcUHpmcHpm(2,2,2,2),              & 
& cplUHpmSdcUHpmcSd(2,6,2,6),cplUHpmSecUHpmcSe(2,6,2,6),cplUHpmSucUHpmcSu(2,6,2,6),      & 
& cplUHpmSvcUHpmcSv(2,3,2,3),cplUSdcUSdcVWmVWm(6,6),cplUSdcUSdVGVG(6,6),cplUSdcUSdVPVP(6,6),& 
& cplUSdcUSdVZVZ(6,6),cplUSdSdcUSdcSd(6,6,6,6),cplUSdSecUSdcSe(6,6,6,6),cplUSdSucUSdcSu(6,6,6,6),& 
& cplUSdSvcUSdcSv(6,3,6,3),cplUSecUSecVWmVWm(6,6),cplUSecUSeVPVP(6,6),cplUSecUSeVZVZ(6,6),& 
& cplUSeSecUSecSe(6,6,6,6),cplUSeSucUSecSu(6,6,6,6),cplUSeSvcUSecSv(6,3,6,3),            & 
& cplUSucUSucVWmVWm(6,6),cplUSucUSuVGVG(6,6),cplUSucUSuVPVP(6,6),cplUSucUSuVZVZ(6,6)

Complex(dp) :: cplUSuSucUSucSu(6,6,6,6),cplUSuSvcUSucSv(6,3,6,3),cplUSvcUSvcVWmVWm(3,3),              & 
& cplUSvcUSvVZVZ(3,3),cplUSvSvcUSvcSv(3,3,3,3),cplVGVGVG,cplVGVGVGVG1,cplVGVGVGVG2,      & 
& cplVGVGVGVG3

Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1, j2, j3, j4, il, i_count, ierr 
Integer :: i2L, iFin 
Logical :: Convergence2L 
Real(dp) :: Pi2S_EffPot_save(2,2), diff(2,2)
Complex(dp) :: Tad1Loop(2), dmz2  
Real(dp) :: comp(3), tanbQ, vev2, vSM
Real(dp) :: tadpoles_MSSM(2), vevs_MSSM(2), Q2, Pi2S_MSSM(2,2), sinb2, cosb2, sinbcosb 
Real(dp) :: tadpoles_asat(4), vevs_asat(4), Pi2S_asat(4,4),PiP2S_asat(4,4), Mglu_asat(2),mo2_asat, MAps 
complex(dp) :: MDO_asat,BO_asat,MO_asat,lam_asat,LT_asat, ZG_asat(2,2), mu_asat
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopMasses' 
 
kont = 0 
 
! Set Gauge fixing parameters 
RXi= RXiNew 
RXiG = RXi 
RXiP = RXi 
RXiWm = RXi 
RXiZ = RXi 

 ! Running angles 
TanBetaQ = vu/vd

 
Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,MVWm,               & 
& MVWm2,MVZ,MVZ2,pG,TW,UM,UP,v,ZA,ZD,ZDL,ZDR,ZE,ZEL,ZER,ZH,ZN,ZP,ZU,ZUL,ZUR,             & 
& ZV,ZW,ZZ,alphaH,betaH,vd,vu,g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,T2,Tu,Bmu,             & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,GenerationMixing,kont)

mHd2Tree  = mHd2
mHu2Tree  = mHu2
MuTree  = Mu
BmuTree  = Bmu

 
 If (CalculateOneLoopMasses) Then 
 
If ((DecoupleAtRenScale).and.(Abs(1._dp-RXiNew).lt.0.01_dp)) Then 
vSM=vSM_Q 
vd=vSM/Sqrt(1 + TanBetaQ**2) 
vu=TanBetaQ*vd 
Else 
Call CouplingsForVectorBosons(g1,g2,UM,UP,TW,ZP,vd,vu,ZD,ZE,ZU,ZH,ZA,ZN,              & 
& ZDL,ZUL,ZEL,ZV,cplcChaChaVPL,cplcChaChaVPR,cplcFdFdVPL,cplcFdFdVPR,cplcFeFeVPL,        & 
& cplcFeFeVPR,cplcFuFuVPL,cplcFuFuVPR,cplcgWmgWmVP,cplcgWpCgWpCVP,cplHpmcHpmVP,          & 
& cplHpmcVWmVP,cplSdcSdVP,cplSecSeVP,cplSucSuVP,cplcVWmVPVWm,cplHpmcHpmVPVP,             & 
& cplSdcSdVPVP,cplSecSeVPVP,cplSucSuVPVP,cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmVPVPVWm3,& 
& cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,cplChiChiVZR,cplcFdFdVZL,           & 
& cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,cplcFuFuVZR,cplFvFvVZL,cplFvFvVZR,     & 
& cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,cplHpmcHpmVZ,cplHpmcVWmVZ,cplSdcSdVZ,            & 
& cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,cplAhAhVZVZ,cplhhhhVZVZ,cplHpmcHpmVZVZ,  & 
& cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,cplSvcSvVZVZ,cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,   & 
& cplcVWmVWmVZVZ3,cplAhHpmcVWm,cplChiChacVWmL,cplChiChacVWmR,cplcFuFdcVWmL,              & 
& cplcFuFdcVWmR,cplFvFecVWmL,cplFvFecVWmR,cplcgWpCgAcVWm,cplcgAgWmcVWm,cplcgZgWmcVWm,    & 
& cplcgWpCgZcVWm,cplhhHpmcVWm,cplhhcVWmVWm,cplSdcSucVWm,cplSecSvcVWm,cplAhAhcVWmVWm,     & 
& cplhhhhcVWmVWm,cplHpmcHpmcVWmVWm,cplSdcSdcVWmVWm,cplSecSecVWmVWm,cplSucSucVWmVWm,      & 
& cplSvcSvcVWmVWm,cplcVWmcVWmVWmVWm1,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,              & 
& cplcHpmVWmVZ,cplcHpmVPVWm,cplHpmcHpmVPVZ,cplSdcSdVPVZ,cplSecSeVPVZ,cplSucSuVPVZ,       & 
& cplcVWmVPVWmVZ1,cplcVWmVPVWmVZ2,cplcVWmVPVWmVZ3)

Call Pi1LoopVZ(mZ2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,              & 
& MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,cplChiChiVZR,              & 
& cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplFvFvVZL,cplFvFvVZR,cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,cplHpmcHpmVZ,              & 
& cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,cplAhAhVZVZ,     & 
& cplhhhhVZVZ,cplHpmcHpmVZVZ,cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,cplSvcSvVZVZ,        & 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,dmZ2)

vev2=4._dp*Real(mZ2+dmz2,dp)/(g1**2+g2**2) -0 
vSM=sqrt(vev2) 
vd=vSM/Sqrt(1 + TanBetaQ**2) 
vu=TanBetaQ*vd 
End if 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,T2,Tu,Bmu,             & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,vd,vu,(/ ZeroC, ZeroC /))

Call TreeMasses(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,            & 
& MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,MVWm,               & 
& MVWm2,MVZ,MVZ2,pG,TW,UM,UP,v,ZA,ZD,ZDL,ZDR,ZE,ZEL,ZER,ZH,ZN,ZP,ZU,ZUL,ZUR,             & 
& ZV,ZW,ZZ,alphaH,betaH,vd,vu,g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,T2,Tu,Bmu,             & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,GenerationMixing,kont)

Call CouplingsForLoopMasses(Mu,Yd,Td,ZD,ZA,g2,Yu,UM,UP,ZUL,ZUR,g1,ZN,ZDL,             & 
& ZDR,L2,g3,pG,ZEL,vd,vu,ZH,Tu,ZU,ZP,T2,ZV,TW,Ye,ZE,L1,ZER,Te,T1,cplAhSdcUSd,            & 
& cplChaFucUSdL,cplChaFucUSdR,cplChiFdcUSdL,cplChiFdcUSdR,cplFvFdcUSdL,cplFvFdcUSdR,     & 
& cplGluFdcUSdL,cplGluFdcUSdR,cplFeFucUSdL,cplFeFucUSdR,cplhhSdcUSd,cplHpmSucUSd,        & 
& cplSdcUSdcSv,cplSdcUSdVG,cplSdcUSdVP,cplSdcUSdVZ,cplSeSucUSd,cplSucUSdVWm,             & 
& cplAhAhUSdcUSd,cplhhhhUSdcUSd,cplHpmUSdcHpmcUSd,cplUSdSdcUSdcSd,cplUSdSecUSdcSe,       & 
& cplUSdSucUSdcSu,cplUSdSvcUSdcSv,cplUSdcUSdVGVG,cplUSdcUSdVPVP,cplUSdcUSdcVWmVWm,       & 
& cplUSdcUSdVZVZ,cplChiFvcUSvL,cplChiFvcUSvR,cplcFdFdcUSvL,cplcFdFdcUSvR,cplcChaFecUSvL, & 
& cplcChaFecUSvR,cplcFeFecUSvL,cplcFeFecUSvR,cplhhSvcUSv,cplSdcSdcUSv,cplSecHpmcUSv,     & 
& cplSecSecUSv,cplSecUSvcVWm,cplSvcUSvVZ,cplAhAhUSvcUSv,cplhhhhUSvcUSv,cplHpmUSvcHpmcUSv,& 
& cplSdUSvcSdcUSv,cplSeUSvcSecUSv,cplSuUSvcSucUSv,cplUSvSvcUSvcSv,cplUSvcUSvcVWmVWm,     & 
& cplUSvcUSvVZVZ,cplAhSucUSu,cplChiFucUSuL,cplChiFucUSuR,cplcChaFdcUSuL,cplcChaFdcUSuR,  & 
& cplcFeFdcUSuL,cplcFeFdcUSuR,cplGluFucUSuL,cplGluFucUSuR,cplhhSucUSu,cplSdcHpmcUSu,     & 
& cplSdcSecUSu,cplSdcUSucVWm,cplSucUSuVG,cplSucUSuVP,cplSucUSuVZ,cplAhAhUSucUSu,         & 
& cplhhhhUSucUSu,cplHpmUSucHpmcUSu,cplSdUSucSdcUSu,cplSeUSucSecUSu,cplUSuSucUSucSu,      & 
& cplUSuSvcUSucSv,cplUSucUSuVGVG,cplUSucUSuVPVP,cplUSucUSucVWmVWm,cplUSucUSuVZVZ,        & 
& cplAhSecUSe,cplFvChacUSeL,cplFvChacUSeR,cplChiFecUSeL,cplChiFecUSeR,cplcFuFdcUSeL,     & 
& cplcFuFdcUSeR,cplFvFecUSeL,cplFvFecUSeR,cplhhSecUSe,cplHpmSvcUSe,cplSdcUSecSu,         & 
& cplSecUSecSv,cplSecUSeVP,cplSecUSeVZ,cplSvcUSeVWm,cplAhAhUSecUSe,cplhhhhUSecUSe,       & 
& cplHpmUSecHpmcUSe,cplSdUSecSdcUSe,cplUSeSecUSecSe,cplUSeSucUSecSu,cplUSeSvcUSecSv,     & 
& cplUSecUSeVPVP,cplUSecUSecVWmVWm,cplUSecUSeVZVZ,cplAhAhUhh,cplAhUhhVZ,cplcChaChaUhhL,  & 
& cplcChaChaUhhR,cplChiChiUhhL,cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,     & 
& cplcFeFeUhhR,cplcFuFuUhhL,cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,      & 
& cplUhhhhhh,cplUhhHpmcHpm,cplUhhHpmcVWm,cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,            & 
& cplUhhSvcSv,cplUhhcVWmVWm,cplUhhVZVZ,cplAhAhUhhUhh,cplUhhUhhhhhh,cplUhhUhhHpmcHpm,     & 
& cplUhhUhhSdcSd,cplUhhUhhSecSe,cplUhhUhhSucSu,cplUhhUhhSvcSv,cplUhhUhhcVWmVWm,          & 
& cplUhhUhhVZVZ,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,cplChiChiUAhL,cplChiChiUAhR,    & 
& cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,cplcFuFuUAhL,cplcFuFuUAhR,         & 
& cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,cplUAhHpmcVWm,cplUAhSdcSd,      & 
& cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,   & 
& cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,           & 
& cplAhHpmcUHpm,cplAhcUHpmVWm,cplChiChacUHpmL,cplChiChacUHpmR,cplcFuFdcUHpmL,            & 
& cplcFuFdcUHpmR,cplFvFecUHpmL,cplFvFecUHpmR,cplcgZgWmcUHpm,cplcgWmgZUHpm,               & 
& cplcgWpCgZcUHpm,cplcgZgWpCUHpm,cplhhHpmcUHpm,cplhhcUHpmVWm,cplHpmcUHpmVP,              & 
& cplHpmcUHpmVZ,cplSdcUHpmcSu,cplSecUHpmcSv,cplcUHpmVPVWm,cplcUHpmVWmVZ,cplAhAhUHpmcUHpm,& 
& cplhhhhUHpmcUHpm,cplUHpmHpmcUHpmcHpm,cplUHpmSdcUHpmcSd,cplUHpmSecUHpmcSe,              & 
& cplUHpmSucUHpmcSu,cplUHpmSvcUHpmcSv,cplUHpmcUHpmVPVP,cplUHpmcUHpmcVWmVWm,              & 
& cplUHpmcUHpmVZVZ,cplUChiChiAhL,cplUChiChiAhR,cplUChiChacHpmL,cplUChiChacHpmR,          & 
& cplUChiChacVWmL,cplUChiChacVWmR,cplUChiChihhL,cplUChiChihhR,cplUChiChiVZL,             & 
& cplUChiChiVZR,cplUChiFdcSdL,cplUChiFdcSdR,cplUChiFecSeL,cplUChiFecSeR,cplUChiFucSuL,   & 
& cplUChiFucSuR,cplUChiFvcSvL,cplUChiFvcSvR,cplUChiFvSvL,cplUChiFvSvR,cplcChaUChiHpmL,   & 
& cplcChaUChiHpmR,cplcFdUChiSdL,cplcFdUChiSdR,cplcFeUChiSeL,cplcFeUChiSeR,               & 
& cplcFuUChiSuL,cplcFuUChiSuR,cplcChaUChiVWmL,cplcChaUChiVWmR,cplcUChaChaAhL,            & 
& cplcUChaChaAhR,cplcUChaChahhL,cplcUChaChahhR,cplcUChaChaVPL,cplcUChaChaVPR,            & 
& cplcUChaChaVZL,cplcUChaChaVZR,cplcUChaChiHpmL,cplcUChaChiHpmR,cplcUChaChiVWmL,         & 
& cplcUChaChiVWmR,cplcUChaFdcSuL,cplcUChaFdcSuR,cplcUChaFecSvL,cplcUChaFecSvR,           & 
& cplcUChaFvSeL,cplcUChaFvSeR,cplcUChacFuSdL,cplcUChacFuSdR,cplcUFeFeAhL,cplcUFeFeAhR,   & 
& cplcUFeChaSvL,cplcUFeChaSvR,cplcUFeChiSeL,cplcUFeChiSeR,cplcUFeFdcSuL,cplcUFeFdcSuR,   & 
& cplcUFeFehhL,cplcUFeFehhR,cplcUFeFecSvL,cplcUFeFecSvR,cplcUFeFeVPL,cplcUFeFeVPR,       & 
& cplcUFeFeVZL,cplcUFeFeVZR,cplcUFeFvHpmL,cplcUFeFvHpmR,cplcUFeFvSeL,cplcUFeFvSeR,       & 
& cplcUFeFvVWmL,cplcUFeFvVWmR,cplcUFecFuSdL,cplcUFecFuSdR,cplcUFdFdAhL,cplcUFdFdAhR,     & 
& cplcUFdChaSuL,cplcUFdChaSuR,cplcUFdChiSdL,cplcUFdChiSdR,cplcUFdFdhhL,cplcUFdFdhhR,     & 
& cplcUFdFdcSvL,cplcUFdFdcSvR,cplcUFdFdVGL,cplcUFdFdVGR,cplcUFdFdVPL,cplcUFdFdVPR,       & 
& cplcUFdFdVZL,cplcUFdFdVZR,cplcUFdFeSuL,cplcUFdFeSuR,cplcUFdFuHpmL,cplcUFdFuHpmR,       & 
& cplcUFdFuSeL,cplcUFdFuSeR,cplcUFdFuVWmL,cplcUFdFuVWmR,cplcUFdFvSdL,cplcUFdFvSdR,       & 
& cplcUFdGluSdL,cplcUFdGluSdR,cplcUFuFuAhL,cplcUFuFuAhR,cplcUFuChiSuL,cplcUFuChiSuR,     & 
& cplcUFuFdcHpmL,cplcUFuFdcHpmR,cplcUFuFdcSeL,cplcUFuFdcSeR,cplcUFuFdcVWmL,              & 
& cplcUFuFdcVWmR,cplcUFuFuhhL,cplcUFuFuhhR,cplcUFuFuVGL,cplcUFuFuVGR,cplcUFuFuVPL,       & 
& cplcUFuFuVPR,cplcUFuFuVZL,cplcUFuFuVZR,cplcUFuGluSuL,cplcUFuGluSuR,cplcChacUFuSdL,     & 
& cplcChacUFuSdR,cplcFecUFuSdL,cplcFecUFuSdR,cplGluFdcSdL,cplGluFdcSdR,cplGluFucSuL,     & 
& cplGluFucSuR,cplGluGluVGL,cplGluGluVGR,cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,         & 
& cplcFuGluSuR,cplUFvChacSeL,cplUFvChacSeR,cplChiUFvcSvL,cplChiUFvcSvR,cplChiUFvSvL,     & 
& cplChiUFvSvR,cplUFvFdcSdL,cplUFvFdcSdR,cplUFvFecHpmL,cplUFvFecHpmR,cplUFvFecSeL,       & 
& cplUFvFecSeR,cplUFvFecVWmL,cplUFvFecVWmR,cplUFvFvVZL,cplUFvFvVZR,cplcFeUFvHpmL,        & 
& cplcFeUFvHpmR,cplcFdUFvSdL,cplcFdUFvSdR,cplcChaUFvSeL,cplcChaUFvSeR,cplcFeUFvSeL,      & 
& cplcFeUFvSeR,cplcFeUFvVWmL,cplcFeUFvVWmR,cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,          & 
& cplcFuFuVGR,cplcgGgGVG,cplSdcSdVG,cplSucSuVG,cplVGVGVG,cplSdcSdVGVG,cplSucSuVGVG,      & 
& cplVGVGVGVG1,cplVGVGVGVG2,cplVGVGVGVG3,cplcChaChaVPL,cplcChaChaVPR,cplcFdFdVPL,        & 
& cplcFdFdVPR,cplcFeFeVPL,cplcFeFeVPR,cplcFuFuVPL,cplcFuFuVPR,cplcgWmgWmVP,              & 
& cplcgWpCgWpCVP,cplHpmcHpmVP,cplHpmcVWmVP,cplSdcSdVP,cplSecSeVP,cplSucSuVP,             & 
& cplcVWmVPVWm,cplHpmcHpmVPVP,cplSdcSdVPVP,cplSecSeVPVP,cplSucSuVPVP,cplcVWmVPVPVWm1,    & 
& cplcVWmVPVPVWm2,cplcVWmVPVPVWm3,cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,    & 
& cplChiChiVZR,cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,              & 
& cplcFuFuVZR,cplFvFvVZL,cplFvFvVZR,cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,               & 
& cplHpmcHpmVZ,cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,    & 
& cplAhAhVZVZ,cplhhhhVZVZ,cplHpmcHpmVZVZ,cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,         & 
& cplSvcSvVZVZ,cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,cplAhHpmcVWm,             & 
& cplChiChacVWmL,cplChiChacVWmR,cplcFuFdcVWmL,cplcFuFdcVWmR,cplFvFecVWmL,cplFvFecVWmR,   & 
& cplcgWpCgAcVWm,cplcgAgWmcVWm,cplcgZgWmcVWm,cplcgWpCgZcVWm,cplhhHpmcVWm,cplhhcVWmVWm,   & 
& cplSdcSucVWm,cplSecSvcVWm,cplAhAhcVWmVWm,cplhhhhcVWmVWm,cplHpmcHpmcVWmVWm,             & 
& cplSdcSdcVWmVWm,cplSecSecVWmVWm,cplSucSucVWmVWm,cplSvcSvcVWmVWm,cplcVWmcVWmVWmVWm1,    & 
& cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcHpmVWmVZ,cplcHpmVPVWm,cplHpmcHpmVPVZ,        & 
& cplSdcSdVPVZ,cplSecSeVPVZ,cplSucSuVPVZ,cplcVWmVPVWmVZ1,cplcVWmVPVWmVZ2,cplcVWmVPVWmVZ3)

Call OneLoopTadpoleshh(vd,vu,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,             & 
& MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,MVWm,            & 
& MVWm2,MVZ,MVZ2,cplAhAhUhh,cplcChaChaUhhL,cplcChaChaUhhR,cplChiChiUhhL,cplChiChiUhhR,   & 
& cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,cplcFeFeUhhR,cplcFuFuUhhL,cplcFuFuUhhR,         & 
& cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,cplUhhhhhh,cplUhhHpmcHpm,cplUhhSdcSd,        & 
& cplUhhSecSe,cplUhhSucSu,cplUhhSvcSv,cplUhhcVWmVWm,cplUhhVZVZ,Tad1Loop(1:2))

mHd2Tree  = mHd2
mHu2Tree  = mHu2
MuTree  = Mu
BmuTree  = Bmu
If (CalculateTwoLoopHiggsMasses) Then 
    If(GaugelessLimit) Then 
  vdFix = 0._dp 
  vuFix = 0._dp 
   g1_saveEP =g1
   g1 = 0._dp 
   g2_saveEP =g2
   g2 = 0._dp 
     Else 
  vdFix = vd 
  vuFix = vu 
     End if 

SELECT CASE (TwoLoopMethod) 
CASE ( 1 , 2 ) 
  ! Make sure that there are no exactly degenerated masses! 
   Yd_saveEP =Yd
   where (aint(Abs(Yd)).eq.Yd) Yd=Yd*(1 + 1*1.0E-12_dp)
   Ye_saveEP =Ye
   where (aint(Abs(Ye)).eq.Ye) Ye=Ye*(1 + 2*1.0E-12_dp)
   Yu_saveEP =Yu
   where (aint(Abs(Yu)).eq.Yu) Yu=Yu*(1 + 3*1.0E-12_dp)
   Td_saveEP =Td
   where (aint(Abs(Td)).eq.Td) Td=Td*(1 + 4*1.0E-12_dp)
   Te_saveEP =Te
   where (aint(Abs(Te)).eq.Te) Te=Te*(1 + 5*1.0E-12_dp)
   Tu_saveEP =Tu
   where (aint(Abs(Tu)).eq.Tu) Tu=Tu*(1 + 6*1.0E-12_dp)
   mq2_saveEP =mq2
   where (aint(Abs(mq2)).eq.mq2) mq2=mq2*(1 + 7*1.0E-12_dp)
   ml2_saveEP =ml2
   where (aint(Abs(ml2)).eq.ml2) ml2=ml2*(1 + 8*1.0E-12_dp)
   md2_saveEP =md2
   where (aint(Abs(md2)).eq.md2) md2=md2*(1 + 9*1.0E-12_dp)
   mu2_saveEP =mu2
   where (aint(Abs(mu2)).eq.mu2) mu2=mu2*(1 + 10*1.0E-12_dp)
   me2_saveEP =me2
   where (aint(Abs(me2)).eq.me2) me2=me2*(1 + 11*1.0E-12_dp)

If (TwoLoopSafeMode) Then 
  iFin = 12 
  Convergence2L = .false. 
  hstep_pn = 2.0_dp 
  hstep_sa = 2.0_dp 
Else 
  iFin = 1 
  Convergence2L = .true. 
End if 

Pi2S_EffPot_save = 0._dp 
Pi2S_EffPot = 0._dp 

Do i2L = 1, iFin 
Call CalculateCorrectionsEffPot(ti_ep2L,pi_ep2L,vd,vu,g1,g2,g3,Yd,Ye,L1,              & 
& L2,Yu,Mu,Td,Te,T1,T2,Tu,Bmu,mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,kont)

Pi2S_EffPot(1,1) = pi_ep2L(1,1)!-ti_ep2L(1)/vd
Pi2S_EffPot(1,2) = pi_ep2L(1,2)
Pi2S_EffPot(2,1) = pi_ep2L(2,1)
Pi2S_EffPot(2,2) = pi_ep2L(2,2)!-ti_ep2L(2)/vu
 diff=(Pi2S_EffPot-Pi2S_EffPot_save)/MaxVal(Abs(Pi2S_EffPot)) 
  If (MaxVal(Abs(diff)).lt.1.0E-4_dp) Then 
    Convergence2L = .True. 
    Exit 
  Else 
    Pi2S_EffPot_save = Pi2S_EffPot 
  hstep_pn = hstep_pn/2._dp 
  hstep_sa = hstep_sa/2._dp 
  End If 
End do 
If (.not.Convergence2L) Then 
 Write(*,*) "WARNING: Two-Loop corrections are numerically unstable! "  
 Call TerminateProgram 
End If  
 Pi2A0 = 0._dp 
   Yd =Yd_saveEP 
   Ye =Ye_saveEP 
   Yu =Yu_saveEP 
   Td =Td_saveEP 
   Te =Te_saveEP 
   Tu =Tu_saveEP 
   mq2 =mq2_saveEP 
   ml2 =ml2_saveEP 
   md2 =md2_saveEP 
   mu2 =mu2_saveEP 
   me2 =me2_saveEP 


 CASE ( 3 ) ! Diagrammatic method 
  ! Make sure that there are no exactly degenerated masses! 
   Yd_saveEP =Yd
   where (aint(Abs(Yd)).eq.Yd) Yd=Yd*(1 + 1*1.0E-12_dp)
   Ye_saveEP =Ye
   where (aint(Abs(Ye)).eq.Ye) Ye=Ye*(1 + 2*1.0E-12_dp)
   Yu_saveEP =Yu
   where (aint(Abs(Yu)).eq.Yu) Yu=Yu*(1 + 3*1.0E-12_dp)
   Td_saveEP =Td
   where (aint(Abs(Td)).eq.Td) Td=Td*(1 + 4*1.0E-12_dp)
   Te_saveEP =Te
   where (aint(Abs(Te)).eq.Te) Te=Te*(1 + 5*1.0E-12_dp)
   Tu_saveEP =Tu
   where (aint(Abs(Tu)).eq.Tu) Tu=Tu*(1 + 6*1.0E-12_dp)
   mq2_saveEP =mq2
   where (aint(Abs(mq2)).eq.mq2) mq2=mq2*(1 + 7*1.0E-12_dp)
   ml2_saveEP =ml2
   where (aint(Abs(ml2)).eq.ml2) ml2=ml2*(1 + 8*1.0E-12_dp)
   md2_saveEP =md2
   where (aint(Abs(md2)).eq.md2) md2=md2*(1 + 9*1.0E-12_dp)
   mu2_saveEP =mu2
   where (aint(Abs(mu2)).eq.mu2) mu2=mu2*(1 + 10*1.0E-12_dp)
   me2_saveEP =me2
   where (aint(Abs(me2)).eq.me2) me2=me2*(1 + 11*1.0E-12_dp)

If (NewGBC) Then 
Call CalculatePi2S(125._dp**2,vd,vu,g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,              & 
& T2,Tu,Bmu,mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,kont,ti_ep2L,Pi2S_EffPot,             & 
& PiP2S_EffPot)

Else 
Call CalculatePi2S(0._dp,vd,vu,g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,T2,Tu,             & 
& Bmu,mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,kont,ti_ep2L,Pi2S_EffPot,PiP2S_EffPot)

End if 
   Yd =Yd_saveEP 
   Ye =Ye_saveEP 
   Yu =Yu_saveEP 
   Td =Td_saveEP 
   Te =Te_saveEP 
   Tu =Tu_saveEP 
   mq2 =mq2_saveEP 
   ml2 =ml2_saveEP 
   md2 =md2_saveEP 
   mu2 =mu2_saveEP 
   me2 =me2_saveEP 


 CASE ( 8 , 9 ) ! Hard-coded routines 
  
! Strong corrections first
vevs_asat = 1._dp
vevs_asat(1) = vd
vevs_asat(2) = vu
Mglu_asat(1) = MGlu
Mglu_asat(2) = 1._dp
ZG_asat(1,1) = 1._dp
ZG_asat(1,2) = 0._dp
ZG_asat(2,1) = 0._dp
ZG_asat(2,2) = 1._dp
MDO_asat = 0._dp
mo2_asat = 1._dp
MO_asat = 0._dp
BO_asat = 0._dp
lam_asat = 0._dp
mu_asat = mu ! Sign changed FS 03/06/15 
LT_asat = 0._dp
Q2 = GetRenormalizationScale() 
call CalculateStrongCorrections2L(Q2,g3,MGlu_asat,ZG_asat, vevs_asat, &
& Real(md2(3,3),dp),Real(mu2(3,3),dp),Real(mq2(3,3),dp),Real(me2(3,3),dp), &
& Real(ml2(3,3),dp), Td(3,3), Tu(3,3), Te(3,3), Yd(3,3), &
& Yu(3,3), Ye(3,3), mu_asat,MDO_asat,mo2_asat,BO_asat,MO_asat,lam_asat,LT_asat,0,pip2s_asat,pi2s_asat,tadpoles_asat,kont )
do i1 = 1,2
ti_ep2L(i1) = tadpoles_asat(i1)
pi2s_effpot(i1,i1) = pi2s_asat(i1,i1)+tadpoles_asat(i1)/vevs_asat(i1)
pip2s_effpot(i1,i1) = pip2s_asat(i1,i1)+tadpoles_asat(i1)/vevs_asat(i1)
do i2 = 1,i1-1
pi2s_effpot(i1,i2) = pi2s_asat(i1,i2)
pi2s_effpot(i2,i1) = pi2s_asat(i2,i1)
pip2s_effpot(i1,i2) = pip2s_asat(i1,i2)
pip2s_effpot(i2,i1) = pip2s_asat(i2,i1)
end do
end do

If ( TwoLoopMethod .eq. 9) then ! Slavich's routines 

vevs_MSSM(1) = vd 
vevs_MSSM(2) = vu 
MAps = Bmu*(vu/vd+vd/vu) 
Q2 = GetRenormalizationScale() 
Call Yukawa2L_Tadpoles_MSSM(MAps,vevs_MSSM, &
& Real(md2(3,3),dp),Real(mu2(3,3),dp),Real(mq2(3,3),dp),Real(me2(3,3),dp), & 
& Real(ml2(3,3),dp), Td(3,3), Tu(3,3), Te(3,3), Yd(3,3), &  
& Yu(3,3), Ye(3,3), mu_asat, tadpoles_MSSM,kont ) 

ti_ep2L(1) = ti_ep2L(1) + tadpoles_MSSM(1)*vevs_MSSM(1) 
ti_ep2L(2) = ti_ep2L(2) + tadpoles_MSSM(2)*vevs_MSSM(2) 
 
vevs_MSSM(1) = vd 
vevs_MSSM(2) = vu 
Call Yukawa2L_Scalar(Q2,MAps,vevs_MSSM, & 
 Real(md2(3,3),dp),Real(mu2(3,3),dp),Real(mq2(3,3),dp),Real(me2(3,3),dp), & 
 Real(ml2(3,3),dp), Td(3,3), Tu(3,3), Te(3,3), Yd(3,3), &  
 Yu(3,3), Ye(3,3), mu_asat,0, Pi2S_MSSM,kont ) 

! Some two loop corrections are absorbed in M_A in Pietro's routines 

 Call Yukawa2L_PseudoScalar(MAps,vevs_MSSM,& 
 & Real(md2(3,3),dp),Real(mu2(3,3),dp),Real(mq2(3,3),dp),Real(me2(3,3),dp),& 
 & Real(ml2(3,3),dp),Td(3,3),Tu(3,3),Te(3,3),Yd(3,3),& 
 & Yu(3,3),Ye(3,3),mu_asat,Pi2A0,kont) 
 tanbQ=vu/vd 
 cosb2=1._dp/(1._dp+tanbQ**2) 
 sinb2=1._dp-cosb2 
 sinbcosb=Sqrt(cosb2*sinb2) 
 Pi2S_MSSM(1,1)=Pi2S_MSSM(1,1)+Pi2A0*sinb2 + tadpoles_MSSM(1) 
 Pi2S_MSSM(1,2)=Pi2S_MSSM(1,2)-Pi2A0*sinbcosb 
 Pi2S_MSSM(2,2)=Pi2S_MSSM(2,2)+Pi2A0*cosb2 + tadpoles_MSSM(2) 
 Pi2S_MSSM(2,1)=Pi2S_MSSM(1,2) 
 
Pi2S_EffPot(1:2,1:2) = Pi2S_EffPot(1:2,1:2) + Pi2S_MSSM 

 PiP2S_effpot(1,1)=pip2s_effpot(1,1)+Pi2A0*sinb2 +  tadpoles_MSSM(1) 
 PiP2S_effpot(1,2)=pip2s_effpot(1,2)+Pi2A0*sinbcosb 
 PiP2S_effpot(2,2)=pip2s_effpot(2,2)+Pi2A0*cosb2 +  tadpoles_MSSM(2) 
 PiP2S_effpot(2,1)=pip2s_effpot(1,2)
 
end if ! Ends slavich routines

 END SELECT
 
   If(GaugelessLimit) Then 
   g1 =g1_saveEP 
   g2 =g2_saveEP 
   End if 

Else ! Two loop turned off 
Pi2S_EffPot = 0._dp 

Pi2A0 = 0._dp 

ti_ep2L = 0._dp 

End if 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,T2,Tu,Bmu,             & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,vd,vu,Tad1Loop)

Mu1L = Mu
Bmu1L = Bmu
mHd21L = mHd2
mHu21L = mHu2
Tad1Loop(1:2) = Tad1Loop(1:2) - ti_ep2L 
Call SolveTadpoleEquations(g1,g2,g3,Yd,Ye,L1,L2,Yu,Mu,Td,Te,T1,T2,Tu,Bmu,             & 
& mq2,ml2,mHd2,mHu2,md2,mu2,me2,M1,M2,M3,vd,vu,Tad1Loop)

Mu2L = Mu
Bmu2L = Bmu
mHd22L = mHd2
mHu22L = mHu2
Call OneLoopSd(g1,g2,Mu1L,Yd,Td,mq2,md2,vd,vu,MSd,MSd2,MAh,MAh2,MFu,MFu2,             & 
& MCha,MCha2,MFd,MFd2,MChi,MChi2,MGlu,MGlu2,MFe,MFe2,Mhh,Mhh2,MSu,MSu2,MHpm,             & 
& MHpm2,MSv,MSv2,MVZ,MVZ2,MSe,MSe2,MVWm,MVWm2,cplAhSdcUSd,cplChaFucUSdL,cplChaFucUSdR,   & 
& cplChiFdcUSdL,cplChiFdcUSdR,cplFvFdcUSdL,cplFvFdcUSdR,cplGluFdcUSdL,cplGluFdcUSdR,     & 
& cplFeFucUSdL,cplFeFucUSdR,cplhhSdcUSd,cplHpmSucUSd,cplSdcUSdcSv,cplSdcUSdVG,           & 
& cplSdcUSdVP,cplSdcUSdVZ,cplSeSucUSd,cplSucUSdVWm,cplAhAhUSdcUSd,cplhhhhUSdcUSd,        & 
& cplHpmUSdcHpmcUSd,cplUSdSdcUSdcSd,cplUSdSecUSdcSe,cplUSdSucUSdcSu,cplUSdSvcUSdcSv,     & 
& cplUSdcUSdVGVG,cplUSdcUSdVPVP,cplUSdcUSdcVWmVWm,cplUSdcUSdVZVZ,0.1_dp*delta_mass,      & 
& MSd_1L,MSd2_1L,ZD_1L,kont)

Call OneLoopSv(g1,g2,ml2,vd,vu,MChi,MChi2,MFd,MFd2,MCha,MCha2,MFe,MFe2,               & 
& MSv,MSv2,Mhh,Mhh2,MSd,MSd2,MHpm,MHpm2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,MAh,MAh2,           & 
& MSu,MSu2,cplChiFvcUSvL,cplChiFvcUSvR,cplcFdFdcUSvL,cplcFdFdcUSvR,cplcChaFecUSvL,       & 
& cplcChaFecUSvR,cplcFeFecUSvL,cplcFeFecUSvR,cplhhSvcUSv,cplSdcSdcUSv,cplSecHpmcUSv,     & 
& cplSecSecUSv,cplSecUSvcVWm,cplSvcUSvVZ,cplAhAhUSvcUSv,cplhhhhUSvcUSv,cplHpmUSvcHpmcUSv,& 
& cplSdUSvcSdcUSv,cplSeUSvcSecUSv,cplSuUSvcSucUSv,cplUSvSvcUSvcSv,cplUSvcUSvcVWmVWm,     & 
& cplUSvcUSvVZVZ,0.1_dp*delta_mass,MSv_1L,MSv2_1L,ZV_1L,kont)

Call OneLoopSu(g1,g2,Mu1L,Yu,Tu,mq2,mu2,vd,vu,MSu,MSu2,MAh,MAh2,MFu,MFu2,             & 
& MChi,MChi2,MCha,MCha2,MFd,MFd2,MFe,MFe2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,            & 
& MSd2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,MSv,MSv2,cplAhSucUSu,cplChiFucUSuL,cplChiFucUSuR,    & 
& cplcChaFdcUSuL,cplcChaFdcUSuR,cplcFeFdcUSuL,cplcFeFdcUSuR,cplGluFucUSuL,               & 
& cplGluFucUSuR,cplhhSucUSu,cplSdcHpmcUSu,cplSdcSecUSu,cplSdcUSucVWm,cplSucUSuVG,        & 
& cplSucUSuVP,cplSucUSuVZ,cplAhAhUSucUSu,cplhhhhUSucUSu,cplHpmUSucHpmcUSu,               & 
& cplSdUSucSdcUSu,cplSeUSucSecUSu,cplUSuSucUSucSu,cplUSuSvcUSucSv,cplUSucUSuVGVG,        & 
& cplUSucUSuVPVP,cplUSucUSucVWmVWm,cplUSucUSuVZVZ,0.1_dp*delta_mass,MSu_1L,              & 
& MSu2_1L,ZU_1L,kont)

Call OneLoopSe(g1,g2,Mu1L,Ye,Te,ml2,me2,vd,vu,MSe,MSe2,MAh,MAh2,MCha,MCha2,           & 
& MFe,MFe2,MChi,MChi2,MFu,MFu2,MFd,MFd2,Mhh,Mhh2,MSv,MSv2,MHpm,MHpm2,MSu,MSu2,           & 
& MSd,MSd2,MVZ,MVZ2,MVWm,MVWm2,cplAhSecUSe,cplFvChacUSeL,cplFvChacUSeR,cplChiFecUSeL,    & 
& cplChiFecUSeR,cplcFuFdcUSeL,cplcFuFdcUSeR,cplFvFecUSeL,cplFvFecUSeR,cplhhSecUSe,       & 
& cplHpmSvcUSe,cplSdcUSecSu,cplSecUSecSv,cplSecUSeVP,cplSecUSeVZ,cplSvcUSeVWm,           & 
& cplAhAhUSecUSe,cplhhhhUSecUSe,cplHpmUSecHpmcUSe,cplSdUSecSdcUSe,cplUSeSecUSecSe,       & 
& cplUSeSucUSecSu,cplUSeSvcUSecSv,cplUSecUSeVPVP,cplUSecUSecVWmVWm,cplUSecUSeVZVZ,       & 
& 0.1_dp*delta_mass,MSe_1L,MSe2_1L,ZE_1L,kont)

Call OneLoophh(g1,g2,Mu2L,Bmu2L,mHd22L,mHu22L,vd,vu,MAh,MAh2,MVZ,MVZ2,MCha,           & 
& MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MVWm,MVWm2,            & 
& MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,cplAhAhUhh,cplAhUhhVZ,cplcChaChaUhhL,              & 
& cplcChaChaUhhR,cplChiChiUhhL,cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,     & 
& cplcFeFeUhhR,cplcFuFuUhhL,cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,      & 
& cplUhhhhhh,cplUhhHpmcHpm,cplUhhHpmcVWm,cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,            & 
& cplUhhSvcSv,cplUhhcVWmVWm,cplUhhVZVZ,cplAhAhUhhUhh,cplUhhUhhhhhh,cplUhhUhhHpmcHpm,     & 
& cplUhhUhhSdcSd,cplUhhUhhSecSe,cplUhhUhhSucSu,cplUhhUhhSvcSv,cplUhhUhhcVWmVWm,          & 
& cplUhhUhhVZVZ,0.1_dp*delta_mass,Mhh_1L,Mhh2_1L,ZH_1L,kont)

If (TwoLoopMethod.gt.2) Then 
Call OneLoopAh(g1,g2,Mu2L,Bmu2L,mHd22L,mHu22L,vd,vu,TW,Mhh,Mhh2,MAh,MAh2,             & 
& MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,             & 
& MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,    & 
& cplChiChiUAhL,cplChiChiUAhR,cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,       & 
& cplcFuFuUAhL,cplcFuFuUAhR,cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,      & 
& cplUAhHpmcVWm,cplUAhSdcSd,cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,         & 
& cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,          & 
& cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,0.1_dp*delta_mass,MAh_1L,MAh2_1L,ZA_1L,kont)

Else 
Call OneLoopAh(g1,g2,Mu1L,Bmu1L,mHd21L,mHu21L,vd,vu,TW,Mhh,Mhh2,MAh,MAh2,             & 
& MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,             & 
& MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,    & 
& cplChiChiUAhL,cplChiChiUAhR,cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,       & 
& cplcFuFuUAhL,cplcFuFuUAhR,cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,      & 
& cplUAhHpmcVWm,cplUAhSdcSd,cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,         & 
& cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,          & 
& cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,0.1_dp*delta_mass,MAh_1L,MAh2_1L,ZA_1L,kont)

End if 
Call OneLoopHpm(g1,g2,Mu1L,Bmu1L,mHd21L,mHu21L,vd,vu,MHpm,MHpm2,MAh,MAh2,             & 
& MVWm,MVWm2,MChi,MChi2,MCha,MCha2,MFu,MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVZ,              & 
& MVZ2,MSu,MSu2,MSd,MSd2,MSv,MSv2,MSe,MSe2,cplAhHpmcUHpm,cplAhcUHpmVWm,cplChiChacUHpmL,  & 
& cplChiChacUHpmR,cplcFuFdcUHpmL,cplcFuFdcUHpmR,cplFvFecUHpmL,cplFvFecUHpmR,             & 
& cplcgZgWmcUHpm,cplcgWmgZUHpm,cplcgWpCgZcUHpm,cplcgZgWpCUHpm,cplhhHpmcUHpm,             & 
& cplhhcUHpmVWm,cplHpmcUHpmVP,cplHpmcUHpmVZ,cplSdcUHpmcSu,cplSecUHpmcSv,cplcUHpmVPVWm,   & 
& cplcUHpmVWmVZ,cplAhAhUHpmcUHpm,cplhhhhUHpmcUHpm,cplUHpmHpmcUHpmcHpm,cplUHpmSdcUHpmcSd, & 
& cplUHpmSecUHpmcSe,cplUHpmSucUHpmcSu,cplUHpmSvcUHpmcSv,cplUHpmcUHpmVPVP,cplUHpmcUHpmcVWmVWm,& 
& cplUHpmcUHpmVZVZ,0.1_dp*delta_mass,MHpm_1L,MHpm2_1L,ZP_1L,kont)

Call OneLoopChi(g1,g2,Mu1L,M1,M2,vd,vu,MChi,MChi2,MAh,MAh2,MHpm,MHpm2,MCha,           & 
& MCha2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MSd,MSd2,MFd,MFd2,MSe,MSe2,MFe,MFe2,MSu,            & 
& MSu2,MFu,MFu2,MSv,MSv2,cplUChiChiAhL,cplUChiChiAhR,cplUChiChacHpmL,cplUChiChacHpmR,    & 
& cplUChiChacVWmL,cplUChiChacVWmR,cplUChiChihhL,cplUChiChihhR,cplUChiChiVZL,             & 
& cplUChiChiVZR,cplUChiFdcSdL,cplUChiFdcSdR,cplUChiFecSeL,cplUChiFecSeR,cplUChiFucSuL,   & 
& cplUChiFucSuR,cplUChiFvcSvL,cplUChiFvcSvR,cplUChiFvSvL,cplUChiFvSvR,cplcChaUChiHpmL,   & 
& cplcChaUChiHpmR,cplcFdUChiSdL,cplcFdUChiSdR,cplcFeUChiSeL,cplcFeUChiSeR,               & 
& cplcFuUChiSuL,cplcFuUChiSuR,cplcChaUChiVWmL,cplcChaUChiVWmR,0.1_dp*delta_mass,         & 
& MChi_1L,MChi2_1L,ZN_1L,kont)

Call OneLoopCha(g2,Mu1L,M2,vd,vu,MCha,MCha2,MAh,MAh2,Mhh,Mhh2,MVZ,MVZ2,               & 
& MHpm,MHpm2,MChi,MChi2,MVWm,MVWm2,MSu,MSu2,MFd,MFd2,MSv,MSv2,MFe,MFe2,MSe,              & 
& MSe2,MFu,MFu2,MSd,MSd2,cplcUChaChaAhL,cplcUChaChaAhR,cplcUChaChahhL,cplcUChaChahhR,    & 
& cplcUChaChaVPL,cplcUChaChaVPR,cplcUChaChaVZL,cplcUChaChaVZR,cplcUChaChiHpmL,           & 
& cplcUChaChiHpmR,cplcUChaChiVWmL,cplcUChaChiVWmR,cplcUChaFdcSuL,cplcUChaFdcSuR,         & 
& cplcUChaFecSvL,cplcUChaFecSvR,cplcUChaFvSeL,cplcUChaFvSeR,cplcUChacFuSdL,              & 
& cplcUChacFuSdR,0.1_dp*delta_mass,MCha_1L,MCha2_1L,UM_1L,UP_1L,kont)

Call OneLoopGlu(M3,MSd,MSd2,MFd,MFd2,MSu,MSu2,MFu,MFu2,MGlu,MGlu2,cplGluFdcSdL,       & 
& cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,cplGluGluVGL,cplGluGluVGR,cplcFdGluSdL,         & 
& cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,0.1_dp*delta_mass,MGlu_1L,MGlu2_1L,kont)

MSd = MSd_1L 
MSd2 = MSd2_1L 
ZD = ZD_1L 
MSv = MSv_1L 
MSv2 = MSv2_1L 
ZV = ZV_1L 
MSu = MSu_1L 
MSu2 = MSu2_1L 
ZU = ZU_1L 
MSe = MSe_1L 
MSe2 = MSe2_1L 
ZE = ZE_1L 
Mhh = Mhh_1L 
Mhh2 = Mhh2_1L 
ZH = ZH_1L 
MAh = MAh_1L 
MAh2 = MAh2_1L 
ZA = ZA_1L 
MHpm = MHpm_1L 
MHpm2 = MHpm2_1L 
ZP = ZP_1L 
MChi = MChi_1L 
MChi2 = MChi2_1L 
ZN = ZN_1L 
MCha = MCha_1L 
MCha2 = MCha2_1L 
UM = UM_1L 
UP = UP_1L 
MGlu = MGlu_1L 
MGlu2 = MGlu2_1L 
End If 
 
Call SortGoldstones(MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,             & 
& MFu2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,               & 
& MVWm,MVWm2,MVZ,MVZ2,pG,TW,UM,UP,v,ZA,ZD,ZDL,ZDR,ZE,ZEL,ZER,ZH,ZN,ZP,ZU,ZUL,            & 
& ZUR,ZV,ZW,ZZ,alphaH,betaH,kont)

! Set pole masses 
MVWm = mW 
MVWm2 = mW2 
MVZ = mZ 
MVZ2 = mZ2 
MFe(1:3) = mf_l 
MFe2(1:3) = mf_l**2 
MFu(1:3) = mf_u 
MFu2(1:3) = mf_u**2 
MFd(1:3) = mf_d 
MFd2(1:3) = mf_d**2 
! Shift Everything to t'Hooft Gauge
RXi=  1._dp 
RXiG = 1._dp 
RXiP = 1._dp 
RXiWm = 1._dp 
RXiZ = 1._dp 
MAh(1)=MVZ
MAh2(1)=MVZ2
MHpm(1)=MVWm
MHpm2(1)=MVWm2
mf_u2 = mf_u**2 
mf_d2 = mf_d**2 
mf_l2 = mf_l**2 
 

 If (WriteTreeLevelTadpoleSolutions) Then 
! Saving tree-level parameters for output
mHd2  = mHd2Tree 
mHu2  = mHu2Tree 
Mu  = MuTree 
Bmu  = BmuTree 
End if 


Iname = Iname -1 
End Subroutine OneLoopMasses 
 
Subroutine OneLoopTadpoleshh(vd,vu,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,           & 
& MFe,MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,             & 
& MVWm,MVWm2,MVZ,MVZ2,cplAhAhUhh,cplcChaChaUhhL,cplcChaChaUhhR,cplChiChiUhhL,            & 
& cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,cplcFeFeUhhR,cplcFuFuUhhL,        & 
& cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,cplUhhhhhh,cplUhhHpmcHpm,       & 
& cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,cplUhhSvcSv,cplUhhcVWmVWm,cplUhhVZVZ,tadpoles)

Implicit None 
Real(dp), Intent(in) :: MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),               & 
& MFe2(3),MFu(3),MFu2(3),Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MSd(6),MSd2(6),MSe(6),          & 
& MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3),MVWm,MVWm2,MVZ,MVZ2

Complex(dp), Intent(in) :: cplAhAhUhh(2,2,2),cplcChaChaUhhL(2,2,2),cplcChaChaUhhR(2,2,2),cplChiChiUhhL(4,4,2),   & 
& cplChiChiUhhR(4,4,2),cplcFdFdUhhL(3,3,2),cplcFdFdUhhR(3,3,2),cplcFeFeUhhL(3,3,2),      & 
& cplcFeFeUhhR(3,3,2),cplcFuFuUhhL(3,3,2),cplcFuFuUhhR(3,3,2),cplcgWmgWmUhh(2),          & 
& cplcgWpCgWpCUhh(2),cplcgZgZUhh(2),cplUhhhhhh(2,2,2),cplUhhHpmcHpm(2,2,2),              & 
& cplUhhSdcSd(2,6,6),cplUhhSecSe(2,6,6),cplUhhSucSu(2,6,6),cplUhhSvcSv(2,3,3),           & 
& cplUhhcVWmVWm(2),cplUhhVZVZ(2)

Real(dp), Intent(in) :: vd,vu

Integer :: i1,i2, gO1, gO2 
Complex(dp) :: coupL, coupR, coup, temp, res, A0m, sumI(2)  
Real(dp) :: m1 
Complex(dp), Intent(out) :: tadpoles(2) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopTadpoleshh'
 
tadpoles = 0._dp 
 
!------------------------ 
! Ah 
!------------------------ 
Do i1 = 1, 2
 A0m = SA_A0(MAh2(i1)) 
  Do gO1 = 1, 2
   coup = cplAhAhUhh(i1,i1,gO1)
   sumI(gO1) = -coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp/2._dp*sumI 
End Do 
 !------------------------ 
! bar[Cha] 
!------------------------ 
Do i1 = 1, 2
 A0m = 2._dp*MCha(i1)*SA_A0(MCha2(i1)) 
  Do gO1 = 1, 2
   coupL = cplcChaChaUhhL(i1,i1,gO1)
   coupR = cplcChaChaUhhR(i1,i1,gO1)
   sumI(gO1) = (coupL+coupR)*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
End Do 
 !------------------------ 
! Chi 
!------------------------ 
Do i1 = 1, 4
 A0m = 2._dp*MChi(i1)*SA_A0(MChi2(i1)) 
  Do gO1 = 1, 2
   coupL = cplChiChiUhhL(i1,i1,gO1)
   coupR = cplChiChiUhhR(i1,i1,gO1)
   sumI(gO1) = (coupL+coupR)*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp/2._dp*sumI 
End Do 
 !------------------------ 
! bar[Fd] 
!------------------------ 
Do i1 = 1, 3
 A0m = 2._dp*MFd(i1)*SA_A0(MFd2(i1)) 
  Do gO1 = 1, 2
   coupL = cplcFdFdUhhL(i1,i1,gO1)
   coupR = cplcFdFdUhhR(i1,i1,gO1)
   sumI(gO1) = (coupL+coupR)*A0m 
  End Do 
 
tadpoles =  tadpoles + 3._dp*sumI 
End Do 
 !------------------------ 
! bar[Fe] 
!------------------------ 
Do i1 = 1, 3
 A0m = 2._dp*MFe(i1)*SA_A0(MFe2(i1)) 
  Do gO1 = 1, 2
   coupL = cplcFeFeUhhL(i1,i1,gO1)
   coupR = cplcFeFeUhhR(i1,i1,gO1)
   sumI(gO1) = (coupL+coupR)*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
End Do 
 !------------------------ 
! bar[Fu] 
!------------------------ 
Do i1 = 1, 3
 A0m = 2._dp*MFu(i1)*SA_A0(MFu2(i1)) 
  Do gO1 = 1, 2
   coupL = cplcFuFuUhhL(i1,i1,gO1)
   coupR = cplcFuFuUhhR(i1,i1,gO1)
   sumI(gO1) = (coupL+coupR)*A0m 
  End Do 
 
tadpoles =  tadpoles + 3._dp*sumI 
End Do 
 !------------------------ 
! bar[gWm] 
!------------------------ 
A0m = 1._dp*SA_A0(MVWm2*RXi) 
  Do gO1 = 1, 2
    coup = cplcgWmgWmUhh(gO1)
    sumI(gO1) = coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
!------------------------ 
! bar[gWmC] 
!------------------------ 
A0m = 1._dp*SA_A0(MVWm2*RXi) 
  Do gO1 = 1, 2
    coup = cplcgWpCgWpCUhh(gO1)
    sumI(gO1) = coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
!------------------------ 
! bar[gZ] 
!------------------------ 
A0m = 1._dp*SA_A0(MVZ2*RXi) 
  Do gO1 = 1, 2
    coup = cplcgZgZUhh(gO1)
    sumI(gO1) = coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
!------------------------ 
! hh 
!------------------------ 
Do i1 = 1, 2
 A0m = SA_A0(Mhh2(i1)) 
  Do gO1 = 1, 2
   coup = cplUhhhhhh(gO1,i1,i1)
   sumI(gO1) = -coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp/2._dp*sumI 
End Do 
 !------------------------ 
! conj[Hpm] 
!------------------------ 
Do i1 = 1, 2
 A0m = SA_A0(MHpm2(i1)) 
  Do gO1 = 1, 2
   coup = cplUhhHpmcHpm(gO1,i1,i1)
   sumI(gO1) = -coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
Do i1 = 1, 6
 A0m = SA_A0(MSd2(i1)) 
  Do gO1 = 1, 2
   coup = cplUhhSdcSd(gO1,i1,i1)
   sumI(gO1) = -coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 3._dp*sumI 
End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
Do i1 = 1, 6
 A0m = SA_A0(MSe2(i1)) 
  Do gO1 = 1, 2
   coup = cplUhhSecSe(gO1,i1,i1)
   sumI(gO1) = -coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
Do i1 = 1, 6
 A0m = SA_A0(MSu2(i1)) 
  Do gO1 = 1, 2
   coup = cplUhhSucSu(gO1,i1,i1)
   sumI(gO1) = -coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 3._dp*sumI 
End Do 
 !------------------------ 
! conj[Sv] 
!------------------------ 
Do i1 = 1, 3
 A0m = SA_A0(MSv2(i1)) 
  Do gO1 = 1, 2
   coup = cplUhhSvcSv(gO1,i1,i1)
   sumI(gO1) = -coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
End Do 
 !------------------------ 
! conj[VWm] 
!------------------------ 
A0m = 3._dp*SA_A0(MVWm2)+RXi*SA_A0(MVWm2*RXi) - 2._dp*MVWm2*rMS 
  Do gO1 = 1, 2
    coup = cplUhhcVWmVWm(gO1)
    sumI(gO1) = coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp*sumI 
!------------------------ 
! VZ 
!------------------------ 
A0m = 3._dp*SA_A0(MVZ2)+RXi*SA_A0(MVZ2*RXi) - 2._dp*MVZ2*rMS 
  Do gO1 = 1, 2
    coup = cplUhhVZVZ(gO1)
    sumI(gO1) = coup*A0m 
  End Do 
 
tadpoles =  tadpoles + 1._dp/2._dp*sumI 



tadpoles = oo16pi2*tadpoles 
Iname = Iname - 1 
End Subroutine OneLoopTadpoleshh 
 
Subroutine OneLoopSd(g1,g2,Mu,Yd,Td,mq2,md2,vd,vu,MSd,MSd2,MAh,MAh2,MFu,              & 
& MFu2,MCha,MCha2,MFd,MFd2,MChi,MChi2,MGlu,MGlu2,MFe,MFe2,Mhh,Mhh2,MSu,MSu2,             & 
& MHpm,MHpm2,MSv,MSv2,MVZ,MVZ2,MSe,MSe2,MVWm,MVWm2,cplAhSdcUSd,cplChaFucUSdL,            & 
& cplChaFucUSdR,cplChiFdcUSdL,cplChiFdcUSdR,cplFvFdcUSdL,cplFvFdcUSdR,cplGluFdcUSdL,     & 
& cplGluFdcUSdR,cplFeFucUSdL,cplFeFucUSdR,cplhhSdcUSd,cplHpmSucUSd,cplSdcUSdcSv,         & 
& cplSdcUSdVG,cplSdcUSdVP,cplSdcUSdVZ,cplSeSucUSd,cplSucUSdVWm,cplAhAhUSdcUSd,           & 
& cplhhhhUSdcUSd,cplHpmUSdcHpmcUSd,cplUSdSdcUSdcSd,cplUSdSecUSdcSe,cplUSdSucUSdcSu,      & 
& cplUSdSvcUSdcSv,cplUSdcUSdVGVG,cplUSdcUSdVPVP,cplUSdcUSdcVWmVWm,cplUSdcUSdVZVZ,        & 
& delta,mass,mass2,RS,kont)

Implicit None 
Real(dp), Intent(in) :: MSd(6),MSd2(6),MAh(2),MAh2(2),MFu(3),MFu2(3),MCha(2),MCha2(2),MFd(3),MFd2(3),         & 
& MChi(4),MChi2(4),MGlu,MGlu2,MFe(3),MFe2(3),Mhh(2),Mhh2(2),MSu(6),MSu2(6),              & 
& MHpm(2),MHpm2(2),MSv(3),MSv2(3),MVZ,MVZ2,MSe(6),MSe2(6),MVWm,MVWm2

Real(dp), Intent(in) :: g1,g2,vd,vu

Complex(dp), Intent(in) :: Mu,Yd(3,3),Td(3,3),mq2(3,3),md2(3,3)

Complex(dp), Intent(in) :: cplAhSdcUSd(2,6,6),cplChaFucUSdL(2,3,6),cplChaFucUSdR(2,3,6),cplChiFdcUSdL(4,3,6),    & 
& cplChiFdcUSdR(4,3,6),cplFvFdcUSdL(3,3,6),cplFvFdcUSdR(3,3,6),cplGluFdcUSdL(3,6),       & 
& cplGluFdcUSdR(3,6),cplFeFucUSdL(3,3,6),cplFeFucUSdR(3,3,6),cplhhSdcUSd(2,6,6),         & 
& cplHpmSucUSd(2,6,6),cplSdcUSdcSv(6,6,3),cplSdcUSdVG(6,6),cplSdcUSdVP(6,6),             & 
& cplSdcUSdVZ(6,6),cplSeSucUSd(6,6,6),cplSucUSdVWm(6,6),cplAhAhUSdcUSd(2,2,6,6),         & 
& cplhhhhUSdcUSd(2,2,6,6),cplHpmUSdcHpmcUSd(2,6,2,6),cplUSdSdcUSdcSd(6,6,6,6),           & 
& cplUSdSecUSdcSe(6,6,6,6),cplUSdSucUSdcSu(6,6,6,6),cplUSdSvcUSdcSv(6,3,6,3),            & 
& cplUSdcUSdVGVG(6,6),cplUSdcUSdVPVP(6,6),cplUSdcUSdcVWmVWm(6,6),cplUSdcUSdVZVZ(6,6)

Complex(dp) :: mat2a(6,6), mat2(6,6),  PiSf(6,6,6)
Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(6), test_m2(6),p2, test(6) 
Real(dp), Intent(out) :: mass(6), mass2(6) 
Complex(dp), Intent(out) ::  RS(6,6) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopSd'
 
mat2a(1,1) = 0._dp 
mat2a(1,1) = mat2a(1,1)-(g1**2*vd**2)/24._dp
mat2a(1,1) = mat2a(1,1)-(g2**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(g1**2*vu**2)/24._dp
mat2a(1,1) = mat2a(1,1)+(g2**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+mq2(1,1)
Do j1 = 1,3
mat2a(1,1) = mat2a(1,1)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,1))/2._dp
End Do 
mat2a(1,2) = 0._dp 
mat2a(1,2) = mat2a(1,2)+mq2(1,2)
Do j1 = 1,3
mat2a(1,2) = mat2a(1,2)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,2))/2._dp
End Do 
mat2a(1,3) = 0._dp 
mat2a(1,3) = mat2a(1,3)+mq2(1,3)
Do j1 = 1,3
mat2a(1,3) = mat2a(1,3)+(vd**2*Conjg(Yd(j1,1))*Yd(j1,3))/2._dp
End Do 
mat2a(1,4) = 0._dp 
mat2a(1,4) = mat2a(1,4)-((vu*Mu*Conjg(Yd(1,1)))/sqrt(2._dp))
mat2a(1,4) = mat2a(1,4)+(vd*Conjg(Td(1,1)))/sqrt(2._dp)
mat2a(1,5) = 0._dp 
mat2a(1,5) = mat2a(1,5)-((vu*Mu*Conjg(Yd(2,1)))/sqrt(2._dp))
mat2a(1,5) = mat2a(1,5)+(vd*Conjg(Td(2,1)))/sqrt(2._dp)
mat2a(1,6) = 0._dp 
mat2a(1,6) = mat2a(1,6)-((vu*Mu*Conjg(Yd(3,1)))/sqrt(2._dp))
mat2a(1,6) = mat2a(1,6)+(vd*Conjg(Td(3,1)))/sqrt(2._dp)
mat2a(2,2) = 0._dp 
mat2a(2,2) = mat2a(2,2)-(g1**2*vd**2)/24._dp
mat2a(2,2) = mat2a(2,2)-(g2**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g1**2*vu**2)/24._dp
mat2a(2,2) = mat2a(2,2)+(g2**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+mq2(2,2)
Do j1 = 1,3
mat2a(2,2) = mat2a(2,2)+(vd**2*Conjg(Yd(j1,2))*Yd(j1,2))/2._dp
End Do 
mat2a(2,3) = 0._dp 
mat2a(2,3) = mat2a(2,3)+mq2(2,3)
Do j1 = 1,3
mat2a(2,3) = mat2a(2,3)+(vd**2*Conjg(Yd(j1,2))*Yd(j1,3))/2._dp
End Do 
mat2a(2,4) = 0._dp 
mat2a(2,4) = mat2a(2,4)-((vu*Mu*Conjg(Yd(1,2)))/sqrt(2._dp))
mat2a(2,4) = mat2a(2,4)+(vd*Conjg(Td(1,2)))/sqrt(2._dp)
mat2a(2,5) = 0._dp 
mat2a(2,5) = mat2a(2,5)-((vu*Mu*Conjg(Yd(2,2)))/sqrt(2._dp))
mat2a(2,5) = mat2a(2,5)+(vd*Conjg(Td(2,2)))/sqrt(2._dp)
mat2a(2,6) = 0._dp 
mat2a(2,6) = mat2a(2,6)-((vu*Mu*Conjg(Yd(3,2)))/sqrt(2._dp))
mat2a(2,6) = mat2a(2,6)+(vd*Conjg(Td(3,2)))/sqrt(2._dp)
mat2a(3,3) = 0._dp 
mat2a(3,3) = mat2a(3,3)-(g1**2*vd**2)/24._dp
mat2a(3,3) = mat2a(3,3)-(g2**2*vd**2)/8._dp
mat2a(3,3) = mat2a(3,3)+(g1**2*vu**2)/24._dp
mat2a(3,3) = mat2a(3,3)+(g2**2*vu**2)/8._dp
mat2a(3,3) = mat2a(3,3)+mq2(3,3)
Do j1 = 1,3
mat2a(3,3) = mat2a(3,3)+(vd**2*Conjg(Yd(j1,3))*Yd(j1,3))/2._dp
End Do 
mat2a(3,4) = 0._dp 
mat2a(3,4) = mat2a(3,4)-((vu*Mu*Conjg(Yd(1,3)))/sqrt(2._dp))
mat2a(3,4) = mat2a(3,4)+(vd*Conjg(Td(1,3)))/sqrt(2._dp)
mat2a(3,5) = 0._dp 
mat2a(3,5) = mat2a(3,5)-((vu*Mu*Conjg(Yd(2,3)))/sqrt(2._dp))
mat2a(3,5) = mat2a(3,5)+(vd*Conjg(Td(2,3)))/sqrt(2._dp)
mat2a(3,6) = 0._dp 
mat2a(3,6) = mat2a(3,6)-((vu*Mu*Conjg(Yd(3,3)))/sqrt(2._dp))
mat2a(3,6) = mat2a(3,6)+(vd*Conjg(Td(3,3)))/sqrt(2._dp)
mat2a(4,4) = 0._dp 
mat2a(4,4) = mat2a(4,4)-(g1**2*vd**2)/12._dp
mat2a(4,4) = mat2a(4,4)+(g1**2*vu**2)/12._dp
mat2a(4,4) = mat2a(4,4)+md2(1,1)
Do j1 = 1,3
mat2a(4,4) = mat2a(4,4)+(vd**2*Conjg(Yd(1,j1))*Yd(1,j1))/2._dp
End Do 
mat2a(4,5) = 0._dp 
mat2a(4,5) = mat2a(4,5)+md2(1,2)
Do j1 = 1,3
mat2a(4,5) = mat2a(4,5)+(vd**2*Conjg(Yd(2,j1))*Yd(1,j1))/2._dp
End Do 
mat2a(4,6) = 0._dp 
mat2a(4,6) = mat2a(4,6)+md2(1,3)
Do j1 = 1,3
mat2a(4,6) = mat2a(4,6)+(vd**2*Conjg(Yd(3,j1))*Yd(1,j1))/2._dp
End Do 
mat2a(5,5) = 0._dp 
mat2a(5,5) = mat2a(5,5)-(g1**2*vd**2)/12._dp
mat2a(5,5) = mat2a(5,5)+(g1**2*vu**2)/12._dp
mat2a(5,5) = mat2a(5,5)+md2(2,2)
Do j1 = 1,3
mat2a(5,5) = mat2a(5,5)+(vd**2*Conjg(Yd(2,j1))*Yd(2,j1))/2._dp
End Do 
mat2a(5,6) = 0._dp 
mat2a(5,6) = mat2a(5,6)+md2(2,3)
Do j1 = 1,3
mat2a(5,6) = mat2a(5,6)+(vd**2*Conjg(Yd(3,j1))*Yd(2,j1))/2._dp
End Do 
mat2a(6,6) = 0._dp 
mat2a(6,6) = mat2a(6,6)-(g1**2*vd**2)/12._dp
mat2a(6,6) = mat2a(6,6)+(g1**2*vu**2)/12._dp
mat2a(6,6) = mat2a(6,6)+md2(3,3)
Do j1 = 1,3
mat2a(6,6) = mat2a(6,6)+(vd**2*Conjg(Yd(3,j1))*Yd(3,j1))/2._dp
End Do 

 
 Do i1=2,6
  Do i2 = 1, i1-1 
  mat2a(i1,i2) = Conjg(mat2a(i2,i1)) 
  End do 
End do 

 
! Rotation matrix for p2=0 
PiSf(1,:,:) = ZeroC 
p2 = 0._dp 
Call Pi1LoopSd(p2,MSd,MSd2,MAh,MAh2,MFu,MFu2,MCha,MCha2,MFd,MFd2,MChi,MChi2,          & 
& MGlu,MGlu2,MFe,MFe2,Mhh,Mhh2,MSu,MSu2,MHpm,MHpm2,MSv,MSv2,MVZ,MVZ2,MSe,MSe2,           & 
& MVWm,MVWm2,cplAhSdcUSd,cplChaFucUSdL,cplChaFucUSdR,cplChiFdcUSdL,cplChiFdcUSdR,        & 
& cplFvFdcUSdL,cplFvFdcUSdR,cplGluFdcUSdL,cplGluFdcUSdR,cplFeFucUSdL,cplFeFucUSdR,       & 
& cplhhSdcUSd,cplHpmSucUSd,cplSdcUSdcSv,cplSdcUSdVG,cplSdcUSdVP,cplSdcUSdVZ,             & 
& cplSeSucUSd,cplSucUSdVWm,cplAhAhUSdcUSd,cplhhhhUSdcUSd,cplHpmUSdcHpmcUSd,              & 
& cplUSdSdcUSdcSd,cplUSdSecUSdcSe,cplUSdSucUSdcSu,cplUSdSvcUSdcSv,cplUSdcUSdVGVG,        & 
& cplUSdcUSdVPVP,cplUSdcUSdcVWmVWm,cplUSdcUSdVZVZ,kont,PiSf(1,:,:))

mat2 = mat2a - Real(PiSf(1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZDOS_0 = RS 
 
 
! Now with momenta 
Do i1=1,6
PiSf(i1,:,:) = ZeroC 
p2 = MSd2(i1)
Call Pi1LoopSd(p2,MSd,MSd2,MAh,MAh2,MFu,MFu2,MCha,MCha2,MFd,MFd2,MChi,MChi2,          & 
& MGlu,MGlu2,MFe,MFe2,Mhh,Mhh2,MSu,MSu2,MHpm,MHpm2,MSv,MSv2,MVZ,MVZ2,MSe,MSe2,           & 
& MVWm,MVWm2,cplAhSdcUSd,cplChaFucUSdL,cplChaFucUSdR,cplChiFdcUSdL,cplChiFdcUSdR,        & 
& cplFvFdcUSdL,cplFvFdcUSdR,cplGluFdcUSdL,cplGluFdcUSdR,cplFeFucUSdL,cplFeFucUSdR,       & 
& cplhhSdcUSd,cplHpmSucUSd,cplSdcUSdcSv,cplSdcUSdVG,cplSdcUSdVP,cplSdcUSdVZ,             & 
& cplSeSucUSd,cplSucUSdVWm,cplAhAhUSdcUSd,cplhhhhUSdcUSd,cplHpmUSdcHpmcUSd,              & 
& cplUSdSdcUSdcSd,cplUSdSecUSdcSe,cplUSdSucUSdcSu,cplUSdSvcUSdcSv,cplUSdcUSdVGVG,        & 
& cplUSdcUSdVPVP,cplUSdcUSdcVWmVWm,cplUSdcUSdVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=6,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
 
Do i1=1,6
  If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
  If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = Sqrt(mass2(i1)) 
  Else 
   If (ErrorLevel.Ge.0) Then 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
   Call TerminateProgram 
   End If 
   kont = -301 
   mass(i1) = 0._dp 
  End If 
End Do 
 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
Do i1=1,6
PiSf(i1,:,:) = ZeroC 
p2 =  mass2(i1) 
Call Pi1LoopSd(p2,MSd,MSd2,MAh,MAh2,MFu,MFu2,MCha,MCha2,MFd,MFd2,MChi,MChi2,          & 
& MGlu,MGlu2,MFe,MFe2,Mhh,Mhh2,MSu,MSu2,MHpm,MHpm2,MSv,MSv2,MVZ,MVZ2,MSe,MSe2,           & 
& MVWm,MVWm2,cplAhSdcUSd,cplChaFucUSdL,cplChaFucUSdR,cplChiFdcUSdL,cplChiFdcUSdR,        & 
& cplFvFdcUSdL,cplFvFdcUSdR,cplGluFdcUSdL,cplGluFdcUSdR,cplFeFucUSdL,cplFeFucUSdR,       & 
& cplhhSdcUSd,cplHpmSucUSd,cplSdcUSdcSv,cplSdcUSdVG,cplSdcUSdVP,cplSdcUSdVZ,             & 
& cplSeSucUSd,cplSucUSdVWm,cplAhAhUSdcUSd,cplhhhhUSdcUSd,cplHpmUSdcHpmcUSd,              & 
& cplUSdSdcUSdcSd,cplUSdSecUSdcSe,cplUSdSucUSdcSu,cplUSdSvcUSdcSv,cplUSdcUSdVGVG,        & 
& cplUSdcUSdVPVP,cplUSdcUSdcVWmVWm,cplUSdcUSdVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=6,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZDOS_p2(i1,:) = RS(i1,:) 
 If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
Do i1=1,6
 If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
 If (test_m2(i1).Ne.0._dp) Then 
    test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1) 
 Else 
    test_m2(i1) = Abs(mass2(i1)) 
 End If 
 If (Abs(mass2(i1)).lt.1.0E-30_dp) test_m2(i1) = 0._dp 
 If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = sqrt(mass2(i1)) 
  Else 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses occurred a negative mass squared!' 
     Write(*,*) 'generation: ',i1 
     Write(*,*) 'mass: ',mass2(i1) 
   SignOfMassChanged = .True. 
   mass(i1) = 0._dp 
  End If 
End Do 
 
If (Maxval(test_m2).LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopSd
 
 
Subroutine Pi1LoopSd(p2,MSd,MSd2,MAh,MAh2,MFu,MFu2,MCha,MCha2,MFd,MFd2,               & 
& MChi,MChi2,MGlu,MGlu2,MFe,MFe2,Mhh,Mhh2,MSu,MSu2,MHpm,MHpm2,MSv,MSv2,MVZ,              & 
& MVZ2,MSe,MSe2,MVWm,MVWm2,cplAhSdcUSd,cplChaFucUSdL,cplChaFucUSdR,cplChiFdcUSdL,        & 
& cplChiFdcUSdR,cplFvFdcUSdL,cplFvFdcUSdR,cplGluFdcUSdL,cplGluFdcUSdR,cplFeFucUSdL,      & 
& cplFeFucUSdR,cplhhSdcUSd,cplHpmSucUSd,cplSdcUSdcSv,cplSdcUSdVG,cplSdcUSdVP,            & 
& cplSdcUSdVZ,cplSeSucUSd,cplSucUSdVWm,cplAhAhUSdcUSd,cplhhhhUSdcUSd,cplHpmUSdcHpmcUSd,  & 
& cplUSdSdcUSdcSd,cplUSdSecUSdcSe,cplUSdSucUSdcSu,cplUSdSvcUSdcSv,cplUSdcUSdVGVG,        & 
& cplUSdcUSdVPVP,cplUSdcUSdcVWmVWm,cplUSdcUSdVZVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MSd(6),MSd2(6),MAh(2),MAh2(2),MFu(3),MFu2(3),MCha(2),MCha2(2),MFd(3),MFd2(3),         & 
& MChi(4),MChi2(4),MGlu,MGlu2,MFe(3),MFe2(3),Mhh(2),Mhh2(2),MSu(6),MSu2(6),              & 
& MHpm(2),MHpm2(2),MSv(3),MSv2(3),MVZ,MVZ2,MSe(6),MSe2(6),MVWm,MVWm2

Complex(dp), Intent(in) :: cplAhSdcUSd(2,6,6),cplChaFucUSdL(2,3,6),cplChaFucUSdR(2,3,6),cplChiFdcUSdL(4,3,6),    & 
& cplChiFdcUSdR(4,3,6),cplFvFdcUSdL(3,3,6),cplFvFdcUSdR(3,3,6),cplGluFdcUSdL(3,6),       & 
& cplGluFdcUSdR(3,6),cplFeFucUSdL(3,3,6),cplFeFucUSdR(3,3,6),cplhhSdcUSd(2,6,6),         & 
& cplHpmSucUSd(2,6,6),cplSdcUSdcSv(6,6,3),cplSdcUSdVG(6,6),cplSdcUSdVP(6,6),             & 
& cplSdcUSdVZ(6,6),cplSeSucUSd(6,6,6),cplSucUSdVWm(6,6),cplAhAhUSdcUSd(2,2,6,6),         & 
& cplhhhhUSdcUSd(2,2,6,6),cplHpmUSdcHpmcUSd(2,6,2,6),cplUSdSdcUSdcSd(6,6,6,6),           & 
& cplUSdSecUSdcSe(6,6,6,6),cplUSdSucUSdcSu(6,6,6,6),cplUSdSvcUSdcSv(6,3,6,3),            & 
& cplUSdcUSdVGVG(6,6),cplUSdcUSdVPVP(6,6),cplUSdcUSdcVWmVWm(6,6),cplUSdcUSdVZVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(6,6) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumI(6,6) 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! Sd, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSd2(i1),MAh2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplAhSdcUSd(i2,i1,gO1)
coup2 = Conjg(cplAhSdcUSd(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fu, Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 2
 G0m2 = Real(SA_Gloop(p2,MFu2(i1),MCha2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*MCha(i2)*Real(SA_B0(p2,MFu2(i1),MCha2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplChaFucUSdL(i2,i1,gO1)
coupR1 = cplChaFucUSdR(i2,i1,gO1)
coupL2 =  Conjg(cplChaFucUSdL(i2,i1,gO2))
coupR2 =  Conjg(cplChaFucUSdR(i2,i1,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fd, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 4
 G0m2 = Real(SA_Gloop(p2,MFd2(i1),MChi2(i2)),dp) 
B0m2 = -2._dp*MFd(i1)*MChi(i2)*Real(SA_B0(p2,MFd2(i1),MChi2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplChiFdcUSdL(i2,i1,gO1)
coupR1 = cplChiFdcUSdR(i2,i1,gO1)
coupL2 =  Conjg(cplChiFdcUSdL(i2,i1,gO2))
coupR2 =  Conjg(cplChiFdcUSdR(i2,i1,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,0._dp,MFd2(i2)),dp) 
B0m2 = -2._dp*0.*MFd(i2)*Real(SA_B0(p2,0._dp,MFd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplFvFdcUSdL(i1,i2,gO1)
coupR1 = cplFvFdcUSdR(i1,i2,gO1)
coupL2 =  Conjg(cplFvFdcUSdL(i1,i2,gO2))
coupR2 =  Conjg(cplFvFdcUSdR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Glu, Fd 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MGlu2,MFd2(i2)),dp) 
B0m2 = -2._dp*MGlu*MFd(i2)*Real(SA_B0(p2,MGlu2,MFd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplGluFdcUSdL(i2,gO1)
coupR1 = cplGluFdcUSdR(i2,gO1)
coupL2 =  Conjg(cplGluFdcUSdL(i2,gO2))
coupR2 =  Conjg(cplGluFdcUSdR(i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +4._dp/3._dp* SumI  
    End Do 
 !------------------------ 
! Fu, Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFu2(i1),MFe2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*MFe(i2)*Real(SA_B0(p2,MFu2(i1),MFe2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplFeFucUSdL(i2,i1,gO1)
coupR1 = cplFeFucUSdR(i2,i1,gO1)
coupL2 =  Conjg(cplFeFucUSdL(i2,i1,gO2))
coupR2 =  Conjg(cplFeFucUSdR(i2,i1,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Sd, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSd2(i1),Mhh2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplhhSdcUSd(i2,i1,gO1)
coup2 = Conjg(cplhhSdcUSd(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Su, Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSu2(i1),MHpm2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplHpmSucUSd(i2,i1,gO1)
coup2 = Conjg(cplHpmSucUSd(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSv2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcUSdcSv(i2,gO1,i1)
coup2 = Conjg(cplSdcUSdcSv(i2,gO2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VG, Sd 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSd2(i2),0._dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcUSdVG(i2,gO1)
coup2 =  Conjg(cplSdcUSdVG(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +4._dp/3._dp* SumI  
    End Do 
 !------------------------ 
! VP, Sd 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSd2(i2),0._dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcUSdVP(i2,gO1)
coup2 =  Conjg(cplSdcUSdVP(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VZ, Sd 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSd2(i2),MVZ2) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcUSdVZ(i2,gO1)
coup2 =  Conjg(cplSdcUSdVZ(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! Su, Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSu2(i1),MSe2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSeSucUSd(i2,i1,gO1)
coup2 = Conjg(cplSeSucUSd(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, Su 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSu2(i2),MVWm2) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSucUSdVWm(i2,gO1)
coup2 =  Conjg(cplSucUSdVWm(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MAh2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplAhAhUSdcUSd(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(Mhh2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplhhhhUSdcUSd(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MHpm2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplHpmUSdcHpmcUSd(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSd2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSdSdcUSdcSd(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSe2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSdSecUSdcSe(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSu2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSdSucUSdcSu(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
 A0m2 = SA_A0(MSv2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSdSvcUSdcSv(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! VG, VG 
!------------------------ 
sumI = 0._dp 
 
!------------------------ 
! VP, VP 
!------------------------ 
sumI = 0._dp 
 
!------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVWm2) + 0.25_dp*RXi*SA_A0(MVWm2*RXi) - 0.5_dp*MVWm2*rMS 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSdcUSdcVWmVWm(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +4._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVZ2) + 0.25_dp*RXi*SA_A0(MVZ2*RXi) - 0.5_dp*MVZ2*rMS 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSdcUSdVZVZ(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +2._dp* SumI  


Do gO2 = 1, 6
  Do gO1 = gO2+1, 6
     res(gO1,gO2) = Conjg(res(gO2,gO1))  
   End Do 
End Do 
 
res = oo16pi2*res 
 
End Subroutine Pi1LoopSd 
 
Subroutine OneLoopSv(g1,g2,ml2,vd,vu,MChi,MChi2,MFd,MFd2,MCha,MCha2,MFe,              & 
& MFe2,MSv,MSv2,Mhh,Mhh2,MSd,MSd2,MHpm,MHpm2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,               & 
& MAh,MAh2,MSu,MSu2,cplChiFvcUSvL,cplChiFvcUSvR,cplcFdFdcUSvL,cplcFdFdcUSvR,             & 
& cplcChaFecUSvL,cplcChaFecUSvR,cplcFeFecUSvL,cplcFeFecUSvR,cplhhSvcUSv,cplSdcSdcUSv,    & 
& cplSecHpmcUSv,cplSecSecUSv,cplSecUSvcVWm,cplSvcUSvVZ,cplAhAhUSvcUSv,cplhhhhUSvcUSv,    & 
& cplHpmUSvcHpmcUSv,cplSdUSvcSdcUSv,cplSeUSvcSecUSv,cplSuUSvcSucUSv,cplUSvSvcUSvcSv,     & 
& cplUSvcUSvcVWmVWm,cplUSvcUSvVZVZ,delta,mass,mass2,RS,kont)

Implicit None 
Real(dp), Intent(in) :: MChi(4),MChi2(4),MFd(3),MFd2(3),MCha(2),MCha2(2),MFe(3),MFe2(3),MSv(3),               & 
& MSv2(3),Mhh(2),Mhh2(2),MSd(6),MSd2(6),MHpm(2),MHpm2(2),MSe(6),MSe2(6),MVWm,            & 
& MVWm2,MVZ,MVZ2,MAh(2),MAh2(2),MSu(6),MSu2(6)

Real(dp), Intent(in) :: g1,g2,vd,vu

Complex(dp), Intent(in) :: ml2(3,3)

Complex(dp), Intent(in) :: cplChiFvcUSvL(4,3,3),cplChiFvcUSvR(4,3,3),cplcFdFdcUSvL(3,3,3),cplcFdFdcUSvR(3,3,3),  & 
& cplcChaFecUSvL(2,3,3),cplcChaFecUSvR(2,3,3),cplcFeFecUSvL(3,3,3),cplcFeFecUSvR(3,3,3), & 
& cplhhSvcUSv(2,3,3),cplSdcSdcUSv(6,6,3),cplSecHpmcUSv(6,2,3),cplSecSecUSv(6,6,3),       & 
& cplSecUSvcVWm(6,3),cplSvcUSvVZ(3,3),cplAhAhUSvcUSv(2,2,3,3),cplhhhhUSvcUSv(2,2,3,3),   & 
& cplHpmUSvcHpmcUSv(2,3,2,3),cplSdUSvcSdcUSv(6,3,6,3),cplSeUSvcSecUSv(6,3,6,3),          & 
& cplSuUSvcSucUSv(6,3,6,3),cplUSvSvcUSvcSv(3,3,3,3),cplUSvcUSvcVWmVWm(3,3),              & 
& cplUSvcUSvVZVZ(3,3)

Complex(dp) :: mat2a(3,3), mat2(3,3),  PiSf(3,3,3)
Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(3), test_m2(3),p2, test(3) 
Real(dp), Intent(out) :: mass(3), mass2(3) 
Complex(dp), Intent(out) ::  RS(3,3) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopSv'
 
mat2a(1,1) = 0._dp 
mat2a(1,1) = mat2a(1,1)+(g1**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(g2**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g1**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g2**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+ml2(1,1)
mat2a(1,2) = 0._dp 
mat2a(1,2) = mat2a(1,2)+ml2(1,2)
mat2a(1,3) = 0._dp 
mat2a(1,3) = mat2a(1,3)+ml2(1,3)
mat2a(2,2) = 0._dp 
mat2a(2,2) = mat2a(2,2)+(g1**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g2**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)-(g1**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)-(g2**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+ml2(2,2)
mat2a(2,3) = 0._dp 
mat2a(2,3) = mat2a(2,3)+ml2(2,3)
mat2a(3,3) = 0._dp 
mat2a(3,3) = mat2a(3,3)+(g1**2*vd**2)/8._dp
mat2a(3,3) = mat2a(3,3)+(g2**2*vd**2)/8._dp
mat2a(3,3) = mat2a(3,3)-(g1**2*vu**2)/8._dp
mat2a(3,3) = mat2a(3,3)-(g2**2*vu**2)/8._dp
mat2a(3,3) = mat2a(3,3)+ml2(3,3)

 
 Do i1=2,3
  Do i2 = 1, i1-1 
  mat2a(i1,i2) = Conjg(mat2a(i2,i1)) 
  End do 
End do 

 
! Rotation matrix for p2=0 
PiSf(1,:,:) = ZeroC 
p2 = 0._dp 
Call Pi1LoopSv(p2,MChi,MChi2,MFd,MFd2,MCha,MCha2,MFe,MFe2,MSv,MSv2,Mhh,               & 
& Mhh2,MSd,MSd2,MHpm,MHpm2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,MAh,MAh2,MSu,MSu2,               & 
& cplChiFvcUSvL,cplChiFvcUSvR,cplcFdFdcUSvL,cplcFdFdcUSvR,cplcChaFecUSvL,cplcChaFecUSvR, & 
& cplcFeFecUSvL,cplcFeFecUSvR,cplhhSvcUSv,cplSdcSdcUSv,cplSecHpmcUSv,cplSecSecUSv,       & 
& cplSecUSvcVWm,cplSvcUSvVZ,cplAhAhUSvcUSv,cplhhhhUSvcUSv,cplHpmUSvcHpmcUSv,             & 
& cplSdUSvcSdcUSv,cplSeUSvcSecUSv,cplSuUSvcSucUSv,cplUSvSvcUSvcSv,cplUSvcUSvcVWmVWm,     & 
& cplUSvcUSvVZVZ,kont,PiSf(1,:,:))

mat2 = mat2a - Real(PiSf(1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZVOS_0 = RS 
 
 
! Now with momenta 
Do i1=1,3
PiSf(i1,:,:) = ZeroC 
p2 = MSv2(i1)
Call Pi1LoopSv(p2,MChi,MChi2,MFd,MFd2,MCha,MCha2,MFe,MFe2,MSv,MSv2,Mhh,               & 
& Mhh2,MSd,MSd2,MHpm,MHpm2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,MAh,MAh2,MSu,MSu2,               & 
& cplChiFvcUSvL,cplChiFvcUSvR,cplcFdFdcUSvL,cplcFdFdcUSvR,cplcChaFecUSvL,cplcChaFecUSvR, & 
& cplcFeFecUSvL,cplcFeFecUSvR,cplhhSvcUSv,cplSdcSdcUSv,cplSecHpmcUSv,cplSecSecUSv,       & 
& cplSecUSvcVWm,cplSvcUSvVZ,cplAhAhUSvcUSv,cplhhhhUSvcUSv,cplHpmUSvcHpmcUSv,             & 
& cplSdUSvcSdcUSv,cplSeUSvcSecUSv,cplSuUSvcSucUSv,cplUSvSvcUSvcSv,cplUSvcUSvcVWmVWm,     & 
& cplUSvcUSvVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=3,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
 
Do i1=1,3
  If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
  If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = Sqrt(mass2(i1)) 
  Else 
   If (ErrorLevel.Ge.0) Then 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
   Call TerminateProgram 
   End If 
   kont = -301 
   mass(i1) = 0._dp 
  End If 
End Do 
 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
Do i1=1,3
PiSf(i1,:,:) = ZeroC 
p2 =  mass2(i1) 
Call Pi1LoopSv(p2,MChi,MChi2,MFd,MFd2,MCha,MCha2,MFe,MFe2,MSv,MSv2,Mhh,               & 
& Mhh2,MSd,MSd2,MHpm,MHpm2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,MAh,MAh2,MSu,MSu2,               & 
& cplChiFvcUSvL,cplChiFvcUSvR,cplcFdFdcUSvL,cplcFdFdcUSvR,cplcChaFecUSvL,cplcChaFecUSvR, & 
& cplcFeFecUSvL,cplcFeFecUSvR,cplhhSvcUSv,cplSdcSdcUSv,cplSecHpmcUSv,cplSecSecUSv,       & 
& cplSecUSvcVWm,cplSvcUSvVZ,cplAhAhUSvcUSv,cplhhhhUSvcUSv,cplHpmUSvcHpmcUSv,             & 
& cplSdUSvcSdcUSv,cplSeUSvcSecUSv,cplSuUSvcSucUSv,cplUSvSvcUSvcSv,cplUSvcUSvcVWmVWm,     & 
& cplUSvcUSvVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=3,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZVOS_p2(i1,:) = RS(i1,:) 
 If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
Do i1=1,3
 If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
 If (test_m2(i1).Ne.0._dp) Then 
    test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1) 
 Else 
    test_m2(i1) = Abs(mass2(i1)) 
 End If 
 If (Abs(mass2(i1)).lt.1.0E-30_dp) test_m2(i1) = 0._dp 
 If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = sqrt(mass2(i1)) 
  Else 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses occurred a negative mass squared!' 
     Write(*,*) 'generation: ',i1 
     Write(*,*) 'mass: ',mass2(i1) 
   SignOfMassChanged = .True. 
   mass(i1) = 0._dp 
  End If 
End Do 
 
If (Maxval(test_m2).LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopSv
 
 
Subroutine Pi1LoopSv(p2,MChi,MChi2,MFd,MFd2,MCha,MCha2,MFe,MFe2,MSv,MSv2,             & 
& Mhh,Mhh2,MSd,MSd2,MHpm,MHpm2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,MAh,MAh2,MSu,MSu2,           & 
& cplChiFvcUSvL,cplChiFvcUSvR,cplcFdFdcUSvL,cplcFdFdcUSvR,cplcChaFecUSvL,cplcChaFecUSvR, & 
& cplcFeFecUSvL,cplcFeFecUSvR,cplhhSvcUSv,cplSdcSdcUSv,cplSecHpmcUSv,cplSecSecUSv,       & 
& cplSecUSvcVWm,cplSvcUSvVZ,cplAhAhUSvcUSv,cplhhhhUSvcUSv,cplHpmUSvcHpmcUSv,             & 
& cplSdUSvcSdcUSv,cplSeUSvcSecUSv,cplSuUSvcSucUSv,cplUSvSvcUSvcSv,cplUSvcUSvcVWmVWm,     & 
& cplUSvcUSvVZVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MChi(4),MChi2(4),MFd(3),MFd2(3),MCha(2),MCha2(2),MFe(3),MFe2(3),MSv(3),               & 
& MSv2(3),Mhh(2),Mhh2(2),MSd(6),MSd2(6),MHpm(2),MHpm2(2),MSe(6),MSe2(6),MVWm,            & 
& MVWm2,MVZ,MVZ2,MAh(2),MAh2(2),MSu(6),MSu2(6)

Complex(dp), Intent(in) :: cplChiFvcUSvL(4,3,3),cplChiFvcUSvR(4,3,3),cplcFdFdcUSvL(3,3,3),cplcFdFdcUSvR(3,3,3),  & 
& cplcChaFecUSvL(2,3,3),cplcChaFecUSvR(2,3,3),cplcFeFecUSvL(3,3,3),cplcFeFecUSvR(3,3,3), & 
& cplhhSvcUSv(2,3,3),cplSdcSdcUSv(6,6,3),cplSecHpmcUSv(6,2,3),cplSecSecUSv(6,6,3),       & 
& cplSecUSvcVWm(6,3),cplSvcUSvVZ(3,3),cplAhAhUSvcUSv(2,2,3,3),cplhhhhUSvcUSv(2,2,3,3),   & 
& cplHpmUSvcHpmcUSv(2,3,2,3),cplSdUSvcSdcUSv(6,3,6,3),cplSeUSvcSecUSv(6,3,6,3),          & 
& cplSuUSvcSucUSv(6,3,6,3),cplUSvSvcUSvcSv(3,3,3,3),cplUSvcUSvcVWmVWm(3,3),              & 
& cplUSvcUSvVZVZ(3,3)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumI(3,3) 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! Fv, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 4
 G0m2 = Real(SA_Gloop(p2,0._dp,MChi2(i2)),dp) 
B0m2 = -2._dp*0.*MChi(i2)*Real(SA_B0(p2,0._dp,MChi2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coupL1 = cplChiFvcUSvL(i2,i1,gO1)
coupR1 = cplChiFvcUSvR(i2,i1,gO1)
coupL2 =  Conjg(cplChiFvcUSvL(i2,i1,gO2))
coupR2 =  Conjg(cplChiFvcUSvR(i2,i1,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = -2._dp*MFd(i1)*MFd(i2)*Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coupL1 = cplcFdFdcUSvL(i1,i2,gO1)
coupR1 = cplcFdFdcUSvR(i1,i2,gO1)
coupL2 =  Conjg(cplcFdFdcUSvL(i1,i2,gO2))
coupR2 =  Conjg(cplcFdFdcUSvR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +3._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MCha2(i1),MFe2(i2)),dp) 
B0m2 = -2._dp*MCha(i1)*MFe(i2)*Real(SA_B0(p2,MCha2(i1),MFe2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coupL1 = cplcChaFecUSvL(i1,i2,gO1)
coupR1 = cplcChaFecUSvR(i1,i2,gO1)
coupL2 =  Conjg(cplcChaFecUSvL(i1,i2,gO2))
coupR2 =  Conjg(cplcChaFecUSvR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = -2._dp*MFe(i1)*MFe(i2)*Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coupL1 = cplcFeFecUSvL(i1,i2,gO1)
coupR1 = cplcFeFecUSvR(i1,i2,gO1)
coupL2 =  Conjg(cplcFeFecUSvL(i1,i2,gO2))
coupR2 =  Conjg(cplcFeFecUSvR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Sv, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSv2(i1),Mhh2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplhhSvcUSv(i2,i1,gO1)
coup2 = Conjg(cplhhSvcUSv(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSd2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSdcSdcUSv(i2,i1,gO1)
coup2 = Conjg(cplSdcSdcUSv(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +3._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MHpm2(i1),MSe2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSecHpmcUSv(i2,i1,gO1)
coup2 = Conjg(cplSecHpmcUSv(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSe2(i1),MSe2(i2)),dp) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSecSecUSv(i2,i1,gO1)
coup2 = Conjg(cplSecSecUSv(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Se 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSe2(i2),MVWm2) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSecUSvcVWm(i2,gO1)
coup2 =  Conjg(cplSecUSvcVWm(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VZ, Sv 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 3
 F0m2 = FloopRXi(p2,MSv2(i2),MVZ2) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSvcUSvVZ(i2,gO1)
coup2 =  Conjg(cplSvcUSvVZ(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MAh2(i1)) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplAhAhUSvcUSv(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(Mhh2(i1)) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplhhhhUSvcUSv(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MHpm2(i1)) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplHpmUSvcHpmcUSv(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSd2(i1)) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSdUSvcSdcUSv(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSe2(i1)) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSeUSvcSecUSv(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSu2(i1)) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplSuUSvcSucUSv(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
 A0m2 = SA_A0(MSv2(i1)) 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplUSvSvcUSvcSv(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVWm2) + 0.25_dp*RXi*SA_A0(MVWm2*RXi) - 0.5_dp*MVWm2*rMS 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplUSvcUSvcVWmVWm(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +4._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVZ2) + 0.25_dp*RXi*SA_A0(MVZ2*RXi) - 0.5_dp*MVZ2*rMS 
Do gO1 = 1, 3
  Do gO2 = gO1, 3
coup1 = cplUSvcUSvVZVZ(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +2._dp* SumI  


Do gO2 = 1, 3
  Do gO1 = gO2+1, 3
     res(gO1,gO2) = Conjg(res(gO2,gO1))  
   End Do 
End Do 
 
res = oo16pi2*res 
 
End Subroutine Pi1LoopSv 
 
Subroutine OneLoopSu(g1,g2,Mu,Yu,Tu,mq2,mu2,vd,vu,MSu,MSu2,MAh,MAh2,MFu,              & 
& MFu2,MChi,MChi2,MCha,MCha2,MFd,MFd2,MFe,MFe2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,           & 
& MSd,MSd2,MSe,MSe2,MVWm,MVWm2,MVZ,MVZ2,MSv,MSv2,cplAhSucUSu,cplChiFucUSuL,              & 
& cplChiFucUSuR,cplcChaFdcUSuL,cplcChaFdcUSuR,cplcFeFdcUSuL,cplcFeFdcUSuR,               & 
& cplGluFucUSuL,cplGluFucUSuR,cplhhSucUSu,cplSdcHpmcUSu,cplSdcSecUSu,cplSdcUSucVWm,      & 
& cplSucUSuVG,cplSucUSuVP,cplSucUSuVZ,cplAhAhUSucUSu,cplhhhhUSucUSu,cplHpmUSucHpmcUSu,   & 
& cplSdUSucSdcUSu,cplSeUSucSecUSu,cplUSuSucUSucSu,cplUSuSvcUSucSv,cplUSucUSuVGVG,        & 
& cplUSucUSuVPVP,cplUSucUSucVWmVWm,cplUSucUSuVZVZ,delta,mass,mass2,RS,kont)

Implicit None 
Real(dp), Intent(in) :: MSu(6),MSu2(6),MAh(2),MAh2(2),MFu(3),MFu2(3),MChi(4),MChi2(4),MCha(2),MCha2(2),       & 
& MFd(3),MFd2(3),MFe(3),MFe2(3),MGlu,MGlu2,Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),              & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MVWm,MVWm2,MVZ,MVZ2,MSv(3),MSv2(3)

Real(dp), Intent(in) :: g1,g2,vd,vu

Complex(dp), Intent(in) :: Mu,Yu(3,3),Tu(3,3),mq2(3,3),mu2(3,3)

Complex(dp), Intent(in) :: cplAhSucUSu(2,6,6),cplChiFucUSuL(4,3,6),cplChiFucUSuR(4,3,6),cplcChaFdcUSuL(2,3,6),   & 
& cplcChaFdcUSuR(2,3,6),cplcFeFdcUSuL(3,3,6),cplcFeFdcUSuR(3,3,6),cplGluFucUSuL(3,6),    & 
& cplGluFucUSuR(3,6),cplhhSucUSu(2,6,6),cplSdcHpmcUSu(6,2,6),cplSdcSecUSu(6,6,6),        & 
& cplSdcUSucVWm(6,6),cplSucUSuVG(6,6),cplSucUSuVP(6,6),cplSucUSuVZ(6,6),cplAhAhUSucUSu(2,2,6,6),& 
& cplhhhhUSucUSu(2,2,6,6),cplHpmUSucHpmcUSu(2,6,2,6),cplSdUSucSdcUSu(6,6,6,6),           & 
& cplSeUSucSecUSu(6,6,6,6),cplUSuSucUSucSu(6,6,6,6),cplUSuSvcUSucSv(6,3,6,3),            & 
& cplUSucUSuVGVG(6,6),cplUSucUSuVPVP(6,6),cplUSucUSucVWmVWm(6,6),cplUSucUSuVZVZ(6,6)

Complex(dp) :: mat2a(6,6), mat2(6,6),  PiSf(6,6,6)
Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(6), test_m2(6),p2, test(6) 
Real(dp), Intent(out) :: mass(6), mass2(6) 
Complex(dp), Intent(out) ::  RS(6,6) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopSu'
 
mat2a(1,1) = 0._dp 
mat2a(1,1) = mat2a(1,1)-(g1**2*vd**2)/24._dp
mat2a(1,1) = mat2a(1,1)+(g2**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(g1**2*vu**2)/24._dp
mat2a(1,1) = mat2a(1,1)-(g2**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+mq2(1,1)
Do j1 = 1,3
mat2a(1,1) = mat2a(1,1)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,1))/2._dp
End Do 
mat2a(1,2) = 0._dp 
mat2a(1,2) = mat2a(1,2)+mq2(1,2)
Do j1 = 1,3
mat2a(1,2) = mat2a(1,2)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,2))/2._dp
End Do 
mat2a(1,3) = 0._dp 
mat2a(1,3) = mat2a(1,3)+mq2(1,3)
Do j1 = 1,3
mat2a(1,3) = mat2a(1,3)+(vu**2*Conjg(Yu(j1,1))*Yu(j1,3))/2._dp
End Do 
mat2a(1,4) = 0._dp 
mat2a(1,4) = mat2a(1,4)-((vd*Mu*Conjg(Yu(1,1)))/sqrt(2._dp))
mat2a(1,4) = mat2a(1,4)+(vu*Conjg(Tu(1,1)))/sqrt(2._dp)
mat2a(1,5) = 0._dp 
mat2a(1,5) = mat2a(1,5)-((vd*Mu*Conjg(Yu(2,1)))/sqrt(2._dp))
mat2a(1,5) = mat2a(1,5)+(vu*Conjg(Tu(2,1)))/sqrt(2._dp)
mat2a(1,6) = 0._dp 
mat2a(1,6) = mat2a(1,6)-((vd*Mu*Conjg(Yu(3,1)))/sqrt(2._dp))
mat2a(1,6) = mat2a(1,6)+(vu*Conjg(Tu(3,1)))/sqrt(2._dp)
mat2a(2,2) = 0._dp 
mat2a(2,2) = mat2a(2,2)-(g1**2*vd**2)/24._dp
mat2a(2,2) = mat2a(2,2)+(g2**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g1**2*vu**2)/24._dp
mat2a(2,2) = mat2a(2,2)-(g2**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+mq2(2,2)
Do j1 = 1,3
mat2a(2,2) = mat2a(2,2)+(vu**2*Conjg(Yu(j1,2))*Yu(j1,2))/2._dp
End Do 
mat2a(2,3) = 0._dp 
mat2a(2,3) = mat2a(2,3)+mq2(2,3)
Do j1 = 1,3
mat2a(2,3) = mat2a(2,3)+(vu**2*Conjg(Yu(j1,2))*Yu(j1,3))/2._dp
End Do 
mat2a(2,4) = 0._dp 
mat2a(2,4) = mat2a(2,4)-((vd*Mu*Conjg(Yu(1,2)))/sqrt(2._dp))
mat2a(2,4) = mat2a(2,4)+(vu*Conjg(Tu(1,2)))/sqrt(2._dp)
mat2a(2,5) = 0._dp 
mat2a(2,5) = mat2a(2,5)-((vd*Mu*Conjg(Yu(2,2)))/sqrt(2._dp))
mat2a(2,5) = mat2a(2,5)+(vu*Conjg(Tu(2,2)))/sqrt(2._dp)
mat2a(2,6) = 0._dp 
mat2a(2,6) = mat2a(2,6)-((vd*Mu*Conjg(Yu(3,2)))/sqrt(2._dp))
mat2a(2,6) = mat2a(2,6)+(vu*Conjg(Tu(3,2)))/sqrt(2._dp)
mat2a(3,3) = 0._dp 
mat2a(3,3) = mat2a(3,3)-(g1**2*vd**2)/24._dp
mat2a(3,3) = mat2a(3,3)+(g2**2*vd**2)/8._dp
mat2a(3,3) = mat2a(3,3)+(g1**2*vu**2)/24._dp
mat2a(3,3) = mat2a(3,3)-(g2**2*vu**2)/8._dp
mat2a(3,3) = mat2a(3,3)+mq2(3,3)
Do j1 = 1,3
mat2a(3,3) = mat2a(3,3)+(vu**2*Conjg(Yu(j1,3))*Yu(j1,3))/2._dp
End Do 
mat2a(3,4) = 0._dp 
mat2a(3,4) = mat2a(3,4)-((vd*Mu*Conjg(Yu(1,3)))/sqrt(2._dp))
mat2a(3,4) = mat2a(3,4)+(vu*Conjg(Tu(1,3)))/sqrt(2._dp)
mat2a(3,5) = 0._dp 
mat2a(3,5) = mat2a(3,5)-((vd*Mu*Conjg(Yu(2,3)))/sqrt(2._dp))
mat2a(3,5) = mat2a(3,5)+(vu*Conjg(Tu(2,3)))/sqrt(2._dp)
mat2a(3,6) = 0._dp 
mat2a(3,6) = mat2a(3,6)-((vd*Mu*Conjg(Yu(3,3)))/sqrt(2._dp))
mat2a(3,6) = mat2a(3,6)+(vu*Conjg(Tu(3,3)))/sqrt(2._dp)
mat2a(4,4) = 0._dp 
mat2a(4,4) = mat2a(4,4)+(g1**2*vd**2)/6._dp
mat2a(4,4) = mat2a(4,4)-(g1**2*vu**2)/6._dp
mat2a(4,4) = mat2a(4,4)+mu2(1,1)
Do j1 = 1,3
mat2a(4,4) = mat2a(4,4)+(vu**2*Conjg(Yu(1,j1))*Yu(1,j1))/2._dp
End Do 
mat2a(4,5) = 0._dp 
mat2a(4,5) = mat2a(4,5)+mu2(1,2)
Do j1 = 1,3
mat2a(4,5) = mat2a(4,5)+(vu**2*Conjg(Yu(2,j1))*Yu(1,j1))/2._dp
End Do 
mat2a(4,6) = 0._dp 
mat2a(4,6) = mat2a(4,6)+mu2(1,3)
Do j1 = 1,3
mat2a(4,6) = mat2a(4,6)+(vu**2*Conjg(Yu(3,j1))*Yu(1,j1))/2._dp
End Do 
mat2a(5,5) = 0._dp 
mat2a(5,5) = mat2a(5,5)+(g1**2*vd**2)/6._dp
mat2a(5,5) = mat2a(5,5)-(g1**2*vu**2)/6._dp
mat2a(5,5) = mat2a(5,5)+mu2(2,2)
Do j1 = 1,3
mat2a(5,5) = mat2a(5,5)+(vu**2*Conjg(Yu(2,j1))*Yu(2,j1))/2._dp
End Do 
mat2a(5,6) = 0._dp 
mat2a(5,6) = mat2a(5,6)+mu2(2,3)
Do j1 = 1,3
mat2a(5,6) = mat2a(5,6)+(vu**2*Conjg(Yu(3,j1))*Yu(2,j1))/2._dp
End Do 
mat2a(6,6) = 0._dp 
mat2a(6,6) = mat2a(6,6)+(g1**2*vd**2)/6._dp
mat2a(6,6) = mat2a(6,6)-(g1**2*vu**2)/6._dp
mat2a(6,6) = mat2a(6,6)+mu2(3,3)
Do j1 = 1,3
mat2a(6,6) = mat2a(6,6)+(vu**2*Conjg(Yu(3,j1))*Yu(3,j1))/2._dp
End Do 

 
 Do i1=2,6
  Do i2 = 1, i1-1 
  mat2a(i1,i2) = Conjg(mat2a(i2,i1)) 
  End do 
End do 

 
! Rotation matrix for p2=0 
PiSf(1,:,:) = ZeroC 
p2 = 0._dp 
Call Pi1LoopSu(p2,MSu,MSu2,MAh,MAh2,MFu,MFu2,MChi,MChi2,MCha,MCha2,MFd,               & 
& MFd2,MFe,MFe2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MVWm,MVWm2,             & 
& MVZ,MVZ2,MSv,MSv2,cplAhSucUSu,cplChiFucUSuL,cplChiFucUSuR,cplcChaFdcUSuL,              & 
& cplcChaFdcUSuR,cplcFeFdcUSuL,cplcFeFdcUSuR,cplGluFucUSuL,cplGluFucUSuR,cplhhSucUSu,    & 
& cplSdcHpmcUSu,cplSdcSecUSu,cplSdcUSucVWm,cplSucUSuVG,cplSucUSuVP,cplSucUSuVZ,          & 
& cplAhAhUSucUSu,cplhhhhUSucUSu,cplHpmUSucHpmcUSu,cplSdUSucSdcUSu,cplSeUSucSecUSu,       & 
& cplUSuSucUSucSu,cplUSuSvcUSucSv,cplUSucUSuVGVG,cplUSucUSuVPVP,cplUSucUSucVWmVWm,       & 
& cplUSucUSuVZVZ,kont,PiSf(1,:,:))

mat2 = mat2a - Real(PiSf(1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZUOS_0 = RS 
 
 
! Now with momenta 
Do i1=1,6
PiSf(i1,:,:) = ZeroC 
p2 = MSu2(i1)
Call Pi1LoopSu(p2,MSu,MSu2,MAh,MAh2,MFu,MFu2,MChi,MChi2,MCha,MCha2,MFd,               & 
& MFd2,MFe,MFe2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MVWm,MVWm2,             & 
& MVZ,MVZ2,MSv,MSv2,cplAhSucUSu,cplChiFucUSuL,cplChiFucUSuR,cplcChaFdcUSuL,              & 
& cplcChaFdcUSuR,cplcFeFdcUSuL,cplcFeFdcUSuR,cplGluFucUSuL,cplGluFucUSuR,cplhhSucUSu,    & 
& cplSdcHpmcUSu,cplSdcSecUSu,cplSdcUSucVWm,cplSucUSuVG,cplSucUSuVP,cplSucUSuVZ,          & 
& cplAhAhUSucUSu,cplhhhhUSucUSu,cplHpmUSucHpmcUSu,cplSdUSucSdcUSu,cplSeUSucSecUSu,       & 
& cplUSuSucUSucSu,cplUSuSvcUSucSv,cplUSucUSuVGVG,cplUSucUSuVPVP,cplUSucUSucVWmVWm,       & 
& cplUSucUSuVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=6,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
 
Do i1=1,6
  If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
  If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = Sqrt(mass2(i1)) 
  Else 
   If (ErrorLevel.Ge.0) Then 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
   Call TerminateProgram 
   End If 
   kont = -301 
   mass(i1) = 0._dp 
  End If 
End Do 
 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
Do i1=1,6
PiSf(i1,:,:) = ZeroC 
p2 =  mass2(i1) 
Call Pi1LoopSu(p2,MSu,MSu2,MAh,MAh2,MFu,MFu2,MChi,MChi2,MCha,MCha2,MFd,               & 
& MFd2,MFe,MFe2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MVWm,MVWm2,             & 
& MVZ,MVZ2,MSv,MSv2,cplAhSucUSu,cplChiFucUSuL,cplChiFucUSuR,cplcChaFdcUSuL,              & 
& cplcChaFdcUSuR,cplcFeFdcUSuL,cplcFeFdcUSuR,cplGluFucUSuL,cplGluFucUSuR,cplhhSucUSu,    & 
& cplSdcHpmcUSu,cplSdcSecUSu,cplSdcUSucVWm,cplSucUSuVG,cplSucUSuVP,cplSucUSuVZ,          & 
& cplAhAhUSucUSu,cplhhhhUSucUSu,cplHpmUSucHpmcUSu,cplSdUSucSdcUSu,cplSeUSucSecUSu,       & 
& cplUSuSucUSucSu,cplUSuSvcUSucSv,cplUSucUSuVGVG,cplUSucUSuVPVP,cplUSucUSucVWmVWm,       & 
& cplUSucUSuVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=6,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZUOS_p2(i1,:) = RS(i1,:) 
 If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
Do i1=1,6
 If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
 If (test_m2(i1).Ne.0._dp) Then 
    test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1) 
 Else 
    test_m2(i1) = Abs(mass2(i1)) 
 End If 
 If (Abs(mass2(i1)).lt.1.0E-30_dp) test_m2(i1) = 0._dp 
 If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = sqrt(mass2(i1)) 
  Else 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses occurred a negative mass squared!' 
     Write(*,*) 'generation: ',i1 
     Write(*,*) 'mass: ',mass2(i1) 
   SignOfMassChanged = .True. 
   mass(i1) = 0._dp 
  End If 
End Do 
 
If (Maxval(test_m2).LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopSu
 
 
Subroutine Pi1LoopSu(p2,MSu,MSu2,MAh,MAh2,MFu,MFu2,MChi,MChi2,MCha,MCha2,             & 
& MFd,MFd2,MFe,MFe2,MGlu,MGlu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MVWm,               & 
& MVWm2,MVZ,MVZ2,MSv,MSv2,cplAhSucUSu,cplChiFucUSuL,cplChiFucUSuR,cplcChaFdcUSuL,        & 
& cplcChaFdcUSuR,cplcFeFdcUSuL,cplcFeFdcUSuR,cplGluFucUSuL,cplGluFucUSuR,cplhhSucUSu,    & 
& cplSdcHpmcUSu,cplSdcSecUSu,cplSdcUSucVWm,cplSucUSuVG,cplSucUSuVP,cplSucUSuVZ,          & 
& cplAhAhUSucUSu,cplhhhhUSucUSu,cplHpmUSucHpmcUSu,cplSdUSucSdcUSu,cplSeUSucSecUSu,       & 
& cplUSuSucUSucSu,cplUSuSvcUSucSv,cplUSucUSuVGVG,cplUSucUSuVPVP,cplUSucUSucVWmVWm,       & 
& cplUSucUSuVZVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MSu(6),MSu2(6),MAh(2),MAh2(2),MFu(3),MFu2(3),MChi(4),MChi2(4),MCha(2),MCha2(2),       & 
& MFd(3),MFd2(3),MFe(3),MFe2(3),MGlu,MGlu2,Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),              & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MVWm,MVWm2,MVZ,MVZ2,MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplAhSucUSu(2,6,6),cplChiFucUSuL(4,3,6),cplChiFucUSuR(4,3,6),cplcChaFdcUSuL(2,3,6),   & 
& cplcChaFdcUSuR(2,3,6),cplcFeFdcUSuL(3,3,6),cplcFeFdcUSuR(3,3,6),cplGluFucUSuL(3,6),    & 
& cplGluFucUSuR(3,6),cplhhSucUSu(2,6,6),cplSdcHpmcUSu(6,2,6),cplSdcSecUSu(6,6,6),        & 
& cplSdcUSucVWm(6,6),cplSucUSuVG(6,6),cplSucUSuVP(6,6),cplSucUSuVZ(6,6),cplAhAhUSucUSu(2,2,6,6),& 
& cplhhhhUSucUSu(2,2,6,6),cplHpmUSucHpmcUSu(2,6,2,6),cplSdUSucSdcUSu(6,6,6,6),           & 
& cplSeUSucSecUSu(6,6,6,6),cplUSuSucUSucSu(6,6,6,6),cplUSuSvcUSucSv(6,3,6,3),            & 
& cplUSucUSuVGVG(6,6),cplUSucUSuVPVP(6,6),cplUSucUSucVWmVWm(6,6),cplUSucUSuVZVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(6,6) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumI(6,6) 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! Su, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSu2(i1),MAh2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplAhSucUSu(i2,i1,gO1)
coup2 = Conjg(cplAhSucUSu(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fu, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 4
 G0m2 = Real(SA_Gloop(p2,MFu2(i1),MChi2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*MChi(i2)*Real(SA_B0(p2,MFu2(i1),MChi2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplChiFucUSuL(i2,i1,gO1)
coupR1 = cplChiFucUSuR(i2,i1,gO1)
coupL2 =  Conjg(cplChiFucUSuL(i2,i1,gO2))
coupR2 =  Conjg(cplChiFucUSuR(i2,i1,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MCha2(i1),MFd2(i2)),dp) 
B0m2 = -2._dp*MCha(i1)*MFd(i2)*Real(SA_B0(p2,MCha2(i1),MFd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplcChaFdcUSuL(i1,i2,gO1)
coupR1 = cplcChaFdcUSuR(i1,i2,gO1)
coupL2 =  Conjg(cplcChaFdcUSuL(i1,i2,gO2))
coupR2 =  Conjg(cplcChaFdcUSuR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFe2(i1),MFd2(i2)),dp) 
B0m2 = -2._dp*MFe(i1)*MFd(i2)*Real(SA_B0(p2,MFe2(i1),MFd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplcFeFdcUSuL(i1,i2,gO1)
coupR1 = cplcFeFdcUSuR(i1,i2,gO1)
coupL2 =  Conjg(cplcFeFdcUSuL(i1,i2,gO2))
coupR2 =  Conjg(cplcFeFdcUSuR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Glu, Fu 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MGlu2,MFu2(i2)),dp) 
B0m2 = -2._dp*MGlu*MFu(i2)*Real(SA_B0(p2,MGlu2,MFu2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplGluFucUSuL(i2,gO1)
coupR1 = cplGluFucUSuR(i2,gO1)
coupL2 =  Conjg(cplGluFucUSuL(i2,gO2))
coupR2 =  Conjg(cplGluFucUSuR(i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +4._dp/3._dp* SumI  
    End Do 
 !------------------------ 
! Su, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSu2(i1),Mhh2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplhhSucUSu(i2,i1,gO1)
coup2 = Conjg(cplhhSucUSu(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MHpm2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcHpmcUSu(i2,i1,gO1)
coup2 = Conjg(cplSdcHpmcUSu(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSe2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcSecUSu(i2,i1,gO1)
coup2 = Conjg(cplSdcSecUSu(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Sd 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSd2(i2),MVWm2) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcUSucVWm(i2,gO1)
coup2 =  Conjg(cplSdcUSucVWm(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VG, Su 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSu2(i2),0._dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSucUSuVG(i2,gO1)
coup2 =  Conjg(cplSucUSuVG(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +4._dp/3._dp* SumI  
    End Do 
 !------------------------ 
! VP, Su 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSu2(i2),0._dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSucUSuVP(i2,gO1)
coup2 =  Conjg(cplSucUSuVP(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VZ, Su 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSu2(i2),MVZ2) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSucUSuVZ(i2,gO1)
coup2 =  Conjg(cplSucUSuVZ(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MAh2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplAhAhUSucUSu(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(Mhh2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplhhhhUSucUSu(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MHpm2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplHpmUSucHpmcUSu(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSd2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdUSucSdcUSu(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSe2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSeUSucSecUSu(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSu2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSuSucUSucSu(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
 A0m2 = SA_A0(MSv2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSuSvcUSucSv(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! VG, VG 
!------------------------ 
sumI = 0._dp 
 
!------------------------ 
! VP, VP 
!------------------------ 
sumI = 0._dp 
 
!------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVWm2) + 0.25_dp*RXi*SA_A0(MVWm2*RXi) - 0.5_dp*MVWm2*rMS 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSucUSucVWmVWm(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +4._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVZ2) + 0.25_dp*RXi*SA_A0(MVZ2*RXi) - 0.5_dp*MVZ2*rMS 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSucUSuVZVZ(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +2._dp* SumI  


Do gO2 = 1, 6
  Do gO1 = gO2+1, 6
     res(gO1,gO2) = Conjg(res(gO2,gO1))  
   End Do 
End Do 
 
res = oo16pi2*res 
 
End Subroutine Pi1LoopSu 
 
Subroutine OneLoopSe(g1,g2,Mu,Ye,Te,ml2,me2,vd,vu,MSe,MSe2,MAh,MAh2,MCha,             & 
& MCha2,MFe,MFe2,MChi,MChi2,MFu,MFu2,MFd,MFd2,Mhh,Mhh2,MSv,MSv2,MHpm,MHpm2,              & 
& MSu,MSu2,MSd,MSd2,MVZ,MVZ2,MVWm,MVWm2,cplAhSecUSe,cplFvChacUSeL,cplFvChacUSeR,         & 
& cplChiFecUSeL,cplChiFecUSeR,cplcFuFdcUSeL,cplcFuFdcUSeR,cplFvFecUSeL,cplFvFecUSeR,     & 
& cplhhSecUSe,cplHpmSvcUSe,cplSdcUSecSu,cplSecUSecSv,cplSecUSeVP,cplSecUSeVZ,            & 
& cplSvcUSeVWm,cplAhAhUSecUSe,cplhhhhUSecUSe,cplHpmUSecHpmcUSe,cplSdUSecSdcUSe,          & 
& cplUSeSecUSecSe,cplUSeSucUSecSu,cplUSeSvcUSecSv,cplUSecUSeVPVP,cplUSecUSecVWmVWm,      & 
& cplUSecUSeVZVZ,delta,mass,mass2,RS,kont)

Implicit None 
Real(dp), Intent(in) :: MSe(6),MSe2(6),MAh(2),MAh2(2),MCha(2),MCha2(2),MFe(3),MFe2(3),MChi(4),MChi2(4),       & 
& MFu(3),MFu2(3),MFd(3),MFd2(3),Mhh(2),Mhh2(2),MSv(3),MSv2(3),MHpm(2),MHpm2(2),          & 
& MSu(6),MSu2(6),MSd(6),MSd2(6),MVZ,MVZ2,MVWm,MVWm2

Real(dp), Intent(in) :: g1,g2,vd,vu

Complex(dp), Intent(in) :: Mu,Ye(3,3),Te(3,3),ml2(3,3),me2(3,3)

Complex(dp), Intent(in) :: cplAhSecUSe(2,6,6),cplFvChacUSeL(3,2,6),cplFvChacUSeR(3,2,6),cplChiFecUSeL(4,3,6),    & 
& cplChiFecUSeR(4,3,6),cplcFuFdcUSeL(3,3,6),cplcFuFdcUSeR(3,3,6),cplFvFecUSeL(3,3,6),    & 
& cplFvFecUSeR(3,3,6),cplhhSecUSe(2,6,6),cplHpmSvcUSe(2,3,6),cplSdcUSecSu(6,6,6),        & 
& cplSecUSecSv(6,6,3),cplSecUSeVP(6,6),cplSecUSeVZ(6,6),cplSvcUSeVWm(3,6),               & 
& cplAhAhUSecUSe(2,2,6,6),cplhhhhUSecUSe(2,2,6,6),cplHpmUSecHpmcUSe(2,6,2,6),            & 
& cplSdUSecSdcUSe(6,6,6,6),cplUSeSecUSecSe(6,6,6,6),cplUSeSucUSecSu(6,6,6,6),            & 
& cplUSeSvcUSecSv(6,3,6,3),cplUSecUSeVPVP(6,6),cplUSecUSecVWmVWm(6,6),cplUSecUSeVZVZ(6,6)

Complex(dp) :: mat2a(6,6), mat2(6,6),  PiSf(6,6,6)
Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(6), test_m2(6),p2, test(6) 
Real(dp), Intent(out) :: mass(6), mass2(6) 
Complex(dp), Intent(out) ::  RS(6,6) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopSe'
 
mat2a(1,1) = 0._dp 
mat2a(1,1) = mat2a(1,1)+(g1**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g2**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g1**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(g2**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+ml2(1,1)
Do j1 = 1,3
mat2a(1,1) = mat2a(1,1)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,1))/2._dp
End Do 
mat2a(1,2) = 0._dp 
mat2a(1,2) = mat2a(1,2)+ml2(1,2)
Do j1 = 1,3
mat2a(1,2) = mat2a(1,2)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,2))/2._dp
End Do 
mat2a(1,3) = 0._dp 
mat2a(1,3) = mat2a(1,3)+ml2(1,3)
Do j1 = 1,3
mat2a(1,3) = mat2a(1,3)+(vd**2*Conjg(Ye(j1,1))*Ye(j1,3))/2._dp
End Do 
mat2a(1,4) = 0._dp 
mat2a(1,4) = mat2a(1,4)-((vu*Mu*Conjg(Ye(1,1)))/sqrt(2._dp))
mat2a(1,4) = mat2a(1,4)+(vd*Conjg(Te(1,1)))/sqrt(2._dp)
mat2a(1,5) = 0._dp 
mat2a(1,5) = mat2a(1,5)-((vu*Mu*Conjg(Ye(2,1)))/sqrt(2._dp))
mat2a(1,5) = mat2a(1,5)+(vd*Conjg(Te(2,1)))/sqrt(2._dp)
mat2a(1,6) = 0._dp 
mat2a(1,6) = mat2a(1,6)-((vu*Mu*Conjg(Ye(3,1)))/sqrt(2._dp))
mat2a(1,6) = mat2a(1,6)+(vd*Conjg(Te(3,1)))/sqrt(2._dp)
mat2a(2,2) = 0._dp 
mat2a(2,2) = mat2a(2,2)+(g1**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)-(g2**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)-(g1**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g2**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+ml2(2,2)
Do j1 = 1,3
mat2a(2,2) = mat2a(2,2)+(vd**2*Conjg(Ye(j1,2))*Ye(j1,2))/2._dp
End Do 
mat2a(2,3) = 0._dp 
mat2a(2,3) = mat2a(2,3)+ml2(2,3)
Do j1 = 1,3
mat2a(2,3) = mat2a(2,3)+(vd**2*Conjg(Ye(j1,2))*Ye(j1,3))/2._dp
End Do 
mat2a(2,4) = 0._dp 
mat2a(2,4) = mat2a(2,4)-((vu*Mu*Conjg(Ye(1,2)))/sqrt(2._dp))
mat2a(2,4) = mat2a(2,4)+(vd*Conjg(Te(1,2)))/sqrt(2._dp)
mat2a(2,5) = 0._dp 
mat2a(2,5) = mat2a(2,5)-((vu*Mu*Conjg(Ye(2,2)))/sqrt(2._dp))
mat2a(2,5) = mat2a(2,5)+(vd*Conjg(Te(2,2)))/sqrt(2._dp)
mat2a(2,6) = 0._dp 
mat2a(2,6) = mat2a(2,6)-((vu*Mu*Conjg(Ye(3,2)))/sqrt(2._dp))
mat2a(2,6) = mat2a(2,6)+(vd*Conjg(Te(3,2)))/sqrt(2._dp)
mat2a(3,3) = 0._dp 
mat2a(3,3) = mat2a(3,3)+(g1**2*vd**2)/8._dp
mat2a(3,3) = mat2a(3,3)-(g2**2*vd**2)/8._dp
mat2a(3,3) = mat2a(3,3)-(g1**2*vu**2)/8._dp
mat2a(3,3) = mat2a(3,3)+(g2**2*vu**2)/8._dp
mat2a(3,3) = mat2a(3,3)+ml2(3,3)
Do j1 = 1,3
mat2a(3,3) = mat2a(3,3)+(vd**2*Conjg(Ye(j1,3))*Ye(j1,3))/2._dp
End Do 
mat2a(3,4) = 0._dp 
mat2a(3,4) = mat2a(3,4)-((vu*Mu*Conjg(Ye(1,3)))/sqrt(2._dp))
mat2a(3,4) = mat2a(3,4)+(vd*Conjg(Te(1,3)))/sqrt(2._dp)
mat2a(3,5) = 0._dp 
mat2a(3,5) = mat2a(3,5)-((vu*Mu*Conjg(Ye(2,3)))/sqrt(2._dp))
mat2a(3,5) = mat2a(3,5)+(vd*Conjg(Te(2,3)))/sqrt(2._dp)
mat2a(3,6) = 0._dp 
mat2a(3,6) = mat2a(3,6)-((vu*Mu*Conjg(Ye(3,3)))/sqrt(2._dp))
mat2a(3,6) = mat2a(3,6)+(vd*Conjg(Te(3,3)))/sqrt(2._dp)
mat2a(4,4) = 0._dp 
mat2a(4,4) = mat2a(4,4)-(g1**2*vd**2)/4._dp
mat2a(4,4) = mat2a(4,4)+(g1**2*vu**2)/4._dp
mat2a(4,4) = mat2a(4,4)+me2(1,1)
Do j1 = 1,3
mat2a(4,4) = mat2a(4,4)+(vd**2*Conjg(Ye(1,j1))*Ye(1,j1))/2._dp
End Do 
mat2a(4,5) = 0._dp 
mat2a(4,5) = mat2a(4,5)+me2(1,2)
Do j1 = 1,3
mat2a(4,5) = mat2a(4,5)+(vd**2*Conjg(Ye(2,j1))*Ye(1,j1))/2._dp
End Do 
mat2a(4,6) = 0._dp 
mat2a(4,6) = mat2a(4,6)+me2(1,3)
Do j1 = 1,3
mat2a(4,6) = mat2a(4,6)+(vd**2*Conjg(Ye(3,j1))*Ye(1,j1))/2._dp
End Do 
mat2a(5,5) = 0._dp 
mat2a(5,5) = mat2a(5,5)-(g1**2*vd**2)/4._dp
mat2a(5,5) = mat2a(5,5)+(g1**2*vu**2)/4._dp
mat2a(5,5) = mat2a(5,5)+me2(2,2)
Do j1 = 1,3
mat2a(5,5) = mat2a(5,5)+(vd**2*Conjg(Ye(2,j1))*Ye(2,j1))/2._dp
End Do 
mat2a(5,6) = 0._dp 
mat2a(5,6) = mat2a(5,6)+me2(2,3)
Do j1 = 1,3
mat2a(5,6) = mat2a(5,6)+(vd**2*Conjg(Ye(3,j1))*Ye(2,j1))/2._dp
End Do 
mat2a(6,6) = 0._dp 
mat2a(6,6) = mat2a(6,6)-(g1**2*vd**2)/4._dp
mat2a(6,6) = mat2a(6,6)+(g1**2*vu**2)/4._dp
mat2a(6,6) = mat2a(6,6)+me2(3,3)
Do j1 = 1,3
mat2a(6,6) = mat2a(6,6)+(vd**2*Conjg(Ye(3,j1))*Ye(3,j1))/2._dp
End Do 

 
 Do i1=2,6
  Do i2 = 1, i1-1 
  mat2a(i1,i2) = Conjg(mat2a(i2,i1)) 
  End do 
End do 

 
! Rotation matrix for p2=0 
PiSf(1,:,:) = ZeroC 
p2 = 0._dp 
Call Pi1LoopSe(p2,MSe,MSe2,MAh,MAh2,MCha,MCha2,MFe,MFe2,MChi,MChi2,MFu,               & 
& MFu2,MFd,MFd2,Mhh,Mhh2,MSv,MSv2,MHpm,MHpm2,MSu,MSu2,MSd,MSd2,MVZ,MVZ2,MVWm,            & 
& MVWm2,cplAhSecUSe,cplFvChacUSeL,cplFvChacUSeR,cplChiFecUSeL,cplChiFecUSeR,             & 
& cplcFuFdcUSeL,cplcFuFdcUSeR,cplFvFecUSeL,cplFvFecUSeR,cplhhSecUSe,cplHpmSvcUSe,        & 
& cplSdcUSecSu,cplSecUSecSv,cplSecUSeVP,cplSecUSeVZ,cplSvcUSeVWm,cplAhAhUSecUSe,         & 
& cplhhhhUSecUSe,cplHpmUSecHpmcUSe,cplSdUSecSdcUSe,cplUSeSecUSecSe,cplUSeSucUSecSu,      & 
& cplUSeSvcUSecSv,cplUSecUSeVPVP,cplUSecUSecVWmVWm,cplUSecUSeVZVZ,kont,PiSf(1,:,:))

mat2 = mat2a - Real(PiSf(1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZEOS_0 = RS 
 
 
! Now with momenta 
Do i1=1,6
PiSf(i1,:,:) = ZeroC 
p2 = MSe2(i1)
Call Pi1LoopSe(p2,MSe,MSe2,MAh,MAh2,MCha,MCha2,MFe,MFe2,MChi,MChi2,MFu,               & 
& MFu2,MFd,MFd2,Mhh,Mhh2,MSv,MSv2,MHpm,MHpm2,MSu,MSu2,MSd,MSd2,MVZ,MVZ2,MVWm,            & 
& MVWm2,cplAhSecUSe,cplFvChacUSeL,cplFvChacUSeR,cplChiFecUSeL,cplChiFecUSeR,             & 
& cplcFuFdcUSeL,cplcFuFdcUSeR,cplFvFecUSeL,cplFvFecUSeR,cplhhSecUSe,cplHpmSvcUSe,        & 
& cplSdcUSecSu,cplSecUSecSv,cplSecUSeVP,cplSecUSeVZ,cplSvcUSeVWm,cplAhAhUSecUSe,         & 
& cplhhhhUSecUSe,cplHpmUSecHpmcUSe,cplSdUSecSdcUSe,cplUSeSecUSecSe,cplUSeSucUSecSu,      & 
& cplUSeSvcUSecSv,cplUSecUSeVPVP,cplUSecUSecVWmVWm,cplUSecUSeVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=6,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
 
Do i1=1,6
  If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
  If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = Sqrt(mass2(i1)) 
  Else 
   If (ErrorLevel.Ge.0) Then 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
   Call TerminateProgram 
   End If 
   kont = -301 
   mass(i1) = 0._dp 
  End If 
End Do 
 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
Do i1=1,6
PiSf(i1,:,:) = ZeroC 
p2 =  mass2(i1) 
Call Pi1LoopSe(p2,MSe,MSe2,MAh,MAh2,MCha,MCha2,MFe,MFe2,MChi,MChi2,MFu,               & 
& MFu2,MFd,MFd2,Mhh,Mhh2,MSv,MSv2,MHpm,MHpm2,MSu,MSu2,MSd,MSd2,MVZ,MVZ2,MVWm,            & 
& MVWm2,cplAhSecUSe,cplFvChacUSeL,cplFvChacUSeR,cplChiFecUSeL,cplChiFecUSeR,             & 
& cplcFuFdcUSeL,cplcFuFdcUSeR,cplFvFecUSeL,cplFvFecUSeR,cplhhSecUSe,cplHpmSvcUSe,        & 
& cplSdcUSecSu,cplSecUSecSv,cplSecUSeVP,cplSecUSeVZ,cplSvcUSeVWm,cplAhAhUSecUSe,         & 
& cplhhhhUSecUSe,cplHpmUSecHpmcUSe,cplSdUSecSdcUSe,cplUSeSecUSecSe,cplUSeSucUSecSu,      & 
& cplUSeSvcUSecSv,cplUSecUSeVPVP,cplUSecUSecVWmVWm,cplUSecUSeVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=6,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZEOS_p2(i1,:) = RS(i1,:) 
 If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
Do i1=1,6
 If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
 If (test_m2(i1).Ne.0._dp) Then 
    test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1) 
 Else 
    test_m2(i1) = Abs(mass2(i1)) 
 End If 
 If (Abs(mass2(i1)).lt.1.0E-30_dp) test_m2(i1) = 0._dp 
 If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = sqrt(mass2(i1)) 
  Else 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses occurred a negative mass squared!' 
     Write(*,*) 'generation: ',i1 
     Write(*,*) 'mass: ',mass2(i1) 
   SignOfMassChanged = .True. 
   mass(i1) = 0._dp 
  End If 
End Do 
 
If (Maxval(test_m2).LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopSe
 
 
Subroutine Pi1LoopSe(p2,MSe,MSe2,MAh,MAh2,MCha,MCha2,MFe,MFe2,MChi,MChi2,             & 
& MFu,MFu2,MFd,MFd2,Mhh,Mhh2,MSv,MSv2,MHpm,MHpm2,MSu,MSu2,MSd,MSd2,MVZ,MVZ2,             & 
& MVWm,MVWm2,cplAhSecUSe,cplFvChacUSeL,cplFvChacUSeR,cplChiFecUSeL,cplChiFecUSeR,        & 
& cplcFuFdcUSeL,cplcFuFdcUSeR,cplFvFecUSeL,cplFvFecUSeR,cplhhSecUSe,cplHpmSvcUSe,        & 
& cplSdcUSecSu,cplSecUSecSv,cplSecUSeVP,cplSecUSeVZ,cplSvcUSeVWm,cplAhAhUSecUSe,         & 
& cplhhhhUSecUSe,cplHpmUSecHpmcUSe,cplSdUSecSdcUSe,cplUSeSecUSecSe,cplUSeSucUSecSu,      & 
& cplUSeSvcUSecSv,cplUSecUSeVPVP,cplUSecUSecVWmVWm,cplUSecUSeVZVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MSe(6),MSe2(6),MAh(2),MAh2(2),MCha(2),MCha2(2),MFe(3),MFe2(3),MChi(4),MChi2(4),       & 
& MFu(3),MFu2(3),MFd(3),MFd2(3),Mhh(2),Mhh2(2),MSv(3),MSv2(3),MHpm(2),MHpm2(2),          & 
& MSu(6),MSu2(6),MSd(6),MSd2(6),MVZ,MVZ2,MVWm,MVWm2

Complex(dp), Intent(in) :: cplAhSecUSe(2,6,6),cplFvChacUSeL(3,2,6),cplFvChacUSeR(3,2,6),cplChiFecUSeL(4,3,6),    & 
& cplChiFecUSeR(4,3,6),cplcFuFdcUSeL(3,3,6),cplcFuFdcUSeR(3,3,6),cplFvFecUSeL(3,3,6),    & 
& cplFvFecUSeR(3,3,6),cplhhSecUSe(2,6,6),cplHpmSvcUSe(2,3,6),cplSdcUSecSu(6,6,6),        & 
& cplSecUSecSv(6,6,3),cplSecUSeVP(6,6),cplSecUSeVZ(6,6),cplSvcUSeVWm(3,6),               & 
& cplAhAhUSecUSe(2,2,6,6),cplhhhhUSecUSe(2,2,6,6),cplHpmUSecHpmcUSe(2,6,2,6),            & 
& cplSdUSecSdcUSe(6,6,6,6),cplUSeSecUSecSe(6,6,6,6),cplUSeSucUSecSu(6,6,6,6),            & 
& cplUSeSvcUSecSv(6,3,6,3),cplUSecUSeVPVP(6,6),cplUSecUSecVWmVWm(6,6),cplUSecUSeVZVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(6,6) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumI(6,6) 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! Se, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSe2(i1),MAh2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplAhSecUSe(i2,i1,gO1)
coup2 = Conjg(cplAhSecUSe(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 2
 G0m2 = Real(SA_Gloop(p2,0._dp,MCha2(i2)),dp) 
B0m2 = -2._dp*0.*MCha(i2)*Real(SA_B0(p2,0._dp,MCha2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplFvChacUSeL(i1,i2,gO1)
coupR1 = cplFvChacUSeR(i1,i2,gO1)
coupL2 =  Conjg(cplFvChacUSeL(i1,i2,gO2))
coupR2 =  Conjg(cplFvChacUSeR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fe, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 4
 G0m2 = Real(SA_Gloop(p2,MFe2(i1),MChi2(i2)),dp) 
B0m2 = -2._dp*MFe(i1)*MChi(i2)*Real(SA_B0(p2,MFe2(i1),MChi2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplChiFecUSeL(i2,i1,gO1)
coupR1 = cplChiFecUSeR(i2,i1,gO1)
coupL2 =  Conjg(cplChiFecUSeL(i2,i1,gO2))
coupR2 =  Conjg(cplChiFecUSeR(i2,i1,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFu2(i1),MFd2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*MFd(i2)*Real(SA_B0(p2,MFu2(i1),MFd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplcFuFdcUSeL(i1,i2,gO1)
coupR1 = cplcFuFdcUSeR(i1,i2,gO1)
coupL2 =  Conjg(cplcFuFdcUSeL(i1,i2,gO2))
coupR2 =  Conjg(cplcFuFdcUSeR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,0._dp,MFe2(i2)),dp) 
B0m2 = -2._dp*0.*MFe(i2)*Real(SA_B0(p2,0._dp,MFe2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coupL1 = cplFvFecUSeL(i1,i2,gO1)
coupR1 = cplFvFecUSeR(i1,i2,gO1)
coupL2 =  Conjg(cplFvFecUSeL(i1,i2,gO2))
coupR2 =  Conjg(cplFvFecUSeR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Se, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSe2(i1),Mhh2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplhhSecUSe(i2,i1,gO1)
coup2 = Conjg(cplhhSecUSe(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Sv, Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MSv2(i1),MHpm2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplHpmSvcUSe(i2,i1,gO1)
coup2 = Conjg(cplHpmSvcUSe(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSu2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdcUSecSu(i2,gO1,i1)
coup2 = Conjg(cplSdcUSecSu(i2,gO2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSv2(i1),MSe2(i2)),dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSecUSecSv(i2,gO1,i1)
coup2 = Conjg(cplSecUSecSv(i2,gO2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VP, Se 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSe2(i2),0._dp) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSecUSeVP(i2,gO1)
coup2 =  Conjg(cplSecUSeVP(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VZ, Se 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 6
 F0m2 = FloopRXi(p2,MSe2(i2),MVZ2) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSecUSeVZ(i2,gO1)
coup2 =  Conjg(cplSecUSeVZ(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VWm, Sv 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 3
 F0m2 = FloopRXi(p2,MSv2(i2),MVWm2) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSvcUSeVWm(i2,gO1)
coup2 =  Conjg(cplSvcUSeVWm(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MAh2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplAhAhUSecUSe(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(Mhh2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplhhhhUSecUSe(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MHpm2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplHpmUSecHpmcUSe(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSd2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplSdUSecSdcUSe(i1,gO2,i1,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSe2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSeSecUSecSe(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSu2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSeSucUSecSu(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
 A0m2 = SA_A0(MSv2(i1)) 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSeSvcUSecSv(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! VP, VP 
!------------------------ 
sumI = 0._dp 
 
!------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVWm2) + 0.25_dp*RXi*SA_A0(MVWm2*RXi) - 0.5_dp*MVWm2*rMS 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSecUSecVWmVWm(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +4._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVZ2) + 0.25_dp*RXi*SA_A0(MVZ2*RXi) - 0.5_dp*MVZ2*rMS 
Do gO1 = 1, 6
  Do gO2 = gO1, 6
coup1 = cplUSecUSeVZVZ(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +2._dp* SumI  


Do gO2 = 1, 6
  Do gO1 = gO2+1, 6
     res(gO1,gO2) = Conjg(res(gO2,gO1))  
   End Do 
End Do 
 
res = oo16pi2*res 
 
End Subroutine Pi1LoopSe 
 
Subroutine OneLoophh(g1,g2,Mu,Bmu,mHd2,mHu2,vd,vu,MAh,MAh2,MVZ,MVZ2,MCha,             & 
& MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MVWm,MVWm2,            & 
& MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,cplAhAhUhh,cplAhUhhVZ,cplcChaChaUhhL,              & 
& cplcChaChaUhhR,cplChiChiUhhL,cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,     & 
& cplcFeFeUhhR,cplcFuFuUhhL,cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,      & 
& cplUhhhhhh,cplUhhHpmcHpm,cplUhhHpmcVWm,cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,            & 
& cplUhhSvcSv,cplUhhcVWmVWm,cplUhhVZVZ,cplAhAhUhhUhh,cplUhhUhhhhhh,cplUhhUhhHpmcHpm,     & 
& cplUhhUhhSdcSd,cplUhhUhhSecSe,cplUhhUhhSucSu,cplUhhUhhSvcSv,cplUhhUhhcVWmVWm,          & 
& cplUhhUhhVZVZ,delta,mass,mass2,RS,kont)

Implicit None 
Real(dp), Intent(in) :: MAh(2),MAh2(2),MVZ,MVZ2,MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),             & 
& MFe(3),MFe2(3),MFu(3),MFu2(3),Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MVWm,MVWm2,              & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3)

Real(dp), Intent(in) :: g1,g2,mHd2,mHu2,vd,vu

Complex(dp), Intent(in) :: Mu,Bmu

Complex(dp), Intent(in) :: cplAhAhUhh(2,2,2),cplAhUhhVZ(2,2),cplcChaChaUhhL(2,2,2),cplcChaChaUhhR(2,2,2),        & 
& cplChiChiUhhL(4,4,2),cplChiChiUhhR(4,4,2),cplcFdFdUhhL(3,3,2),cplcFdFdUhhR(3,3,2),     & 
& cplcFeFeUhhL(3,3,2),cplcFeFeUhhR(3,3,2),cplcFuFuUhhL(3,3,2),cplcFuFuUhhR(3,3,2),       & 
& cplcgWmgWmUhh(2),cplcgWpCgWpCUhh(2),cplcgZgZUhh(2),cplUhhhhhh(2,2,2),cplUhhHpmcHpm(2,2,2),& 
& cplUhhHpmcVWm(2,2),cplUhhSdcSd(2,6,6),cplUhhSecSe(2,6,6),cplUhhSucSu(2,6,6),           & 
& cplUhhSvcSv(2,3,3),cplUhhcVWmVWm(2),cplUhhVZVZ(2),cplAhAhUhhUhh(2,2,2,2),              & 
& cplUhhUhhhhhh(2,2,2,2),cplUhhUhhHpmcHpm(2,2,2,2),cplUhhUhhSdcSd(2,2,6,6),              & 
& cplUhhUhhSecSe(2,2,6,6),cplUhhUhhSucSu(2,2,6,6),cplUhhUhhSvcSv(2,2,3,3),               & 
& cplUhhUhhcVWmVWm(2,2),cplUhhUhhVZVZ(2,2)

Complex(dp) :: mat2a(2,2), mat2(2,2),  PiSf(2,2,2)
Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(2), test_m2(2),p2, test(2) 
Real(dp), Intent(out) :: mass(2), mass2(2) 
Complex(dp), Intent(out) ::  RS(2,2) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoophh'
 
mat2a(1,1) = 0._dp 
mat2a(1,1) = mat2a(1,1)+mHd2
mat2a(1,1) = mat2a(1,1)+(3*g1**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(3*g2**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g1**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g2**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+Mu*Conjg(Mu)
mat2a(1,2) = 0._dp 
mat2a(1,2) = mat2a(1,2)-(g1**2*vd*vu)/4._dp
mat2a(1,2) = mat2a(1,2)-(g2**2*vd*vu)/4._dp
mat2a(1,2) = mat2a(1,2)-Bmu/2._dp
mat2a(1,2) = mat2a(1,2)-Conjg(Bmu)/2._dp
mat2a(2,2) = 0._dp 
mat2a(2,2) = mat2a(2,2)+mHu2
mat2a(2,2) = mat2a(2,2)-(g1**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)-(g2**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(3*g1**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(3*g2**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+Mu*Conjg(Mu)

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat2a(i1,i2) = (mat2a(i2,i1)) 
  End do 
End do 

 
! Rotation matrix for p2=0 
PiSf(1,:,:) = ZeroC 
p2 = 0._dp 
Call Pi1Loophh(p2,MAh,MAh2,MVZ,MVZ2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplAhAhUhh,cplAhUhhVZ,cplcChaChaUhhL,cplcChaChaUhhR,cplChiChiUhhL,            & 
& cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,cplcFeFeUhhR,cplcFuFuUhhL,        & 
& cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,cplUhhhhhh,cplUhhHpmcHpm,       & 
& cplUhhHpmcVWm,cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,cplUhhSvcSv,cplUhhcVWmVWm,           & 
& cplUhhVZVZ,cplAhAhUhhUhh,cplUhhUhhhhhh,cplUhhUhhHpmcHpm,cplUhhUhhSdcSd,cplUhhUhhSecSe, & 
& cplUhhUhhSucSu,cplUhhUhhSvcSv,cplUhhUhhcVWmVWm,cplUhhUhhVZVZ,kont,PiSf(1,:,:))

PiSf(1,:,:) = PiSf(1,:,:) - Pi2S_EffPot 
mat2 = mat2a - Real(PiSf(1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZHOS_0 = RS 
 
 
! Now with momenta 
Do i1=1,2
PiSf(i1,:,:) = ZeroC 
p2 = Mhh2(i1)
Call Pi1Loophh(p2,MAh,MAh2,MVZ,MVZ2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplAhAhUhh,cplAhUhhVZ,cplcChaChaUhhL,cplcChaChaUhhR,cplChiChiUhhL,            & 
& cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,cplcFeFeUhhR,cplcFuFuUhhL,        & 
& cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,cplUhhhhhh,cplUhhHpmcHpm,       & 
& cplUhhHpmcVWm,cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,cplUhhSvcSv,cplUhhcVWmVWm,           & 
& cplUhhVZVZ,cplAhAhUhhUhh,cplUhhUhhhhhh,cplUhhUhhHpmcHpm,cplUhhUhhSdcSd,cplUhhUhhSecSe, & 
& cplUhhUhhSucSu,cplUhhUhhSvcSv,cplUhhUhhcVWmVWm,cplUhhUhhVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=2,1,-1 
PiSf(i1,:,:) = PiSf(i1,:,:) - Pi2S_EffPot 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
 
Do i1=1,2
  If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
  If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = Sqrt(mass2(i1)) 
  Else 
   If (ErrorLevel.Ge.0) Then 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
   Call TerminateProgram 
   End If 
   kont = -301 
   mass(i1) = 0._dp 
  End If 
End Do 
 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
Do i1=1,2
PiSf(i1,:,:) = ZeroC 
p2 =  mass2(i1) 
Call Pi1Loophh(p2,MAh,MAh2,MVZ,MVZ2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplAhAhUhh,cplAhUhhVZ,cplcChaChaUhhL,cplcChaChaUhhR,cplChiChiUhhL,            & 
& cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,cplcFeFeUhhR,cplcFuFuUhhL,        & 
& cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,cplUhhhhhh,cplUhhHpmcHpm,       & 
& cplUhhHpmcVWm,cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,cplUhhSvcSv,cplUhhcVWmVWm,           & 
& cplUhhVZVZ,cplAhAhUhhUhh,cplUhhUhhhhhh,cplUhhUhhHpmcHpm,cplUhhUhhSdcSd,cplUhhUhhSecSe, & 
& cplUhhUhhSucSu,cplUhhUhhSvcSv,cplUhhUhhcVWmVWm,cplUhhUhhVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=2,1,-1 
PiSf(i1,:,:) = PiSf(i1,:,:) - Pi2S_EffPot 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZHOS_p2(i1,:) = RS(i1,:) 
 If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
Do i1=1,2
 If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
 If (test_m2(i1).Ne.0._dp) Then 
    test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1) 
 Else 
    test_m2(i1) = Abs(mass2(i1)) 
 End If 
 If (Abs(mass2(i1)).lt.1.0E-30_dp) test_m2(i1) = 0._dp 
 If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = sqrt(mass2(i1)) 
  Else 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses occurred a negative mass squared!' 
     Write(*,*) 'generation: ',i1 
     Write(*,*) 'mass: ',mass2(i1) 
   SignOfMassChanged = .True. 
   mass(i1) = 0._dp 
  End If 
End Do 
 
If (Maxval(test_m2).LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoophh
 
 
Subroutine Pi1Loophh(p2,MAh,MAh2,MVZ,MVZ2,MCha,MCha2,MChi,MChi2,MFd,MFd2,             & 
& MFe,MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,           & 
& MSv,MSv2,cplAhAhUhh,cplAhUhhVZ,cplcChaChaUhhL,cplcChaChaUhhR,cplChiChiUhhL,            & 
& cplChiChiUhhR,cplcFdFdUhhL,cplcFdFdUhhR,cplcFeFeUhhL,cplcFeFeUhhR,cplcFuFuUhhL,        & 
& cplcFuFuUhhR,cplcgWmgWmUhh,cplcgWpCgWpCUhh,cplcgZgZUhh,cplUhhhhhh,cplUhhHpmcHpm,       & 
& cplUhhHpmcVWm,cplUhhSdcSd,cplUhhSecSe,cplUhhSucSu,cplUhhSvcSv,cplUhhcVWmVWm,           & 
& cplUhhVZVZ,cplAhAhUhhUhh,cplUhhUhhhhhh,cplUhhUhhHpmcHpm,cplUhhUhhSdcSd,cplUhhUhhSecSe, & 
& cplUhhUhhSucSu,cplUhhUhhSvcSv,cplUhhUhhcVWmVWm,cplUhhUhhVZVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MAh(2),MAh2(2),MVZ,MVZ2,MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),             & 
& MFe(3),MFe2(3),MFu(3),MFu2(3),Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MVWm,MVWm2,              & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplAhAhUhh(2,2,2),cplAhUhhVZ(2,2),cplcChaChaUhhL(2,2,2),cplcChaChaUhhR(2,2,2),        & 
& cplChiChiUhhL(4,4,2),cplChiChiUhhR(4,4,2),cplcFdFdUhhL(3,3,2),cplcFdFdUhhR(3,3,2),     & 
& cplcFeFeUhhL(3,3,2),cplcFeFeUhhR(3,3,2),cplcFuFuUhhL(3,3,2),cplcFuFuUhhR(3,3,2),       & 
& cplcgWmgWmUhh(2),cplcgWpCgWpCUhh(2),cplcgZgZUhh(2),cplUhhhhhh(2,2,2),cplUhhHpmcHpm(2,2,2),& 
& cplUhhHpmcVWm(2,2),cplUhhSdcSd(2,6,6),cplUhhSecSe(2,6,6),cplUhhSucSu(2,6,6),           & 
& cplUhhSvcSv(2,3,3),cplUhhcVWmVWm(2),cplUhhVZVZ(2),cplAhAhUhhUhh(2,2,2,2),              & 
& cplUhhUhhhhhh(2,2,2,2),cplUhhUhhHpmcHpm(2,2,2,2),cplUhhUhhSdcSd(2,2,6,6),              & 
& cplUhhUhhSecSe(2,2,6,6),cplUhhUhhSucSu(2,2,6,6),cplUhhUhhSvcSv(2,2,3,3),               & 
& cplUhhUhhcVWmVWm(2,2),cplUhhUhhVZVZ(2,2)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2,2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumI(2,2) 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MAh2(i1),MAh2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplAhAhUhh(i1,i2,gO1)
coup2 = Conjg(cplAhAhUhh(i1,i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VZ, Ah 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,MAh2(i2),MVZ2) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplAhUhhVZ(i2,gO1)
coup2 =  Conjg(cplAhUhhVZ(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 G0m2 = Real(SA_Gloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = -2._dp*MCha(i1)*MCha(i2)*Real(SA_B0(p2,MCha2(i1),MCha2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcChaChaUhhL(i1,i2,gO1)
coupR1 = cplcChaChaUhhR(i1,i2,gO1)
coupL2 =  Conjg(cplcChaChaUhhL(i1,i2,gO2))
coupR2 =  Conjg(cplcChaChaUhhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 G0m2 = Real(SA_Gloop(p2,MChi2(i1),MChi2(i2)),dp) 
B0m2 = -2._dp*MChi(i1)*MChi(i2)*Real(SA_B0(p2,MChi2(i1),MChi2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplChiChiUhhL(i1,i2,gO1)
coupR1 = cplChiChiUhhR(i1,i2,gO1)
coupL2 =  Conjg(cplChiChiUhhL(i1,i2,gO2))
coupR2 =  Conjg(cplChiChiUhhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = -2._dp*MFd(i1)*MFd(i2)*Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcFdFdUhhL(i1,i2,gO1)
coupR1 = cplcFdFdUhhR(i1,i2,gO1)
coupL2 =  Conjg(cplcFdFdUhhL(i1,i2,gO2))
coupR2 =  Conjg(cplcFdFdUhhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = -2._dp*MFe(i1)*MFe(i2)*Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcFeFeUhhL(i1,i2,gO1)
coupR1 = cplcFeFeUhhR(i1,i2,gO1)
coupL2 =  Conjg(cplcFeFeUhhL(i1,i2,gO2))
coupR2 =  Conjg(cplcFeFeUhhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*MFu(i2)*Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcFuFuUhhL(i1,i2,gO1)
coupR1 = cplcFuFuUhhR(i1,i2,gO1)
coupL2 =  Conjg(cplcFuFuUhhL(i1,i2,gO2))
coupR2 =  Conjg(cplcFuFuUhhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
F0m2 = -Real(SA_B0(p2,MVWm2*RXi,MVWm2*RXi),dp) 
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcgWmgWmUhh(gO1)
coup2 =  cplcgWmgWmUhh(gO2) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
F0m2 = -Real(SA_B0(p2,MVWm2*RXi,MVWm2*RXi),dp) 
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcgWpCgWpCUhh(gO1)
coup2 =  cplcgWpCgWpCUhh(gO2) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gZ], gZ 
!------------------------ 
sumI = 0._dp 
 
F0m2 = -Real(SA_B0(p2,MVZ2*RXi,MVZ2*RXi),dp) 
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcgZgZUhh(gO1)
coup2 =  cplcgZgZUhh(gO2) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,Mhh2(i1),Mhh2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhhhhh(gO1,i1,i2)
coup2 = Conjg(cplUhhhhhh(gO2,i1,i2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MHpm2(i1),MHpm2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhHpmcHpm(gO1,i2,i1)
coup2 = Conjg(cplUhhHpmcHpm(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,MHpm2(i2),MVWm2) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhHpmcVWm(gO1,i2)
coup2 =  Conjg(cplUhhHpmcVWm(gO2,i2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +2._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSd2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhSdcSd(gO1,i2,i1)
coup2 = Conjg(cplUhhSdcSd(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSe2(i1),MSe2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhSecSe(gO1,i2,i1)
coup2 = Conjg(cplUhhSecSe(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSu2(i1),MSu2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhSucSu(gO1,i2,i1)
coup2 = Conjg(cplUhhSucSu(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 B0m2 = Real(SA_B0(p2,MSv2(i1),MSv2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhSvcSv(gO1,i2,i1)
coup2 = Conjg(cplUhhSvcSv(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
F0m2 = Real(SVVloop(p2,MVWm2,MVWm2),dp)   
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhcVWmVWm(gO1)
coup2 =  Conjg(cplUhhcVWmVWm(gO2)) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
F0m2 = Real(SVVloop(p2,MVZ2,MVZ2),dp)   
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhVZVZ(gO1)
coup2 =  Conjg(cplUhhVZVZ(gO2)) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
!------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MAh2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplAhAhUhhUhh(i1,i1,gO1,gO2)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(Mhh2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhhhhh(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MHpm2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhHpmcHpm(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSd2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhSdcSd(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSe2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhSecSe(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSu2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhSucSu(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
 A0m2 = SA_A0(MSv2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhSvcSv(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVWm2) + 0.25_dp*RXi*SA_A0(MVWm2*RXi) - 0.5_dp*MVWm2*rMS 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhcVWmVWm(gO1,gO2)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +4._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVZ2) + 0.25_dp*RXi*SA_A0(MVZ2*RXi) - 0.5_dp*MVZ2*rMS 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUhhUhhVZVZ(gO1,gO2)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +2._dp* SumI  


Do gO2 = 1, 2
  Do gO1 = gO2+1, 2
     res(gO1,gO2) = (res(gO2,gO1))  
   End Do 
End Do 
 
res = oo16pi2*res 
 
End Subroutine Pi1Loophh 
 
Subroutine OneLoopAh(g1,g2,Mu,Bmu,mHd2,mHu2,vd,vu,TW,Mhh,Mhh2,MAh,MAh2,               & 
& MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,             & 
& MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,    & 
& cplChiChiUAhL,cplChiChiUAhR,cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,       & 
& cplcFuFuUAhL,cplcFuFuUAhR,cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,      & 
& cplUAhHpmcVWm,cplUAhSdcSd,cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,         & 
& cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,          & 
& cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,delta,mass,mass2,RS,kont)

Implicit None 
Real(dp), Intent(in) :: Mhh(2),Mhh2(2),MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),               & 
& MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MVZ,MVZ2,MHpm(2),MHpm2(2),MVWm,MVWm2,            & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3)

Real(dp), Intent(in) :: g1,g2,mHd2,mHu2,vd,vu,TW

Complex(dp), Intent(in) :: Mu,Bmu

Complex(dp), Intent(in) :: cplUAhAhhh(2,2,2),cplcChaChaUAhL(2,2,2),cplcChaChaUAhR(2,2,2),cplChiChiUAhL(4,4,2),   & 
& cplChiChiUAhR(4,4,2),cplcFdFdUAhL(3,3,2),cplcFdFdUAhR(3,3,2),cplcFeFeUAhL(3,3,2),      & 
& cplcFeFeUAhR(3,3,2),cplcFuFuUAhL(3,3,2),cplcFuFuUAhR(3,3,2),cplcgWmgWmUAh(2),          & 
& cplcgWpCgWpCUAh(2),cplUAhhhVZ(2,2),cplUAhHpmcHpm(2,2,2),cplUAhHpmcVWm(2,2),            & 
& cplUAhSdcSd(2,6,6),cplUAhSecSe(2,6,6),cplUAhSucSu(2,6,6),cplUAhUAhAhAh(2,2,2,2),       & 
& cplUAhUAhhhhh(2,2,2,2),cplUAhUAhHpmcHpm(2,2,2,2),cplUAhUAhSdcSd(2,2,6,6),              & 
& cplUAhUAhSecSe(2,2,6,6),cplUAhUAhSucSu(2,2,6,6),cplUAhUAhSvcSv(2,2,3,3),               & 
& cplUAhUAhcVWmVWm(2,2),cplUAhUAhVZVZ(2,2)

Complex(dp) :: mat2a(2,2), mat2(2,2),  PiSf(2,2,2)
Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(2), test_m2(2),p2, test(2) 
Real(dp), Intent(out) :: mass(2), mass2(2) 
Complex(dp), Intent(out) ::  RS(2,2) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopAh'
 
mat2a(1,1) = 0._dp 
mat2a(1,1) = mat2a(1,1)+mHd2
mat2a(1,1) = mat2a(1,1)+(g1**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(g2**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g1**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g2**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+Mu*Conjg(Mu)
mat2a(1,1) = mat2a(1,1)+(g2**2*vd**2*Cos(TW)**2*RXiZ)/4._dp
mat2a(1,1) = mat2a(1,1)+(g1*g2*vd**2*Cos(TW)*RXiZ*Sin(TW))/2._dp
mat2a(1,1) = mat2a(1,1)+(g1**2*vd**2*RXiZ*Sin(TW)**2)/4._dp
mat2a(1,2) = 0._dp 
mat2a(1,2) = mat2a(1,2)+Bmu/2._dp
mat2a(1,2) = mat2a(1,2)+Conjg(Bmu)/2._dp
mat2a(1,2) = mat2a(1,2)-(g2**2*vd*vu*Cos(TW)**2*RXiZ)/4._dp
mat2a(1,2) = mat2a(1,2)-(g1*g2*vd*vu*Cos(TW)*RXiZ*Sin(TW))/2._dp
mat2a(1,2) = mat2a(1,2)-(g1**2*vd*vu*RXiZ*Sin(TW)**2)/4._dp
mat2a(2,2) = 0._dp 
mat2a(2,2) = mat2a(2,2)+mHu2
mat2a(2,2) = mat2a(2,2)-(g1**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)-(g2**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g1**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g2**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+Mu*Conjg(Mu)
mat2a(2,2) = mat2a(2,2)+(g2**2*vu**2*Cos(TW)**2*RXiZ)/4._dp
mat2a(2,2) = mat2a(2,2)+(g1*g2*vu**2*Cos(TW)*RXiZ*Sin(TW))/2._dp
mat2a(2,2) = mat2a(2,2)+(g1**2*vu**2*RXiZ*Sin(TW)**2)/4._dp

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat2a(i1,i2) = (mat2a(i2,i1)) 
  End do 
End do 

 
! Rotation matrix for p2=0 
PiSf(1,:,:) = ZeroC 
p2 = 0._dp 
If (i1.eq.1) p2 = 0._dp 
Call Pi1LoopAh(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,cplChiChiUAhL,cplChiChiUAhR,         & 
& cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,cplcFuFuUAhL,cplcFuFuUAhR,         & 
& cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,cplUAhHpmcVWm,cplUAhSdcSd,      & 
& cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,   & 
& cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,           & 
& kont,PiSf(1,:,:))

PiSf(1,:,:) = PiSf(1,:,:) - PiP2S_EffPot 
mat2 = mat2a - Real(PiSf(1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZAOS_0 = RS 
 
 
! Now with momenta 
Do i1=1,2
PiSf(i1,:,:) = ZeroC 
p2 = MAh2(i1)
If (i1.eq.1) p2 = 0._dp 
Call Pi1LoopAh(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,cplChiChiUAhL,cplChiChiUAhR,         & 
& cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,cplcFuFuUAhL,cplcFuFuUAhR,         & 
& cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,cplUAhHpmcVWm,cplUAhSdcSd,      & 
& cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,   & 
& cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,           & 
& kont,PiSf(i1,:,:))

End Do 
Do i1=2,1,-1 
PiSf(i1,:,:) = PiSf(i1,:,:) - PiP2S_EffPot 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
 
Do i1=1,2
  If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
  If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = Sqrt(mass2(i1)) 
  Else 
   If (ErrorLevel.Ge.0) Then 
   If ((i1.Gt.1).or.(Abs(mass2(i1)).gt.MaxVal(Abs(mass2)))) Then 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
   Call TerminateProgram 
   End If 
   End If 
   kont = -301 
   mass(i1) = 0._dp 
  End If 
End Do 
 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
Do i1=1,2
PiSf(i1,:,:) = ZeroC 
p2 =  mass2(i1) 
If (i1.eq.1) p2 = 0._dp 
Call Pi1LoopAh(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,cplChiChiUAhL,cplChiChiUAhR,         & 
& cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,cplcFuFuUAhL,cplcFuFuUAhR,         & 
& cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,cplUAhHpmcVWm,cplUAhSdcSd,      & 
& cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,   & 
& cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,           & 
& kont,PiSf(i1,:,:))

End Do 
Do i1=2,1,-1 
PiSf(i1,:,:) = PiSf(i1,:,:) - PiP2S_EffPot 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZAOS_p2(i1,:) = RS(i1,:) 
 If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
Do i1=1,2
 If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
 If (test_m2(i1).Ne.0._dp) Then 
    test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1) 
 Else 
    test_m2(i1) = Abs(mass2(i1)) 
 End If 
 If (Abs(mass2(i1)).lt.1.0E-30_dp) test_m2(i1) = 0._dp 
 If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = sqrt(mass2(i1)) 
  Else 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses occurred a negative mass squared!' 
     Write(*,*) 'generation: ',i1 
     Write(*,*) 'mass: ',mass2(i1) 
   SignOfMassChanged = .True. 
   mass(i1) = 0._dp 
  End If 
End Do 
 
If (Maxval(test_m2).LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopAh
 
 
Subroutine Pi1LoopAh(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,             & 
& MFe,MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,           & 
& MSv,MSv2,cplUAhAhhh,cplcChaChaUAhL,cplcChaChaUAhR,cplChiChiUAhL,cplChiChiUAhR,         & 
& cplcFdFdUAhL,cplcFdFdUAhR,cplcFeFeUAhL,cplcFeFeUAhR,cplcFuFuUAhL,cplcFuFuUAhR,         & 
& cplcgWmgWmUAh,cplcgWpCgWpCUAh,cplUAhhhVZ,cplUAhHpmcHpm,cplUAhHpmcVWm,cplUAhSdcSd,      & 
& cplUAhSecSe,cplUAhSucSu,cplUAhUAhAhAh,cplUAhUAhhhhh,cplUAhUAhHpmcHpm,cplUAhUAhSdcSd,   & 
& cplUAhUAhSecSe,cplUAhUAhSucSu,cplUAhUAhSvcSv,cplUAhUAhcVWmVWm,cplUAhUAhVZVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: Mhh(2),Mhh2(2),MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),               & 
& MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MVZ,MVZ2,MHpm(2),MHpm2(2),MVWm,MVWm2,            & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplUAhAhhh(2,2,2),cplcChaChaUAhL(2,2,2),cplcChaChaUAhR(2,2,2),cplChiChiUAhL(4,4,2),   & 
& cplChiChiUAhR(4,4,2),cplcFdFdUAhL(3,3,2),cplcFdFdUAhR(3,3,2),cplcFeFeUAhL(3,3,2),      & 
& cplcFeFeUAhR(3,3,2),cplcFuFuUAhL(3,3,2),cplcFuFuUAhR(3,3,2),cplcgWmgWmUAh(2),          & 
& cplcgWpCgWpCUAh(2),cplUAhhhVZ(2,2),cplUAhHpmcHpm(2,2,2),cplUAhHpmcVWm(2,2),            & 
& cplUAhSdcSd(2,6,6),cplUAhSecSe(2,6,6),cplUAhSucSu(2,6,6),cplUAhUAhAhAh(2,2,2,2),       & 
& cplUAhUAhhhhh(2,2,2,2),cplUAhUAhHpmcHpm(2,2,2,2),cplUAhUAhSdcSd(2,2,6,6),              & 
& cplUAhUAhSecSe(2,2,6,6),cplUAhUAhSucSu(2,2,6,6),cplUAhUAhSvcSv(2,2,3,3),               & 
& cplUAhUAhcVWmVWm(2,2),cplUAhUAhVZVZ(2,2)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2,2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumI(2,2) 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! hh, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,Mhh2(i1),MAh2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhAhhh(gO1,i2,i1)
coup2 = Conjg(cplUAhAhhh(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 G0m2 = Real(SA_Gloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = -2._dp*MCha(i1)*MCha(i2)*Real(SA_B0(p2,MCha2(i1),MCha2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcChaChaUAhL(i1,i2,gO1)
coupR1 = cplcChaChaUAhR(i1,i2,gO1)
coupL2 =  Conjg(cplcChaChaUAhL(i1,i2,gO2))
coupR2 =  Conjg(cplcChaChaUAhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 G0m2 = Real(SA_Gloop(p2,MChi2(i1),MChi2(i2)),dp) 
B0m2 = -2._dp*MChi(i1)*MChi(i2)*Real(SA_B0(p2,MChi2(i1),MChi2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplChiChiUAhL(i1,i2,gO1)
coupR1 = cplChiChiUAhR(i1,i2,gO1)
coupL2 =  Conjg(cplChiChiUAhL(i1,i2,gO2))
coupR2 =  Conjg(cplChiChiUAhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = -2._dp*MFd(i1)*MFd(i2)*Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcFdFdUAhL(i1,i2,gO1)
coupR1 = cplcFdFdUAhR(i1,i2,gO1)
coupL2 =  Conjg(cplcFdFdUAhL(i1,i2,gO2))
coupR2 =  Conjg(cplcFdFdUAhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = -2._dp*MFe(i1)*MFe(i2)*Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcFeFeUAhL(i1,i2,gO1)
coupR1 = cplcFeFeUAhR(i1,i2,gO1)
coupL2 =  Conjg(cplcFeFeUAhL(i1,i2,gO2))
coupR2 =  Conjg(cplcFeFeUAhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*MFu(i2)*Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcFuFuUAhL(i1,i2,gO1)
coupR1 = cplcFuFuUAhR(i1,i2,gO1)
coupL2 =  Conjg(cplcFuFuUAhL(i1,i2,gO2))
coupR2 =  Conjg(cplcFuFuUAhR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
F0m2 = -Real(SA_B0(p2,MVWm2*RXi,MVWm2*RXi),dp) 
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcgWmgWmUAh(gO1)
coup2 =  cplcgWmgWmUAh(gO2) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
F0m2 = -Real(SA_B0(p2,MVWm2*RXi,MVWm2*RXi),dp) 
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcgWpCgWpCUAh(gO1)
coup2 =  cplcgWpCgWpCUAh(gO2) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! VZ, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,Mhh2(i2),MVZ2) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhhhVZ(gO1,i2)
coup2 =  Conjg(cplUAhhhVZ(gO2,i2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MHpm2(i1),MHpm2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhHpmcHpm(gO1,i2,i1)
coup2 = Conjg(cplUAhHpmcHpm(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,MHpm2(i2),MVWm2) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhHpmcVWm(gO1,i2)
coup2 =  Conjg(cplUAhHpmcVWm(gO2,i2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +2._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSd2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhSdcSd(gO1,i2,i1)
coup2 = Conjg(cplUAhSdcSd(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSe2(i1),MSe2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhSecSe(gO1,i2,i1)
coup2 = Conjg(cplUAhSecSe(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSu2(i1),MSu2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhSucSu(gO1,i2,i1)
coup2 = Conjg(cplUAhSucSu(gO2,i2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MAh2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhAhAh(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(Mhh2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhhhhh(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MHpm2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhHpmcHpm(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSd2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhSdcSd(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSe2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhSecSe(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSu2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhSucSu(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
 A0m2 = SA_A0(MSv2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhSvcSv(gO1,gO2,i1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVWm2) + 0.25_dp*RXi*SA_A0(MVWm2*RXi) - 0.5_dp*MVWm2*rMS 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhcVWmVWm(gO1,gO2)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +4._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVZ2) + 0.25_dp*RXi*SA_A0(MVZ2*RXi) - 0.5_dp*MVZ2*rMS 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUAhUAhVZVZ(gO1,gO2)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +2._dp* SumI  


Do gO2 = 1, 2
  Do gO1 = gO2+1, 2
     res(gO1,gO2) = (res(gO2,gO1))  
   End Do 
End Do 
 
res = oo16pi2*res 
 
End Subroutine Pi1LoopAh 
 
Subroutine OneLoopHpm(g1,g2,Mu,Bmu,mHd2,mHu2,vd,vu,MHpm,MHpm2,MAh,MAh2,               & 
& MVWm,MVWm2,MChi,MChi2,MCha,MCha2,MFu,MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVZ,              & 
& MVZ2,MSu,MSu2,MSd,MSd2,MSv,MSv2,MSe,MSe2,cplAhHpmcUHpm,cplAhcUHpmVWm,cplChiChacUHpmL,  & 
& cplChiChacUHpmR,cplcFuFdcUHpmL,cplcFuFdcUHpmR,cplFvFecUHpmL,cplFvFecUHpmR,             & 
& cplcgZgWmcUHpm,cplcgWmgZUHpm,cplcgWpCgZcUHpm,cplcgZgWpCUHpm,cplhhHpmcUHpm,             & 
& cplhhcUHpmVWm,cplHpmcUHpmVP,cplHpmcUHpmVZ,cplSdcUHpmcSu,cplSecUHpmcSv,cplcUHpmVPVWm,   & 
& cplcUHpmVWmVZ,cplAhAhUHpmcUHpm,cplhhhhUHpmcUHpm,cplUHpmHpmcUHpmcHpm,cplUHpmSdcUHpmcSd, & 
& cplUHpmSecUHpmcSe,cplUHpmSucUHpmcSu,cplUHpmSvcUHpmcSv,cplUHpmcUHpmVPVP,cplUHpmcUHpmcVWmVWm,& 
& cplUHpmcUHpmVZVZ,delta,mass,mass2,RS,kont)

Implicit None 
Real(dp), Intent(in) :: MHpm(2),MHpm2(2),MAh(2),MAh2(2),MVWm,MVWm2,MChi(4),MChi2(4),MCha(2),MCha2(2),         & 
& MFu(3),MFu2(3),MFd(3),MFd2(3),MFe(3),MFe2(3),Mhh(2),Mhh2(2),MVZ,MVZ2,MSu(6),           & 
& MSu2(6),MSd(6),MSd2(6),MSv(3),MSv2(3),MSe(6),MSe2(6)

Real(dp), Intent(in) :: g1,g2,mHd2,mHu2,vd,vu

Complex(dp), Intent(in) :: Mu,Bmu

Complex(dp), Intent(in) :: cplAhHpmcUHpm(2,2,2),cplAhcUHpmVWm(2,2),cplChiChacUHpmL(4,2,2),cplChiChacUHpmR(4,2,2),& 
& cplcFuFdcUHpmL(3,3,2),cplcFuFdcUHpmR(3,3,2),cplFvFecUHpmL(3,3,2),cplFvFecUHpmR(3,3,2), & 
& cplcgZgWmcUHpm(2),cplcgWmgZUHpm(2),cplcgWpCgZcUHpm(2),cplcgZgWpCUHpm(2),               & 
& cplhhHpmcUHpm(2,2,2),cplhhcUHpmVWm(2,2),cplHpmcUHpmVP(2,2),cplHpmcUHpmVZ(2,2),         & 
& cplSdcUHpmcSu(6,2,6),cplSecUHpmcSv(6,2,3),cplcUHpmVPVWm(2),cplcUHpmVWmVZ(2),           & 
& cplAhAhUHpmcUHpm(2,2,2,2),cplhhhhUHpmcUHpm(2,2,2,2),cplUHpmHpmcUHpmcHpm(2,2,2,2),      & 
& cplUHpmSdcUHpmcSd(2,6,2,6),cplUHpmSecUHpmcSe(2,6,2,6),cplUHpmSucUHpmcSu(2,6,2,6),      & 
& cplUHpmSvcUHpmcSv(2,3,2,3),cplUHpmcUHpmVPVP(2,2),cplUHpmcUHpmcVWmVWm(2,2),             & 
& cplUHpmcUHpmVZVZ(2,2)

Complex(dp) :: mat2a(2,2), mat2(2,2),  PiSf(2,2,2)
Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(2), test_m2(2),p2, test(2) 
Real(dp), Intent(out) :: mass(2), mass2(2) 
Complex(dp), Intent(out) ::  RS(2,2) 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopHpm'
 
mat2a(1,1) = 0._dp 
mat2a(1,1) = mat2a(1,1)+mHd2
mat2a(1,1) = mat2a(1,1)+(g1**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(g2**2*vd**2)/8._dp
mat2a(1,1) = mat2a(1,1)-(g1**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+(g2**2*vu**2)/8._dp
mat2a(1,1) = mat2a(1,1)+Mu*Conjg(Mu)
mat2a(1,1) = mat2a(1,1)+(g2**2*vd**2*RXiWm)/4._dp
mat2a(1,2) = 0._dp 
mat2a(1,2) = mat2a(1,2)+(g2**2*vd*vu)/4._dp
mat2a(1,2) = mat2a(1,2)+Conjg(Bmu)
mat2a(1,2) = mat2a(1,2)-(g2**2*vd*vu*RXiWm)/4._dp
mat2a(2,2) = 0._dp 
mat2a(2,2) = mat2a(2,2)+mHu2
mat2a(2,2) = mat2a(2,2)-(g1**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g2**2*vd**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g1**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+(g2**2*vu**2)/8._dp
mat2a(2,2) = mat2a(2,2)+Mu*Conjg(Mu)
mat2a(2,2) = mat2a(2,2)+(g2**2*vu**2*RXiWm)/4._dp

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat2a(i1,i2) = Conjg(mat2a(i2,i1)) 
  End do 
End do 

 
! Rotation matrix for p2=0 
PiSf(1,:,:) = ZeroC 
p2 = 0._dp 
If (i1.eq.1) p2 = 0._dp 
Call Pi1LoopHpm(p2,MHpm,MHpm2,MAh,MAh2,MVWm,MVWm2,MChi,MChi2,MCha,MCha2,              & 
& MFu,MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,MSv2,               & 
& MSe,MSe2,cplAhHpmcUHpm,cplAhcUHpmVWm,cplChiChacUHpmL,cplChiChacUHpmR,cplcFuFdcUHpmL,   & 
& cplcFuFdcUHpmR,cplFvFecUHpmL,cplFvFecUHpmR,cplcgZgWmcUHpm,cplcgWmgZUHpm,               & 
& cplcgWpCgZcUHpm,cplcgZgWpCUHpm,cplhhHpmcUHpm,cplhhcUHpmVWm,cplHpmcUHpmVP,              & 
& cplHpmcUHpmVZ,cplSdcUHpmcSu,cplSecUHpmcSv,cplcUHpmVPVWm,cplcUHpmVWmVZ,cplAhAhUHpmcUHpm,& 
& cplhhhhUHpmcUHpm,cplUHpmHpmcUHpmcHpm,cplUHpmSdcUHpmcSd,cplUHpmSecUHpmcSe,              & 
& cplUHpmSucUHpmcSu,cplUHpmSvcUHpmcSv,cplUHpmcUHpmVPVP,cplUHpmcUHpmcVWmVWm,              & 
& cplUHpmcUHpmVZVZ,kont,PiSf(1,:,:))

mat2 = mat2a - Real(PiSf(1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZPOS_0 = RS 
 
 
! Now with momenta 
Do i1=1,2
PiSf(i1,:,:) = ZeroC 
p2 = MHpm2(i1)
If (i1.eq.1) p2 = 0._dp 
Call Pi1LoopHpm(p2,MHpm,MHpm2,MAh,MAh2,MVWm,MVWm2,MChi,MChi2,MCha,MCha2,              & 
& MFu,MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,MSv2,               & 
& MSe,MSe2,cplAhHpmcUHpm,cplAhcUHpmVWm,cplChiChacUHpmL,cplChiChacUHpmR,cplcFuFdcUHpmL,   & 
& cplcFuFdcUHpmR,cplFvFecUHpmL,cplFvFecUHpmR,cplcgZgWmcUHpm,cplcgWmgZUHpm,               & 
& cplcgWpCgZcUHpm,cplcgZgWpCUHpm,cplhhHpmcUHpm,cplhhcUHpmVWm,cplHpmcUHpmVP,              & 
& cplHpmcUHpmVZ,cplSdcUHpmcSu,cplSecUHpmcSv,cplcUHpmVPVWm,cplcUHpmVWmVZ,cplAhAhUHpmcUHpm,& 
& cplhhhhUHpmcUHpm,cplUHpmHpmcUHpmcHpm,cplUHpmSdcUHpmcSd,cplUHpmSecUHpmcSe,              & 
& cplUHpmSucUHpmcSu,cplUHpmSvcUHpmcSv,cplUHpmcUHpmVPVP,cplUHpmcUHpmcVWmVWm,              & 
& cplUHpmcUHpmVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=2,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
 
Do i1=1,2
  If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
  If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = Sqrt(mass2(i1)) 
  Else 
   If (ErrorLevel.Ge.0) Then 
   If ((i1.Gt.1).or.(Abs(mass2(i1)).gt.MaxVal(Abs(mass2)))) Then 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
   Call TerminateProgram 
   End If 
   End If 
   kont = -301 
   mass(i1) = 0._dp 
  End If 
End Do 
 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
Do i1=1,2
PiSf(i1,:,:) = ZeroC 
p2 =  mass2(i1) 
If (i1.eq.1) p2 = 0._dp 
Call Pi1LoopHpm(p2,MHpm,MHpm2,MAh,MAh2,MVWm,MVWm2,MChi,MChi2,MCha,MCha2,              & 
& MFu,MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,MSv2,               & 
& MSe,MSe2,cplAhHpmcUHpm,cplAhcUHpmVWm,cplChiChacUHpmL,cplChiChacUHpmR,cplcFuFdcUHpmL,   & 
& cplcFuFdcUHpmR,cplFvFecUHpmL,cplFvFecUHpmR,cplcgZgWmcUHpm,cplcgWmgZUHpm,               & 
& cplcgWpCgZcUHpm,cplcgZgWpCUHpm,cplhhHpmcUHpm,cplhhcUHpmVWm,cplHpmcUHpmVP,              & 
& cplHpmcUHpmVZ,cplSdcUHpmcSu,cplSecUHpmcSv,cplcUHpmVPVWm,cplcUHpmVWmVZ,cplAhAhUHpmcUHpm,& 
& cplhhhhUHpmcUHpm,cplUHpmHpmcUHpmcHpm,cplUHpmSdcUHpmcSd,cplUHpmSecUHpmcSe,              & 
& cplUHpmSucUHpmcSu,cplUHpmSvcUHpmcSv,cplUHpmcUHpmVPVP,cplUHpmcUHpmcVWmVWm,              & 
& cplUHpmcUHpmVZVZ,kont,PiSf(i1,:,:))

End Do 
Do i1=2,1,-1 
mat2 = mat2a - Real(PiSf(i1,:,:),dp) 
Call Chop(mat2) 
Call Eigensystem(mat2,mi2,RS,kont,test) 
ZPOS_p2(i1,:) = RS(i1,:) 
 If ((kont.Eq.-8).Or.(kont.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  kont = 0 
End If 
If ((kont.Ne.0).And.(ErrorLevel.Ge.0)) Then 
  Write(ErrCan,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  Write(*,*) "Diagonalization did not work in routine "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
End If 
mass2(i1) = mi2(i1) 
End do 
Do i1=1,2
 If (Abs(mass2(i1)).Le.MaxMassNumericalZero**2) mass2(i1) = 0._dp 
 If (test_m2(i1).Ne.0._dp) Then 
    test_m2(i1) = Abs(test_m2(i1) - mass2(i1)) / test_m2(i1) 
 Else 
    test_m2(i1) = Abs(mass2(i1)) 
 End If 
 If (Abs(mass2(i1)).lt.1.0E-30_dp) test_m2(i1) = 0._dp 
 If (mass2(i1).Ge.0._dp) Then 
    mass(i1) = sqrt(mass2(i1)) 
  Else 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses occurred a negative mass squared!' 
     Write(*,*) 'generation: ',i1 
     Write(*,*) 'mass: ',mass2(i1) 
   SignOfMassChanged = .True. 
   mass(i1) = 0._dp 
  End If 
End Do 
 
If (Maxval(test_m2).LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopHpm
 
 
Subroutine Pi1LoopHpm(p2,MHpm,MHpm2,MAh,MAh2,MVWm,MVWm2,MChi,MChi2,MCha,              & 
& MCha2,MFu,MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,              & 
& MSv2,MSe,MSe2,cplAhHpmcUHpm,cplAhcUHpmVWm,cplChiChacUHpmL,cplChiChacUHpmR,             & 
& cplcFuFdcUHpmL,cplcFuFdcUHpmR,cplFvFecUHpmL,cplFvFecUHpmR,cplcgZgWmcUHpm,              & 
& cplcgWmgZUHpm,cplcgWpCgZcUHpm,cplcgZgWpCUHpm,cplhhHpmcUHpm,cplhhcUHpmVWm,              & 
& cplHpmcUHpmVP,cplHpmcUHpmVZ,cplSdcUHpmcSu,cplSecUHpmcSv,cplcUHpmVPVWm,cplcUHpmVWmVZ,   & 
& cplAhAhUHpmcUHpm,cplhhhhUHpmcUHpm,cplUHpmHpmcUHpmcHpm,cplUHpmSdcUHpmcSd,               & 
& cplUHpmSecUHpmcSe,cplUHpmSucUHpmcSu,cplUHpmSvcUHpmcSv,cplUHpmcUHpmVPVP,cplUHpmcUHpmcVWmVWm,& 
& cplUHpmcUHpmVZVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MHpm(2),MHpm2(2),MAh(2),MAh2(2),MVWm,MVWm2,MChi(4),MChi2(4),MCha(2),MCha2(2),         & 
& MFu(3),MFu2(3),MFd(3),MFd2(3),MFe(3),MFe2(3),Mhh(2),Mhh2(2),MVZ,MVZ2,MSu(6),           & 
& MSu2(6),MSd(6),MSd2(6),MSv(3),MSv2(3),MSe(6),MSe2(6)

Complex(dp), Intent(in) :: cplAhHpmcUHpm(2,2,2),cplAhcUHpmVWm(2,2),cplChiChacUHpmL(4,2,2),cplChiChacUHpmR(4,2,2),& 
& cplcFuFdcUHpmL(3,3,2),cplcFuFdcUHpmR(3,3,2),cplFvFecUHpmL(3,3,2),cplFvFecUHpmR(3,3,2), & 
& cplcgZgWmcUHpm(2),cplcgWmgZUHpm(2),cplcgWpCgZcUHpm(2),cplcgZgWpCUHpm(2),               & 
& cplhhHpmcUHpm(2,2,2),cplhhcUHpmVWm(2,2),cplHpmcUHpmVP(2,2),cplHpmcUHpmVZ(2,2),         & 
& cplSdcUHpmcSu(6,2,6),cplSecUHpmcSv(6,2,3),cplcUHpmVPVWm(2),cplcUHpmVWmVZ(2),           & 
& cplAhAhUHpmcUHpm(2,2,2,2),cplhhhhUHpmcUHpm(2,2,2,2),cplUHpmHpmcUHpmcHpm(2,2,2,2),      & 
& cplUHpmSdcUHpmcSd(2,6,2,6),cplUHpmSecUHpmcSe(2,6,2,6),cplUHpmSucUHpmcSu(2,6,2,6),      & 
& cplUHpmSvcUHpmcSv(2,3,2,3),cplUHpmcUHpmVPVP(2,2),cplUHpmcUHpmcVWmVWm(2,2),             & 
& cplUHpmcUHpmVZVZ(2,2)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2,2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumI(2,2) 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! Hpm, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MHpm2(i1),MAh2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplAhHpmcUHpm(i2,i1,gO1)
coup2 = Conjg(cplAhHpmcUHpm(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, Ah 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,MAh2(i2),MVWm2) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplAhcUHpmVWm(i2,gO1)
coup2 =  Conjg(cplAhcUHpmVWm(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! Chi, Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 2
 G0m2 = Real(SA_Gloop(p2,MChi2(i1),MCha2(i2)),dp) 
B0m2 = -2._dp*MChi(i1)*MCha(i2)*Real(SA_B0(p2,MChi2(i1),MCha2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplChiChacUHpmL(i1,i2,gO1)
coupR1 = cplChiChacUHpmR(i1,i2,gO1)
coupL2 =  Conjg(cplChiChacUHpmL(i1,i2,gO2))
coupR2 =  Conjg(cplChiChacUHpmR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,MFu2(i1),MFd2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*MFd(i2)*Real(SA_B0(p2,MFu2(i1),MFd2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplcFuFdcUHpmL(i1,i2,gO1)
coupR1 = cplcFuFdcUHpmR(i1,i2,gO1)
coupL2 =  Conjg(cplcFuFdcUHpmL(i1,i2,gO2))
coupR2 =  Conjg(cplcFuFdcUHpmR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 G0m2 = Real(SA_Gloop(p2,0._dp,MFe2(i2)),dp) 
B0m2 = -2._dp*0.*MFe(i2)*Real(SA_B0(p2,0._dp,MFe2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coupL1 = cplFvFecUHpmL(i1,i2,gO1)
coupR1 = cplFvFecUHpmR(i1,i2,gO1)
coupL2 =  Conjg(cplFvFecUHpmL(i1,i2,gO2))
coupR2 =  Conjg(cplFvFecUHpmR(i1,i2,gO2))
    SumI(gO1,gO2) = (coupL1*coupL2+coupR1*coupR2)*G0m2 & 
                & + (coupL1*coupR2+coupR1*coupL2)*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gZ], gWm 
!------------------------ 
sumI = 0._dp 
 
F0m2 = -Real(SA_B0(p2,MVWm2*RXi,MVZ2*RXi),dp) 
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcgZgWmcUHpm(gO1)
coup2 =  cplcgWmgZUHpm(gO2) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gZ 
!------------------------ 
sumI = 0._dp 
 
F0m2 = -Real(SA_B0(p2,MVZ2*RXi,MVWm2*RXi),dp) 
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcgWpCgZcUHpm(gO1)
coup2 =  cplcgZgWpCUHpm(gO2) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! Hpm, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B0m2 = Real(SA_B0(p2,MHpm2(i1),Mhh2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplhhHpmcUHpm(i2,i1,gO1)
coup2 = Conjg(cplhhHpmcUHpm(i2,i1,gO2))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,Mhh2(i2),MVWm2) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplhhcUHpmVWm(i2,gO1)
coup2 =  Conjg(cplhhcUHpmVWm(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VP, Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,MHpm2(i2),0._dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplHpmcUHpmVP(i2,gO1)
coup2 =  Conjg(cplHpmcUHpmVP(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VZ, Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 F0m2 = FloopRXi(p2,MHpm2(i2),MVZ2) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplHpmcUHpmVZ(i2,gO1)
coup2 =  Conjg(cplHpmcUHpmVZ(i2,gO2))
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Su], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSu2(i1),MSd2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplSdcUHpmcSu(i2,gO1,i1)
coup2 = Conjg(cplSdcUHpmcSu(i2,gO2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 6
 B0m2 = Real(SA_B0(p2,MSv2(i1),MSe2(i2)),dp) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplSecUHpmcSv(i2,gO1,i1)
coup2 = Conjg(cplSecUHpmcSv(i2,gO2,i1))
    SumI(gO1,gO2) = coup1*coup2*B0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, VP 
!------------------------ 
sumI = 0._dp 
 
F0m2 = Real(SVVloop(p2,0._dp,MVWm2),dp)   
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcUHpmVPVWm(gO1)
coup2 =  Conjg(cplcUHpmVPVWm(gO2)) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! VZ, VWm 
!------------------------ 
sumI = 0._dp 
 
F0m2 = Real(SVVloop(p2,MVWm2,MVZ2),dp)   
 Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplcUHpmVWmVZ(gO1)
coup2 =  Conjg(cplcUHpmVWmVZ(gO2)) 
    SumI(gO1,gO2) = coup1*coup2*F0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
!------------------------ 
! Ah, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MAh2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplAhAhUHpmcUHpm(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(Mhh2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplhhhhUHpmcUHpm(i1,i1,gO2,gO1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 A0m2 = SA_A0(MHpm2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUHpmHpmcUHpmcHpm(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSd2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUHpmSdcUHpmcSd(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSe2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUHpmSecUHpmcSe(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
 A0m2 = SA_A0(MSu2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUHpmSucUHpmcSu(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
 A0m2 = SA_A0(MSv2(i1)) 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUHpmSvcUHpmcSv(gO2,i1,gO1,i1)
    SumI(gO1,gO2) = -coup1*A0m2 
   End Do 
End Do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! VP, VP 
!------------------------ 
sumI = 0._dp 
 
!------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVWm2) + 0.25_dp*RXi*SA_A0(MVWm2*RXi) - 0.5_dp*MVWm2*rMS 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUHpmcUHpmcVWmVWm(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +4._dp* SumI  
!------------------------ 
! VZ, VZ 
!------------------------ 
sumI = 0._dp 
 
A0m2 = 0.75_dp*SA_A0(MVZ2) + 0.25_dp*RXi*SA_A0(MVZ2*RXi) - 0.5_dp*MVZ2*rMS 
Do gO1 = 1, 2
  Do gO2 = gO1, 2
coup1 = cplUHpmcUHpmVZVZ(gO2,gO1)
    SumI(gO1,gO2) = coup1*A0m2 
   End Do 
End Do 
res = res +2._dp* SumI  


Do gO2 = 1, 2
  Do gO1 = gO2+1, 2
     res(gO1,gO2) = Conjg(res(gO2,gO1))  
   End Do 
End Do 
 
res = oo16pi2*res 
 
End Subroutine Pi1LoopHpm 
 
Subroutine OneLoopChi(g1,g2,Mu,M1,M2,vd,vu,MChi,MChi2,MAh,MAh2,MHpm,MHpm2,            & 
& MCha,MCha2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MSd,MSd2,MFd,MFd2,MSe,MSe2,MFe,MFe2,           & 
& MSu,MSu2,MFu,MFu2,MSv,MSv2,cplUChiChiAhL,cplUChiChiAhR,cplUChiChacHpmL,cplUChiChacHpmR,& 
& cplUChiChacVWmL,cplUChiChacVWmR,cplUChiChihhL,cplUChiChihhR,cplUChiChiVZL,             & 
& cplUChiChiVZR,cplUChiFdcSdL,cplUChiFdcSdR,cplUChiFecSeL,cplUChiFecSeR,cplUChiFucSuL,   & 
& cplUChiFucSuR,cplUChiFvcSvL,cplUChiFvcSvR,cplUChiFvSvL,cplUChiFvSvR,cplcChaUChiHpmL,   & 
& cplcChaUChiHpmR,cplcFdUChiSdL,cplcFdUChiSdR,cplcFeUChiSeL,cplcFeUChiSeR,               & 
& cplcFuUChiSuL,cplcFuUChiSuR,cplcChaUChiVWmL,cplcChaUChiVWmR,delta,MChi_1L,             & 
& MChi2_1L,ZN_1L,ierr)

Implicit None 
Real(dp), Intent(in) :: MChi(4),MChi2(4),MAh(2),MAh2(2),MHpm(2),MHpm2(2),MCha(2),MCha2(2),MVWm,               & 
& MVWm2,Mhh(2),Mhh2(2),MVZ,MVZ2,MSd(6),MSd2(6),MFd(3),MFd2(3),MSe(6),MSe2(6),            & 
& MFe(3),MFe2(3),MSu(6),MSu2(6),MFu(3),MFu2(3),MSv(3),MSv2(3)

Real(dp), Intent(in) :: g1,g2,vd,vu

Complex(dp), Intent(in) :: Mu,M1,M2

Complex(dp), Intent(in) :: cplUChiChiAhL(4,4,2),cplUChiChiAhR(4,4,2),cplUChiChacHpmL(4,2,2),cplUChiChacHpmR(4,2,2),& 
& cplUChiChacVWmL(4,2),cplUChiChacVWmR(4,2),cplUChiChihhL(4,4,2),cplUChiChihhR(4,4,2),   & 
& cplUChiChiVZL(4,4),cplUChiChiVZR(4,4),cplUChiFdcSdL(4,3,6),cplUChiFdcSdR(4,3,6),       & 
& cplUChiFecSeL(4,3,6),cplUChiFecSeR(4,3,6),cplUChiFucSuL(4,3,6),cplUChiFucSuR(4,3,6),   & 
& cplUChiFvcSvL(4,3,3),cplUChiFvcSvR(4,3,3),cplUChiFvSvL(4,3,3),cplUChiFvSvR(4,3,3),     & 
& cplcChaUChiHpmL(2,4,2),cplcChaUChiHpmR(2,4,2),cplcFdUChiSdL(3,4,6),cplcFdUChiSdR(3,4,6),& 
& cplcFeUChiSeL(3,4,6),cplcFeUChiSeR(3,4,6),cplcFuUChiSuL(3,4,6),cplcFuUChiSuR(3,4,6),   & 
& cplcChaUChiVWmL(2,4),cplcChaUChiVWmR(2,4)

Complex(dp) :: mat1a(4,4), mat1(4,4), mat2(4,4) 
Integer , Intent(inout):: ierr 
Integer :: i1,i2,i3,i4,j1, j2,j3,j4,il,i_count 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(4), test_m2(4),p2 
Real(dp), Intent(out) :: MChi_1L(4),MChi2_1L(4) 
Complex(dp), Intent(out) ::  ZN_1L(4,4) 
Real(dp) :: MChi_t(4),MChi2_t(4) 
Complex(dp) ::  ZN_t(4,4) 
Complex(dp) ::  phaseM, E4(4), sigL(4,4), sigR(4,4), sigSL(4,4), sigSR(4,4) 
Real(dp) :: ZNa(4,4), test(2), eig(4) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopMChi'
 
mat1a(1,1) = 0._dp 
mat1a(1,1) = mat1a(1,1)+M1
mat1a(1,2) = 0._dp 
mat1a(1,3) = 0._dp 
mat1a(1,3) = mat1a(1,3)-(g1*vd)/2._dp
mat1a(1,4) = 0._dp 
mat1a(1,4) = mat1a(1,4)+(g1*vu)/2._dp
mat1a(2,1) = 0._dp 
mat1a(2,2) = 0._dp 
mat1a(2,2) = mat1a(2,2)+M2
mat1a(2,3) = 0._dp 
mat1a(2,3) = mat1a(2,3)+(g2*vd)/2._dp
mat1a(2,4) = 0._dp 
mat1a(2,4) = mat1a(2,4)-(g2*vu)/2._dp
mat1a(3,1) = 0._dp 
mat1a(3,1) = mat1a(3,1)-(g1*vd)/2._dp
mat1a(3,2) = 0._dp 
mat1a(3,2) = mat1a(3,2)+(g2*vd)/2._dp
mat1a(3,3) = 0._dp 
mat1a(3,4) = 0._dp 
mat1a(3,4) = mat1a(3,4)-1._dp*(Mu)
mat1a(4,1) = 0._dp 
mat1a(4,1) = mat1a(4,1)+(g1*vu)/2._dp
mat1a(4,2) = 0._dp 
mat1a(4,2) = mat1a(4,2)-(g2*vu)/2._dp
mat1a(4,3) = 0._dp 
mat1a(4,3) = mat1a(4,3)-1._dp*(Mu)
mat1a(4,4) = 0._dp 

 
 !---------------------------------------- 
! Rotation matrix for p2=0 
!----------------------------------------- 
 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = 0._dp 
Call Sigma1LoopChi(p2,MChi,MChi2,MAh,MAh2,MHpm,MHpm2,MCha,MCha2,MVWm,MVWm2,           & 
& Mhh,Mhh2,MVZ,MVZ2,MSd,MSd2,MFd,MFd2,MSe,MSe2,MFe,MFe2,MSu,MSu2,MFu,MFu2,               & 
& MSv,MSv2,cplUChiChiAhL,cplUChiChiAhR,cplUChiChacHpmL,cplUChiChacHpmR,cplUChiChacVWmL,  & 
& cplUChiChacVWmR,cplUChiChihhL,cplUChiChihhR,cplUChiChiVZL,cplUChiChiVZR,               & 
& cplUChiFdcSdL,cplUChiFdcSdR,cplUChiFecSeL,cplUChiFecSeR,cplUChiFucSuL,cplUChiFucSuR,   & 
& cplUChiFvcSvL,cplUChiFvcSvR,cplUChiFvSvL,cplUChiFvSvR,cplcChaUChiHpmL,cplcChaUChiHpmR, & 
& cplcFdUChiSdL,cplcFdUChiSdR,cplcFeUChiSeL,cplcFeUChiSeR,cplcFuUChiSuL,cplcFuUChiSuR,   & 
& cplcChaUChiVWmL,cplcChaUChiVWmR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - 0.5_dp*(Conjg(SigSL) + SigSR + & 
      & 0.5_dp*(MatMul(Transpose(SigL),mat1a) + MatMul(SigR,mat1a) + & 
      & MatMul(mat1a,Transpose(SigR)) + MatMul(mat1a,SigL))) 
 
If (ForceRealMatrices) mat1 = Real(mat1,dp) 
If (Maxval(Abs(Aimag(mat1))).Eq.0._dp) Then 
Call EigenSystem(Real(mat1,dp),Eig,ZNa,ierr,test) 
 
   Do i1=1,4
   If (Eig(i1).Lt.0._dp) Then 
    MChi_t(i1) = - Eig(i1) 
    ZN_1L(i1,:) = (0._dp,1._dp)*ZNa(i1,:) 
   Else 
    MChi_t(i1) = Eig(i1) 
    ZN_1L(i1,:) = ZNa(i1,:)
    End If 
   End Do 
 
Do i1=1,3
  Do i2=i1+1,4
    If (Abs(MChi_t(i1)).Gt.Abs(MChi_t(i2))) Then 
      Eig(1) = MChi_t(i1) 
      MChi_t(i1) = MChi_t(i2) 
      MChi_t(i2) = Eig(1) 
      E4 = ZN_1L(i1,:) 
      ZN_1L(i1,:) = ZN_1L(i2,:) 
      ZN_1L(i2,:) = E4
    End If 
   End Do 
End Do 
 
ZNOS_0 = ZN_1L 
 Else 
 
mat2 = Matmul( Transpose(Conjg( mat1) ), mat1 ) 
Call Eigensystem(mat2, Eig, ZN_1L, ierr, test) 
mat2 = Matmul( Conjg(ZN_1L), Matmul( mat1, Transpose( Conjg(ZN_1L)))) 
Do i1=1,4
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
  phaseM = Sqrt( mat2(i1,i1) / Abs(mat2(i1,i1))) 
  ZN_1L(i1,:)= phaseM * ZN_1L(i1,:) 
End if 
ZNOS_0 = ZN_1L 
 End Do 
End If 
 
!---------------------------------------- 
! Now, with momenta
!----------------------------------------- 
 
Do il=4,1,-1
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MChi2(il)
Call Sigma1LoopChi(p2,MChi,MChi2,MAh,MAh2,MHpm,MHpm2,MCha,MCha2,MVWm,MVWm2,           & 
& Mhh,Mhh2,MVZ,MVZ2,MSd,MSd2,MFd,MFd2,MSe,MSe2,MFe,MFe2,MSu,MSu2,MFu,MFu2,               & 
& MSv,MSv2,cplUChiChiAhL,cplUChiChiAhR,cplUChiChacHpmL,cplUChiChacHpmR,cplUChiChacVWmL,  & 
& cplUChiChacVWmR,cplUChiChihhL,cplUChiChihhR,cplUChiChiVZL,cplUChiChiVZR,               & 
& cplUChiFdcSdL,cplUChiFdcSdR,cplUChiFecSeL,cplUChiFecSeR,cplUChiFucSuL,cplUChiFucSuR,   & 
& cplUChiFvcSvL,cplUChiFvcSvR,cplUChiFvSvL,cplUChiFvSvR,cplcChaUChiHpmL,cplcChaUChiHpmR, & 
& cplcFdUChiSdL,cplcFdUChiSdR,cplcFeUChiSeL,cplcFeUChiSeR,cplcFuUChiSuL,cplcFuUChiSuR,   & 
& cplcChaUChiVWmL,cplcChaUChiVWmR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - 0.5_dp*(Conjg(SigSL) + SigSR + & 
      & 0.5_dp*(MatMul(Transpose(SigL),mat1a) + MatMul(SigR,mat1a) + & 
      & MatMul(mat1a,Transpose(SigR)) + MatMul(mat1a,SigL))) 
 
If (ForceRealMatrices) mat1 = Real(mat1,dp) 
If (Maxval(Abs(Aimag(mat1))).Eq.0._dp) Then 
Call EigenSystem(Real(mat1,dp),Eig,ZNa,ierr,test) 
 
   Do i1=1,4
   If (Eig(i1).Lt.0._dp) Then 
    MChi_t(i1) = - Eig(i1) 
    ZN_1L(i1,:) = (0._dp,1._dp)*ZNa(i1,:) 
   Else 
    MChi_t(i1) = Eig(i1) 
    ZN_1L(i1,:) = ZNa(i1,:)
    End If 
   End Do 
 
Do i1=1,3
  Do i2=i1+1,4
    If (Abs(MChi_t(i1)).Gt.Abs(MChi_t(i2))) Then 
      Eig(1) = MChi_t(i1) 
      MChi_t(i1) = MChi_t(i2) 
      MChi_t(i2) = Eig(1) 
      E4 = ZN_1L(i1,:) 
      ZN_1L(i1,:) = ZN_1L(i2,:) 
      ZN_1L(i2,:) = E4
    End If 
   End Do 
End Do 
 
MChi_1L(iL) = MChi_t(iL) 
MChi2_1L(iL) = MChi_t(iL)**2 
ZNOS_p2(il,:) = ZN_1L(il,:) 
 Else 
 
mat2 = Matmul( Transpose(Conjg( mat1) ), mat1 ) 
Call Eigensystem(mat2, Eig, ZN_1L, ierr, test) 
mat2 = Matmul( Conjg(ZN_1L), Matmul( mat1, Transpose( Conjg(ZN_1L)))) 
Do i1=1,4
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
  phaseM = Sqrt( mat2(i1,i1) / Abs(mat2(i1,i1))) 
  ZN_1L(i1,:)= phaseM * ZN_1L(i1,:) 
End if 
ZNOS_p2(il,:) = ZN_1L(il,:) 
   If (Eig(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,Eig(i1) 
  Eig(i1) = 1._dp 
   SignOfMassChanged = .True. 
 End if 
End Do 
MChi_1L = Sqrt( Eig ) 
 
MChi2_1L = Eig 
 
End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
!---------------------------------------- 
! Redoing Calculation using redefined p2 
!----------------------------------------- 
 
i_count = 0 
p2_loop: Do  
i_count = i_count + 1 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MChi2_1L(iL)
Call Sigma1LoopChi(p2,MChi,MChi2,MAh,MAh2,MHpm,MHpm2,MCha,MCha2,MVWm,MVWm2,           & 
& Mhh,Mhh2,MVZ,MVZ2,MSd,MSd2,MFd,MFd2,MSe,MSe2,MFe,MFe2,MSu,MSu2,MFu,MFu2,               & 
& MSv,MSv2,cplUChiChiAhL,cplUChiChiAhR,cplUChiChacHpmL,cplUChiChacHpmR,cplUChiChacVWmL,  & 
& cplUChiChacVWmR,cplUChiChihhL,cplUChiChihhR,cplUChiChiVZL,cplUChiChiVZR,               & 
& cplUChiFdcSdL,cplUChiFdcSdR,cplUChiFecSeL,cplUChiFecSeR,cplUChiFucSuL,cplUChiFucSuR,   & 
& cplUChiFvcSvL,cplUChiFvcSvR,cplUChiFvSvL,cplUChiFvSvR,cplcChaUChiHpmL,cplcChaUChiHpmR, & 
& cplcFdUChiSdL,cplcFdUChiSdR,cplcFeUChiSeL,cplcFeUChiSeR,cplcFuUChiSuL,cplcFuUChiSuR,   & 
& cplcChaUChiVWmL,cplcChaUChiVWmR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - 0.5_dp*(Conjg(SigSL) + SigSR + & 
      & 0.5_dp*(MatMul(Transpose(SigL),mat1a) + MatMul(SigR,mat1a) + & 
      & MatMul(mat1a,Transpose(SigR)) + MatMul(mat1a,SigL))) 
 
If (ForceRealMatrices) mat1 = Real(mat1,dp) 
If (Maxval(Abs(Aimag(mat1))).Eq.0._dp) Then 
Call EigenSystem(Real(mat1,dp),Eig,ZNa,ierr,test) 
 
   Do i1=1,4
   If (Eig(i1).Lt.0._dp) Then 
    MChi_t(i1) = - Eig(i1) 
    ZN_1L(i1,:) = (0._dp,1._dp)*ZNa(i1,:) 
   Else 
    MChi_t(i1) = Eig(i1) 
    ZN_1L(i1,:) = ZNa(i1,:)
    End If 
   End Do 
 
Do i1=1,3
  Do i2=i1+1,4
    If (Abs(MChi_t(i1)).Gt.Abs(MChi_t(i2))) Then 
      Eig(1) = MChi_t(i1) 
      MChi_t(i1) = MChi_t(i2) 
      MChi_t(i2) = Eig(1) 
      E4 = ZN_1L(i1,:) 
      ZN_1L(i1,:) = ZN_1L(i2,:) 
      ZN_1L(i2,:) = E4
    End If 
   End Do 
End Do 
 
MChi_1L(iL) = MChi_t(iL) 
MChi2_1L(iL) = MChi_t(iL)**2 
ZNOS_p2(il,:) = ZN_1L(il,:) 
 Else 
 
mat2 = Matmul( Transpose(Conjg( mat1) ), mat1 ) 
Call Eigensystem(mat2, Eig, ZN_1L, ierr, test) 
mat2 = Matmul( Conjg(ZN_1L), Matmul( mat1, Transpose( Conjg(ZN_1L)))) 
Do i1=1,4
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
  phaseM = Sqrt( mat2(i1,i1) / Abs(mat2(i1,i1))) 
  ZN_1L(i1,:)= phaseM * ZN_1L(i1,:) 
End if 
ZNOS_p2(il,:) = ZN_1L(il,:) 
   If (Eig(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,Eig(i1) 
  Eig(i1) = 1._dp 
   SignOfMassChanged = .True. 
 End if 
End Do 
MChi_1L = Sqrt( Eig ) 
 
MChi2_1L = Eig 
 
End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
If (p2.Ne.0._dp) Then 
  test(1) = Abs(MChi2_1L(il)-p2)/p2
Else 
  test(2) = Abs(MChi2_1L(il))
End If 
If (Abs(MChi2_1L(il)).lt.1.0E-30_dp) Exit p2_loop 
If (test(1).lt.0.1_dp*delta) Exit p2_loop 
If(i_count.gt.30) then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Exit p2_loop 
End if
End Do p2_loop 
End Do 
 
Iname = Iname -1 
End Subroutine OneLoopChi
 
 
Subroutine Sigma1LoopChi(p2,MChi,MChi2,MAh,MAh2,MHpm,MHpm2,MCha,MCha2,MVWm,           & 
& MVWm2,Mhh,Mhh2,MVZ,MVZ2,MSd,MSd2,MFd,MFd2,MSe,MSe2,MFe,MFe2,MSu,MSu2,MFu,              & 
& MFu2,MSv,MSv2,cplUChiChiAhL,cplUChiChiAhR,cplUChiChacHpmL,cplUChiChacHpmR,             & 
& cplUChiChacVWmL,cplUChiChacVWmR,cplUChiChihhL,cplUChiChihhR,cplUChiChiVZL,             & 
& cplUChiChiVZR,cplUChiFdcSdL,cplUChiFdcSdR,cplUChiFecSeL,cplUChiFecSeR,cplUChiFucSuL,   & 
& cplUChiFucSuR,cplUChiFvcSvL,cplUChiFvcSvR,cplUChiFvSvL,cplUChiFvSvR,cplcChaUChiHpmL,   & 
& cplcChaUChiHpmR,cplcFdUChiSdL,cplcFdUChiSdR,cplcFeUChiSeL,cplcFeUChiSeR,               & 
& cplcFuUChiSuL,cplcFuUChiSuR,cplcChaUChiVWmL,cplcChaUChiVWmR,sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MChi(4),MChi2(4),MAh(2),MAh2(2),MHpm(2),MHpm2(2),MCha(2),MCha2(2),MVWm,               & 
& MVWm2,Mhh(2),Mhh2(2),MVZ,MVZ2,MSd(6),MSd2(6),MFd(3),MFd2(3),MSe(6),MSe2(6),            & 
& MFe(3),MFe2(3),MSu(6),MSu2(6),MFu(3),MFu2(3),MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplUChiChiAhL(4,4,2),cplUChiChiAhR(4,4,2),cplUChiChacHpmL(4,2,2),cplUChiChacHpmR(4,2,2),& 
& cplUChiChacVWmL(4,2),cplUChiChacVWmR(4,2),cplUChiChihhL(4,4,2),cplUChiChihhR(4,4,2),   & 
& cplUChiChiVZL(4,4),cplUChiChiVZR(4,4),cplUChiFdcSdL(4,3,6),cplUChiFdcSdR(4,3,6),       & 
& cplUChiFecSeL(4,3,6),cplUChiFecSeR(4,3,6),cplUChiFucSuL(4,3,6),cplUChiFucSuR(4,3,6),   & 
& cplUChiFvcSvL(4,3,3),cplUChiFvcSvR(4,3,3),cplUChiFvSvL(4,3,3),cplUChiFvSvR(4,3,3),     & 
& cplcChaUChiHpmL(2,4,2),cplcChaUChiHpmR(2,4,2),cplcFdUChiSdL(3,4,6),cplcFdUChiSdR(3,4,6),& 
& cplcFeUChiSeL(3,4,6),cplcFeUChiSeR(3,4,6),cplcFuUChiSuL(3,4,6),cplcFuUChiSuR(3,4,6),   & 
& cplcChaUChiVWmL(2,4),cplcChaUChiVWmR(2,4)

Complex(dp), Intent(out) :: SigL(4,4),SigR(4,4), SigSL(4,4), SigSR(4,4) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(4,4), sumR(4,4), sumSL(4,4), sumSR(4,4) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Chi, Ah 
!------------------------ 
    Do i1 = 1, 4
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MChi2(i1),MAh2(i2)),dp) 
B0m2 = 2._dp*MChi(i1)*Real(SA_B0(p2,MChi2(i1),MAh2(i2)),dp) 
coupL1 = cplUChiChiAhL(gO1,i1,i2)
coupR1 = cplUChiChiAhR(gO1,i1,i2)
coupL2 =  Conjg(cplUChiChiAhL(gO2,i1,i2))
coupR2 =  Conjg(cplUChiChiAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], Cha 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MCha2(i2),MHpm2(i1)),dp) 
B0m2 = 2._dp*MCha(i2)*Real(SA_B0(p2,MCha2(i2),MHpm2(i1)),dp) 
coupL1 = cplUChiChacHpmL(gO1,i2,i1)
coupR1 = cplUChiChacHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplUChiChacHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplUChiChacHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Cha 
!------------------------ 
      Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -4._dp*(Real(SA_B1(p2,MCha2(i2),MVWm2),dp) + 0.5_dp*rMS) 
B0m2 = -8._dp*MCha(i2)*(Real(SA_B0(p2,MCha2(i2),MVWm2),dp) - 0.5_dp*rMS) 
coupL1 = cplUChiChacVWmL(gO1,i2)
coupR1 = cplUChiChacVWmR(gO1,i2)
coupL2 =  Conjg(cplUChiChacVWmL(gO2,i2))
coupR2 =  Conjg(cplUChiChacVWmR(gO2,i2))
SumSR(gO1,gO2) = coupL2*coupR1*B0m2 
SumSL(gO1,gO2) = coupR2*coupL1*B0m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
    End Do 
 !------------------------ 
! hh, Chi 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MChi2(i2),Mhh2(i1)),dp) 
B0m2 = 2._dp*MChi(i2)*Real(SA_B0(p2,MChi2(i2),Mhh2(i1)),dp) 
coupL1 = cplUChiChihhL(gO1,i2,i1)
coupR1 = cplUChiChihhR(gO1,i2,i1)
coupL2 =  Conjg(cplUChiChihhL(gO2,i2,i1))
coupR2 =  Conjg(cplUChiChihhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VZ, Chi 
!------------------------ 
      Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -4._dp*(Real(SA_B1(p2,MChi2(i2),MVZ2),dp) + 0.5_dp*rMS) 
B0m2 = -8._dp*MChi(i2)*(Real(SA_B0(p2,MChi2(i2),MVZ2),dp) - 0.5_dp*rMS) 
coupL1 = cplUChiChiVZL(gO1,i2)
coupR1 = cplUChiChiVZR(gO1,i2)
coupL2 =  Conjg(cplUChiChiVZL(gO2,i2))
coupR2 =  Conjg(cplUChiChiVZR(gO2,i2))
SumSR(gO1,gO2) = coupL2*coupR1*B0m2 
SumSL(gO1,gO2) = coupR2*coupL1*B0m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
    End Do 
 !------------------------ 
! conj[Sd], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),MSd2(i1)),dp) 
B0m2 = 2._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSd2(i1)),dp) 
coupL1 = cplUChiFdcSdL(gO1,i2,i1)
coupR1 = cplUChiFdcSdR(gO1,i2,i1)
coupL2 =  Conjg(cplUChiFdcSdL(gO2,i2,i1))
coupR2 =  Conjg(cplUChiFdcSdR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp/2._dp* sumL
SigR = SigR +3._dp/2._dp* sumR 
SigSL = SigSL +3._dp/2._dp* sumSL 
SigSR = SigSR +3._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Fe 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i2),MSe2(i1)),dp) 
B0m2 = 2._dp*MFe(i2)*Real(SA_B0(p2,MFe2(i2),MSe2(i1)),dp) 
coupL1 = cplUChiFecSeL(gO1,i2,i1)
coupR1 = cplUChiFecSeR(gO1,i2,i1)
coupL2 =  Conjg(cplUChiFecSeL(gO2,i2,i1))
coupR2 =  Conjg(cplUChiFecSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Fu 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),MSu2(i1)),dp) 
B0m2 = 2._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),MSu2(i1)),dp) 
coupL1 = cplUChiFucSuL(gO1,i2,i1)
coupR1 = cplUChiFucSuR(gO1,i2,i1)
coupL2 =  Conjg(cplUChiFucSuL(gO2,i2,i1))
coupR2 =  Conjg(cplUChiFucSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp/2._dp* sumL
SigR = SigR +3._dp/2._dp* sumR 
SigSL = SigSL +3._dp/2._dp* sumSL 
SigSR = SigSR +3._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Fv 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,0._dp,MSv2(i1)),dp) 
B0m2 = 2._dp*0.*Real(SA_B0(p2,0._dp,MSv2(i1)),dp) 
coupL1 = cplUChiFvcSvL(gO1,i2,i1)
coupR1 = cplUChiFvcSvR(gO1,i2,i1)
coupL2 =  Conjg(cplUChiFvcSvL(gO2,i2,i1))
coupR2 =  Conjg(cplUChiFvcSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sv, Fv 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,0._dp,MSv2(i1)),dp) 
B0m2 = 2._dp*0.*Real(SA_B0(p2,0._dp,MSv2(i1)),dp) 
coupL1 = cplUChiFvSvL(gO1,i2,i1)
coupR1 = cplUChiFvSvR(gO1,i2,i1)
coupL2 =  Conjg(cplUChiFvSvL(gO2,i2,i1))
coupR2 =  Conjg(cplUChiFvSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Hpm 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MCha2(i1),MHpm2(i2)),dp) 
B0m2 = 2._dp*MCha(i1)*Real(SA_B0(p2,MCha2(i1),MHpm2(i2)),dp) 
coupL1 = cplcChaUChiHpmL(i1,gO1,i2)
coupR1 = cplcChaUChiHpmR(i1,gO1,i2)
coupL2 =  Conjg(cplcChaUChiHpmL(i1,gO2,i2))
coupR2 =  Conjg(cplcChaUChiHpmR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i1),MSd2(i2)),dp) 
B0m2 = 2._dp*MFd(i1)*Real(SA_B0(p2,MFd2(i1),MSd2(i2)),dp) 
coupL1 = cplcFdUChiSdL(i1,gO1,i2)
coupR1 = cplcFdUChiSdR(i1,gO1,i2)
coupL2 =  Conjg(cplcFdUChiSdL(i1,gO2,i2))
coupR2 =  Conjg(cplcFdUChiSdR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp/2._dp* sumL
SigR = SigR +3._dp/2._dp* sumR 
SigSL = SigSL +3._dp/2._dp* sumSL 
SigSR = SigSR +3._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Se 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i1),MSe2(i2)),dp) 
B0m2 = 2._dp*MFe(i1)*Real(SA_B0(p2,MFe2(i1),MSe2(i2)),dp) 
coupL1 = cplcFeUChiSeL(i1,gO1,i2)
coupR1 = cplcFeUChiSeR(i1,gO1,i2)
coupL2 =  Conjg(cplcFeUChiSeL(i1,gO2,i2))
coupR2 =  Conjg(cplcFeUChiSeR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Su 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i1),MSu2(i2)),dp) 
B0m2 = 2._dp*MFu(i1)*Real(SA_B0(p2,MFu2(i1),MSu2(i2)),dp) 
coupL1 = cplcFuUChiSuL(i1,gO1,i2)
coupR1 = cplcFuUChiSuR(i1,gO1,i2)
coupL2 =  Conjg(cplcFuUChiSuL(i1,gO2,i2))
coupR2 =  Conjg(cplcFuUChiSuR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp/2._dp* sumL
SigR = SigR +3._dp/2._dp* sumR 
SigSL = SigSL +3._dp/2._dp* sumSL 
SigSR = SigSR +3._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], VWm 
!------------------------ 
    Do i1 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 4
  Do gO2 = 1, 4
B1m2 = -4._dp*(Real(SA_B1(p2,MCha2(i1),MVWm2),dp) + 0.5_dp*rMS) 
B0m2 = -8._dp*MCha(i1)*(Real(SA_B0(p2,MCha2(i1),MVWm2),dp) - 0.5_dp*rMS) 
coupL1 = cplcChaUChiVWmL(i1,gO1)
coupR1 = cplcChaUChiVWmR(i1,gO1)
coupL2 =  Conjg(cplcChaUChiVWmL(i1,gO2))
coupR2 =  Conjg(cplcChaUChiVWmR(i1,gO2))
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
 

SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopChi 
 
Subroutine OneLoopCha(g2,Mu,M2,vd,vu,MCha,MCha2,MAh,MAh2,Mhh,Mhh2,MVZ,MVZ2,           & 
& MHpm,MHpm2,MChi,MChi2,MVWm,MVWm2,MSu,MSu2,MFd,MFd2,MSv,MSv2,MFe,MFe2,MSe,              & 
& MSe2,MFu,MFu2,MSd,MSd2,cplcUChaChaAhL,cplcUChaChaAhR,cplcUChaChahhL,cplcUChaChahhR,    & 
& cplcUChaChaVPL,cplcUChaChaVPR,cplcUChaChaVZL,cplcUChaChaVZR,cplcUChaChiHpmL,           & 
& cplcUChaChiHpmR,cplcUChaChiVWmL,cplcUChaChiVWmR,cplcUChaFdcSuL,cplcUChaFdcSuR,         & 
& cplcUChaFecSvL,cplcUChaFecSvR,cplcUChaFvSeL,cplcUChaFvSeR,cplcUChacFuSdL,              & 
& cplcUChacFuSdR,delta,MCha_1L,MCha2_1L,UM_1L,UP_1L,ierr)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MAh(2),MAh2(2),Mhh(2),Mhh2(2),MVZ,MVZ2,MHpm(2),MHpm2(2),             & 
& MChi(4),MChi2(4),MVWm,MVWm2,MSu(6),MSu2(6),MFd(3),MFd2(3),MSv(3),MSv2(3),              & 
& MFe(3),MFe2(3),MSe(6),MSe2(6),MFu(3),MFu2(3),MSd(6),MSd2(6)

Real(dp), Intent(in) :: g2,vd,vu

Complex(dp), Intent(in) :: Mu,M2

Complex(dp), Intent(in) :: cplcUChaChaAhL(2,2,2),cplcUChaChaAhR(2,2,2),cplcUChaChahhL(2,2,2),cplcUChaChahhR(2,2,2),& 
& cplcUChaChaVPL(2,2),cplcUChaChaVPR(2,2),cplcUChaChaVZL(2,2),cplcUChaChaVZR(2,2),       & 
& cplcUChaChiHpmL(2,4,2),cplcUChaChiHpmR(2,4,2),cplcUChaChiVWmL(2,4),cplcUChaChiVWmR(2,4),& 
& cplcUChaFdcSuL(2,3,6),cplcUChaFdcSuR(2,3,6),cplcUChaFecSvL(2,3,3),cplcUChaFecSvR(2,3,3),& 
& cplcUChaFvSeL(2,3,6),cplcUChaFvSeR(2,3,6),cplcUChacFuSdL(2,3,6),cplcUChacFuSdR(2,3,6)

Complex(dp) :: mat1a(2,2), mat1(2,2) 
Integer , Intent(inout):: ierr 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(2), test_m2(2), p2 
Real(dp), Intent(out) :: MCha_1L(2),MCha2_1L(2) 
 Complex(dp), Intent(out) :: UM_1L(2,2), UP_1L(2,2) 
 
 Real(dp) :: MCha_t(2),MCha2_t(2) 
 Complex(dp) :: UM_t(2,2), UP_t(2,2), sigL(2,2), sigR(2,2), sigSL(2,2), sigSR(2,2) 
 
 Complex(dp) :: mat(2,2)=0._dp, mat2(2,2)=0._dp, phaseM 

Complex(dp) :: UM2(2,2), UP2(2,2) 
 
 Real(dp) :: UM1(2,2), UP1(2,2), test(2) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopMCha'
 
mat1a(1,1) = 0._dp 
mat1a(1,1) = mat1a(1,1)+M2
mat1a(1,2) = 0._dp 
mat1a(1,2) = mat1a(1,2)+(g2*vu)/sqrt(2._dp)
mat1a(2,1) = 0._dp 
mat1a(2,1) = mat1a(2,1)+(g2*vd)/sqrt(2._dp)
mat1a(2,2) = 0._dp 
mat1a(2,2) = mat1a(2,2)+Mu

 
 !---------------------------------------- 
! Rotation matrix for p2=0 
!----------------------------------------- 
 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = 0._dp 
Call Sigma1LoopCha(p2,MCha,MCha2,MAh,MAh2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,               & 
& MChi,MChi2,MVWm,MVWm2,MSu,MSu2,MFd,MFd2,MSv,MSv2,MFe,MFe2,MSe,MSe2,MFu,MFu2,           & 
& MSd,MSd2,cplcUChaChaAhL,cplcUChaChaAhR,cplcUChaChahhL,cplcUChaChahhR,cplcUChaChaVPL,   & 
& cplcUChaChaVPR,cplcUChaChaVZL,cplcUChaChaVZR,cplcUChaChiHpmL,cplcUChaChiHpmR,          & 
& cplcUChaChiVWmL,cplcUChaChiVWmR,cplcUChaFdcSuL,cplcUChaFdcSuR,cplcUChaFecSvL,          & 
& cplcUChaFecSvR,cplcUChaFvSeL,cplcUChaFvSeR,cplcUChacFuSdL,cplcUChacFuSdR,              & 
& sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MCha2_t,UP1,ierr,test) 
UP2 = UP1 
Else 
Call EigenSystem(mat2,MCha2_t,UP2,ierr,test) 
 End If 
 
UPOS_0 = UP2 
 mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MCha2_t,UM1,ierr,test) 
 
 
UM2 = UM1 
Else 
Call EigenSystem(mat2,MCha2_t,UM2,ierr,test) 
 
 
End If 
UM2 = Conjg(UM2) 
UMOS_0 = UM2 
 
!---------------------------------------- 
! Now, with momenta
!----------------------------------------- 
 
Do il=2,1,-1
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MCha2(il) 
Call Sigma1LoopCha(p2,MCha,MCha2,MAh,MAh2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,               & 
& MChi,MChi2,MVWm,MVWm2,MSu,MSu2,MFd,MFd2,MSv,MSv2,MFe,MFe2,MSe,MSe2,MFu,MFu2,           & 
& MSd,MSd2,cplcUChaChaAhL,cplcUChaChaAhR,cplcUChaChahhL,cplcUChaChahhR,cplcUChaChaVPL,   & 
& cplcUChaChaVPR,cplcUChaChaVZL,cplcUChaChaVZR,cplcUChaChiHpmL,cplcUChaChiHpmR,          & 
& cplcUChaChiVWmL,cplcUChaChiVWmR,cplcUChaFdcSuL,cplcUChaFdcSuR,cplcUChaFecSvL,          & 
& cplcUChaFecSvR,cplcUChaFvSeL,cplcUChaFvSeR,cplcUChacFuSdL,cplcUChacFuSdR,              & 
& sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MCha2_t,UP1,ierr,test) 
UP2 = UP1 
Else 
Call EigenSystem(mat2,MCha2_t,UP2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
!---------------------------------------- 
! Redoing Calculation using redefined p2 
!----------------------------------------- 
 
i_count = 0 
p2_loop: Do  
i_count = i_count + 1 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MCha2_t(iL)
Call Sigma1LoopCha(p2,MCha,MCha2,MAh,MAh2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,               & 
& MChi,MChi2,MVWm,MVWm2,MSu,MSu2,MFd,MFd2,MSv,MSv2,MFe,MFe2,MSe,MSe2,MFu,MFu2,           & 
& MSd,MSd2,cplcUChaChaAhL,cplcUChaChaAhR,cplcUChaChahhL,cplcUChaChahhR,cplcUChaChaVPL,   & 
& cplcUChaChaVPR,cplcUChaChaVZL,cplcUChaChaVZR,cplcUChaChiHpmL,cplcUChaChiHpmR,          & 
& cplcUChaChiVWmL,cplcUChaChiVWmR,cplcUChaFdcSuL,cplcUChaFdcSuR,cplcUChaFecSvL,          & 
& cplcUChaFecSvR,cplcUChaFvSeL,cplcUChaFvSeR,cplcUChacFuSdL,cplcUChacFuSdR,              & 
& sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MCha2_t,UP1,ierr,test) 
UP2 = UP1 
Else 
Call EigenSystem(mat2,MCha2_t,UP2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
MCha2_1L(il) = MCha2_t(il) 
MCha_1L(il) = Sqrt(MCha2_1L(il)) 
 
If (p2.Ne.0._dp) Then 
  test(1) = Abs(MCha2_1L(il)-p2)/p2
Else 
  test(2) = Abs(MCha2_1L(il))
End If 
If (Abs(MCha2_1L(il)).lt.1.0E-30_dp) Exit p2_loop 
If (test(1).lt.0.1_dp*delta) Exit p2_loop 
If(i_count.gt.30) then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Exit p2_loop 
End if
End Do p2_loop 
mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MCha2_t,UM1,ierr,test) 
 
 
UM2 = UM1 
Else 
Call EigenSystem(mat2,MCha2_t,UM2,ierr,test) 
 
 
End If 
UM2 = Conjg(UM2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(UM2),mat1),Transpose( Conjg(UP2))) 
Do i1=1,2
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
UP2(i1,:) = phaseM *UP2(i1,:) 
 End if 
End Do 
 
UMOS_p2(il,:) = UM2(il,:) 
 UPOS_p2(il,:) = UP2(il,:) 
 UM_1L = UM2 
 UP_1L = UP2 
 End Do  
 
Iname = Iname -1 
End Subroutine OneLoopCha
 
 
Subroutine Sigma1LoopCha(p2,MCha,MCha2,MAh,MAh2,Mhh,Mhh2,MVZ,MVZ2,MHpm,               & 
& MHpm2,MChi,MChi2,MVWm,MVWm2,MSu,MSu2,MFd,MFd2,MSv,MSv2,MFe,MFe2,MSe,MSe2,              & 
& MFu,MFu2,MSd,MSd2,cplcUChaChaAhL,cplcUChaChaAhR,cplcUChaChahhL,cplcUChaChahhR,         & 
& cplcUChaChaVPL,cplcUChaChaVPR,cplcUChaChaVZL,cplcUChaChaVZR,cplcUChaChiHpmL,           & 
& cplcUChaChiHpmR,cplcUChaChiVWmL,cplcUChaChiVWmR,cplcUChaFdcSuL,cplcUChaFdcSuR,         & 
& cplcUChaFecSvL,cplcUChaFecSvR,cplcUChaFvSeL,cplcUChaFvSeR,cplcUChacFuSdL,              & 
& cplcUChacFuSdR,sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MAh(2),MAh2(2),Mhh(2),Mhh2(2),MVZ,MVZ2,MHpm(2),MHpm2(2),             & 
& MChi(4),MChi2(4),MVWm,MVWm2,MSu(6),MSu2(6),MFd(3),MFd2(3),MSv(3),MSv2(3),              & 
& MFe(3),MFe2(3),MSe(6),MSe2(6),MFu(3),MFu2(3),MSd(6),MSd2(6)

Complex(dp), Intent(in) :: cplcUChaChaAhL(2,2,2),cplcUChaChaAhR(2,2,2),cplcUChaChahhL(2,2,2),cplcUChaChahhR(2,2,2),& 
& cplcUChaChaVPL(2,2),cplcUChaChaVPR(2,2),cplcUChaChaVZL(2,2),cplcUChaChaVZR(2,2),       & 
& cplcUChaChiHpmL(2,4,2),cplcUChaChiHpmR(2,4,2),cplcUChaChiVWmL(2,4),cplcUChaChiVWmR(2,4),& 
& cplcUChaFdcSuL(2,3,6),cplcUChaFdcSuR(2,3,6),cplcUChaFecSvL(2,3,3),cplcUChaFecSvR(2,3,3),& 
& cplcUChaFvSeL(2,3,6),cplcUChaFvSeR(2,3,6),cplcUChacFuSdL(2,3,6),cplcUChacFuSdR(2,3,6)

Complex(dp), Intent(out) :: SigL(2,2),SigR(2,2), SigSL(2,2), SigSR(2,2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(2,2), sumR(2,2), sumSL(2,2), sumSR(2,2) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
 SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Cha, Ah 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -Real(SA_B1(p2,MCha2(i1),MAh2(i2)),dp) 
B0m2 = MCha(i1)*Real(SA_B0(p2,MCha2(i1),MAh2(i2)),dp) 
coupL1 = cplcUChaChaAhL(gO1,i1,i2)
coupR1 = cplcUChaChaAhR(gO1,i1,i2)
coupL2 =  Conjg(cplcUChaChaAhL(gO2,i1,i2))
coupR2 =  Conjg(cplcUChaChaAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! hh, Cha 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -Real(SA_B1(p2,MCha2(i2),Mhh2(i1)),dp) 
B0m2 = MCha(i2)*Real(SA_B0(p2,MCha2(i2),Mhh2(i1)),dp) 
coupL1 = cplcUChaChahhL(gO1,i2,i1)
coupR1 = cplcUChaChahhR(gO1,i2,i1)
coupL2 =  Conjg(cplcUChaChahhL(gO2,i2,i1))
coupR2 =  Conjg(cplcUChaChahhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VP, Cha 
!------------------------ 
      Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -2._dp*Real(SA_B1(p2,MCha2(i2),0._dp)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MCha(i2)*Real(SA_B0(p2,MCha2(i2),0._dp)-0.5_dp*rMS,dp) 
coupL1 = cplcUChaChaVPL(gO1,i2)
coupR1 = cplcUChaChaVPR(gO1,i2)
coupL2 =  Conjg(cplcUChaChaVPL(gO2,i2))
coupR2 =  Conjg(cplcUChaChaVPR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! VZ, Cha 
!------------------------ 
      Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -2._dp*Real(SA_B1(p2,MCha2(i2),MVZ2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MCha(i2)*Real(SA_B0(p2,MCha2(i2),MVZ2)-0.5_dp*rMS,dp) 
coupL1 = cplcUChaChaVZL(gO1,i2)
coupR1 = cplcUChaChaVZR(gO1,i2)
coupL2 =  Conjg(cplcUChaChaVZL(gO2,i2))
coupR2 =  Conjg(cplcUChaChaVZR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Hpm, Chi 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -Real(SA_B1(p2,MChi2(i2),MHpm2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(p2,MChi2(i2),MHpm2(i1)),dp) 
coupL1 = cplcUChaChiHpmL(gO1,i2,i1)
coupR1 = cplcUChaChiHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplcUChaChiHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplcUChaChiHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VWm, Chi 
!------------------------ 
      Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -2._dp*Real(SA_B1(p2,MChi2(i2),MVWm2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MChi(i2)*Real(SA_B0(p2,MChi2(i2),MVWm2)-0.5_dp*rMS,dp) 
coupL1 = cplcUChaChiVWmL(gO1,i2)
coupR1 = cplcUChaChiVWmR(gO1,i2)
coupL2 =  Conjg(cplcUChaChiVWmL(gO2,i2))
coupR2 =  Conjg(cplcUChaChiVWmR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! conj[Su], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -Real(SA_B1(p2,MFd2(i2),MSu2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSu2(i1)),dp) 
coupL1 = cplcUChaFdcSuL(gO1,i2,i1)
coupR1 = cplcUChaFdcSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUChaFdcSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUChaFdcSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp* sumL
SigR = SigR +3._dp* sumR 
SigSL = SigSL +3._dp* sumSL 
SigSR = SigSR +3._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Fe 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -Real(SA_B1(p2,MFe2(i2),MSv2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(p2,MFe2(i2),MSv2(i1)),dp) 
coupL1 = cplcUChaFecSvL(gO1,i2,i1)
coupR1 = cplcUChaFecSvR(gO1,i2,i1)
coupL2 =  Conjg(cplcUChaFecSvL(gO2,i2,i1))
coupR2 =  Conjg(cplcUChaFecSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Se, Fv 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -Real(SA_B1(p2,0._dp,MSe2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(p2,0._dp,MSe2(i1)),dp) 
coupL1 = cplcUChaFvSeL(gO1,i2,i1)
coupR1 = cplcUChaFvSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUChaFvSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUChaFvSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 2
  Do gO2 = 1, 2
B1m2 = -Real(SA_B1(p2,MFu2(i1),MSd2(i2)),dp) 
B0m2 = MFu(i1)*Real(SA_B0(p2,MFu2(i1),MSd2(i2)),dp) 
coupL1 = cplcUChacFuSdL(gO1,i1,i2)
coupR1 = cplcUChacFuSdR(gO1,i1,i2)
coupL2 =  Conjg(cplcUChacFuSdL(gO2,i1,i2))
coupR2 =  Conjg(cplcUChacFuSdR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp* sumL
SigR = SigR +3._dp* sumR 
SigSL = SigSL +3._dp* sumSL 
SigSR = SigSR +3._dp* sumSR 
      End Do 
     End Do 
 SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopCha 
 
Subroutine OneLoopFe(Ye,vd,MFe,MFe2,MAh,MAh2,MSv,MSv2,MCha,MCha2,MSe,MSe2,            & 
& MChi,MChi2,MSu,MSu2,MFd,MFd2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MFu,              & 
& MFu2,MSd,MSd2,cplcUFeFeAhL,cplcUFeFeAhR,cplcUFeChaSvL,cplcUFeChaSvR,cplcUFeChiSeL,     & 
& cplcUFeChiSeR,cplcUFeFdcSuL,cplcUFeFdcSuR,cplcUFeFehhL,cplcUFeFehhR,cplcUFeFecSvL,     & 
& cplcUFeFecSvR,cplcUFeFeVPL,cplcUFeFeVPR,cplcUFeFeVZL,cplcUFeFeVZR,cplcUFeFvHpmL,       & 
& cplcUFeFvHpmR,cplcUFeFvSeL,cplcUFeFvSeR,cplcUFeFvVWmL,cplcUFeFvVWmR,cplcUFecFuSdL,     & 
& cplcUFecFuSdR,delta,MFe_1L,MFe2_1L,ZEL_1L,ZER_1L,ierr)

Implicit None 
Real(dp), Intent(in) :: MFe(3),MFe2(3),MAh(2),MAh2(2),MSv(3),MSv2(3),MCha(2),MCha2(2),MSe(6),MSe2(6),         & 
& MChi(4),MChi2(4),MSu(6),MSu2(6),MFd(3),MFd2(3),Mhh(2),Mhh2(2),MVZ,MVZ2,MHpm(2),        & 
& MHpm2(2),MVWm,MVWm2,MFu(3),MFu2(3),MSd(6),MSd2(6)

Real(dp), Intent(in) :: vd

Complex(dp), Intent(in) :: Ye(3,3)

Complex(dp), Intent(in) :: cplcUFeFeAhL(3,3,2),cplcUFeFeAhR(3,3,2),cplcUFeChaSvL(3,2,3),cplcUFeChaSvR(3,2,3),    & 
& cplcUFeChiSeL(3,4,6),cplcUFeChiSeR(3,4,6),cplcUFeFdcSuL(3,3,6),cplcUFeFdcSuR(3,3,6),   & 
& cplcUFeFehhL(3,3,2),cplcUFeFehhR(3,3,2),cplcUFeFecSvL(3,3,3),cplcUFeFecSvR(3,3,3),     & 
& cplcUFeFeVPL(3,3),cplcUFeFeVPR(3,3),cplcUFeFeVZL(3,3),cplcUFeFeVZR(3,3),               & 
& cplcUFeFvHpmL(3,3,2),cplcUFeFvHpmR(3,3,2),cplcUFeFvSeL(3,3,6),cplcUFeFvSeR(3,3,6),     & 
& cplcUFeFvVWmL(3,3),cplcUFeFvVWmR(3,3),cplcUFecFuSdL(3,3,6),cplcUFecFuSdR(3,3,6)

Complex(dp) :: mat1a(3,3), mat1(3,3) 
Integer , Intent(inout):: ierr 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(3), test_m2(3), p2 
Real(dp), Intent(out) :: MFe_1L(3),MFe2_1L(3) 
 Complex(dp), Intent(out) :: ZEL_1L(3,3), ZER_1L(3,3) 
 
 Real(dp) :: MFe_t(3),MFe2_t(3) 
 Complex(dp) :: ZEL_t(3,3), ZER_t(3,3), sigL(3,3), sigR(3,3), sigSL(3,3), sigSR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZEL2(3,3), ZER2(3,3) 
 
 Real(dp) :: ZEL1(3,3), ZER1(3,3), test(2) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopMFe'
 
mat1a(1,1) = 0._dp 
mat1a(1,1) = mat1a(1,1)+(vd*Ye(1,1))/sqrt(2._dp)
mat1a(1,2) = 0._dp 
mat1a(1,2) = mat1a(1,2)+(vd*Ye(2,1))/sqrt(2._dp)
mat1a(1,3) = 0._dp 
mat1a(1,3) = mat1a(1,3)+(vd*Ye(3,1))/sqrt(2._dp)
mat1a(2,1) = 0._dp 
mat1a(2,1) = mat1a(2,1)+(vd*Ye(1,2))/sqrt(2._dp)
mat1a(2,2) = 0._dp 
mat1a(2,2) = mat1a(2,2)+(vd*Ye(2,2))/sqrt(2._dp)
mat1a(2,3) = 0._dp 
mat1a(2,3) = mat1a(2,3)+(vd*Ye(3,2))/sqrt(2._dp)
mat1a(3,1) = 0._dp 
mat1a(3,1) = mat1a(3,1)+(vd*Ye(1,3))/sqrt(2._dp)
mat1a(3,2) = 0._dp 
mat1a(3,2) = mat1a(3,2)+(vd*Ye(2,3))/sqrt(2._dp)
mat1a(3,3) = 0._dp 
mat1a(3,3) = mat1a(3,3)+(vd*Ye(3,3))/sqrt(2._dp)

 
 !---------------------------------------- 
! Rotation matrix for p2=0 
!----------------------------------------- 
 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = 0._dp 
Call Sigma1LoopFe(p2,MFe,MFe2,MAh,MAh2,MSv,MSv2,MCha,MCha2,MSe,MSe2,MChi,             & 
& MChi2,MSu,MSu2,MFd,MFd2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MFu,MFu2,              & 
& MSd,MSd2,cplcUFeFeAhL,cplcUFeFeAhR,cplcUFeChaSvL,cplcUFeChaSvR,cplcUFeChiSeL,          & 
& cplcUFeChiSeR,cplcUFeFdcSuL,cplcUFeFdcSuR,cplcUFeFehhL,cplcUFeFehhR,cplcUFeFecSvL,     & 
& cplcUFeFecSvR,cplcUFeFeVPL,cplcUFeFeVPR,cplcUFeFeVZL,cplcUFeFeVZR,cplcUFeFvHpmL,       & 
& cplcUFeFvHpmR,cplcUFeFvSeL,cplcUFeFvSeR,cplcUFeFvVWmL,cplcUFeFvVWmR,cplcUFecFuSdL,     & 
& cplcUFecFuSdR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFe2_t,ZER1,ierr,test) 
ZER2 = ZER1 
Else 
Call EigenSystem(mat2,MFe2_t,ZER2,ierr,test) 
 End If 
 
ZEROS_0 = ZER2 
 mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFe2_t,ZEL1,ierr,test) 
 
 
ZEL2 = ZEL1 
Else 
Call EigenSystem(mat2,MFe2_t,ZEL2,ierr,test) 
 
 
End If 
ZEL2 = Conjg(ZEL2) 
ZELOS_0 = ZEL2 
 
!---------------------------------------- 
! Now, with momenta
!----------------------------------------- 
 
Do il=3,1,-1
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MFe2(il) 
Call Sigma1LoopFe(p2,MFe,MFe2,MAh,MAh2,MSv,MSv2,MCha,MCha2,MSe,MSe2,MChi,             & 
& MChi2,MSu,MSu2,MFd,MFd2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MFu,MFu2,              & 
& MSd,MSd2,cplcUFeFeAhL,cplcUFeFeAhR,cplcUFeChaSvL,cplcUFeChaSvR,cplcUFeChiSeL,          & 
& cplcUFeChiSeR,cplcUFeFdcSuL,cplcUFeFdcSuR,cplcUFeFehhL,cplcUFeFehhR,cplcUFeFecSvL,     & 
& cplcUFeFecSvR,cplcUFeFeVPL,cplcUFeFeVPR,cplcUFeFeVZL,cplcUFeFeVZR,cplcUFeFvHpmL,       & 
& cplcUFeFvHpmR,cplcUFeFvSeL,cplcUFeFvSeR,cplcUFeFvVWmL,cplcUFeFvVWmR,cplcUFecFuSdL,     & 
& cplcUFecFuSdR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFe2_t,ZER1,ierr,test) 
ZER2 = ZER1 
Else 
Call EigenSystem(mat2,MFe2_t,ZER2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
!---------------------------------------- 
! Redoing Calculation using redefined p2 
!----------------------------------------- 
 
i_count = 0 
p2_loop: Do  
i_count = i_count + 1 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MFe2_t(iL)
Call Sigma1LoopFe(p2,MFe,MFe2,MAh,MAh2,MSv,MSv2,MCha,MCha2,MSe,MSe2,MChi,             & 
& MChi2,MSu,MSu2,MFd,MFd2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MFu,MFu2,              & 
& MSd,MSd2,cplcUFeFeAhL,cplcUFeFeAhR,cplcUFeChaSvL,cplcUFeChaSvR,cplcUFeChiSeL,          & 
& cplcUFeChiSeR,cplcUFeFdcSuL,cplcUFeFdcSuR,cplcUFeFehhL,cplcUFeFehhR,cplcUFeFecSvL,     & 
& cplcUFeFecSvR,cplcUFeFeVPL,cplcUFeFeVPR,cplcUFeFeVZL,cplcUFeFeVZR,cplcUFeFvHpmL,       & 
& cplcUFeFvHpmR,cplcUFeFvSeL,cplcUFeFvSeR,cplcUFeFvVWmL,cplcUFeFvVWmR,cplcUFecFuSdL,     & 
& cplcUFecFuSdR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFe2_t,ZER1,ierr,test) 
ZER2 = ZER1 
Else 
Call EigenSystem(mat2,MFe2_t,ZER2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
MFe2_1L(il) = MFe2_t(il) 
MFe_1L(il) = Sqrt(MFe2_1L(il)) 
 
If (p2.Ne.0._dp) Then 
  test(1) = Abs(MFe2_1L(il)-p2)/p2
Else 
  test(2) = Abs(MFe2_1L(il))
End If 
If (Abs(MFe2_1L(il)).lt.1.0E-30_dp) Exit p2_loop 
If (test(1).lt.0.1_dp*delta) Exit p2_loop 
If(i_count.gt.30) then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Exit p2_loop 
End if
End Do p2_loop 
mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFe2_t,ZEL1,ierr,test) 
 
 
ZEL2 = ZEL1 
Else 
Call EigenSystem(mat2,MFe2_t,ZEL2,ierr,test) 
 
 
End If 
ZEL2 = Conjg(ZEL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZEL2),mat1),Transpose( Conjg(ZER2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZER2(i1,:) = phaseM *ZER2(i1,:) 
 End if 
End Do 
 
ZELOS_p2(il,:) = ZEL2(il,:) 
 ZEROS_p2(il,:) = ZER2(il,:) 
 ZEL_1L = ZEL2 
 ZER_1L = ZER2 
 End Do  
 
Iname = Iname -1 
End Subroutine OneLoopFe
 
 
Subroutine Sigma1LoopFe(p2,MFe,MFe2,MAh,MAh2,MSv,MSv2,MCha,MCha2,MSe,MSe2,            & 
& MChi,MChi2,MSu,MSu2,MFd,MFd2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MFu,              & 
& MFu2,MSd,MSd2,cplcUFeFeAhL,cplcUFeFeAhR,cplcUFeChaSvL,cplcUFeChaSvR,cplcUFeChiSeL,     & 
& cplcUFeChiSeR,cplcUFeFdcSuL,cplcUFeFdcSuR,cplcUFeFehhL,cplcUFeFehhR,cplcUFeFecSvL,     & 
& cplcUFeFecSvR,cplcUFeFeVPL,cplcUFeFeVPR,cplcUFeFeVZL,cplcUFeFeVZR,cplcUFeFvHpmL,       & 
& cplcUFeFvHpmR,cplcUFeFvSeL,cplcUFeFvSeR,cplcUFeFvVWmL,cplcUFeFvVWmR,cplcUFecFuSdL,     & 
& cplcUFecFuSdR,sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MFe(3),MFe2(3),MAh(2),MAh2(2),MSv(3),MSv2(3),MCha(2),MCha2(2),MSe(6),MSe2(6),         & 
& MChi(4),MChi2(4),MSu(6),MSu2(6),MFd(3),MFd2(3),Mhh(2),Mhh2(2),MVZ,MVZ2,MHpm(2),        & 
& MHpm2(2),MVWm,MVWm2,MFu(3),MFu2(3),MSd(6),MSd2(6)

Complex(dp), Intent(in) :: cplcUFeFeAhL(3,3,2),cplcUFeFeAhR(3,3,2),cplcUFeChaSvL(3,2,3),cplcUFeChaSvR(3,2,3),    & 
& cplcUFeChiSeL(3,4,6),cplcUFeChiSeR(3,4,6),cplcUFeFdcSuL(3,3,6),cplcUFeFdcSuR(3,3,6),   & 
& cplcUFeFehhL(3,3,2),cplcUFeFehhR(3,3,2),cplcUFeFecSvL(3,3,3),cplcUFeFecSvR(3,3,3),     & 
& cplcUFeFeVPL(3,3),cplcUFeFeVPR(3,3),cplcUFeFeVZL(3,3),cplcUFeFeVZR(3,3),               & 
& cplcUFeFvHpmL(3,3,2),cplcUFeFvHpmR(3,3,2),cplcUFeFvSeL(3,3,6),cplcUFeFvSeR(3,3,6),     & 
& cplcUFeFvVWmL(3,3),cplcUFeFvVWmR(3,3),cplcUFecFuSdL(3,3,6),cplcUFecFuSdR(3,3,6)

Complex(dp), Intent(out) :: SigL(3,3),SigR(3,3), SigSL(3,3), SigSR(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(3,3), sumR(3,3), sumSL(3,3), sumSR(3,3) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
 SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Fe, Ah 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFe2(i1),MAh2(i2)),dp) 
B0m2 = MFe(i1)*Real(SA_B0(p2,MFe2(i1),MAh2(i2)),dp) 
coupL1 = cplcUFeFeAhL(gO1,i1,i2)
coupR1 = cplcUFeFeAhR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFeFeAhL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFeFeAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sv, Cha 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MCha2(i2),MSv2(i1)),dp) 
B0m2 = MCha(i2)*Real(SA_B0(p2,MCha2(i2),MSv2(i1)),dp) 
coupL1 = cplcUFeChaSvL(gO1,i2,i1)
coupR1 = cplcUFeChaSvR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeChaSvL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeChaSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Se, Chi 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MChi2(i2),MSe2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSe2(i1)),dp) 
coupL1 = cplcUFeChiSeL(gO1,i2,i1)
coupR1 = cplcUFeChiSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeChiSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeChiSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFd2(i2),MSu2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSu2(i1)),dp) 
coupL1 = cplcUFeFdcSuL(gO1,i2,i1)
coupR1 = cplcUFeFdcSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFdcSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFdcSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp* sumL
SigR = SigR +3._dp* sumR 
SigSL = SigSL +3._dp* sumSL 
SigSR = SigSR +3._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! hh, Fe 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFe2(i2),Mhh2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(p2,MFe2(i2),Mhh2(i1)),dp) 
coupL1 = cplcUFeFehhL(gO1,i2,i1)
coupR1 = cplcUFeFehhR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFehhL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFehhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Fe 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFe2(i2),MSv2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(p2,MFe2(i2),MSv2(i1)),dp) 
coupL1 = cplcUFeFecSvL(gO1,i2,i1)
coupR1 = cplcUFeFecSvR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFecSvL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFecSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VP, Fe 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i2),0._dp)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFe(i2)*Real(SA_B0(p2,MFe2(i2),0._dp)-0.5_dp*rMS,dp) 
coupL1 = cplcUFeFeVPL(gO1,i2)
coupR1 = cplcUFeFeVPR(gO1,i2)
coupL2 =  Conjg(cplcUFeFeVPL(gO2,i2))
coupR2 =  Conjg(cplcUFeFeVPR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! VZ, Fe 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i2),MVZ2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFe(i2)*Real(SA_B0(p2,MFe2(i2),MVZ2)-0.5_dp*rMS,dp) 
coupL1 = cplcUFeFeVZL(gO1,i2)
coupR1 = cplcUFeFeVZR(gO1,i2)
coupL2 =  Conjg(cplcUFeFeVZL(gO2,i2))
coupR2 =  Conjg(cplcUFeFeVZR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Hpm, Fv 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,0._dp,MHpm2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(p2,0._dp,MHpm2(i1)),dp) 
coupL1 = cplcUFeFvHpmL(gO1,i2,i1)
coupR1 = cplcUFeFvHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFvHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFvHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Se, Fv 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,0._dp,MSe2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(p2,0._dp,MSe2(i1)),dp) 
coupL1 = cplcUFeFvSeL(gO1,i2,i1)
coupR1 = cplcUFeFvSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFvSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFvSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VWm, Fv 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,0._dp,MVWm2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*0.*Real(SA_B0(p2,0._dp,MVWm2)-0.5_dp*rMS,dp) 
coupL1 = cplcUFeFvVWmL(gO1,i2)
coupR1 = cplcUFeFvVWmR(gO1,i2)
coupL2 =  Conjg(cplcUFeFvVWmL(gO2,i2))
coupR2 =  Conjg(cplcUFeFvVWmR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! bar[Fu], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFu2(i1),MSd2(i2)),dp) 
B0m2 = MFu(i1)*Real(SA_B0(p2,MFu2(i1),MSd2(i2)),dp) 
coupL1 = cplcUFecFuSdL(gO1,i1,i2)
coupR1 = cplcUFecFuSdR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFecFuSdL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFecFuSdR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp* sumL
SigR = SigR +3._dp* sumR 
SigSL = SigSL +3._dp* sumSL 
SigSR = SigSR +3._dp* sumSR 
      End Do 
     End Do 
 SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopFe 
 
Subroutine OneLoopFd(Yd,vd,MFd,MFd2,MAh,MAh2,MSu,MSu2,MCha,MCha2,MSd,MSd2,            & 
& MChi,MChi2,Mhh,Mhh2,MSv,MSv2,MVZ,MVZ2,MFe,MFe2,MHpm,MHpm2,MFu,MFu2,MSe,MSe2,           & 
& MVWm,MVWm2,MGlu,MGlu2,cplcUFdFdAhL,cplcUFdFdAhR,cplcUFdChaSuL,cplcUFdChaSuR,           & 
& cplcUFdChiSdL,cplcUFdChiSdR,cplcUFdFdhhL,cplcUFdFdhhR,cplcUFdFdcSvL,cplcUFdFdcSvR,     & 
& cplcUFdFdVGL,cplcUFdFdVGR,cplcUFdFdVPL,cplcUFdFdVPR,cplcUFdFdVZL,cplcUFdFdVZR,         & 
& cplcUFdFeSuL,cplcUFdFeSuR,cplcUFdFuHpmL,cplcUFdFuHpmR,cplcUFdFuSeL,cplcUFdFuSeR,       & 
& cplcUFdFuVWmL,cplcUFdFuVWmR,cplcUFdFvSdL,cplcUFdFvSdR,cplcUFdGluSdL,cplcUFdGluSdR,     & 
& delta,MFd_1L,MFd2_1L,ZDL_1L,ZDR_1L,ierr)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MAh(2),MAh2(2),MSu(6),MSu2(6),MCha(2),MCha2(2),MSd(6),MSd2(6),         & 
& MChi(4),MChi2(4),Mhh(2),Mhh2(2),MSv(3),MSv2(3),MVZ,MVZ2,MFe(3),MFe2(3),MHpm(2),        & 
& MHpm2(2),MFu(3),MFu2(3),MSe(6),MSe2(6),MVWm,MVWm2,MGlu,MGlu2

Real(dp), Intent(in) :: vd

Complex(dp), Intent(in) :: Yd(3,3)

Complex(dp), Intent(in) :: cplcUFdFdAhL(3,3,2),cplcUFdFdAhR(3,3,2),cplcUFdChaSuL(3,2,6),cplcUFdChaSuR(3,2,6),    & 
& cplcUFdChiSdL(3,4,6),cplcUFdChiSdR(3,4,6),cplcUFdFdhhL(3,3,2),cplcUFdFdhhR(3,3,2),     & 
& cplcUFdFdcSvL(3,3,3),cplcUFdFdcSvR(3,3,3),cplcUFdFdVGL(3,3),cplcUFdFdVGR(3,3),         & 
& cplcUFdFdVPL(3,3),cplcUFdFdVPR(3,3),cplcUFdFdVZL(3,3),cplcUFdFdVZR(3,3),               & 
& cplcUFdFeSuL(3,3,6),cplcUFdFeSuR(3,3,6),cplcUFdFuHpmL(3,3,2),cplcUFdFuHpmR(3,3,2),     & 
& cplcUFdFuSeL(3,3,6),cplcUFdFuSeR(3,3,6),cplcUFdFuVWmL(3,3),cplcUFdFuVWmR(3,3),         & 
& cplcUFdFvSdL(3,3,6),cplcUFdFvSdR(3,3,6),cplcUFdGluSdL(3,6),cplcUFdGluSdR(3,6)

Complex(dp) :: mat1a(3,3), mat1(3,3) 
Integer , Intent(inout):: ierr 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(3), test_m2(3), p2 
Real(dp), Intent(out) :: MFd_1L(3),MFd2_1L(3) 
 Complex(dp), Intent(out) :: ZDL_1L(3,3), ZDR_1L(3,3) 
 
 Real(dp) :: MFd_t(3),MFd2_t(3) 
 Complex(dp) :: ZDL_t(3,3), ZDR_t(3,3), sigL(3,3), sigR(3,3), sigSL(3,3), sigSR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZDL2(3,3), ZDR2(3,3) 
 
 Real(dp) :: ZDL1(3,3), ZDR1(3,3), test(2) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopMFd'
 
mat1a(1,1) = 0._dp 
mat1a(1,1) = mat1a(1,1)+(vd*Yd(1,1))/sqrt(2._dp)
mat1a(1,2) = 0._dp 
mat1a(1,2) = mat1a(1,2)+(vd*Yd(2,1))/sqrt(2._dp)
mat1a(1,3) = 0._dp 
mat1a(1,3) = mat1a(1,3)+(vd*Yd(3,1))/sqrt(2._dp)
mat1a(2,1) = 0._dp 
mat1a(2,1) = mat1a(2,1)+(vd*Yd(1,2))/sqrt(2._dp)
mat1a(2,2) = 0._dp 
mat1a(2,2) = mat1a(2,2)+(vd*Yd(2,2))/sqrt(2._dp)
mat1a(2,3) = 0._dp 
mat1a(2,3) = mat1a(2,3)+(vd*Yd(3,2))/sqrt(2._dp)
mat1a(3,1) = 0._dp 
mat1a(3,1) = mat1a(3,1)+(vd*Yd(1,3))/sqrt(2._dp)
mat1a(3,2) = 0._dp 
mat1a(3,2) = mat1a(3,2)+(vd*Yd(2,3))/sqrt(2._dp)
mat1a(3,3) = 0._dp 
mat1a(3,3) = mat1a(3,3)+(vd*Yd(3,3))/sqrt(2._dp)

 
 !---------------------------------------- 
! Rotation matrix for p2=0 
!----------------------------------------- 
 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = 0._dp 
Call Sigma1LoopFd(p2,MFd,MFd2,MAh,MAh2,MSu,MSu2,MCha,MCha2,MSd,MSd2,MChi,             & 
& MChi2,Mhh,Mhh2,MSv,MSv2,MVZ,MVZ2,MFe,MFe2,MHpm,MHpm2,MFu,MFu2,MSe,MSe2,MVWm,           & 
& MVWm2,MGlu,MGlu2,cplcUFdFdAhL,cplcUFdFdAhR,cplcUFdChaSuL,cplcUFdChaSuR,cplcUFdChiSdL,  & 
& cplcUFdChiSdR,cplcUFdFdhhL,cplcUFdFdhhR,cplcUFdFdcSvL,cplcUFdFdcSvR,cplcUFdFdVGL,      & 
& cplcUFdFdVGR,cplcUFdFdVPL,cplcUFdFdVPR,cplcUFdFdVZL,cplcUFdFdVZR,cplcUFdFeSuL,         & 
& cplcUFdFeSuR,cplcUFdFuHpmL,cplcUFdFuHpmR,cplcUFdFuSeL,cplcUFdFuSeR,cplcUFdFuVWmL,      & 
& cplcUFdFuVWmR,cplcUFdFvSdL,cplcUFdFvSdR,cplcUFdGluSdL,cplcUFdGluSdR,sigL,              & 
& sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2_t,ZDR1,ierr,test) 
ZDR2 = ZDR1 
Else 
Call EigenSystem(mat2,MFd2_t,ZDR2,ierr,test) 
 End If 
 
ZDROS_0 = ZDR2 
 mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2_t,ZDL1,ierr,test) 
 
 
ZDL2 = ZDL1 
Else 
Call EigenSystem(mat2,MFd2_t,ZDL2,ierr,test) 
 
 
End If 
ZDL2 = Conjg(ZDL2) 
ZDLOS_0 = ZDL2 
 
!---------------------------------------- 
! Now, with momenta
!----------------------------------------- 
 
Do il=3,1,-1
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MFd2(il) 
Call Sigma1LoopFd(p2,MFd,MFd2,MAh,MAh2,MSu,MSu2,MCha,MCha2,MSd,MSd2,MChi,             & 
& MChi2,Mhh,Mhh2,MSv,MSv2,MVZ,MVZ2,MFe,MFe2,MHpm,MHpm2,MFu,MFu2,MSe,MSe2,MVWm,           & 
& MVWm2,MGlu,MGlu2,cplcUFdFdAhL,cplcUFdFdAhR,cplcUFdChaSuL,cplcUFdChaSuR,cplcUFdChiSdL,  & 
& cplcUFdChiSdR,cplcUFdFdhhL,cplcUFdFdhhR,cplcUFdFdcSvL,cplcUFdFdcSvR,cplcUFdFdVGL,      & 
& cplcUFdFdVGR,cplcUFdFdVPL,cplcUFdFdVPR,cplcUFdFdVZL,cplcUFdFdVZR,cplcUFdFeSuL,         & 
& cplcUFdFeSuR,cplcUFdFuHpmL,cplcUFdFuHpmR,cplcUFdFuSeL,cplcUFdFuSeR,cplcUFdFuVWmL,      & 
& cplcUFdFuVWmR,cplcUFdFvSdL,cplcUFdFvSdR,cplcUFdGluSdL,cplcUFdGluSdR,sigL,              & 
& sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2_t,ZDR1,ierr,test) 
ZDR2 = ZDR1 
Else 
Call EigenSystem(mat2,MFd2_t,ZDR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
!---------------------------------------- 
! Redoing Calculation using redefined p2 
!----------------------------------------- 
 
i_count = 0 
p2_loop: Do  
i_count = i_count + 1 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MFd2_t(iL)
Call Sigma1LoopFd(p2,MFd,MFd2,MAh,MAh2,MSu,MSu2,MCha,MCha2,MSd,MSd2,MChi,             & 
& MChi2,Mhh,Mhh2,MSv,MSv2,MVZ,MVZ2,MFe,MFe2,MHpm,MHpm2,MFu,MFu2,MSe,MSe2,MVWm,           & 
& MVWm2,MGlu,MGlu2,cplcUFdFdAhL,cplcUFdFdAhR,cplcUFdChaSuL,cplcUFdChaSuR,cplcUFdChiSdL,  & 
& cplcUFdChiSdR,cplcUFdFdhhL,cplcUFdFdhhR,cplcUFdFdcSvL,cplcUFdFdcSvR,cplcUFdFdVGL,      & 
& cplcUFdFdVGR,cplcUFdFdVPL,cplcUFdFdVPR,cplcUFdFdVZL,cplcUFdFdVZR,cplcUFdFeSuL,         & 
& cplcUFdFeSuR,cplcUFdFuHpmL,cplcUFdFuHpmR,cplcUFdFuSeL,cplcUFdFuSeR,cplcUFdFuVWmL,      & 
& cplcUFdFuVWmR,cplcUFdFvSdL,cplcUFdFvSdR,cplcUFdGluSdL,cplcUFdGluSdR,sigL,              & 
& sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2_t,ZDR1,ierr,test) 
ZDR2 = ZDR1 
Else 
Call EigenSystem(mat2,MFd2_t,ZDR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
MFd2_1L(il) = MFd2_t(il) 
MFd_1L(il) = Sqrt(MFd2_1L(il)) 
 
If (p2.Ne.0._dp) Then 
  test(1) = Abs(MFd2_1L(il)-p2)/p2
Else 
  test(2) = Abs(MFd2_1L(il))
End If 
If (Abs(MFd2_1L(il)).lt.1.0E-30_dp) Exit p2_loop 
If (test(1).lt.0.1_dp*delta) Exit p2_loop 
If(i_count.gt.30) then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Exit p2_loop 
End if
End Do p2_loop 
mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2_t,ZDL1,ierr,test) 
 
 
ZDL2 = ZDL1 
Else 
Call EigenSystem(mat2,MFd2_t,ZDL2,ierr,test) 
 
 
End If 
ZDL2 = Conjg(ZDL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZDL2),mat1),Transpose( Conjg(ZDR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZDR2(i1,:) = phaseM *ZDR2(i1,:) 
 End if 
End Do 
 
ZDLOS_p2(il,:) = ZDL2(il,:) 
 ZDROS_p2(il,:) = ZDR2(il,:) 
 ZDL_1L = ZDL2 
 ZDR_1L = ZDR2 
 End Do  
 
Iname = Iname -1 
End Subroutine OneLoopFd
 
 
Subroutine Sigma1LoopFd(p2,MFd,MFd2,MAh,MAh2,MSu,MSu2,MCha,MCha2,MSd,MSd2,            & 
& MChi,MChi2,Mhh,Mhh2,MSv,MSv2,MVZ,MVZ2,MFe,MFe2,MHpm,MHpm2,MFu,MFu2,MSe,MSe2,           & 
& MVWm,MVWm2,MGlu,MGlu2,cplcUFdFdAhL,cplcUFdFdAhR,cplcUFdChaSuL,cplcUFdChaSuR,           & 
& cplcUFdChiSdL,cplcUFdChiSdR,cplcUFdFdhhL,cplcUFdFdhhR,cplcUFdFdcSvL,cplcUFdFdcSvR,     & 
& cplcUFdFdVGL,cplcUFdFdVGR,cplcUFdFdVPL,cplcUFdFdVPR,cplcUFdFdVZL,cplcUFdFdVZR,         & 
& cplcUFdFeSuL,cplcUFdFeSuR,cplcUFdFuHpmL,cplcUFdFuHpmR,cplcUFdFuSeL,cplcUFdFuSeR,       & 
& cplcUFdFuVWmL,cplcUFdFuVWmR,cplcUFdFvSdL,cplcUFdFvSdR,cplcUFdGluSdL,cplcUFdGluSdR,     & 
& sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MAh(2),MAh2(2),MSu(6),MSu2(6),MCha(2),MCha2(2),MSd(6),MSd2(6),         & 
& MChi(4),MChi2(4),Mhh(2),Mhh2(2),MSv(3),MSv2(3),MVZ,MVZ2,MFe(3),MFe2(3),MHpm(2),        & 
& MHpm2(2),MFu(3),MFu2(3),MSe(6),MSe2(6),MVWm,MVWm2,MGlu,MGlu2

Complex(dp), Intent(in) :: cplcUFdFdAhL(3,3,2),cplcUFdFdAhR(3,3,2),cplcUFdChaSuL(3,2,6),cplcUFdChaSuR(3,2,6),    & 
& cplcUFdChiSdL(3,4,6),cplcUFdChiSdR(3,4,6),cplcUFdFdhhL(3,3,2),cplcUFdFdhhR(3,3,2),     & 
& cplcUFdFdcSvL(3,3,3),cplcUFdFdcSvR(3,3,3),cplcUFdFdVGL(3,3),cplcUFdFdVGR(3,3),         & 
& cplcUFdFdVPL(3,3),cplcUFdFdVPR(3,3),cplcUFdFdVZL(3,3),cplcUFdFdVZR(3,3),               & 
& cplcUFdFeSuL(3,3,6),cplcUFdFeSuR(3,3,6),cplcUFdFuHpmL(3,3,2),cplcUFdFuHpmR(3,3,2),     & 
& cplcUFdFuSeL(3,3,6),cplcUFdFuSeR(3,3,6),cplcUFdFuVWmL(3,3),cplcUFdFuVWmR(3,3),         & 
& cplcUFdFvSdL(3,3,6),cplcUFdFvSdR(3,3,6),cplcUFdGluSdL(3,6),cplcUFdGluSdR(3,6)

Complex(dp), Intent(out) :: SigL(3,3),SigR(3,3), SigSL(3,3), SigSR(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(3,3), sumR(3,3), sumSL(3,3), sumSR(3,3) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
 SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Fd, Ah 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFd2(i1),MAh2(i2)),dp) 
B0m2 = MFd(i1)*Real(SA_B0(p2,MFd2(i1),MAh2(i2)),dp) 
coupL1 = cplcUFdFdAhL(gO1,i1,i2)
coupR1 = cplcUFdFdAhR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFdFdAhL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFdFdAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Su, Cha 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MCha2(i2),MSu2(i1)),dp) 
B0m2 = MCha(i2)*Real(SA_B0(p2,MCha2(i2),MSu2(i1)),dp) 
coupL1 = cplcUFdChaSuL(gO1,i2,i1)
coupR1 = cplcUFdChaSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdChaSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdChaSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sd, Chi 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MChi2(i2),MSd2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSd2(i1)),dp) 
coupL1 = cplcUFdChiSdL(gO1,i2,i1)
coupR1 = cplcUFdChiSdR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdChiSdL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdChiSdR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! hh, Fd 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFd2(i2),Mhh2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),Mhh2(i1)),dp) 
coupL1 = cplcUFdFdhhL(gO1,i2,i1)
coupR1 = cplcUFdFdhhR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFdhhL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFdhhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Fd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFd2(i2),MSv2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSv2(i1)),dp) 
coupL1 = cplcUFdFdcSvL(gO1,i2,i1)
coupR1 = cplcUFdFdcSvR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFdcSvL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFdcSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VG, Fd 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),0._dp)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),0._dp)-0.5_dp*rMS,dp) 
coupL1 = cplcUFdFdVGL(gO1,i2)
coupR1 = cplcUFdFdVGR(gO1,i2)
coupL2 =  Conjg(cplcUFdFdVGL(gO2,i2))
coupR2 =  Conjg(cplcUFdFdVGR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +4._dp/3._dp* sumL
SigR = SigR +4._dp/3._dp* sumR 
SigSL = SigSL +4._dp/3._dp* sumSL 
SigSR = SigSR +4._dp/3._dp* sumSR 
    End Do 
 !------------------------ 
! VP, Fd 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),0._dp)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),0._dp)-0.5_dp*rMS,dp) 
coupL1 = cplcUFdFdVPL(gO1,i2)
coupR1 = cplcUFdFdVPR(gO1,i2)
coupL2 =  Conjg(cplcUFdFdVPL(gO2,i2))
coupR2 =  Conjg(cplcUFdFdVPR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! VZ, Fd 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),MVZ2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),MVZ2)-0.5_dp*rMS,dp) 
coupL1 = cplcUFdFdVZL(gO1,i2)
coupR1 = cplcUFdFdVZR(gO1,i2)
coupL2 =  Conjg(cplcUFdFdVZL(gO2,i2))
coupR2 =  Conjg(cplcUFdFdVZR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Su, Fe 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFe2(i2),MSu2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(p2,MFe2(i2),MSu2(i1)),dp) 
coupL1 = cplcUFdFeSuL(gO1,i2,i1)
coupR1 = cplcUFdFeSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFeSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFeSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Hpm, Fu 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFu2(i2),MHpm2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(p2,MFu2(i2),MHpm2(i1)),dp) 
coupL1 = cplcUFdFuHpmL(gO1,i2,i1)
coupR1 = cplcUFdFuHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFuHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFuHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Se, Fu 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFu2(i2),MSe2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(p2,MFu2(i2),MSe2(i1)),dp) 
coupL1 = cplcUFdFuSeL(gO1,i2,i1)
coupR1 = cplcUFdFuSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFuSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFuSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VWm, Fu 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),MVWm2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),MVWm2)-0.5_dp*rMS,dp) 
coupL1 = cplcUFdFuVWmL(gO1,i2)
coupR1 = cplcUFdFuVWmR(gO1,i2)
coupL2 =  Conjg(cplcUFdFuVWmL(gO2,i2))
coupR2 =  Conjg(cplcUFdFuVWmR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Sd, Fv 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,0._dp,MSd2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(p2,0._dp,MSd2(i1)),dp) 
coupL1 = cplcUFdFvSdL(gO1,i2,i1)
coupR1 = cplcUFdFvSdR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFvSdL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFvSdR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sd, Glu 
!------------------------ 
    Do i1 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MGlu2,MSd2(i1)),dp) 
B0m2 = MGlu*Real(SA_B0(p2,MGlu2,MSd2(i1)),dp) 
coupL1 = cplcUFdGluSdL(gO1,i1)
coupR1 = cplcUFdGluSdR(gO1,i1)
coupL2 =  Conjg(cplcUFdGluSdL(gO2,i1))
coupR2 =  Conjg(cplcUFdGluSdR(gO2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +4._dp/3._dp* sumL
SigR = SigR +4._dp/3._dp* sumR 
SigSL = SigSL +4._dp/3._dp* sumSL 
SigSR = SigSR +4._dp/3._dp* sumSR 
      End Do 
 SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopFd 
 
Subroutine OneLoopFu(Yu,vu,MFu,MFu2,MAh,MAh2,MSu,MSu2,MChi,MChi2,MHpm,MHpm2,          & 
& MFd,MFd2,MSe,MSe2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MGlu,MGlu2,MCha,MCha2,MSd,              & 
& MSd2,MFe,MFe2,cplcUFuFuAhL,cplcUFuFuAhR,cplcUFuChiSuL,cplcUFuChiSuR,cplcUFuFdcHpmL,    & 
& cplcUFuFdcHpmR,cplcUFuFdcSeL,cplcUFuFdcSeR,cplcUFuFdcVWmL,cplcUFuFdcVWmR,              & 
& cplcUFuFuhhL,cplcUFuFuhhR,cplcUFuFuVGL,cplcUFuFuVGR,cplcUFuFuVPL,cplcUFuFuVPR,         & 
& cplcUFuFuVZL,cplcUFuFuVZR,cplcUFuGluSuL,cplcUFuGluSuR,cplcChacUFuSdL,cplcChacUFuSdR,   & 
& cplcFecUFuSdL,cplcFecUFuSdR,delta,MFu_1L,MFu2_1L,ZUL_1L,ZUR_1L,ierr)

Implicit None 
Real(dp), Intent(in) :: MFu(3),MFu2(3),MAh(2),MAh2(2),MSu(6),MSu2(6),MChi(4),MChi2(4),MHpm(2),MHpm2(2),       & 
& MFd(3),MFd2(3),MSe(6),MSe2(6),MVWm,MVWm2,Mhh(2),Mhh2(2),MVZ,MVZ2,MGlu,MGlu2,           & 
& MCha(2),MCha2(2),MSd(6),MSd2(6),MFe(3),MFe2(3)

Real(dp), Intent(in) :: vu

Complex(dp), Intent(in) :: Yu(3,3)

Complex(dp), Intent(in) :: cplcUFuFuAhL(3,3,2),cplcUFuFuAhR(3,3,2),cplcUFuChiSuL(3,4,6),cplcUFuChiSuR(3,4,6),    & 
& cplcUFuFdcHpmL(3,3,2),cplcUFuFdcHpmR(3,3,2),cplcUFuFdcSeL(3,3,6),cplcUFuFdcSeR(3,3,6), & 
& cplcUFuFdcVWmL(3,3),cplcUFuFdcVWmR(3,3),cplcUFuFuhhL(3,3,2),cplcUFuFuhhR(3,3,2),       & 
& cplcUFuFuVGL(3,3),cplcUFuFuVGR(3,3),cplcUFuFuVPL(3,3),cplcUFuFuVPR(3,3),               & 
& cplcUFuFuVZL(3,3),cplcUFuFuVZR(3,3),cplcUFuGluSuL(3,6),cplcUFuGluSuR(3,6),             & 
& cplcChacUFuSdL(2,3,6),cplcChacUFuSdR(2,3,6),cplcFecUFuSdL(3,3,6),cplcFecUFuSdR(3,3,6)

Complex(dp) :: mat1a(3,3), mat1(3,3) 
Integer , Intent(inout):: ierr 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi2(3), test_m2(3), p2 
Real(dp), Intent(out) :: MFu_1L(3),MFu2_1L(3) 
 Complex(dp), Intent(out) :: ZUL_1L(3,3), ZUR_1L(3,3) 
 
 Real(dp) :: MFu_t(3),MFu2_t(3) 
 Complex(dp) :: ZUL_t(3,3), ZUR_t(3,3), sigL(3,3), sigR(3,3), sigSL(3,3), sigSR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZUL2(3,3), ZUR2(3,3) 
 
 Real(dp) :: ZUL1(3,3), ZUR1(3,3), test(2) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopMFu'
 
mat1a(1,1) = 0._dp 
mat1a(1,1) = mat1a(1,1)+(vu*Yu(1,1))/sqrt(2._dp)
mat1a(1,2) = 0._dp 
mat1a(1,2) = mat1a(1,2)+(vu*Yu(2,1))/sqrt(2._dp)
mat1a(1,3) = 0._dp 
mat1a(1,3) = mat1a(1,3)+(vu*Yu(3,1))/sqrt(2._dp)
mat1a(2,1) = 0._dp 
mat1a(2,1) = mat1a(2,1)+(vu*Yu(1,2))/sqrt(2._dp)
mat1a(2,2) = 0._dp 
mat1a(2,2) = mat1a(2,2)+(vu*Yu(2,2))/sqrt(2._dp)
mat1a(2,3) = 0._dp 
mat1a(2,3) = mat1a(2,3)+(vu*Yu(3,2))/sqrt(2._dp)
mat1a(3,1) = 0._dp 
mat1a(3,1) = mat1a(3,1)+(vu*Yu(1,3))/sqrt(2._dp)
mat1a(3,2) = 0._dp 
mat1a(3,2) = mat1a(3,2)+(vu*Yu(2,3))/sqrt(2._dp)
mat1a(3,3) = 0._dp 
mat1a(3,3) = mat1a(3,3)+(vu*Yu(3,3))/sqrt(2._dp)

 
 !---------------------------------------- 
! Rotation matrix for p2=0 
!----------------------------------------- 
 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = 0._dp 
Call Sigma1LoopFu(p2,MFu,MFu2,MAh,MAh2,MSu,MSu2,MChi,MChi2,MHpm,MHpm2,MFd,            & 
& MFd2,MSe,MSe2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MGlu,MGlu2,MCha,MCha2,MSd,MSd2,             & 
& MFe,MFe2,cplcUFuFuAhL,cplcUFuFuAhR,cplcUFuChiSuL,cplcUFuChiSuR,cplcUFuFdcHpmL,         & 
& cplcUFuFdcHpmR,cplcUFuFdcSeL,cplcUFuFdcSeR,cplcUFuFdcVWmL,cplcUFuFdcVWmR,              & 
& cplcUFuFuhhL,cplcUFuFuhhR,cplcUFuFuVGL,cplcUFuFuVGR,cplcUFuFuVPL,cplcUFuFuVPR,         & 
& cplcUFuFuVZL,cplcUFuFuVZR,cplcUFuGluSuL,cplcUFuGluSuR,cplcChacUFuSdL,cplcChacUFuSdR,   & 
& cplcFecUFuSdL,cplcFecUFuSdR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2_t,ZUR1,ierr,test) 
ZUR2 = ZUR1 
Else 
Call EigenSystem(mat2,MFu2_t,ZUR2,ierr,test) 
 End If 
 
ZUROS_0 = ZUR2 
 mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2_t,ZUL1,ierr,test) 
 
 
ZUL2 = ZUL1 
Else 
Call EigenSystem(mat2,MFu2_t,ZUL2,ierr,test) 
 
 
End If 
ZUL2 = Conjg(ZUL2) 
ZULOS_0 = ZUL2 
 
!---------------------------------------- 
! Now, with momenta
!----------------------------------------- 
 
Do il=3,1,-1
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MFu2(il) 
Call Sigma1LoopFu(p2,MFu,MFu2,MAh,MAh2,MSu,MSu2,MChi,MChi2,MHpm,MHpm2,MFd,            & 
& MFd2,MSe,MSe2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MGlu,MGlu2,MCha,MCha2,MSd,MSd2,             & 
& MFe,MFe2,cplcUFuFuAhL,cplcUFuFuAhR,cplcUFuChiSuL,cplcUFuChiSuR,cplcUFuFdcHpmL,         & 
& cplcUFuFdcHpmR,cplcUFuFdcSeL,cplcUFuFdcSeR,cplcUFuFdcVWmL,cplcUFuFdcVWmR,              & 
& cplcUFuFuhhL,cplcUFuFuhhR,cplcUFuFuVGL,cplcUFuFuVGR,cplcUFuFuVPL,cplcUFuFuVPR,         & 
& cplcUFuFuVZL,cplcUFuFuVZR,cplcUFuGluSuL,cplcUFuGluSuR,cplcChacUFuSdL,cplcChacUFuSdR,   & 
& cplcFecUFuSdL,cplcFecUFuSdR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2_t,ZUR1,ierr,test) 
ZUR2 = ZUR1 
Else 
Call EigenSystem(mat2,MFu2_t,ZUR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
!---------------------------------------- 
! Redoing Calculation using redefined p2 
!----------------------------------------- 
 
i_count = 0 
p2_loop: Do  
i_count = i_count + 1 
sigL=0._dp 
sigR=0._dp 
sigSL=0._dp 
sigSR=0._dp 
p2 = MFu2_t(iL)
Call Sigma1LoopFu(p2,MFu,MFu2,MAh,MAh2,MSu,MSu2,MChi,MChi2,MHpm,MHpm2,MFd,            & 
& MFd2,MSe,MSe2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MGlu,MGlu2,MCha,MCha2,MSd,MSd2,             & 
& MFe,MFe2,cplcUFuFuAhL,cplcUFuFuAhR,cplcUFuChiSuL,cplcUFuChiSuR,cplcUFuFdcHpmL,         & 
& cplcUFuFdcHpmR,cplcUFuFdcSeL,cplcUFuFdcSeR,cplcUFuFdcVWmL,cplcUFuFdcVWmR,              & 
& cplcUFuFuhhL,cplcUFuFuhhR,cplcUFuFuVGL,cplcUFuFuVGR,cplcUFuFuVPL,cplcUFuFuVPR,         & 
& cplcUFuFuVZL,cplcUFuFuVZR,cplcUFuGluSuL,cplcUFuGluSuR,cplcChacUFuSdL,cplcChacUFuSdR,   & 
& cplcFecUFuSdL,cplcFecUFuSdR,sigL,sigR,sigSL,sigSR)

mat1 = mat1a - SigSL - 0.5_dp*(MatMul(SigR,mat1a) + MatMul(mat1a,SigL)) 
 
mat2 = Matmul(Transpose(Conjg(mat1)),mat1) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2_t,ZUR1,ierr,test) 
ZUR2 = ZUR1 
Else 
Call EigenSystem(mat2,MFu2_t,ZUR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
MFu2_1L(il) = MFu2_t(il) 
MFu_1L(il) = Sqrt(MFu2_1L(il)) 
 
If (p2.Ne.0._dp) Then 
  test(1) = Abs(MFu2_1L(il)-p2)/p2
Else 
  test(2) = Abs(MFu2_1L(il))
End If 
If (Abs(MFu2_1L(il)).lt.1.0E-30_dp) Exit p2_loop 
If (test(1).lt.0.1_dp*delta) Exit p2_loop 
If(i_count.gt.30) then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Exit p2_loop 
End if
End Do p2_loop 
mat2 = Matmul(mat1,Transpose(Conjg(mat1))) 
If (ForceRealMatrices) mat2 = Real(mat2,dp) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2_t,ZUL1,ierr,test) 
 
 
ZUL2 = ZUL1 
Else 
Call EigenSystem(mat2,MFu2_t,ZUL2,ierr,test) 
 
 
End If 
ZUL2 = Conjg(ZUL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Call TerminateProgram 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZUL2),mat1),Transpose( Conjg(ZUR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZUR2(i1,:) = phaseM *ZUR2(i1,:) 
 End if 
End Do 
 
ZULOS_p2(il,:) = ZUL2(il,:) 
 ZUROS_p2(il,:) = ZUR2(il,:) 
 ZUL_1L = ZUL2 
 ZUR_1L = ZUR2 
 End Do  
 
Iname = Iname -1 
End Subroutine OneLoopFu
 
 
Subroutine Sigma1LoopFu(p2,MFu,MFu2,MAh,MAh2,MSu,MSu2,MChi,MChi2,MHpm,MHpm2,          & 
& MFd,MFd2,MSe,MSe2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MGlu,MGlu2,MCha,MCha2,MSd,              & 
& MSd2,MFe,MFe2,cplcUFuFuAhL,cplcUFuFuAhR,cplcUFuChiSuL,cplcUFuChiSuR,cplcUFuFdcHpmL,    & 
& cplcUFuFdcHpmR,cplcUFuFdcSeL,cplcUFuFdcSeR,cplcUFuFdcVWmL,cplcUFuFdcVWmR,              & 
& cplcUFuFuhhL,cplcUFuFuhhR,cplcUFuFuVGL,cplcUFuFuVGR,cplcUFuFuVPL,cplcUFuFuVPR,         & 
& cplcUFuFuVZL,cplcUFuFuVZR,cplcUFuGluSuL,cplcUFuGluSuR,cplcChacUFuSdL,cplcChacUFuSdR,   & 
& cplcFecUFuSdL,cplcFecUFuSdR,sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MFu(3),MFu2(3),MAh(2),MAh2(2),MSu(6),MSu2(6),MChi(4),MChi2(4),MHpm(2),MHpm2(2),       & 
& MFd(3),MFd2(3),MSe(6),MSe2(6),MVWm,MVWm2,Mhh(2),Mhh2(2),MVZ,MVZ2,MGlu,MGlu2,           & 
& MCha(2),MCha2(2),MSd(6),MSd2(6),MFe(3),MFe2(3)

Complex(dp), Intent(in) :: cplcUFuFuAhL(3,3,2),cplcUFuFuAhR(3,3,2),cplcUFuChiSuL(3,4,6),cplcUFuChiSuR(3,4,6),    & 
& cplcUFuFdcHpmL(3,3,2),cplcUFuFdcHpmR(3,3,2),cplcUFuFdcSeL(3,3,6),cplcUFuFdcSeR(3,3,6), & 
& cplcUFuFdcVWmL(3,3),cplcUFuFdcVWmR(3,3),cplcUFuFuhhL(3,3,2),cplcUFuFuhhR(3,3,2),       & 
& cplcUFuFuVGL(3,3),cplcUFuFuVGR(3,3),cplcUFuFuVPL(3,3),cplcUFuFuVPR(3,3),               & 
& cplcUFuFuVZL(3,3),cplcUFuFuVZR(3,3),cplcUFuGluSuL(3,6),cplcUFuGluSuR(3,6),             & 
& cplcChacUFuSdL(2,3,6),cplcChacUFuSdR(2,3,6),cplcFecUFuSdL(3,3,6),cplcFecUFuSdR(3,3,6)

Complex(dp), Intent(out) :: SigL(3,3),SigR(3,3), SigSL(3,3), SigSR(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(3,3), sumR(3,3), sumSL(3,3), sumSR(3,3) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
 SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Fu, Ah 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFu2(i1),MAh2(i2)),dp) 
B0m2 = MFu(i1)*Real(SA_B0(p2,MFu2(i1),MAh2(i2)),dp) 
coupL1 = cplcUFuFuAhL(gO1,i1,i2)
coupR1 = cplcUFuFuAhR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFuFuAhL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFuFuAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Su, Chi 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MChi2(i2),MSu2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSu2(i1)),dp) 
coupL1 = cplcUFuChiSuL(gO1,i2,i1)
coupR1 = cplcUFuChiSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuChiSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuChiSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], Fd 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFd2(i2),MHpm2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MHpm2(i1)),dp) 
coupL1 = cplcUFuFdcHpmL(gO1,i2,i1)
coupR1 = cplcUFuFdcHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuFdcHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuFdcHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFd2(i2),MSe2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSe2(i1)),dp) 
coupL1 = cplcUFuFdcSeL(gO1,i2,i1)
coupR1 = cplcUFuFdcSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuFdcSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuFdcSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Fd 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),MVWm2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),MVWm2)-0.5_dp*rMS,dp) 
coupL1 = cplcUFuFdcVWmL(gO1,i2)
coupR1 = cplcUFuFdcVWmR(gO1,i2)
coupL2 =  Conjg(cplcUFuFdcVWmL(gO2,i2))
coupR2 =  Conjg(cplcUFuFdcVWmR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! hh, Fu 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFu2(i2),Mhh2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(p2,MFu2(i2),Mhh2(i1)),dp) 
coupL1 = cplcUFuFuhhL(gO1,i2,i1)
coupR1 = cplcUFuFuhhR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuFuhhL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuFuhhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VG, Fu 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),0._dp)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),0._dp)-0.5_dp*rMS,dp) 
coupL1 = cplcUFuFuVGL(gO1,i2)
coupR1 = cplcUFuFuVGR(gO1,i2)
coupL2 =  Conjg(cplcUFuFuVGL(gO2,i2))
coupR2 =  Conjg(cplcUFuFuVGR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +4._dp/3._dp* sumL
SigR = SigR +4._dp/3._dp* sumR 
SigSL = SigSL +4._dp/3._dp* sumSL 
SigSR = SigSR +4._dp/3._dp* sumSR 
    End Do 
 !------------------------ 
! VP, Fu 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),0._dp)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),0._dp)-0.5_dp*rMS,dp) 
coupL1 = cplcUFuFuVPL(gO1,i2)
coupR1 = cplcUFuFuVPR(gO1,i2)
coupL2 =  Conjg(cplcUFuFuVPL(gO2,i2))
coupR2 =  Conjg(cplcUFuFuVPR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! VZ, Fu 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),MVZ2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),MVZ2)-0.5_dp*rMS,dp) 
coupL1 = cplcUFuFuVZL(gO1,i2)
coupR1 = cplcUFuFuVZR(gO1,i2)
coupL2 =  Conjg(cplcUFuFuVZL(gO2,i2))
coupR2 =  Conjg(cplcUFuFuVZR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Su, Glu 
!------------------------ 
    Do i1 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MGlu2,MSu2(i1)),dp) 
B0m2 = MGlu*Real(SA_B0(p2,MGlu2,MSu2(i1)),dp) 
coupL1 = cplcUFuGluSuL(gO1,i1)
coupR1 = cplcUFuGluSuR(gO1,i1)
coupL2 =  Conjg(cplcUFuGluSuL(gO2,i1))
coupR2 =  Conjg(cplcUFuGluSuR(gO2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +4._dp/3._dp* sumL
SigR = SigR +4._dp/3._dp* sumR 
SigSL = SigSL +4._dp/3._dp* sumSL 
SigSR = SigSR +4._dp/3._dp* sumSR 
      End Do 
 !------------------------ 
! bar[Cha], Sd 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MCha2(i1),MSd2(i2)),dp) 
B0m2 = MCha(i1)*Real(SA_B0(p2,MCha2(i1),MSd2(i2)),dp) 
coupL1 = cplcChacUFuSdL(i1,gO1,i2)
coupR1 = cplcChacUFuSdR(i1,gO1,i2)
coupL2 =  Conjg(cplcChacUFuSdL(i1,gO2,i2))
coupR2 =  Conjg(cplcChacUFuSdR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -Real(SA_B1(p2,MFe2(i1),MSd2(i2)),dp) 
B0m2 = MFe(i1)*Real(SA_B0(p2,MFe2(i1),MSd2(i2)),dp) 
coupL1 = cplcFecUFuSdL(i1,gO1,i2)
coupR1 = cplcFecUFuSdR(i1,gO1,i2)
coupL2 =  Conjg(cplcFecUFuSdL(i1,gO2,i2))
coupR2 =  Conjg(cplcFecUFuSdR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopFu 
 
Subroutine OneLoopGlu(M3,MSd,MSd2,MFd,MFd2,MSu,MSu2,MFu,MFu2,MGlu,MGlu2,              & 
& cplGluFdcSdL,cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,cplGluGluVGL,cplGluGluVGR,         & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,delta,mass,mass2,kont)

Real(dp), Intent(in) :: MSd(6),MSd2(6),MFd(3),MFd2(3),MSu(6),MSu2(6),MFu(3),MFu2(3),MGlu,MGlu2

Complex(dp), Intent(in) :: M3

Complex(dp), Intent(in) :: cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6),              & 
& cplGluGluVGL,cplGluGluVGR,cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),       & 
& cplcFuGluSuR(3,6)

Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi, mi2, p2, test_m2 
Complex(dp) :: PiSf, SigL, SigR, SigSL, SigSR 
Real(dp), Intent(out) :: mass, mass2 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopGlu'
 
mi = MGlu 

 
p2 = MGlu2
sigL = ZeroC 
sigR = ZeroC 
sigSL = ZeroC 
sigSR = ZeroC 
Call Sigma1LoopGlu(p2,MSd,MSd2,MFd,MFd2,MSu,MSu2,MFu,MFu2,MGlu,MGlu2,cplGluFdcSdL,    & 
& cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,cplGluGluVGL,cplGluGluVGR,cplcFdGluSdL,         & 
& cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,sigL,sigR,sigSL,sigSR)

mass = mi - 0.5_dp*(sigSL + sigSR)- 0.5_dp*MGlu*(SigR+SigL) 
mass2= mass**2 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
p2 =  mass2 
sig = ZeroC 
Call Sigma1LoopGlu(p2,MSd,MSd2,MFd,MFd2,MSu,MSu2,MFu,MFu2,MGlu,MGlu2,cplGluFdcSdL,    & 
& cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,cplGluGluVGL,cplGluGluVGR,cplcFdGluSdL,         & 
& cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,sigL,sigR,sigSL,sigSR)

mass = mi - SigSR - 0.5_dp*MGlu*(SigR+SigL) 
mass2= mass**2 
 If (test_m2.Ne.0._dp) Then 
    test_m2 = Abs(test_m2 - mass2) / test_m2 
 Else 
    test_m2 = Abs(mass2) 
 End If 
 If (mass2.Ge.0._dp) Then 
   If (RotateNegativeFermionMasses) Then 
    mass = sqrt(mass2) 
   End if 
  Else 
 If (Abs(mass2).lt.1.0E-30_dp) test_m2 = 0._dp 
     Write(ErrCan,*) 'Warning from routine'//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
   SignOfMassChanged = .True. 
   mass = 0._dp 
  End If 
If (test_m2.LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopGlu
 
 
Subroutine Sigma1LoopGlu(p2,MSd,MSd2,MFd,MFd2,MSu,MSu2,MFu,MFu2,MGlu,MGlu2,           & 
& cplGluFdcSdL,cplGluFdcSdR,cplGluFucSuL,cplGluFucSuR,cplGluGluVGL,cplGluGluVGR,         & 
& cplcFdGluSdL,cplcFdGluSdR,cplcFuGluSuL,cplcFuGluSuR,SigL,SigR,SigSL,SigSR)

Implicit None 
Real(dp), Intent(in) :: MSd(6),MSd2(6),MFd(3),MFd2(3),MSu(6),MSu2(6),MFu(3),MFu2(3),MGlu,MGlu2

Complex(dp), Intent(in) :: cplGluFdcSdL(3,6),cplGluFdcSdR(3,6),cplGluFucSuL(3,6),cplGluFucSuR(3,6),              & 
& cplGluGluVGL,cplGluGluVGR,cplcFdGluSdL(3,6),cplcFdGluSdR(3,6),cplcFuGluSuL(3,6),       & 
& cplcFuGluSuR(3,6)

Complex(dp), Intent(out) :: SigL, SigR, SigSL, SigSR 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumSL, sumSR, sumR,sumL 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! conj[Sd], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumR = 0._dp 
SumL = 0._dp 
SumSL = 0._dp 
SumSR = 0._dp 
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),MSd2(i1)),dp) 
B0m2 = -2._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSd2(i1)),dp) 
coupL1 = cplGluFdcSdL(i2,i1)
coupR1 = cplGluFdcSdR(i2,i1)
coupL2 =  Conjg(cplGluFdcSdL(i2,i1))
coupR2 =  Conjg(cplGluFdcSdR(i2,i1))
SumSR = coupL1*coupR2*B0m2 
SumSL = coupR1*coupL2*B0m2 
sumR = coupR1*coupR2*B1m2 
sumL = coupL1*coupL2*B1m2 
SigL = SigL +1._dp/4._dp*sumL
SigR = SigR +1._dp/4._dp*sumR
SigSL = SigSL +1._dp/4._dp*sumSL
SigSR = SigSR +1._dp/4._dp*sumSR
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Fu 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumR = 0._dp 
SumL = 0._dp 
SumSL = 0._dp 
SumSR = 0._dp 
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),MSu2(i1)),dp) 
B0m2 = -2._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),MSu2(i1)),dp) 
coupL1 = cplGluFucSuL(i2,i1)
coupR1 = cplGluFucSuR(i2,i1)
coupL2 =  Conjg(cplGluFucSuL(i2,i1))
coupR2 =  Conjg(cplGluFucSuR(i2,i1))
SumSR = coupL1*coupR2*B0m2 
SumSL = coupR1*coupL2*B0m2 
sumR = coupR1*coupR2*B1m2 
sumL = coupL1*coupL2*B1m2 
SigL = SigL +1._dp/4._dp*sumL
SigR = SigR +1._dp/4._dp*sumR
SigSL = SigSL +1._dp/4._dp*sumSL
SigSR = SigSR +1._dp/4._dp*sumSR
      End Do 
     End Do 
 !------------------------ 
! VG, Glu 
!------------------------ 
SumR = 0._dp 
SumL = 0._dp 
SumSL = 0._dp 
SumSR = 0._dp 
B1m2 = -4._dp*(Real(SA_B1(p2,MGlu2,0._dp),dp)+ 0.5_dp*rMS) 
B0m2 = -8._dp*MGlu*(Real(SA_B0(p2,MGlu2,0._dp),dp) - 0.5_dp*rMS) 
coupL1 = cplGluGluVGL
coupR1 = cplGluGluVGR
coupL2 =  Conjg(cplGluGluVGL)
coupR2 =  Conjg(cplGluGluVGR)
SumSL = coupL1*coupR2*B0m2 
SumSR = coupR1*coupL2*B0m2 
sumR = coupL1*coupL2*B1m2 
sumL = coupR1*coupR2*B1m2 
SigL = SigL +3._dp/2._dp*sumL
SigR = SigR +3._dp/2._dp*sumR
SigSL = SigSL +3._dp/2._dp*sumSL
SigSR = SigSR +3._dp/2._dp*sumSR
!------------------------ 
! bar[Fd], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumR = 0._dp 
SumL = 0._dp 
SumSL = 0._dp 
SumSR = 0._dp 
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i1),MSd2(i2)),dp) 
B0m2 = -2._dp*MFd(i1)*Real(SA_B0(p2,MFd2(i1),MSd2(i2)),dp) 
coupL1 = cplcFdGluSdL(i1,i2)
coupR1 = cplcFdGluSdR(i1,i2)
coupL2 =  Conjg(cplcFdGluSdL(i1,i2))
coupR2 =  Conjg(cplcFdGluSdR(i1,i2))
SumSR = coupL1*coupR2*B0m2 
SumSL = coupR1*coupL2*B0m2 
sumR = coupR1*coupR2*B1m2 
sumL = coupL1*coupL2*B1m2 
SigL = SigL +1._dp/4._dp*sumL
SigR = SigR +1._dp/4._dp*sumR
SigSL = SigSL +1._dp/4._dp*sumSL
SigSR = SigSR +1._dp/4._dp*sumSR
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Su 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumR = 0._dp 
SumL = 0._dp 
SumSL = 0._dp 
SumSR = 0._dp 
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i1),MSu2(i2)),dp) 
B0m2 = -2._dp*MFu(i1)*Real(SA_B0(p2,MFu2(i1),MSu2(i2)),dp) 
coupL1 = cplcFuGluSuL(i1,i2)
coupR1 = cplcFuGluSuR(i1,i2)
coupL2 =  Conjg(cplcFuGluSuL(i1,i2))
coupR2 =  Conjg(cplcFuGluSuR(i1,i2))
SumSR = coupL1*coupR2*B0m2 
SumSL = coupR1*coupL2*B0m2 
sumR = coupR1*coupR2*B1m2 
sumL = coupL1*coupL2*B1m2 
SigL = SigL +1._dp/4._dp*sumL
SigR = SigR +1._dp/4._dp*sumR
SigSL = SigSL +1._dp/4._dp*sumSL
SigSR = SigSR +1._dp/4._dp*sumSR
      End Do 
     End Do 
 

SigL = oo16pi2*SigL 
 SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
 SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopGlu 
 
Subroutine Sigma1LoopFv(p2,MSe,MSe2,MCha,MCha2,MSv,MSv2,MChi,MChi2,MSd,               & 
& MSd2,MFd,MFd2,MHpm,MHpm2,MFe,MFe2,MVWm,MVWm2,MVZ,MVZ2,cplUFvChacSeL,cplUFvChacSeR,     & 
& cplChiUFvcSvL,cplChiUFvcSvR,cplChiUFvSvL,cplChiUFvSvR,cplUFvFdcSdL,cplUFvFdcSdR,       & 
& cplUFvFecHpmL,cplUFvFecHpmR,cplUFvFecSeL,cplUFvFecSeR,cplUFvFecVWmL,cplUFvFecVWmR,     & 
& cplUFvFvVZL,cplUFvFvVZR,cplcFeUFvHpmL,cplcFeUFvHpmR,cplcFdUFvSdL,cplcFdUFvSdR,         & 
& cplcChaUFvSeL,cplcChaUFvSeR,cplcFeUFvSeL,cplcFeUFvSeR,cplcFeUFvVWmL,cplcFeUFvVWmR,     & 
& sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MSe(6),MSe2(6),MCha(2),MCha2(2),MSv(3),MSv2(3),MChi(4),MChi2(4),MSd(6),               & 
& MSd2(6),MFd(3),MFd2(3),MHpm(2),MHpm2(2),MFe(3),MFe2(3),MVWm,MVWm2,MVZ,MVZ2

Complex(dp), Intent(in) :: cplUFvChacSeL(3,2,6),cplUFvChacSeR(3,2,6),cplChiUFvcSvL(4,3,3),cplChiUFvcSvR(4,3,3),  & 
& cplChiUFvSvL(4,3,3),cplChiUFvSvR(4,3,3),cplUFvFdcSdL(3,3,6),cplUFvFdcSdR(3,3,6),       & 
& cplUFvFecHpmL(3,3,2),cplUFvFecHpmR(3,3,2),cplUFvFecSeL(3,3,6),cplUFvFecSeR(3,3,6),     & 
& cplUFvFecVWmL(3,3),cplUFvFecVWmR(3,3),cplUFvFvVZL(3,3),cplUFvFvVZR(3,3),               & 
& cplcFeUFvHpmL(3,3,2),cplcFeUFvHpmR(3,3,2),cplcFdUFvSdL(3,3,6),cplcFdUFvSdR(3,3,6),     & 
& cplcChaUFvSeL(2,3,6),cplcChaUFvSeR(2,3,6),cplcFeUFvSeL(3,3,6),cplcFeUFvSeR(3,3,6),     & 
& cplcFeUFvVWmL(3,3),cplcFeUFvVWmR(3,3)

Complex(dp), Intent(out) :: SigL(3,3),SigR(3,3), SigSL(3,3), SigSR(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(3,3), sumR(3,3), sumSL(3,3), sumSR(3,3) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! conj[Se], Cha 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MCha2(i2),MSe2(i1)),dp) 
B0m2 = 2._dp*MCha(i2)*Real(SA_B0(p2,MCha2(i2),MSe2(i1)),dp) 
coupL1 = cplUFvChacSeL(gO1,i2,i1)
coupR1 = cplUFvChacSeR(gO1,i2,i1)
coupL2 =  Conjg(cplUFvChacSeL(gO2,i2,i1))
coupR2 =  Conjg(cplUFvChacSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Chi 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MChi2(i2),MSv2(i1)),dp) 
B0m2 = 2._dp*MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSv2(i1)),dp) 
coupL1 = cplChiUFvcSvL(i2,gO1,i1)
coupR1 = cplChiUFvcSvR(i2,gO1,i1)
coupL2 =  Conjg(cplChiUFvcSvL(i2,gO2,i1))
coupR2 =  Conjg(cplChiUFvcSvR(i2,gO2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sv, Chi 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MChi2(i2),MSv2(i1)),dp) 
B0m2 = 2._dp*MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSv2(i1)),dp) 
coupL1 = cplChiUFvSvL(i2,gO1,i1)
coupR1 = cplChiUFvSvR(i2,gO1,i1)
coupL2 =  Conjg(cplChiUFvSvL(i2,gO2,i1))
coupR2 =  Conjg(cplChiUFvSvR(i2,gO2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sd], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),MSd2(i1)),dp) 
B0m2 = 2._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSd2(i1)),dp) 
coupL1 = cplUFvFdcSdL(gO1,i2,i1)
coupR1 = cplUFvFdcSdR(gO1,i2,i1)
coupL2 =  Conjg(cplUFvFdcSdL(gO2,i2,i1))
coupR2 =  Conjg(cplUFvFdcSdR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp/2._dp* sumL
SigR = SigR +3._dp/2._dp* sumR 
SigSL = SigSL +3._dp/2._dp* sumSL 
SigSR = SigSR +3._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], Fe 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i2),MHpm2(i1)),dp) 
B0m2 = 2._dp*MFe(i2)*Real(SA_B0(p2,MFe2(i2),MHpm2(i1)),dp) 
coupL1 = cplUFvFecHpmL(gO1,i2,i1)
coupR1 = cplUFvFecHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplUFvFecHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplUFvFecHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Fe 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i2),MSe2(i1)),dp) 
B0m2 = 2._dp*MFe(i2)*Real(SA_B0(p2,MFe2(i2),MSe2(i1)),dp) 
coupL1 = cplUFvFecSeL(gO1,i2,i1)
coupR1 = cplUFvFecSeR(gO1,i2,i1)
coupL2 =  Conjg(cplUFvFecSeL(gO2,i2,i1))
coupR2 =  Conjg(cplUFvFecSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Fe 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -4._dp*(Real(SA_B1(p2,MFe2(i2),MVWm2),dp) + 0.5_dp*rMS) 
B0m2 = -8._dp*MFe(i2)*(Real(SA_B0(p2,MFe2(i2),MVWm2),dp) - 0.5_dp*rMS) 
coupL1 = cplUFvFecVWmL(gO1,i2)
coupR1 = cplUFvFecVWmR(gO1,i2)
coupL2 =  Conjg(cplUFvFecVWmL(gO2,i2))
coupR2 =  Conjg(cplUFvFecVWmR(gO2,i2))
SumSR(gO1,gO2) = coupL2*coupR1*B0m2 
SumSL(gO1,gO2) = coupR2*coupL1*B0m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
    End Do 
 !------------------------ 
! VZ, Fv 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -4._dp*(Real(SA_B1(p2,0._dp,MVZ2),dp) + 0.5_dp*rMS) 
B0m2 = -8._dp*0.*(Real(SA_B0(p2,0._dp,MVZ2),dp) - 0.5_dp*rMS) 
coupL1 = cplUFvFvVZL(gO1,i2)
coupR1 = cplUFvFvVZR(gO1,i2)
coupL2 =  Conjg(cplUFvFvVZL(gO2,i2))
coupR2 =  Conjg(cplUFvFvVZR(gO2,i2))
SumSR(gO1,gO2) = coupL2*coupR1*B0m2 
SumSL(gO1,gO2) = coupR2*coupL1*B0m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
    End Do 
 !------------------------ 
! bar[Fe], Hpm 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i1),MHpm2(i2)),dp) 
B0m2 = 2._dp*MFe(i1)*Real(SA_B0(p2,MFe2(i1),MHpm2(i2)),dp) 
coupL1 = cplcFeUFvHpmL(i1,gO1,i2)
coupR1 = cplcFeUFvHpmR(i1,gO1,i2)
coupL2 =  Conjg(cplcFeUFvHpmL(i1,gO2,i2))
coupR2 =  Conjg(cplcFeUFvHpmR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i1),MSd2(i2)),dp) 
B0m2 = 2._dp*MFd(i1)*Real(SA_B0(p2,MFd2(i1),MSd2(i2)),dp) 
coupL1 = cplcFdUFvSdL(i1,gO1,i2)
coupR1 = cplcFdUFvSdR(i1,gO1,i2)
coupL2 =  Conjg(cplcFdUFvSdL(i1,gO2,i2))
coupR2 =  Conjg(cplcFdUFvSdR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp/2._dp* sumL
SigR = SigR +3._dp/2._dp* sumR 
SigSL = SigSL +3._dp/2._dp* sumSL 
SigSR = SigSR +3._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Se 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MCha2(i1),MSe2(i2)),dp) 
B0m2 = 2._dp*MCha(i1)*Real(SA_B0(p2,MCha2(i1),MSe2(i2)),dp) 
coupL1 = cplcChaUFvSeL(i1,gO1,i2)
coupR1 = cplcChaUFvSeR(i1,gO1,i2)
coupL2 =  Conjg(cplcChaUFvSeL(i1,gO2,i2))
coupR2 =  Conjg(cplcChaUFvSeR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Se 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i1),MSe2(i2)),dp) 
B0m2 = 2._dp*MFe(i1)*Real(SA_B0(p2,MFe2(i1),MSe2(i2)),dp) 
coupL1 = cplcFeUFvSeL(i1,gO1,i2)
coupR1 = cplcFeUFvSeR(i1,gO1,i2)
coupL2 =  Conjg(cplcFeUFvSeL(i1,gO2,i2))
coupR2 =  Conjg(cplcFeUFvSeR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], VWm 
!------------------------ 
    Do i1 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
B1m2 = -4._dp*(Real(SA_B1(p2,MFe2(i1),MVWm2),dp) + 0.5_dp*rMS) 
B0m2 = -8._dp*MFe(i1)*(Real(SA_B0(p2,MFe2(i1),MVWm2),dp) - 0.5_dp*rMS) 
coupL1 = cplcFeUFvVWmL(i1,gO1)
coupR1 = cplcFeUFvVWmR(i1,gO1)
coupL2 =  Conjg(cplcFeUFvVWmL(i1,gO2))
coupR2 =  Conjg(cplcFeUFvVWmR(i1,gO2))
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp/2._dp* sumL
SigR = SigR +1._dp/2._dp* sumR 
SigSL = SigSL +1._dp/2._dp* sumSL 
SigSR = SigSR +1._dp/2._dp* sumSR 
      End Do 
 

SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopFv 
 
Subroutine Pi1LoopVG(p2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MSd,MSd2,MSu,MSu2,               & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplcgGgGVG,cplGluGluVGL,               & 
& cplGluGluVGR,cplSdcSdVG,cplSucSuVG,cplVGVGVG,cplSdcSdVGVG,cplSucSuVGVG,cplVGVGVGVG1,   & 
& cplVGVGVGVG2,cplVGVGVGVG3,kont,res)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,MSd(6),MSd2(6),MSu(6),MSu2(6)

Complex(dp), Intent(in) :: cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplcgGgGVG,       & 
& cplGluGluVGL,cplGluGluVGR,cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplSdcSdVGVG(6,6), & 
& cplSucSuVGVG(6,6),cplVGVGVGVG1,cplVGVGVGVG2,cplVGVGVGVG3

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFd(i1)*MFd(i2)*Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVGL(i1,i2)
coupR1 = cplcFdFdVGR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFu(i2)*Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVGL(i1,i2)
coupR1 = cplcFuFuVGR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gG], gG 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,0._dp,0._dp),dp)
coup1 = cplcgGgGVG
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +3._dp* SumI  
!------------------------ 
! Glu, Glu 
!------------------------ 
sumI = 0._dp 
 
H0m2 = Real(SA_Hloop(p2,MGlu2,MGlu2),dp) 
B0m2 = 4._dp*MGlu*MGlu*Real(SA_B0(p2,MGlu2,MGlu2),dp) 
coupL1 = cplGluGluVGL
coupR1 = cplGluGluVGR
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1.5_dp* SumI  
!------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSd2(i2),MSd2(i1)),dp)  
coup1 = cplSdcSdVG(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSu2(i2),MSu2(i1)),dp)  
coup1 = cplSucSuVG(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VG, VG 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplVGVGVG
coup2 = Conjg(coup1) 
    SumI = -VVVloop(p2,0._dp,0._dp)*coup1*coup2 
res = res +1.5_dp* SumI  
!------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSd2(i1))
 coup1 = cplSdcSdVGVG(i1,i1)
 SumI = coup1*A0m2 
res = res +999* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSu2(i1))
 coup1 = cplSucSuVGVG(i1,i1)
 SumI = coup1*A0m2 
res = res +999* SumI  
      End Do 
 !------------------------ 
! VG 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_A0(0._dp) +RXi/4._dp*SA_A0(0._dp*RXi) 
coup1 = cplVGVGVGVG1
coup2 = cplVGVGVGVG2
coup3 = cplVGVGVGVG3
SumI = ((2._dp*rMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*0._dp-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1._dp/2._dp* SumI  
res = oo16pi2*res 
 
End Subroutine Pi1LoopVG 
 
Subroutine DerPi1LoopVG(p2,MFd,MFd2,MFu,MFu2,MGlu,MGlu2,MSd,MSd2,MSu,MSu2,            & 
& cplcFdFdVGL,cplcFdFdVGR,cplcFuFuVGL,cplcFuFuVGR,cplcgGgGVG,cplGluGluVGL,               & 
& cplGluGluVGR,cplSdcSdVG,cplSucSuVG,cplVGVGVG,cplSdcSdVGVG,cplSucSuVGVG,cplVGVGVGVG1,   & 
& cplVGVGVGVG2,cplVGVGVGVG3,kont,res)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MFu(3),MFu2(3),MGlu,MGlu2,MSd(6),MSd2(6),MSu(6),MSu2(6)

Complex(dp), Intent(in) :: cplcFdFdVGL(3,3),cplcFdFdVGR(3,3),cplcFuFuVGL(3,3),cplcFuFuVGR(3,3),cplcgGgGVG,       & 
& cplGluGluVGL,cplGluGluVGR,cplSdcSdVG(6,6),cplSucSuVG(6,6),cplVGVGVG,cplSdcSdVGVG(6,6), & 
& cplSucSuVGVG(6,6),cplVGVGVGVG1,cplVGVGVGVG2,cplVGVGVGVG3

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFd(i1)*MFd(i2)*Real(SA_DerB0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVGL(i1,i2)
coupR1 = cplcFdFdVGR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFu(i2)*Real(SA_DerB0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVGL(i1,i2)
coupR1 = cplcFuFuVGR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gG], gG 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVG2,MVG2),dp)
coup1 = cplcgGgGVG
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +3._dp* SumI  
!------------------------ 
! Glu, Glu 
!------------------------ 
sumI = 0._dp 
 
H0m2 = Real(SA_DerHloop(p2,MGlu2,MGlu2),dp) 
B0m2 = 4._dp*MGlu*MGlu*Real(SA_DerB0(p2,MGlu2,MGlu2),dp) 
coupL1 = cplGluGluVGL
coupR1 = cplGluGluVGR
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1.5_dp* SumI  
!------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSd2(i2),MSd2(i1)),dp)  
coup1 = cplSdcSdVG(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSu2(i2),MSu2(i1)),dp)  
coup1 = cplSucSuVG(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VG, VG 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplVGVGVG
coup2 = Conjg(coup1) 
    SumI = -DerVVVloop(p2,MVG2,MVG2)*coup1*coup2 
res = res +1.5_dp* SumI  
!------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSd2(i1))
 coup1 = cplSdcSdVGVG(i1,i1)
 SumI = coup1*A0m2 
res = res +999* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSu2(i1))
 coup1 = cplSucSuVGVG(i1,i1)
 SumI = coup1*A0m2 
res = res +999* SumI  
      End Do 
 !------------------------ 
! VG 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_DerA0(MVG2) +RXi/4._dp*SA_DerA0(MVG2*RXi) 
coup1 = cplVGVGVGVG1
coup2 = cplVGVGVGVG2
coup3 = cplVGVGVGVG3
SumI = ((2._dp*DerrMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVG2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1._dp/2._dp* SumI  
res = oo16pi2*res 
 
End Subroutine DerPi1LoopVG 
 
Subroutine Pi1LoopVP(p2,MCha,MCha2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MHpm,MHpm2,             & 
& MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,cplcChaChaVPL,cplcChaChaVPR,cplcFdFdVPL,         & 
& cplcFdFdVPR,cplcFeFeVPL,cplcFeFeVPR,cplcFuFuVPL,cplcFuFuVPR,cplcgWmgWmVP,              & 
& cplcgWpCgWpCVP,cplHpmcHpmVP,cplHpmcVWmVP,cplSdcSdVP,cplSecSeVP,cplSucSuVP,             & 
& cplcVWmVPVWm,cplHpmcHpmVPVP,cplSdcSdVPVP,cplSecSeVPVP,cplSucSuVPVP,cplcVWmVPVPVWm3,    & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,kont,res)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MHpm(2),MHpm2(2),       & 
& MVWm,MVWm2,MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6)

Complex(dp), Intent(in) :: cplcChaChaVPL(2,2),cplcChaChaVPR(2,2),cplcFdFdVPL(3,3),cplcFdFdVPR(3,3),              & 
& cplcFeFeVPL(3,3),cplcFeFeVPR(3,3),cplcFuFuVPL(3,3),cplcFuFuVPR(3,3),cplcgWmgWmVP,      & 
& cplcgWpCgWpCVP,cplHpmcHpmVP(2,2),cplHpmcVWmVP(2),cplSdcSdVP(6,6),cplSecSeVP(6,6),      & 
& cplSucSuVP(6,6),cplcVWmVPVWm,cplHpmcHpmVPVP(2,2),cplSdcSdVPVP(6,6),cplSecSeVPVP(6,6),  & 
& cplSucSuVPVP(6,6),cplcVWmVPVPVWm3,cplcVWmVPVPVWm1,cplcVWmVPVPVWm2

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MCha2(i1).gt.50._dp**2).and.(MCha2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MCha2(i1).lt.50._dp**2).and.(MCha2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_Hloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = 4._dp*MCha(i1)*MCha(i2)*Real(SA_B0(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVPL(i1,i2)
coupR1 = cplcChaChaVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MFd2(i1).gt.50._dp**2).and.(MFd2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MFd2(i1).lt.50._dp**2).and.(MFd2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_Hloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFd(i1)*MFd(i2)*Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVPL(i1,i2)
coupR1 = cplcFdFdVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MFe2(i1).gt.50._dp**2).and.(MFe2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MFe2(i1).lt.50._dp**2).and.(MFe2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_Hloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = 4._dp*MFe(i1)*MFe(i2)*Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVPL(i1,i2)
coupR1 = cplcFeFeVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MFu2(i1).gt.50._dp**2).and.(MFu2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MFu2(i1).lt.50._dp**2).and.(MFu2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_Hloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFu(i2)*Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVPL(i1,i2)
coupR1 = cplcFuFuVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MVWm2.lt.50._dp**2)) )   Then 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWmgWmVP
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
End If 
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MVWm2.lt.50._dp**2)) )   Then 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWpCgWpCVP
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
End If 
!------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MHpm2(i2).gt.50._dp**2).and.(MHpm2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MHpm2(i2).lt.50._dp**2).and.(MHpm2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(VSSloop(p2,MHpm2(i2),MHpm2(i1)),dp)  
coup1 = cplHpmcHpmVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MHpm2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MHpm2(i2).lt.50._dp**2)) )   Then 
B0m2 = Real(VVSloop(p2,MVWm2,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVP(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +2._dp* SumI  
End If 
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSd2(i2).gt.50._dp**2).and.(MSd2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSd2(i2).lt.50._dp**2).and.(MSd2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(VSSloop(p2,MSd2(i2),MSd2(i1)),dp)  
coup1 = cplSdcSdVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSe2(i2).gt.50._dp**2).and.(MSe2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSe2(i2).lt.50._dp**2).and.(MSe2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(VSSloop(p2,MSe2(i2),MSe2(i1)),dp)  
coup1 = cplSecSeVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSu2(i2).gt.50._dp**2).and.(MSu2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSu2(i2).lt.50._dp**2).and.(MSu2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(VSSloop(p2,MSu2(i2),MSu2(i1)),dp)  
coup1 = cplSucSuVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MVWm2.lt.50._dp**2)) )   Then 
coup1 = cplcVWmVPVWm
coup2 = Conjg(coup1) 
    SumI = -VVVloop(p2,MVWm2,MVWm2)*coup1*coup2 
res = res +1._dp* SumI  
End If 
!------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MHpm2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MHpm2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_A0(MHpm2(i1))
 coup1 = cplHpmcHpmVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
End If 
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSd2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSd2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_A0(MSd2(i1))
 coup1 = cplSdcSdVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
End If 
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSe2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSe2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_A0(MSe2(i1))
 coup1 = cplSecSeVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
End If 
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSu2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSu2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_A0(MSu2(i1))
 coup1 = cplSucSuVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
End If 
      End Do 
 !------------------------ 
! conj[VWm] 
!------------------------ 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2)) )   Then 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_A0(MVWm2) +RXi/4._dp*SA_A0(MVWm2*RXi) 
coup1 = cplcVWmVPVPVWm3
coup2 = cplcVWmVPVPVWm1
coup3 = cplcVWmVPVPVWm2
SumI = ((2._dp*rMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
End If 
res = oo16pi2*res 
 
End Subroutine Pi1LoopVP 
 
Subroutine DerPi1LoopVP(p2,MCha,MCha2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MHpm,MHpm2,          & 
& MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,cplcChaChaVPL,cplcChaChaVPR,cplcFdFdVPL,         & 
& cplcFdFdVPR,cplcFeFeVPL,cplcFeFeVPR,cplcFuFuVPL,cplcFuFuVPR,cplcgWmgWmVP,              & 
& cplcgWpCgWpCVP,cplHpmcHpmVP,cplHpmcVWmVP,cplSdcSdVP,cplSecSeVP,cplSucSuVP,             & 
& cplcVWmVPVWm,cplHpmcHpmVPVP,cplSdcSdVPVP,cplSecSeVPVP,cplSucSuVPVP,cplcVWmVPVPVWm3,    & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,kont,res)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MHpm(2),MHpm2(2),       & 
& MVWm,MVWm2,MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6)

Complex(dp), Intent(in) :: cplcChaChaVPL(2,2),cplcChaChaVPR(2,2),cplcFdFdVPL(3,3),cplcFdFdVPR(3,3),              & 
& cplcFeFeVPL(3,3),cplcFeFeVPR(3,3),cplcFuFuVPL(3,3),cplcFuFuVPR(3,3),cplcgWmgWmVP,      & 
& cplcgWpCgWpCVP,cplHpmcHpmVP(2,2),cplHpmcVWmVP(2),cplSdcSdVP(6,6),cplSecSeVP(6,6),      & 
& cplSucSuVP(6,6),cplcVWmVPVWm,cplHpmcHpmVPVP(2,2),cplSdcSdVPVP(6,6),cplSecSeVPVP(6,6),  & 
& cplSucSuVPVP(6,6),cplcVWmVPVPVWm3,cplcVWmVPVPVWm1,cplcVWmVPVPVWm2

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MCha2(i1).gt.50._dp**2).and.(MCha2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MCha2(i1).lt.50._dp**2).and.(MCha2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_DerHloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = 4._dp*MCha(i1)*MCha(i2)*Real(SA_DerB0(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVPL(i1,i2)
coupR1 = cplcChaChaVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MFd2(i1).gt.50._dp**2).and.(MFd2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MFd2(i1).lt.50._dp**2).and.(MFd2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_DerHloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFd(i1)*MFd(i2)*Real(SA_DerB0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVPL(i1,i2)
coupR1 = cplcFdFdVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MFe2(i1).gt.50._dp**2).and.(MFe2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MFe2(i1).lt.50._dp**2).and.(MFe2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_DerHloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = 4._dp*MFe(i1)*MFe(i2)*Real(SA_DerB0(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVPL(i1,i2)
coupR1 = cplcFeFeVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MFu2(i1).gt.50._dp**2).and.(MFu2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MFu2(i1).lt.50._dp**2).and.(MFu2(i2).lt.50._dp**2)) )   Then 
H0m2 = Real(SA_DerHloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFu(i2)*Real(SA_DerB0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVPL(i1,i2)
coupR1 = cplcFuFuVPR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MVWm2.lt.50._dp**2)) )   Then 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWmgWmVP
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
End If 
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MVWm2.lt.50._dp**2)) )   Then 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWpCgWpCVP
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
End If 
!------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MHpm2(i2).gt.50._dp**2).and.(MHpm2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MHpm2(i2).lt.50._dp**2).and.(MHpm2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(DerVSSloop(p2,MHpm2(i2),MHpm2(i1)),dp)  
coup1 = cplHpmcHpmVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MHpm2(i2).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MHpm2(i2).lt.50._dp**2)) )   Then 
B0m2 = Real(DerVVSloop(p2,MVWm2,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVP(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +2._dp* SumI  
End If 
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSd2(i2).gt.50._dp**2).and.(MSd2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSd2(i2).lt.50._dp**2).and.(MSd2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(DerVSSloop(p2,MSd2(i2),MSd2(i1)),dp)  
coup1 = cplSdcSdVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSe2(i2).gt.50._dp**2).and.(MSe2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSe2(i2).lt.50._dp**2).and.(MSe2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(DerVSSloop(p2,MSe2(i2),MSe2(i1)),dp)  
coup1 = cplSecSeVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSu2(i2).gt.50._dp**2).and.(MSu2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSu2(i2).lt.50._dp**2).and.(MSu2(i1).lt.50._dp**2)) )   Then 
B22m2 = Real(DerVSSloop(p2,MSu2(i2),MSu2(i1)),dp)  
coup1 = cplSucSuVP(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
End If 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2).and.(MVWm2.lt.50._dp**2)) )   Then 
coup1 = cplcVWmVPVWm
coup2 = Conjg(coup1) 
    SumI = -DerVVVloop(p2,MVWm2,MVWm2)*coup1*coup2 
res = res +1._dp* SumI  
End If 
!------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MHpm2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MHpm2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_DerA0(MHpm2(i1))
 coup1 = cplHpmcHpmVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
End If 
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSd2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSd2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_DerA0(MSd2(i1))
 coup1 = cplSdcSdVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
End If 
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSe2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSe2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_DerA0(MSe2(i1))
 coup1 = cplSecSeVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
End If 
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MSu2(i1).gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MSu2(i1).lt.50._dp**2)) )   Then 
SumI = 0._dp 
 A0m2 = SA_DerA0(MSu2(i1))
 coup1 = cplSucSuVPVP(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
End If 
      End Do 
 !------------------------ 
! conj[VWm] 
!------------------------ 
If (((.not.OnlyHeavyStates).and.(.not.OnlyLightStates)) & 
  & .or.((OnlyHeavyStates).and.(MVWm2.gt.50._dp**2))   & 
  & .or.((OnlyLightStates).and.(MVWm2.lt.50._dp**2)) )   Then 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_DerA0(MVWm2) +RXi/4._dp*SA_DerA0(MVWm2*RXi) 
coup1 = cplcVWmVPVPVWm3
coup2 = cplcVWmVPVPVWm1
coup3 = cplcVWmVPVPVWm2
SumI = ((2._dp*DerrMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
End If 
res = oo16pi2*res 
 
End Subroutine DerPi1LoopVP 
 
Subroutine OneLoopVZ(g1,g2,vd,vu,TW,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,          & 
& MFd,MFd2,MFe,MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,           & 
& MSu,MSu2,MSv,MSv2,cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,cplChiChiVZR,     & 
& cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplFvFvVZL,cplFvFvVZR,cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,cplHpmcHpmVZ,              & 
& cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,cplAhAhVZVZ,     & 
& cplhhhhVZVZ,cplHpmcHpmVZVZ,cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,cplSvcSvVZVZ,        & 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,delta,mass,mass2,kont)

Real(dp), Intent(in) :: Mhh(2),Mhh2(2),MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),               & 
& MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MVZ,MVZ2,MHpm(2),MHpm2(2),MVWm,MVWm2,            & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3)

Real(dp), Intent(in) :: g1,g2,vd,vu,TW

Complex(dp), Intent(in) :: cplAhhhVZ(2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),cplChiChiVZL(4,4),               & 
& cplChiChiVZR(4,4),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3), & 
& cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplFvFvVZL(3,3),cplFvFvVZR(3,3),cplcgWmgWmVZ,        & 
& cplcgWpCgWpCVZ,cplhhVZVZ(2),cplHpmcHpmVZ(2,2),cplHpmcVWmVZ(2),cplSdcSdVZ(6,6),         & 
& cplSecSeVZ(6,6),cplSucSuVZ(6,6),cplSvcSvVZ(3,3),cplcVWmVWmVZ,cplAhAhVZVZ(2,2),         & 
& cplhhhhVZVZ(2,2),cplHpmcHpmVZVZ(2,2),cplSdcSdVZVZ(6,6),cplSecSeVZVZ(6,6),              & 
& cplSucSuVZVZ(6,6),cplSvcSvVZVZ(3,3),cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3

Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi, mi2, p2, test_m2 
Complex(dp) :: PiSf, SigL, SigR, SigSL, SigSR 
Real(dp), Intent(out) :: mass, mass2 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopVZ'
 
mi2 = MVZ2 

 
p2 = MVZ2
PiSf = ZeroC 
Call Pi1LoopVZ(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,cplChiChiVZR,              & 
& cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplFvFvVZL,cplFvFvVZR,cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,cplHpmcHpmVZ,              & 
& cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,cplAhAhVZVZ,     & 
& cplhhhhVZVZ,cplHpmcHpmVZVZ,cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,cplSvcSvVZVZ,        & 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,PiSf)

mass2 = mi2 + Real(PiSf,dp) 
mass = sqrt(mass2) 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
p2 =  mass2 
PiSf = ZeroC 
Call Pi1LoopVZ(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,               & 
& MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,               & 
& MSv,MSv2,cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,cplChiChiVZR,              & 
& cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplFvFvVZL,cplFvFvVZR,cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,cplHpmcHpmVZ,              & 
& cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,cplAhAhVZVZ,     & 
& cplhhhhVZVZ,cplHpmcHpmVZVZ,cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,cplSvcSvVZVZ,        & 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,PiSf)

mass2 = mi2 + Real(PiSf,dp) 
mass = sqrt(mass2) 
 If (test_m2.Ne.0._dp) Then 
    test_m2 = Abs(test_m2 - mass2) / test_m2 
 Else 
    test_m2 = Abs(mass2) 
 End If 
 If (mass2.Ge.0._dp) Then 
   If (RotateNegativeFermionMasses) Then 
    mass = sqrt(mass2) 
   End if 
  Else 
 If (Abs(mass2).lt.1.0E-30_dp) test_m2 = 0._dp 
     Write(ErrCan,*) 'Warning from routine'//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
   SignOfMassChanged = .True. 
   mass = 0._dp 
  End If 
If (test_m2.LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopVZ
 
 
Subroutine Pi1LoopVZ(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,             & 
& MFe,MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,           & 
& MSv,MSv2,cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,cplChiChiVZR,              & 
& cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplFvFvVZL,cplFvFvVZR,cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,cplHpmcHpmVZ,              & 
& cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,cplAhAhVZVZ,     & 
& cplhhhhVZVZ,cplHpmcHpmVZVZ,cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,cplSvcSvVZVZ,        & 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,res)

Implicit None 
Real(dp), Intent(in) :: Mhh(2),Mhh2(2),MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),               & 
& MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MVZ,MVZ2,MHpm(2),MHpm2(2),MVWm,MVWm2,            & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplAhhhVZ(2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),cplChiChiVZL(4,4),               & 
& cplChiChiVZR(4,4),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3), & 
& cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplFvFvVZL(3,3),cplFvFvVZR(3,3),cplcgWmgWmVZ,        & 
& cplcgWpCgWpCVZ,cplhhVZVZ(2),cplHpmcHpmVZ(2,2),cplHpmcVWmVZ(2),cplSdcSdVZ(6,6),         & 
& cplSecSeVZ(6,6),cplSucSuVZ(6,6),cplSvcSvVZ(3,3),cplcVWmVWmVZ,cplAhAhVZVZ(2,2),         & 
& cplhhhhVZVZ(2,2),cplHpmcHpmVZVZ(2,2),cplSdcSdVZVZ(6,6),cplSecSeVZVZ(6,6),              & 
& cplSucSuVZVZ(6,6),cplSvcSvVZVZ(3,3),cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! hh, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(VSSloop(p2,MAh2(i2),Mhh2(i1)),dp)  
coup1 = cplAhhhVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 H0m2 = Real(SA_Hloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = 4._dp*MCha(i1)*MCha(i2)*Real(SA_B0(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVZL(i1,i2)
coupR1 = cplcChaChaVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 H0m2 = Real(SA_Hloop(p2,MChi2(i1),MChi2(i2)),dp) 
B0m2 = 4._dp*MChi(i1)*MChi(i2)*Real(SA_B0(p2,MChi2(i1),MChi2(i2)),dp) 
coupL1 = cplChiChiVZL(i1,i2)
coupR1 = cplChiChiVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFd(i1)*MFd(i2)*Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = 4._dp*MFe(i1)*MFe(i2)*Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFu(i2)*Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVZL(i1,i2)
coupR1 = cplcFuFuVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Fv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,0._dp,0._dp),dp) 
B0m2 = 4._dp*0.*0.*Real(SA_B0(p2,0._dp,0._dp),dp) 
coupL1 = cplFvFvVZL(i1,i2)
coupR1 = cplFvFvVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWmgWmVZ
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWpCgWpCVZ
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! VZ, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(VVSloop(p2,MVZ2,Mhh2(i2)),dp)
coup1 = cplhhVZVZ(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(VSSloop(p2,MHpm2(i2),MHpm2(i1)),dp)  
coup1 = cplHpmcHpmVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(VVSloop(p2,MVWm2,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVZ(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +2._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSd2(i2),MSd2(i1)),dp)  
coup1 = cplSdcSdVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSe2(i2),MSe2(i1)),dp)  
coup1 = cplSecSeVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSu2(i2),MSu2(i1)),dp)  
coup1 = cplSucSuVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 B22m2 = Real(VSSloop(p2,MSv2(i2),MSv2(i1)),dp)  
coup1 = cplSvcSvVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplcVWmVWmVZ
coup2 = Conjg(coup1) 
    SumI = -VVVloop(p2,MVWm2,MVWm2)*coup1*coup2 
res = res +1._dp* SumI  
!------------------------ 
! Ah 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_A0(MAh2(i1))
 coup1 = cplAhAhVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_A0(Mhh2(i1))
 coup1 = cplhhhhVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_A0(MHpm2(i1))
 coup1 = cplHpmcHpmVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSd2(i1))
 coup1 = cplSdcSdVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSe2(i1))
 coup1 = cplSecSeVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSu2(i1))
 coup1 = cplSucSuVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Sv] 
!------------------------ 
    Do i1 = 1, 3
 SumI = 0._dp 
 A0m2 = SA_A0(MSv2(i1))
 coup1 = cplSvcSvVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[VWm] 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_A0(MVWm2) +RXi/4._dp*SA_A0(MVWm2*RXi) 
coup1 = cplcVWmVWmVZVZ1
coup2 = cplcVWmVWmVZVZ2
coup3 = cplcVWmVWmVZVZ3
SumI = ((2._dp*rMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
res = oo16pi2*res 
 
End Subroutine Pi1LoopVZ 
 
Subroutine DerPi1LoopVZ(p2,Mhh,Mhh2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,               & 
& MFd2,MFe,MFe2,MFu,MFu2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,MSd,MSd2,MSe,MSe2,               & 
& MSu,MSu2,MSv,MSv2,cplAhhhVZ,cplcChaChaVZL,cplcChaChaVZR,cplChiChiVZL,cplChiChiVZR,     & 
& cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplFvFvVZL,cplFvFvVZR,cplcgWmgWmVZ,cplcgWpCgWpCVZ,cplhhVZVZ,cplHpmcHpmVZ,              & 
& cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,cplSvcSvVZ,cplcVWmVWmVZ,cplAhAhVZVZ,     & 
& cplhhhhVZVZ,cplHpmcHpmVZVZ,cplSdcSdVZVZ,cplSecSeVZVZ,cplSucSuVZVZ,cplSvcSvVZVZ,        & 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,res)

Implicit None 
Real(dp), Intent(in) :: Mhh(2),Mhh2(2),MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),               & 
& MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MVZ,MVZ2,MHpm(2),MHpm2(2),MVWm,MVWm2,            & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplAhhhVZ(2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),cplChiChiVZL(4,4),               & 
& cplChiChiVZR(4,4),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3), & 
& cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplFvFvVZL(3,3),cplFvFvVZR(3,3),cplcgWmgWmVZ,        & 
& cplcgWpCgWpCVZ,cplhhVZVZ(2),cplHpmcHpmVZ(2,2),cplHpmcVWmVZ(2),cplSdcSdVZ(6,6),         & 
& cplSecSeVZ(6,6),cplSucSuVZ(6,6),cplSvcSvVZ(3,3),cplcVWmVWmVZ,cplAhAhVZVZ(2,2),         & 
& cplhhhhVZVZ(2,2),cplHpmcHpmVZVZ(2,2),cplSdcSdVZVZ(6,6),cplSecSeVZVZ(6,6),              & 
& cplSucSuVZVZ(6,6),cplSvcSvVZVZ(3,3),cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! hh, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(DerVSSloop(p2,MAh2(i2),Mhh2(i1)),dp)  
coup1 = cplAhhhVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 H0m2 = Real(SA_DerHloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = 4._dp*MCha(i1)*MCha(i2)*Real(SA_DerB0(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVZL(i1,i2)
coupR1 = cplcChaChaVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 H0m2 = Real(SA_DerHloop(p2,MChi2(i1),MChi2(i2)),dp) 
B0m2 = 4._dp*MChi(i1)*MChi(i2)*Real(SA_DerB0(p2,MChi2(i1),MChi2(i2)),dp) 
coupL1 = cplChiChiVZL(i1,i2)
coupR1 = cplChiChiVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFd(i1)*MFd(i2)*Real(SA_DerB0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = 4._dp*MFe(i1)*MFe(i2)*Real(SA_DerB0(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFu(i2)*Real(SA_DerB0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVZL(i1,i2)
coupR1 = cplcFuFuVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Fv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,0._dp,0._dp),dp) 
B0m2 = 4._dp*0.*0.*Real(SA_DerB0(p2,0._dp,0._dp),dp) 
coupL1 = cplFvFvVZL(i1,i2)
coupR1 = cplFvFvVZR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWmgWmVZ
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWpCgWpCVZ
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! VZ, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(DerVVSloop(p2,MVZ2,Mhh2(i2)),dp)
coup1 = cplhhVZVZ(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(DerVSSloop(p2,MHpm2(i2),MHpm2(i1)),dp)  
coup1 = cplHpmcHpmVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(DerVVSloop(p2,MVWm2,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVZ(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +2._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSd2(i2),MSd2(i1)),dp)  
coup1 = cplSdcSdVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSe2(i2),MSe2(i1)),dp)  
coup1 = cplSecSeVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSu2(i2),MSu2(i1)),dp)  
coup1 = cplSucSuVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 B22m2 = Real(DerVSSloop(p2,MSv2(i2),MSv2(i1)),dp)  
coup1 = cplSvcSvVZ(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplcVWmVWmVZ
coup2 = Conjg(coup1) 
    SumI = -DerVVVloop(p2,MVWm2,MVWm2)*coup1*coup2 
res = res +1._dp* SumI  
!------------------------ 
! Ah 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_DerA0(MAh2(i1))
 coup1 = cplAhAhVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_DerA0(Mhh2(i1))
 coup1 = cplhhhhVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_DerA0(MHpm2(i1))
 coup1 = cplHpmcHpmVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSd2(i1))
 coup1 = cplSdcSdVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSe2(i1))
 coup1 = cplSecSeVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSu2(i1))
 coup1 = cplSucSuVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Sv] 
!------------------------ 
    Do i1 = 1, 3
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSv2(i1))
 coup1 = cplSvcSvVZVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[VWm] 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_DerA0(MVWm2) +RXi/4._dp*SA_DerA0(MVWm2*RXi) 
coup1 = cplcVWmVWmVZVZ1
coup2 = cplcVWmVWmVZVZ2
coup3 = cplcVWmVWmVZVZ3
SumI = ((2._dp*DerrMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
res = oo16pi2*res 
 
End Subroutine DerPi1LoopVZ 
 
Subroutine OneLoopVWm(g2,vd,vu,MHpm,MHpm2,MAh,MAh2,MChi,MChi2,MCha,MCha2,             & 
& MFu,MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVWm,MVWm2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,             & 
& MSv,MSv2,MSe,MSe2,cplAhHpmcVWm,cplChiChacVWmL,cplChiChacVWmR,cplcFuFdcVWmL,            & 
& cplcFuFdcVWmR,cplFvFecVWmL,cplFvFecVWmR,cplcgWpCgAcVWm,cplcgAgWmcVWm,cplcgZgWmcVWm,    & 
& cplcgWpCgZcVWm,cplhhHpmcVWm,cplhhcVWmVWm,cplHpmcVWmVP,cplHpmcVWmVZ,cplSdcSucVWm,       & 
& cplSecSvcVWm,cplcVWmVPVWm,cplcVWmVWmVZ,cplAhAhcVWmVWm,cplhhhhcVWmVWm,cplHpmcHpmcVWmVWm,& 
& cplSdcSdcVWmVWm,cplSecSecVWmVWm,cplSucSucVWmVWm,cplSvcSvcVWmVWm,cplcVWmVPVPVWm3,       & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,delta,mass,mass2,kont)

Real(dp), Intent(in) :: MHpm(2),MHpm2(2),MAh(2),MAh2(2),MChi(4),MChi2(4),MCha(2),MCha2(2),MFu(3),             & 
& MFu2(3),MFd(3),MFd2(3),MFe(3),MFe2(3),Mhh(2),Mhh2(2),MVWm,MVWm2,MVZ,MVZ2,              & 
& MSu(6),MSu2(6),MSd(6),MSd2(6),MSv(3),MSv2(3),MSe(6),MSe2(6)

Real(dp), Intent(in) :: g2,vd,vu

Complex(dp), Intent(in) :: cplAhHpmcVWm(2,2),cplChiChacVWmL(4,2),cplChiChacVWmR(4,2),cplcFuFdcVWmL(3,3),         & 
& cplcFuFdcVWmR(3,3),cplFvFecVWmL(3,3),cplFvFecVWmR(3,3),cplcgWpCgAcVWm,cplcgAgWmcVWm,   & 
& cplcgZgWmcVWm,cplcgWpCgZcVWm,cplhhHpmcVWm(2,2),cplhhcVWmVWm(2),cplHpmcVWmVP(2),        & 
& cplHpmcVWmVZ(2),cplSdcSucVWm(6,6),cplSecSvcVWm(6,3),cplcVWmVPVWm,cplcVWmVWmVZ,         & 
& cplAhAhcVWmVWm(2,2),cplhhhhcVWmVWm(2,2),cplHpmcHpmcVWmVWm(2,2),cplSdcSdcVWmVWm(6,6),   & 
& cplSecSecVWmVWm(6,6),cplSucSucVWmVWm(6,6),cplSvcSvcVWmVWm(3,3),cplcVWmVPVPVWm3,        & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3

Integer , Intent(inout):: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,j4,il,i_count, ierr 
Real(dp), Intent(in) :: delta 
Real(dp) :: mi, mi2, p2, test_m2 
Complex(dp) :: PiSf, SigL, SigR, SigSL, SigSR 
Real(dp), Intent(out) :: mass, mass2 
Iname = Iname + 1 
NameOfUnit(Iname) = 'OneLoopVWm'
 
mi2 = MVWm2 

 
p2 = MVWm2
PiSf = ZeroC 
Call Pi1LoopVWm(p2,MHpm,MHpm2,MAh,MAh2,MChi,MChi2,MCha,MCha2,MFu,MFu2,MFd,            & 
& MFd2,MFe,MFe2,Mhh,Mhh2,MVWm,MVWm2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,MSv2,MSe,             & 
& MSe2,cplAhHpmcVWm,cplChiChacVWmL,cplChiChacVWmR,cplcFuFdcVWmL,cplcFuFdcVWmR,           & 
& cplFvFecVWmL,cplFvFecVWmR,cplcgWpCgAcVWm,cplcgAgWmcVWm,cplcgZgWmcVWm,cplcgWpCgZcVWm,   & 
& cplhhHpmcVWm,cplhhcVWmVWm,cplHpmcVWmVP,cplHpmcVWmVZ,cplSdcSucVWm,cplSecSvcVWm,         & 
& cplcVWmVPVWm,cplcVWmVWmVZ,cplAhAhcVWmVWm,cplhhhhcVWmVWm,cplHpmcHpmcVWmVWm,             & 
& cplSdcSdcVWmVWm,cplSecSecVWmVWm,cplSucSucVWmVWm,cplSvcSvcVWmVWm,cplcVWmVPVPVWm3,       & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,PiSf)

mass2 = mi2 + Real(PiSf,dp) 
mass = sqrt(mass2) 
i_count = 0 
Do  
i_count = i_count + 1 
test_m2 = mass2 
p2 =  mass2 
PiSf = ZeroC 
Call Pi1LoopVWm(p2,MHpm,MHpm2,MAh,MAh2,MChi,MChi2,MCha,MCha2,MFu,MFu2,MFd,            & 
& MFd2,MFe,MFe2,Mhh,Mhh2,MVWm,MVWm2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,MSv2,MSe,             & 
& MSe2,cplAhHpmcVWm,cplChiChacVWmL,cplChiChacVWmR,cplcFuFdcVWmL,cplcFuFdcVWmR,           & 
& cplFvFecVWmL,cplFvFecVWmR,cplcgWpCgAcVWm,cplcgAgWmcVWm,cplcgZgWmcVWm,cplcgWpCgZcVWm,   & 
& cplhhHpmcVWm,cplhhcVWmVWm,cplHpmcVWmVP,cplHpmcVWmVZ,cplSdcSucVWm,cplSecSvcVWm,         & 
& cplcVWmVPVWm,cplcVWmVWmVZ,cplAhAhcVWmVWm,cplhhhhcVWmVWm,cplHpmcHpmcVWmVWm,             & 
& cplSdcSdcVWmVWm,cplSecSecVWmVWm,cplSucSucVWmVWm,cplSvcSvcVWmVWm,cplcVWmVPVPVWm3,       & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,PiSf)

mass2 = mi2 + Real(PiSf,dp) 
mass = sqrt(mass2) 
 If (test_m2.Ne.0._dp) Then 
    test_m2 = Abs(test_m2 - mass2) / test_m2 
 Else 
    test_m2 = Abs(mass2) 
 End If 
 If (mass2.Ge.0._dp) Then 
   If (RotateNegativeFermionMasses) Then 
    mass = sqrt(mass2) 
   End if 
  Else 
 If (Abs(mass2).lt.1.0E-30_dp) test_m2 = 0._dp 
     Write(ErrCan,*) 'Warning from routine'//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
   SignOfMassChanged = .True. 
   mass = 0._dp 
  End If 
If (test_m2.LT.0.1_dp*delta) Exit 
If (i_count.Gt.30) Then 
  Write(*,*) "Problem in "//NameOfUnit(Iname), test_m2, mass2 
  kont = -510 
  Call AddError(510) 
 Exit 
End If 
End Do 
 
 
Iname = Iname -1 
End Subroutine OneLoopVWm
 
 
Subroutine Pi1LoopVWm(p2,MHpm,MHpm2,MAh,MAh2,MChi,MChi2,MCha,MCha2,MFu,               & 
& MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVWm,MVWm2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,             & 
& MSv2,MSe,MSe2,cplAhHpmcVWm,cplChiChacVWmL,cplChiChacVWmR,cplcFuFdcVWmL,cplcFuFdcVWmR,  & 
& cplFvFecVWmL,cplFvFecVWmR,cplcgWpCgAcVWm,cplcgAgWmcVWm,cplcgZgWmcVWm,cplcgWpCgZcVWm,   & 
& cplhhHpmcVWm,cplhhcVWmVWm,cplHpmcVWmVP,cplHpmcVWmVZ,cplSdcSucVWm,cplSecSvcVWm,         & 
& cplcVWmVPVWm,cplcVWmVWmVZ,cplAhAhcVWmVWm,cplhhhhcVWmVWm,cplHpmcHpmcVWmVWm,             & 
& cplSdcSdcVWmVWm,cplSecSecVWmVWm,cplSucSucVWmVWm,cplSvcSvcVWmVWm,cplcVWmVPVPVWm3,       & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,res)

Implicit None 
Real(dp), Intent(in) :: MHpm(2),MHpm2(2),MAh(2),MAh2(2),MChi(4),MChi2(4),MCha(2),MCha2(2),MFu(3),             & 
& MFu2(3),MFd(3),MFd2(3),MFe(3),MFe2(3),Mhh(2),Mhh2(2),MVWm,MVWm2,MVZ,MVZ2,              & 
& MSu(6),MSu2(6),MSd(6),MSd2(6),MSv(3),MSv2(3),MSe(6),MSe2(6)

Complex(dp), Intent(in) :: cplAhHpmcVWm(2,2),cplChiChacVWmL(4,2),cplChiChacVWmR(4,2),cplcFuFdcVWmL(3,3),         & 
& cplcFuFdcVWmR(3,3),cplFvFecVWmL(3,3),cplFvFecVWmR(3,3),cplcgWpCgAcVWm,cplcgAgWmcVWm,   & 
& cplcgZgWmcVWm,cplcgWpCgZcVWm,cplhhHpmcVWm(2,2),cplhhcVWmVWm(2),cplHpmcVWmVP(2),        & 
& cplHpmcVWmVZ(2),cplSdcSucVWm(6,6),cplSecSvcVWm(6,3),cplcVWmVPVWm,cplcVWmVWmVZ,         & 
& cplAhAhcVWmVWm(2,2),cplhhhhcVWmVWm(2,2),cplHpmcHpmcVWmVWm(2,2),cplSdcSdcVWmVWm(6,6),   & 
& cplSecSecVWmVWm(6,6),cplSucSucVWmVWm(6,6),cplSvcSvcVWmVWm(3,3),cplcVWmVPVPVWm3,        & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! Hpm, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(VSSloop(p2,MAh2(i2),MHpm2(i1)),dp)  
coup1 = cplAhHpmcVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 2
 H0m2 = Real(SA_Hloop(p2,MChi2(i1),MCha2(i2)),dp) 
B0m2 = 4._dp*MChi(i1)*MCha(i2)*Real(SA_B0(p2,MChi2(i1),MCha2(i2)),dp) 
coupL1 = cplChiChacVWmL(i1,i2)
coupR1 = cplChiChacVWmR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFu2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFd(i2)*Real(SA_B0(p2,MFu2(i1),MFd2(i2)),dp) 
coupL1 = cplcFuFdcVWmL(i1,i2)
coupR1 = cplcFuFdcVWmR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,0._dp,MFe2(i2)),dp) 
B0m2 = 4._dp*0.*MFe(i2)*Real(SA_B0(p2,0._dp,MFe2(i2)),dp) 
coupL1 = cplFvFecVWmL(i1,i2)
coupR1 = cplFvFecVWmR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWmC], gP 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,0._dp,MVWm2),dp)
coup1 = cplcgWpCgAcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gP], gWm 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,MVWm2,0._dp),dp)
coup1 = cplcgAgWmcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gZ], gWm 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,MVWm2,MVZ2),dp)
coup1 = cplcgZgWmcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gZ 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(VGGloop(p2,MVZ2,MVWm2),dp)
coup1 = cplcgWpCgZcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! Hpm, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(VSSloop(p2,Mhh2(i2),MHpm2(i1)),dp)  
coup1 = cplhhHpmcVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(VVSloop(p2,MVWm2,Mhh2(i2)),dp)
coup1 = cplhhcVWmVWm(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VP, Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(VVSloop(p2,0._dp,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVP(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VZ, Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(VVSloop(p2,MVZ2,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVZ(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Su], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSd2(i2),MSu2(i1)),dp)  
coup1 = cplSdcSucVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSe2(i2),MSv2(i1)),dp)  
coup1 = cplSecSvcVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, VP 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplcVWmVPVWm
coup2 = Conjg(coup1) 
    SumI = -VVVloop(p2,MVWm2,0._dp)*coup1*coup2 
res = res +1._dp* SumI  
!------------------------ 
! VZ, VWm 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplcVWmVWmVZ
coup2 = Conjg(coup1) 
    SumI = -VVVloop(p2,MVZ2,MVWm2)*coup1*coup2 
res = res +1._dp* SumI  
!------------------------ 
! Ah 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_A0(MAh2(i1))
 coup1 = cplAhAhcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_A0(Mhh2(i1))
 coup1 = cplhhhhcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_A0(MHpm2(i1))
 coup1 = cplHpmcHpmcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSd2(i1))
 coup1 = cplSdcSdcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSe2(i1))
 coup1 = cplSecSecVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSu2(i1))
 coup1 = cplSucSucVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Sv] 
!------------------------ 
    Do i1 = 1, 3
 SumI = 0._dp 
 A0m2 = SA_A0(MSv2(i1))
 coup1 = cplSvcSvcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! VP 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_A0(0._dp) +RXi/4._dp*SA_A0(0._dp*RXi) 
coup1 = cplcVWmVPVPVWm3
coup2 = cplcVWmVPVPVWm1
coup3 = cplcVWmVPVPVWm2
SumI = ((2._dp*rMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*0._dp-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1._dp/2._dp* SumI  
!------------------------ 
! conj[VWm] 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_A0(MVWm2) +RXi/4._dp*SA_A0(MVWm2*RXi) 
coup1 = cplcVWmcVWmVWmVWm2
coup2 = cplcVWmcVWmVWmVWm3
coup3 = cplcVWmcVWmVWmVWm1
SumI = ((2._dp*rMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
!------------------------ 
! VZ 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_A0(MVZ2) +RXi/4._dp*SA_A0(MVZ2*RXi) 
coup1 = cplcVWmVWmVZVZ1
coup2 = cplcVWmVWmVZVZ2
coup3 = cplcVWmVWmVZVZ3
SumI = ((2._dp*rMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVZ2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1._dp/2._dp* SumI  
res = oo16pi2*res 
 
End Subroutine Pi1LoopVWm 
 
Subroutine DerPi1LoopVWm(p2,MHpm,MHpm2,MAh,MAh2,MChi,MChi2,MCha,MCha2,MFu,            & 
& MFu2,MFd,MFd2,MFe,MFe2,Mhh,Mhh2,MVWm,MVWm2,MVZ,MVZ2,MSu,MSu2,MSd,MSd2,MSv,             & 
& MSv2,MSe,MSe2,cplAhHpmcVWm,cplChiChacVWmL,cplChiChacVWmR,cplcFuFdcVWmL,cplcFuFdcVWmR,  & 
& cplFvFecVWmL,cplFvFecVWmR,cplcgWpCgAcVWm,cplcgAgWmcVWm,cplcgZgWmcVWm,cplcgWpCgZcVWm,   & 
& cplhhHpmcVWm,cplhhcVWmVWm,cplHpmcVWmVP,cplHpmcVWmVZ,cplSdcSucVWm,cplSecSvcVWm,         & 
& cplcVWmVPVWm,cplcVWmVWmVZ,cplAhAhcVWmVWm,cplhhhhcVWmVWm,cplHpmcHpmcVWmVWm,             & 
& cplSdcSdcVWmVWm,cplSecSecVWmVWm,cplSucSucVWmVWm,cplSvcSvcVWmVWm,cplcVWmVPVPVWm3,       & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3,kont,res)

Implicit None 
Real(dp), Intent(in) :: MHpm(2),MHpm2(2),MAh(2),MAh2(2),MChi(4),MChi2(4),MCha(2),MCha2(2),MFu(3),             & 
& MFu2(3),MFd(3),MFd2(3),MFe(3),MFe2(3),Mhh(2),Mhh2(2),MVWm,MVWm2,MVZ,MVZ2,              & 
& MSu(6),MSu2(6),MSd(6),MSd2(6),MSv(3),MSv2(3),MSe(6),MSe2(6)

Complex(dp), Intent(in) :: cplAhHpmcVWm(2,2),cplChiChacVWmL(4,2),cplChiChacVWmR(4,2),cplcFuFdcVWmL(3,3),         & 
& cplcFuFdcVWmR(3,3),cplFvFecVWmL(3,3),cplFvFecVWmR(3,3),cplcgWpCgAcVWm,cplcgAgWmcVWm,   & 
& cplcgZgWmcVWm,cplcgWpCgZcVWm,cplhhHpmcVWm(2,2),cplhhcVWmVWm(2),cplHpmcVWmVP(2),        & 
& cplHpmcVWmVZ(2),cplSdcSucVWm(6,6),cplSecSvcVWm(6,3),cplcVWmVPVWm,cplcVWmVWmVZ,         & 
& cplAhAhcVWmVWm(2,2),cplhhhhcVWmVWm(2,2),cplHpmcHpmcVWmVWm(2,2),cplSdcSdcVWmVWm(6,6),   & 
& cplSecSecVWmVWm(6,6),cplSucSucVWmVWm(6,6),cplSvcSvcVWmVWm(3,3),cplcVWmVPVPVWm3,        & 
& cplcVWmVPVPVWm1,cplcVWmVPVPVWm2,cplcVWmcVWmVWmVWm2,cplcVWmcVWmVWmVWm3,cplcVWmcVWmVWmVWm1,& 
& cplcVWmVWmVZVZ1,cplcVWmVWmVZVZ2,cplcVWmVWmVZVZ3

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2, A0m12, A0m22 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! Hpm, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(DerVSSloop(p2,MAh2(i2),MHpm2(i1)),dp)  
coup1 = cplAhHpmcVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 2
 H0m2 = Real(SA_DerHloop(p2,MChi2(i1),MCha2(i2)),dp) 
B0m2 = 4._dp*MChi(i1)*MCha(i2)*Real(SA_DerB0(p2,MChi2(i1),MCha2(i2)),dp) 
coupL1 = cplChiChacVWmL(i1,i2)
coupR1 = cplChiChacVWmR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFu2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFd(i2)*Real(SA_DerB0(p2,MFu2(i1),MFd2(i2)),dp) 
coupL1 = cplcFuFdcVWmL(i1,i2)
coupR1 = cplcFuFdcVWmR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Fv, Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,0._dp,MFe2(i2)),dp) 
B0m2 = 4._dp*0.*MFe(i2)*Real(SA_DerB0(p2,0._dp,MFe2(i2)),dp) 
coupL1 = cplFvFecVWmL(i1,i2)
coupR1 = cplFvFecVWmR(i1,i2)
    SumI = (Abs(coupL1)**2+Abs(coupR1)**2)*H0m2 & 
                & + (Real(Conjg(coupL1)*coupR1,dp))*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWmC], gP 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVP2,MVWm2),dp)
coup1 = cplcgWpCgAcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gP], gWm 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVP2),dp)
coup1 = cplcgAgWmcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gZ], gWm 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVZ2),dp)
coup1 = cplcgZgWmcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gZ 
!------------------------ 
sumI = 0._dp 
 
SumI = 0._dp 
B0m2 = Real(DerVGGloop(p2,MVZ2,MVWm2),dp)
coup1 = cplcgWpCgZcVWm
coup2 = Conjg(coup1) 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! Hpm, hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(DerVSSloop(p2,Mhh2(i2),MHpm2(i1)),dp)  
coup1 = cplhhHpmcVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(DerVVSloop(p2,MVWm2,Mhh2(i2)),dp)
coup1 = cplhhcVWmVWm(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VP, Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(DerVVSloop(p2,MVP2,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVP(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! VZ, Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(DerVVSloop(p2,MVZ2,MHpm2(i2)),dp)
coup1 = cplHpmcVWmVZ(i2)
    SumI = Abs(coup1)**2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Su], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSd2(i2),MSu2(i1)),dp)  
coup1 = cplSdcSucVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSe2(i2),MSv2(i1)),dp)  
coup1 = cplSecSvcVWm(i2,i1)
    SumI = Abs(coup1)**2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! VWm, VP 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplcVWmVPVWm
coup2 = Conjg(coup1) 
    SumI = -DerVVVloop(p2,MVWm2,MVP2)*coup1*coup2 
res = res +1._dp* SumI  
!------------------------ 
! VZ, VWm 
!------------------------ 
sumI = 0._dp 
 
coup1 = cplcVWmVWmVZ
coup2 = Conjg(coup1) 
    SumI = -DerVVVloop(p2,MVZ2,MVWm2)*coup1*coup2 
res = res +1._dp* SumI  
!------------------------ 
! Ah 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_DerA0(MAh2(i1))
 coup1 = cplAhAhcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! hh 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_DerA0(Mhh2(i1))
 coup1 = cplhhhhcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1._dp/2._dp* SumI  
      End Do 
 !------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_DerA0(MHpm2(i1))
 coup1 = cplHpmcHpmcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSd2(i1))
 coup1 = cplSdcSdcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSe2(i1))
 coup1 = cplSecSecVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSu2(i1))
 coup1 = cplSucSucVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Sv] 
!------------------------ 
    Do i1 = 1, 3
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSv2(i1))
 coup1 = cplSvcSvcVWmVWm(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! VP 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_DerA0(MVP2) +RXi/4._dp*SA_DerA0(MVP2*RXi) 
coup1 = cplcVWmVPVPVWm3
coup2 = cplcVWmVPVPVWm1
coup3 = cplcVWmVPVPVWm2
SumI = ((2._dp*DerrMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVP2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1._dp/2._dp* SumI  
!------------------------ 
! conj[VWm] 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_DerA0(MVWm2) +RXi/4._dp*SA_DerA0(MVWm2*RXi) 
coup1 = cplcVWmcVWmVWmVWm2
coup2 = cplcVWmcVWmVWmVWm3
coup3 = cplcVWmcVWmVWmVWm1
SumI = ((2._dp*DerrMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
!------------------------ 
! VZ 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_DerA0(MVZ2) +RXi/4._dp*SA_DerA0(MVZ2*RXi) 
coup1 = cplcVWmVWmVZVZ1
coup2 = cplcVWmVWmVZVZ2
coup3 = cplcVWmVWmVZVZ3
SumI = ((2._dp*DerrMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVZ2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1._dp/2._dp* SumI  
res = oo16pi2*res 
 
End Subroutine DerPi1LoopVWm 
 
Subroutine Sigma1LoopFeMZ(p2,MFe,MFe2,MAh,MAh2,MSv,MSv2,MCha,MCha2,MSe,               & 
& MSe2,MChi,MChi2,MSu,MSu2,MFd,MFd2,Mhh,Mhh2,MVZ,MVZ2,MHpm,MHpm2,MVWm,MVWm2,             & 
& MFu,MFu2,MSd,MSd2,cplcUFeFeAhL,cplcUFeFeAhR,cplcUFeChaSvL,cplcUFeChaSvR,               & 
& cplcUFeChiSeL,cplcUFeChiSeR,cplcUFeFdcSuL,cplcUFeFdcSuR,cplcUFeFehhL,cplcUFeFehhR,     & 
& cplcUFeFecSvL,cplcUFeFecSvR,cplcUFeFeVPL,cplcUFeFeVPR,cplcUFeFeVZL,cplcUFeFeVZR,       & 
& cplcUFeFvHpmL,cplcUFeFvHpmR,cplcUFeFvSeL,cplcUFeFvSeR,cplcUFeFvVWmL,cplcUFeFvVWmR,     & 
& cplcUFecFuSdL,cplcUFecFuSdR,sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MFe(3),MFe2(3),MAh(2),MAh2(2),MSv(3),MSv2(3),MCha(2),MCha2(2),MSe(6),MSe2(6),         & 
& MChi(4),MChi2(4),MSu(6),MSu2(6),MFd(3),MFd2(3),Mhh(2),Mhh2(2),MVZ,MVZ2,MHpm(2),        & 
& MHpm2(2),MVWm,MVWm2,MFu(3),MFu2(3),MSd(6),MSd2(6)

Complex(dp), Intent(in) :: cplcUFeFeAhL(3,3,2),cplcUFeFeAhR(3,3,2),cplcUFeChaSvL(3,2,3),cplcUFeChaSvR(3,2,3),    & 
& cplcUFeChiSeL(3,4,6),cplcUFeChiSeR(3,4,6),cplcUFeFdcSuL(3,3,6),cplcUFeFdcSuR(3,3,6),   & 
& cplcUFeFehhL(3,3,2),cplcUFeFehhR(3,3,2),cplcUFeFecSvL(3,3,3),cplcUFeFecSvR(3,3,3),     & 
& cplcUFeFeVPL(3,3),cplcUFeFeVPR(3,3),cplcUFeFeVZL(3,3),cplcUFeFeVZR(3,3),               & 
& cplcUFeFvHpmL(3,3,2),cplcUFeFvHpmR(3,3,2),cplcUFeFvSeL(3,3,6),cplcUFeFvSeR(3,3,6),     & 
& cplcUFeFvVWmL(3,3),cplcUFeFvVWmR(3,3),cplcUFecFuSdL(3,3,6),cplcUFecFuSdR(3,3,6)

Complex(dp), Intent(out) :: SigL(3,3),SigR(3,3), SigSL(3,3), SigSR(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(3,3), sumR(3,3), sumSL(3,3), sumSR(3,3) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
 SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Fe, Ah 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),MFe2(i1),MAh2(i2)),dp) 
B0m2 = MFe(i1)*Real(SA_B0(MFe2(gO1),MFe2(i1),MAh2(i2)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFe2(i1),MAh2(i2)),dp) 
B0m2 = MFe(i1)*Real(SA_B0(p2,MFe2(i1),MAh2(i2)),dp) 
End If 
coupL1 = cplcUFeFeAhL(gO1,i1,i2)
coupR1 = cplcUFeFeAhR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFeFeAhL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFeFeAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sv, Cha 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),MCha2(i2),MSv2(i1)),dp) 
B0m2 = MCha(i2)*Real(SA_B0(MFe2(gO1),MCha2(i2),MSv2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MCha2(i2),MSv2(i1)),dp) 
B0m2 = MCha(i2)*Real(SA_B0(p2,MCha2(i2),MSv2(i1)),dp) 
End If 
coupL1 = cplcUFeChaSvL(gO1,i2,i1)
coupR1 = cplcUFeChaSvR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeChaSvL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeChaSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Se, Chi 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),MChi2(i2),MSe2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(MFe2(gO1),MChi2(i2),MSe2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MChi2(i2),MSe2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSe2(i1)),dp) 
End If 
coupL1 = cplcUFeChiSeL(gO1,i2,i1)
coupR1 = cplcUFeChiSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeChiSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeChiSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),MFd2(i2),MSu2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(MFe2(gO1),MFd2(i2),MSu2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFd2(i2),MSu2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSu2(i1)),dp) 
End If 
coupL1 = cplcUFeFdcSuL(gO1,i2,i1)
coupR1 = cplcUFeFdcSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFdcSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFdcSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp* sumL
SigR = SigR +3._dp* sumR 
SigSL = SigSL +3._dp* sumSL 
SigSR = SigSR +3._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! hh, Fe 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),MFe2(i2),Mhh2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(MFe2(gO1),MFe2(i2),Mhh2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFe2(i2),Mhh2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(p2,MFe2(i2),Mhh2(i1)),dp) 
End If 
coupL1 = cplcUFeFehhL(gO1,i2,i1)
coupR1 = cplcUFeFehhR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFehhL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFehhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Fe 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),MFe2(i2),MSv2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(MFe2(gO1),MFe2(i2),MSv2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFe2(i2),MSv2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(p2,MFe2(i2),MSv2(i1)),dp) 
End If 
coupL1 = cplcUFeFecSvL(gO1,i2,i1)
coupR1 = cplcUFeFecSvR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFecSvL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFecSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VZ, Fe 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -2._dp*Real(SA_B1(MFe2(gO1),MFe2(i2),MVZ2),dp) 
B0m2 = -4._dp*MFe(i2)*Real(SA_B0(MFe2(gO1),MFe2(i2),MVZ2)-0.5_dp*rMS,dp) 
Else 
B1m2 = -2._dp*Real(SA_B1(p2,MFe2(i2),MVZ2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFe(i2)*Real(SA_B0(p2,MFe2(i2),MVZ2)-0.5_dp*rMS,dp) 
End If 
coupL1 = cplcUFeFeVZL(gO1,i2)
coupR1 = cplcUFeFeVZR(gO1,i2)
coupL2 =  Conjg(cplcUFeFeVZL(gO2,i2))
coupR2 =  Conjg(cplcUFeFeVZR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Hpm, Fv 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),0._dp,MHpm2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(MFe2(gO1),0._dp,MHpm2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,0._dp,MHpm2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(p2,0._dp,MHpm2(i1)),dp) 
End If 
coupL1 = cplcUFeFvHpmL(gO1,i2,i1)
coupR1 = cplcUFeFvHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFvHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFvHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Se, Fv 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),0._dp,MSe2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(MFe2(gO1),0._dp,MSe2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,0._dp,MSe2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(p2,0._dp,MSe2(i1)),dp) 
End If 
coupL1 = cplcUFeFvSeL(gO1,i2,i1)
coupR1 = cplcUFeFvSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFeFvSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFeFvSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VWm, Fv 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -2._dp*Real(SA_B1(MFe2(gO1),0._dp,MVWm2),dp) 
B0m2 = -4._dp*0.*Real(SA_B0(MFe2(gO1),0._dp,MVWm2)-0.5_dp*rMS,dp) 
Else 
B1m2 = -2._dp*Real(SA_B1(p2,0._dp,MVWm2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*0.*Real(SA_B0(p2,0._dp,MVWm2)-0.5_dp*rMS,dp) 
End If 
coupL1 = cplcUFeFvVWmL(gO1,i2)
coupR1 = cplcUFeFvVWmR(gO1,i2)
coupL2 =  Conjg(cplcUFeFvVWmL(gO2,i2))
coupR2 =  Conjg(cplcUFeFvVWmR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! bar[Fu], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFe2(gO1),MFu2(i1),MSd2(i2)),dp) 
B0m2 = MFu(i1)*Real(SA_B0(MFe2(gO1),MFu2(i1),MSd2(i2)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFu2(i1),MSd2(i2)),dp) 
B0m2 = MFu(i1)*Real(SA_B0(p2,MFu2(i1),MSd2(i2)),dp) 
End If 
coupL1 = cplcUFecFuSdL(gO1,i1,i2)
coupR1 = cplcUFecFuSdR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFecFuSdL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFecFuSdR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +3._dp* sumL
SigR = SigR +3._dp* sumR 
SigSL = SigSL +3._dp* sumSL 
SigSR = SigSR +3._dp* sumSR 
      End Do 
     End Do 
 SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopFeMZ 
 
Subroutine Sigma1LoopFdMZ(p2,MFd,MFd2,MAh,MAh2,MSu,MSu2,MCha,MCha2,MSd,               & 
& MSd2,MChi,MChi2,Mhh,Mhh2,MSv,MSv2,MVZ,MVZ2,MFe,MFe2,MHpm,MHpm2,MFu,MFu2,               & 
& MSe,MSe2,MVWm,MVWm2,MGlu,MGlu2,cplcUFdFdAhL,cplcUFdFdAhR,cplcUFdChaSuL,cplcUFdChaSuR,  & 
& cplcUFdChiSdL,cplcUFdChiSdR,cplcUFdFdhhL,cplcUFdFdhhR,cplcUFdFdcSvL,cplcUFdFdcSvR,     & 
& cplcUFdFdVGL,cplcUFdFdVGR,cplcUFdFdVPL,cplcUFdFdVPR,cplcUFdFdVZL,cplcUFdFdVZR,         & 
& cplcUFdFeSuL,cplcUFdFeSuR,cplcUFdFuHpmL,cplcUFdFuHpmR,cplcUFdFuSeL,cplcUFdFuSeR,       & 
& cplcUFdFuVWmL,cplcUFdFuVWmR,cplcUFdFvSdL,cplcUFdFvSdR,cplcUFdGluSdL,cplcUFdGluSdR,     & 
& sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MAh(2),MAh2(2),MSu(6),MSu2(6),MCha(2),MCha2(2),MSd(6),MSd2(6),         & 
& MChi(4),MChi2(4),Mhh(2),Mhh2(2),MSv(3),MSv2(3),MVZ,MVZ2,MFe(3),MFe2(3),MHpm(2),        & 
& MHpm2(2),MFu(3),MFu2(3),MSe(6),MSe2(6),MVWm,MVWm2,MGlu,MGlu2

Complex(dp), Intent(in) :: cplcUFdFdAhL(3,3,2),cplcUFdFdAhR(3,3,2),cplcUFdChaSuL(3,2,6),cplcUFdChaSuR(3,2,6),    & 
& cplcUFdChiSdL(3,4,6),cplcUFdChiSdR(3,4,6),cplcUFdFdhhL(3,3,2),cplcUFdFdhhR(3,3,2),     & 
& cplcUFdFdcSvL(3,3,3),cplcUFdFdcSvR(3,3,3),cplcUFdFdVGL(3,3),cplcUFdFdVGR(3,3),         & 
& cplcUFdFdVPL(3,3),cplcUFdFdVPR(3,3),cplcUFdFdVZL(3,3),cplcUFdFdVZR(3,3),               & 
& cplcUFdFeSuL(3,3,6),cplcUFdFeSuR(3,3,6),cplcUFdFuHpmL(3,3,2),cplcUFdFuHpmR(3,3,2),     & 
& cplcUFdFuSeL(3,3,6),cplcUFdFuSeR(3,3,6),cplcUFdFuVWmL(3,3),cplcUFdFuVWmR(3,3),         & 
& cplcUFdFvSdL(3,3,6),cplcUFdFvSdR(3,3,6),cplcUFdGluSdL(3,6),cplcUFdGluSdR(3,6)

Complex(dp), Intent(out) :: SigL(3,3),SigR(3,3), SigSL(3,3), SigSR(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(3,3), sumR(3,3), sumSL(3,3), sumSR(3,3) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
 SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Fd, Ah 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MFd2(i1),MAh2(i2)),dp) 
B0m2 = MFd(i1)*Real(SA_B0(MFd2(gO1),MFd2(i1),MAh2(i2)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFd2(i1),MAh2(i2)),dp) 
B0m2 = MFd(i1)*Real(SA_B0(p2,MFd2(i1),MAh2(i2)),dp) 
End If 
coupL1 = cplcUFdFdAhL(gO1,i1,i2)
coupR1 = cplcUFdFdAhR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFdFdAhL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFdFdAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Su, Cha 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MCha2(i2),MSu2(i1)),dp) 
B0m2 = MCha(i2)*Real(SA_B0(MFd2(gO1),MCha2(i2),MSu2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MCha2(i2),MSu2(i1)),dp) 
B0m2 = MCha(i2)*Real(SA_B0(p2,MCha2(i2),MSu2(i1)),dp) 
End If 
coupL1 = cplcUFdChaSuL(gO1,i2,i1)
coupR1 = cplcUFdChaSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdChaSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdChaSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sd, Chi 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MChi2(i2),MSd2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(MFd2(gO1),MChi2(i2),MSd2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MChi2(i2),MSd2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSd2(i1)),dp) 
End If 
coupL1 = cplcUFdChiSdL(gO1,i2,i1)
coupR1 = cplcUFdChiSdR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdChiSdL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdChiSdR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! hh, Fd 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MFd2(i2),Mhh2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(MFd2(gO1),MFd2(i2),Mhh2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFd2(i2),Mhh2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),Mhh2(i1)),dp) 
End If 
coupL1 = cplcUFdFdhhL(gO1,i2,i1)
coupR1 = cplcUFdFdhhR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFdhhL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFdhhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Fd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MFd2(i2),MSv2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(MFd2(gO1),MFd2(i2),MSv2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFd2(i2),MSv2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSv2(i1)),dp) 
End If 
coupL1 = cplcUFdFdcSvL(gO1,i2,i1)
coupR1 = cplcUFdFdcSvR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFdcSvL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFdcSvR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VZ, Fd 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -2._dp*Real(SA_B1(MFd2(gO1),MFd2(i2),MVZ2),dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(MFd2(gO1),MFd2(i2),MVZ2)-0.5_dp*rMS,dp) 
Else 
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),MVZ2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),MVZ2)-0.5_dp*rMS,dp) 
End If 
coupL1 = cplcUFdFdVZL(gO1,i2)
coupR1 = cplcUFdFdVZR(gO1,i2)
coupL2 =  Conjg(cplcUFdFdVZL(gO2,i2))
coupR2 =  Conjg(cplcUFdFdVZR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Su, Fe 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MFe2(i2),MSu2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(MFd2(gO1),MFe2(i2),MSu2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFe2(i2),MSu2(i1)),dp) 
B0m2 = MFe(i2)*Real(SA_B0(p2,MFe2(i2),MSu2(i1)),dp) 
End If 
coupL1 = cplcUFdFeSuL(gO1,i2,i1)
coupR1 = cplcUFdFeSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFeSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFeSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Hpm, Fu 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MFu2(i2),MHpm2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(MFd2(gO1),MFu2(i2),MHpm2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFu2(i2),MHpm2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(p2,MFu2(i2),MHpm2(i1)),dp) 
End If 
coupL1 = cplcUFdFuHpmL(gO1,i2,i1)
coupR1 = cplcUFdFuHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFuHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFuHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Se, Fu 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MFu2(i2),MSe2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(MFd2(gO1),MFu2(i2),MSe2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFu2(i2),MSe2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(p2,MFu2(i2),MSe2(i1)),dp) 
End If 
coupL1 = cplcUFdFuSeL(gO1,i2,i1)
coupR1 = cplcUFdFuSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFuSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFuSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VWm, Fu 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -2._dp*Real(SA_B1(MFd2(gO1),MFu2(i2),MVWm2),dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(MFd2(gO1),MFu2(i2),MVWm2)-0.5_dp*rMS,dp) 
Else 
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),MVWm2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),MVWm2)-0.5_dp*rMS,dp) 
End If 
coupL1 = cplcUFdFuVWmL(gO1,i2)
coupR1 = cplcUFdFuVWmR(gO1,i2)
coupL2 =  Conjg(cplcUFdFuVWmL(gO2,i2))
coupR2 =  Conjg(cplcUFdFuVWmR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Sd, Fv 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),0._dp,MSd2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(MFd2(gO1),0._dp,MSd2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,0._dp,MSd2(i1)),dp) 
B0m2 = 0.*Real(SA_B0(p2,0._dp,MSd2(i1)),dp) 
End If 
coupL1 = cplcUFdFvSdL(gO1,i2,i1)
coupR1 = cplcUFdFvSdR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFdFvSdL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFdFvSdR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Sd, Glu 
!------------------------ 
    Do i1 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFd2(gO1),MGlu2,MSd2(i1)),dp) 
B0m2 = MGlu*Real(SA_B0(MFd2(gO1),MGlu2,MSd2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MGlu2,MSd2(i1)),dp) 
B0m2 = MGlu*Real(SA_B0(p2,MGlu2,MSd2(i1)),dp) 
End If 
coupL1 = cplcUFdGluSdL(gO1,i1)
coupR1 = cplcUFdGluSdR(gO1,i1)
coupL2 =  Conjg(cplcUFdGluSdL(gO2,i1))
coupR2 =  Conjg(cplcUFdGluSdR(gO2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +4._dp/3._dp* sumL
SigR = SigR +4._dp/3._dp* sumR 
SigSL = SigSL +4._dp/3._dp* sumSL 
SigSR = SigSR +4._dp/3._dp* sumSR 
      End Do 
 SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopFdMZ 
 
Subroutine Sigma1LoopFuMZ(p2,MFu,MFu2,MAh,MAh2,MSu,MSu2,MChi,MChi2,MHpm,              & 
& MHpm2,MFd,MFd2,MSe,MSe2,MVWm,MVWm2,Mhh,Mhh2,MVZ,MVZ2,MGlu,MGlu2,MCha,MCha2,            & 
& MSd,MSd2,MFe,MFe2,cplcUFuFuAhL,cplcUFuFuAhR,cplcUFuChiSuL,cplcUFuChiSuR,               & 
& cplcUFuFdcHpmL,cplcUFuFdcHpmR,cplcUFuFdcSeL,cplcUFuFdcSeR,cplcUFuFdcVWmL,              & 
& cplcUFuFdcVWmR,cplcUFuFuhhL,cplcUFuFuhhR,cplcUFuFuVGL,cplcUFuFuVGR,cplcUFuFuVPL,       & 
& cplcUFuFuVPR,cplcUFuFuVZL,cplcUFuFuVZR,cplcUFuGluSuL,cplcUFuGluSuR,cplcChacUFuSdL,     & 
& cplcChacUFuSdR,cplcFecUFuSdL,cplcFecUFuSdR,sigL,sigR,sigSL,sigSR)

Implicit None 
Real(dp), Intent(in) :: MFu(3),MFu2(3),MAh(2),MAh2(2),MSu(6),MSu2(6),MChi(4),MChi2(4),MHpm(2),MHpm2(2),       & 
& MFd(3),MFd2(3),MSe(6),MSe2(6),MVWm,MVWm2,Mhh(2),Mhh2(2),MVZ,MVZ2,MGlu,MGlu2,           & 
& MCha(2),MCha2(2),MSd(6),MSd2(6),MFe(3),MFe2(3)

Complex(dp), Intent(in) :: cplcUFuFuAhL(3,3,2),cplcUFuFuAhR(3,3,2),cplcUFuChiSuL(3,4,6),cplcUFuChiSuR(3,4,6),    & 
& cplcUFuFdcHpmL(3,3,2),cplcUFuFdcHpmR(3,3,2),cplcUFuFdcSeL(3,3,6),cplcUFuFdcSeR(3,3,6), & 
& cplcUFuFdcVWmL(3,3),cplcUFuFdcVWmR(3,3),cplcUFuFuhhL(3,3,2),cplcUFuFuhhR(3,3,2),       & 
& cplcUFuFuVGL(3,3),cplcUFuFuVGR(3,3),cplcUFuFuVPL(3,3),cplcUFuFuVPR(3,3),               & 
& cplcUFuFuVZL(3,3),cplcUFuFuVZR(3,3),cplcUFuGluSuL(3,6),cplcUFuGluSuR(3,6),             & 
& cplcChacUFuSdL(2,3,6),cplcChacUFuSdR(2,3,6),cplcFecUFuSdL(3,3,6),cplcFecUFuSdR(3,3,6)

Complex(dp), Intent(out) :: SigL(3,3),SigR(3,3), SigSL(3,3), SigSR(3,3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2,temp, sumL(3,3), sumR(3,3), sumSL(3,3), sumSR(3,3) 
Real(dp) :: B0m2, F0m2, G0m2,B1m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
SigL = Cmplx(0._dp,0._dp,dp) 
SigR = Cmplx(0._dp,0._dp,dp) 
SigSL = Cmplx(0._dp,0._dp,dp) 
 SigSR = Cmplx(0._dp,0._dp,dp) 
 
!------------------------ 
! Fu, Ah 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 2
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MFu2(i1),MAh2(i2)),dp) 
B0m2 = MFu(i1)*Real(SA_B0(MFu2(gO1),MFu2(i1),MAh2(i2)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFu2(i1),MAh2(i2)),dp) 
B0m2 = MFu(i1)*Real(SA_B0(p2,MFu2(i1),MAh2(i2)),dp) 
End If 
coupL1 = cplcUFuFuAhL(gO1,i1,i2)
coupR1 = cplcUFuFuAhR(gO1,i1,i2)
coupL2 =  Conjg(cplcUFuFuAhL(gO2,i1,i2))
coupR2 =  Conjg(cplcUFuFuAhR(gO2,i1,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! Su, Chi 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 4
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MChi2(i2),MSu2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(MFu2(gO1),MChi2(i2),MSu2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MChi2(i2),MSu2(i1)),dp) 
B0m2 = MChi(i2)*Real(SA_B0(p2,MChi2(i2),MSu2(i1)),dp) 
End If 
coupL1 = cplcUFuChiSuL(gO1,i2,i1)
coupR1 = cplcUFuChiSuR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuChiSuL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuChiSuR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], Fd 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MFd2(i2),MHpm2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(MFu2(gO1),MFd2(i2),MHpm2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFd2(i2),MHpm2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MHpm2(i1)),dp) 
End If 
coupL1 = cplcUFuFdcHpmL(gO1,i2,i1)
coupR1 = cplcUFuFdcHpmR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuFdcHpmL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuFdcHpmR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Fd 
!------------------------ 
    Do i1 = 1, 6
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MFd2(i2),MSe2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(MFu2(gO1),MFd2(i2),MSe2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFd2(i2),MSe2(i1)),dp) 
B0m2 = MFd(i2)*Real(SA_B0(p2,MFd2(i2),MSe2(i1)),dp) 
End If 
coupL1 = cplcUFuFdcSeL(gO1,i2,i1)
coupR1 = cplcUFuFdcSeR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuFdcSeL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuFdcSeR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Fd 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -2._dp*Real(SA_B1(MFu2(gO1),MFd2(i2),MVWm2),dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(MFu2(gO1),MFd2(i2),MVWm2)-0.5_dp*rMS,dp) 
Else 
B1m2 = -2._dp*Real(SA_B1(p2,MFd2(i2),MVWm2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFd(i2)*Real(SA_B0(p2,MFd2(i2),MVWm2)-0.5_dp*rMS,dp) 
End If 
coupL1 = cplcUFuFdcVWmL(gO1,i2)
coupR1 = cplcUFuFdcVWmR(gO1,i2)
coupL2 =  Conjg(cplcUFuFdcVWmL(gO2,i2))
coupR2 =  Conjg(cplcUFuFdcVWmR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! hh, Fu 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MFu2(i2),Mhh2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(MFu2(gO1),MFu2(i2),Mhh2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFu2(i2),Mhh2(i1)),dp) 
B0m2 = MFu(i2)*Real(SA_B0(p2,MFu2(i2),Mhh2(i1)),dp) 
End If 
coupL1 = cplcUFuFuhhL(gO1,i2,i1)
coupR1 = cplcUFuFuhhR(gO1,i2,i1)
coupL2 =  Conjg(cplcUFuFuhhL(gO2,i2,i1))
coupR2 =  Conjg(cplcUFuFuhhR(gO2,i2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! VZ, Fu 
!------------------------ 
      Do i2 = 1, 3
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -2._dp*Real(SA_B1(MFu2(gO1),MFu2(i2),MVZ2),dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(MFu2(gO1),MFu2(i2),MVZ2)-0.5_dp*rMS,dp) 
Else 
B1m2 = -2._dp*Real(SA_B1(p2,MFu2(i2),MVZ2)+ 0.5_dp*rMS,dp) 
B0m2 = -4._dp*MFu(i2)*Real(SA_B0(p2,MFu2(i2),MVZ2)-0.5_dp*rMS,dp) 
End If 
coupL1 = cplcUFuFuVZL(gO1,i2)
coupR1 = cplcUFuFuVZR(gO1,i2)
coupL2 =  Conjg(cplcUFuFuVZL(gO2,i2))
coupR2 =  Conjg(cplcUFuFuVZR(gO2,i2))
SumSL(gO1,gO2) = coupL1*coupR2*B0m2 
SumSR(gO1,gO2) = coupR1*coupL2*B0m2 
sumR(gO1,gO2) = coupL1*coupL2*B1m2 
sumL(gO1,gO2) = coupR1*coupR2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
    End Do 
 !------------------------ 
! Su, Glu 
!------------------------ 
    Do i1 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MGlu2,MSu2(i1)),dp) 
B0m2 = MGlu*Real(SA_B0(MFu2(gO1),MGlu2,MSu2(i1)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MGlu2,MSu2(i1)),dp) 
B0m2 = MGlu*Real(SA_B0(p2,MGlu2,MSu2(i1)),dp) 
End If 
coupL1 = cplcUFuGluSuL(gO1,i1)
coupR1 = cplcUFuGluSuR(gO1,i1)
coupL2 =  Conjg(cplcUFuGluSuL(gO2,i1))
coupR2 =  Conjg(cplcUFuGluSuR(gO2,i1))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +4._dp/3._dp* sumL
SigR = SigR +4._dp/3._dp* sumR 
SigSL = SigSL +4._dp/3._dp* sumSL 
SigSR = SigSR +4._dp/3._dp* sumSR 
      End Do 
 !------------------------ 
! bar[Cha], Sd 
!------------------------ 
    Do i1 = 1, 2
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MCha2(i1),MSd2(i2)),dp) 
B0m2 = MCha(i1)*Real(SA_B0(MFu2(gO1),MCha2(i1),MSd2(i2)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MCha2(i1),MSd2(i2)),dp) 
B0m2 = MCha(i1)*Real(SA_B0(p2,MCha2(i1),MSd2(i2)),dp) 
End If 
coupL1 = cplcChacUFuSdL(i1,gO1,i2)
coupR1 = cplcChacUFuSdR(i1,gO1,i2)
coupL2 =  Conjg(cplcChacUFuSdL(i1,gO2,i2))
coupR2 =  Conjg(cplcChacUFuSdR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Sd 
!------------------------ 
    Do i1 = 1, 3
       Do i2 = 1, 6
 SumSL = 0._dp 
SumSR = 0._dp 
sumR = 0._dp 
sumL = 0._dp 
Do gO1 = 1, 3
  Do gO2 = 1, 3
If(gO1.eq.gO2) Then 
B1m2 = -Real(SA_B1(MFu2(gO1),MFe2(i1),MSd2(i2)),dp) 
B0m2 = MFe(i1)*Real(SA_B0(MFu2(gO1),MFe2(i1),MSd2(i2)),dp) 
Else 
B1m2 = -Real(SA_B1(p2,MFe2(i1),MSd2(i2)),dp) 
B0m2 = MFe(i1)*Real(SA_B0(p2,MFe2(i1),MSd2(i2)),dp) 
End If 
coupL1 = cplcFecUFuSdL(i1,gO1,i2)
coupR1 = cplcFecUFuSdR(i1,gO1,i2)
coupL2 =  Conjg(cplcFecUFuSdL(i1,gO2,i2))
coupR2 =  Conjg(cplcFecUFuSdR(i1,gO2,i2))
SumSL(gO1,gO2) = coupR1*coupL2*B0m2 
SumSR(gO1,gO2) = coupL1*coupR2*B0m2 
sumR(gO1,gO2) = coupR1*coupR2*B1m2 
sumL(gO1,gO2) = coupL1*coupL2*B1m2 
   End Do 
End Do 
SigL = SigL +1._dp* sumL
SigR = SigR +1._dp* sumR 
SigSL = SigSL +1._dp* sumSL 
SigSR = SigSR +1._dp* sumSR 
      End Do 
     End Do 
 SigL = oo16pi2*SigL 
SigR = oo16pi2*SigR 
SigSL = oo16pi2*SigSL 
SigSR = oo16pi2*SigSR 
 
End Subroutine Sigma1LoopFuMZ 
 
Subroutine Pi1LoopVPVZ(p2,MCha,MCha2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MHpm,MHpm2,           & 
& MSd,MSd2,MSe,MSe2,MSu,MSu2,MVWm,MVWm2,cplcChaChaVPL,cplcChaChaVPR,cplcChaChaVZL,       & 
& cplcChaChaVZR,cplcFdFdVPL,cplcFdFdVPR,cplcFdFdVZL,cplcFdFdVZR,cplcFeFeVPL,             & 
& cplcFeFeVPR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVPL,cplcFuFuVPR,cplcFuFuVZL,               & 
& cplcFuFuVZR,cplcgWmgWmVP,cplcgWmgWmVZ,cplcgWpCgWpCVP,cplcgWpCgWpCVZ,cplcHpmVPVWm,      & 
& cplcHpmVWmVZ,cplcVWmVPVWm,cplcVWmVPVWmVZ1,cplcVWmVPVWmVZ2,cplcVWmVPVWmVZ3,             & 
& cplcVWmVWmVZ,cplHpmcHpmVP,cplHpmcHpmVPVZ,cplHpmcHpmVZ,cplHpmcVWmVP,cplHpmcVWmVZ,       & 
& cplSdcSdVP,cplSdcSdVPVZ,cplSdcSdVZ,cplSecSeVP,cplSecSeVPVZ,cplSecSeVZ,cplSucSuVP,      & 
& cplSucSuVPVZ,cplSucSuVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MHpm(2),MHpm2(2),       & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MVWm,MVWm2

Complex(dp), Intent(in) :: cplcChaChaVPL(2,2),cplcChaChaVPR(2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),          & 
& cplcFdFdVPL(3,3),cplcFdFdVPR(3,3),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFeFeVPL(3,3),  & 
& cplcFeFeVPR(3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),cplcFuFuVPL(3,3),cplcFuFuVPR(3,3),  & 
& cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcgWmgWmVP,cplcgWmgWmVZ,cplcgWpCgWpCVP,            & 
& cplcgWpCgWpCVZ,cplcHpmVPVWm(2),cplcHpmVWmVZ(2),cplcVWmVPVWm,cplcVWmVPVWmVZ1,           & 
& cplcVWmVPVWmVZ2,cplcVWmVPVWmVZ3,cplcVWmVWmVZ,cplHpmcHpmVP(2,2),cplHpmcHpmVPVZ(2,2),    & 
& cplHpmcHpmVZ(2,2),cplHpmcVWmVP(2),cplHpmcVWmVZ(2),cplSdcSdVP(6,6),cplSdcSdVPVZ(6,6),   & 
& cplSdcSdVZ(6,6),cplSecSeVP(6,6),cplSecSeVPVZ(6,6),cplSecSeVZ(6,6),cplSucSuVP(6,6),     & 
& cplSucSuVPVZ(6,6),cplSucSuVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 H0m2 = Real(SA_Hloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = 4._dp*MCha(i1)*MCha(i2)*Real(SA_B0(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVPL(i1,i2)
coupR1 = cplcChaChaVPR(i1,i2)
coupL2 = cplcChaChaVZL(i2,i1)
coupR2 = cplcChaChaVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 4._dp*MFd(i1)*MFd(i2)*Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVPL(i1,i2)
coupR1 = cplcFdFdVPR(i1,i2)
coupL2 = cplcFdFdVZL(i2,i1)
coupR2 = cplcFdFdVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = 4._dp*MFe(i1)*MFe(i2)*Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVPL(i1,i2)
coupR1 = cplcFeFeVPR(i1,i2)
coupL2 = cplcFeFeVZL(i2,i1)
coupR2 = cplcFeFeVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_Hloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 4._dp*MFu(i1)*MFu(i2)*Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVPL(i1,i2)
coupR1 = cplcFuFuVPR(i1,i2)
coupL2 = cplcFuFuVZL(i2,i1)
coupR2 = cplcFuFuVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
B0m2 = Real(VGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWmgWmVP
coup2 = cplcgWmgWmVZ 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
B0m2 = Real(VGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWpCgWpCVP
coup2 = cplcgWpCgWpCVZ 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(VSSloop(p2,MHpm2(i2),MHpm2(i1)),dp) 
coup1 = cplHpmcHpmVP(i2,i1)
coup2 = cplHpmcHpmVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(VVSloop(p2,MVWm2,MHpm2(i2)),dp) 
coup1 = cplHpmcVWmVP(i2)
coup2 = cplcHpmVWmVZ(i2)
    SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVP(i2,i1)
coup2 = cplSdcSdVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVP(i2,i1)
coup2 = cplSecSeVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(VSSloop(p2,MSu2(i2),MSu2(i1)),dp) 
coup1 = cplSucSuVP(i2,i1)
coup2 = cplSucSuVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VWm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 B0m2 = Real(VVSloop(p2,MVWm2,MHpm2(i1)),dp) 
coup1 = cplcHpmVPVWm(i1)
coup2 = cplHpmcVWmVZ(i1)
    SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
B0m2 = Real(VVVloop(p2,MVWm2,MVWm2),dp) 
coup1 = cplcVWmVPVWm
coup2 = cplcVWmVWmVZ
    SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_A0(MHpm2(i1))
 coup1 = cplHpmcHpmVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSd2(i1))
 coup1 = cplSdcSdVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSe2(i1))
 coup1 = cplSecSeVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_A0(MSu2(i1))
 coup1 = cplSucSuVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[VWm] 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_A0(MVWm2) +RXi/4._dp*SA_A0(MVWm2*RXi) 
coup1 = cplcVWmVPVWmVZ2
coup2 = cplcVWmVPVWmVZ1
coup3 = cplcVWmVPVWmVZ3
SumI = ((2._dp*rMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
res = oo16pi2*res 
 
End Subroutine Pi1LoopVPVZ 
 
Subroutine DerPi1LoopVPVZ(p2,MCha,MCha2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MHpm,              & 
& MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MVWm,MVWm2,cplcChaChaVPL,cplcChaChaVPR,               & 
& cplcChaChaVZL,cplcChaChaVZR,cplcFdFdVPL,cplcFdFdVPR,cplcFdFdVZL,cplcFdFdVZR,           & 
& cplcFeFeVPL,cplcFeFeVPR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuVPL,cplcFuFuVPR,               & 
& cplcFuFuVZL,cplcFuFuVZR,cplcgWmgWmVP,cplcgWmgWmVZ,cplcgWpCgWpCVP,cplcgWpCgWpCVZ,       & 
& cplcHpmVPVWm,cplcHpmVWmVZ,cplcVWmVPVWm,cplcVWmVPVWmVZ1,cplcVWmVPVWmVZ2,cplcVWmVPVWmVZ3,& 
& cplcVWmVWmVZ,cplHpmcHpmVP,cplHpmcHpmVPVZ,cplHpmcHpmVZ,cplHpmcVWmVP,cplHpmcVWmVZ,       & 
& cplSdcSdVP,cplSdcSdVPVZ,cplSdcSdVZ,cplSecSeVP,cplSecSeVPVZ,cplSecSeVZ,cplSucSuVP,      & 
& cplSucSuVPVZ,cplSucSuVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MHpm(2),MHpm2(2),       & 
& MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MVWm,MVWm2

Complex(dp), Intent(in) :: cplcChaChaVPL(2,2),cplcChaChaVPR(2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),          & 
& cplcFdFdVPL(3,3),cplcFdFdVPR(3,3),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFeFeVPL(3,3),  & 
& cplcFeFeVPR(3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),cplcFuFuVPL(3,3),cplcFuFuVPR(3,3),  & 
& cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcgWmgWmVP,cplcgWmgWmVZ,cplcgWpCgWpCVP,            & 
& cplcgWpCgWpCVZ,cplcHpmVPVWm(2),cplcHpmVWmVZ(2),cplcVWmVPVWm,cplcVWmVPVWmVZ1,           & 
& cplcVWmVPVWmVZ2,cplcVWmVPVWmVZ3,cplcVWmVWmVZ,cplHpmcHpmVP(2,2),cplHpmcHpmVPVZ(2,2),    & 
& cplHpmcHpmVZ(2,2),cplHpmcVWmVP(2),cplHpmcVWmVZ(2),cplSdcSdVP(6,6),cplSdcSdVPVZ(6,6),   & 
& cplSdcSdVZ(6,6),cplSecSeVP(6,6),cplSecSeVPVZ(6,6),cplSecSeVZ(6,6),cplSucSuVP(6,6),     & 
& cplSucSuVPVZ(6,6),cplSucSuVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 H0m2 = Real(SA_DerHloop(p2,MCha2(i1),MCha2(i2)),dp) 
B0m2 = 2._dp*MCha(i1)*MCha(i2)*Real(SA_DerB0(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVPL(i1,i2)
coupR1 = cplcChaChaVPR(i1,i2)
coupL2 = cplcChaChaVZL(i2,i1)
coupR2 = cplcChaChaVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFd2(i1),MFd2(i2)),dp) 
B0m2 = 2._dp*MFd(i1)*MFd(i2)*Real(SA_DerB0(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVPL(i1,i2)
coupR1 = cplcFdFdVPR(i1,i2)
coupL2 = cplcFdFdVZL(i2,i1)
coupR2 = cplcFdFdVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFe2(i1),MFe2(i2)),dp) 
B0m2 = 2._dp*MFe(i1)*MFe(i2)*Real(SA_DerB0(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVPL(i1,i2)
coupR1 = cplcFeFeVPR(i1,i2)
coupL2 = cplcFeFeVZL(i2,i1)
coupR2 = cplcFeFeVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 H0m2 = Real(SA_DerHloop(p2,MFu2(i1),MFu2(i2)),dp) 
B0m2 = 2._dp*MFu(i1)*MFu(i2)*Real(SA_DerB0(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVPL(i1,i2)
coupR1 = cplcFuFuVPR(i1,i2)
coupL2 = cplcFuFuVZL(i2,i1)
coupR2 = cplcFuFuVZR(i2,i1)
    SumI = (coupL1*coupL2+coupR1*coupR2)*H0m2 & 
                & + 0.5_dp*(coupL1*coupR2 + coupL2*coupR1)*B0m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWmgWmVP
coup2 = cplcgWmgWmVZ 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
B0m2 = Real(DerVGGloop(p2,MVWm2,MVWm2),dp)
coup1 = cplcgWpCgWpCVP
coup2 = cplcgWpCgWpCVZ 
   SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 B22m2 = Real(DerVSSloop(p2,MHpm2(i2),MHpm2(i1)),dp) 
coup1 = cplHpmcHpmVP(i2,i1)
coup2 = cplHpmcHpmVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 B0m2 = Real(DerVVSloop(p2,MVWm2,MHpm2(i2)),dp) 
coup1 = cplHpmcVWmVP(i2)
coup2 = cplcHpmVWmVZ(i2)
    SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVP(i2,i1)
coup2 = cplSdcSdVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVP(i2,i1)
coup2 = cplSecSeVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 B22m2 = Real(DerVSSloop(p2,MSu2(i2),MSu2(i1)),dp) 
coup1 = cplSucSuVP(i2,i1)
coup2 = cplSucSuVZ(i1,i2)
    SumI = coup1*coup2*B22m2 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VWm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 B0m2 = Real(DerVVSloop(p2,MVWm2,MHpm2(i1)),dp) 
coup1 = cplcHpmVPVWm(i1)
coup2 = cplHpmcVWmVZ(i1)
    SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
B0m2 = Real(DerVVVloop(p2,MVWm2,MVWm2),dp) 
coup1 = cplcVWmVPVWm
coup2 = cplcVWmVWmVZ
    SumI = coup1*coup2*B0m2 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm] 
!------------------------ 
    Do i1 = 1, 2
 SumI = 0._dp 
 A0m2 = SA_DerA0(MHpm2(i1))
 coup1 = cplHpmcHpmVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Sd] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSd2(i1))
 coup1 = cplSdcSdVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[Se] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSe2(i1))
 coup1 = cplSecSeVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +1* SumI  
      End Do 
 !------------------------ 
! conj[Su] 
!------------------------ 
    Do i1 = 1, 6
 SumI = 0._dp 
 A0m2 = SA_DerA0(MSu2(i1))
 coup1 = cplSucSuVPVZ(i1,i1)
 SumI = coup1*A0m2 
res = res +3* SumI  
      End Do 
 !------------------------ 
! conj[VWm] 
!------------------------ 
SumI = 0._dp 
A0m2 = 3._dp/4._dp*SA_DerA0(MVWm2) +RXi/4._dp*SA_DerA0(MVWm2*RXi) 
coup1 = cplcVWmVPVWmVZ2
coup2 = cplcVWmVPVWmVZ1
coup3 = cplcVWmVPVWmVZ3
SumI = ((2._dp*DerrMS*coup1+(1-RXi**2)/8._dp*(coup2+coup3))*MVWm2-(4._dp*coup1+coup2+coup3)*A0m2)
res = res +1* SumI  
res = oo16pi2*res 
 
End Subroutine DerPi1LoopVPVZ 
 
Subroutine Pi1LoopVZSv(p2,MFd,MFd2,MFe,MFe2,MSd,MSd2,MSe,MSe2,cplcFdFdcSvL,           & 
& cplcFdFdcSvR,cplcFdFdVZL,cplcFdFdVZR,cplcFeFecSvL,cplcFeFecSvR,cplcFeFeVZL,            & 
& cplcFeFeVZR,cplSdcSdcSv,cplSdcSdVZ,cplSecSecSv,cplSecSeVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MFe(3),MFe2(3),MSd(6),MSd2(6),MSe(6),MSe2(6)

Complex(dp), Intent(in) :: cplcFdFdcSvL(3,3,3),cplcFdFdcSvR(3,3,3),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),            & 
& cplcFeFecSvL(3,3,3),cplcFeFecSvR(3,3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),             & 
& cplSdcSdcSv(6,6,3),cplSdcSdVZ(6,6),cplSecSecSv(6,6,3),cplSecSeVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,3 
B0m2 = Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
coupL2 = cplcFdFdcSvL(i2,i1,gO2)
coupR2 = cplcFdFdcSvR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFd(i2)*B1m2  
End do 
res = res +1.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,3 
B0m2 = Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
coupL2 = cplcFeFecSvL(i2,i1,gO2)
coupR2 = cplcFeFecSvR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFe(i2)*B1m2  
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,3 
B0m2 = Real(SA_B0(p2,MSd2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVZ(i2,i1)
coup2 = cplSdcSdcSv(i1,i2,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,3 
B0m2 = Real(SA_B0(p2,MSe2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVZ(i2,i1)
coup2 = cplSecSecSv(i1,i2,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 res = oo16pi2*res 
 
End Subroutine Pi1LoopVZSv 
 
Subroutine DerPi1LoopVZSv(p2,MFd,MFd2,MFe,MFe2,MSd,MSd2,MSe,MSe2,cplcFdFdcSvL,        & 
& cplcFdFdcSvR,cplcFdFdVZL,cplcFdFdVZR,cplcFeFecSvL,cplcFeFecSvR,cplcFeFeVZL,            & 
& cplcFeFeVZR,cplSdcSdcSv,cplSdcSdVZ,cplSecSecSv,cplSecSeVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MFe(3),MFe2(3),MSd(6),MSd2(6),MSe(6),MSe2(6)

Complex(dp), Intent(in) :: cplcFdFdcSvL(3,3,3),cplcFdFdcSvR(3,3,3),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),            & 
& cplcFeFecSvL(3,3,3),cplcFeFecSvR(3,3,3),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),             & 
& cplSdcSdcSv(6,6,3),cplSdcSdVZ(6,6),cplSecSecSv(6,6,3),cplSecSeVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(3) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,3 
B0m2 = Real(SA_DerB0(p2,MFd2(i1),MFd2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
coupL2 = cplcFdFdcSvL(i2,i1,gO2)
coupR2 = cplcFdFdcSvR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFd(i2)*B1m2  
End do 
res = res +1.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,3 
B0m2 = Real(SA_DerB0(p2,MFe2(i1),MFe2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
coupL2 = cplcFeFecSvL(i2,i1,gO2)
coupR2 = cplcFeFecSvR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFe(i2)*B1m2  
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,3 
B0m2 = Real(SA_DerB0(p2,MSd2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVZ(i2,i1)
coup2 = cplSdcSdcSv(i1,i2,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,3 
B0m2 = Real(SA_DerB0(p2,MSe2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVZ(i2,i1)
coup2 = cplSecSecSv(i1,i2,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 res = oo16pi2*res 
 
End Subroutine DerPi1LoopVZSv 
 
Subroutine Pi1LoopVZhh(p2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,MFu2,           & 
& MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,MVWm,MVWm2,cplcChaChahhL,               & 
& cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcFdFdhhL,cplcFdFdhhR,cplcFdFdVZL,         & 
& cplcFdFdVZR,cplcFeFehhL,cplcFeFehhR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuhhL,               & 
& cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,cplcgWmgWmhh,cplcgWmgWmVZ,cplcgWpCgWpChh,          & 
& cplcgWpCgWpCVZ,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,cplcHpmVWmVZ,       & 
& cplcVWmVWmVZ,cplhhcHpmVWm,cplhhcVWmVWm,cplhhHpmcHpm,cplhhHpmcVWm,cplhhSdcSd,           & 
& cplhhSecSe,cplhhSucSu,cplhhSvcSv,cplHpmcHpmVZ,cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,      & 
& cplSucSuVZ,cplSvcSvVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),               & 
& MFu2(3),MHpm(2),MHpm2(2),MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),          & 
& MSv2(3),MVWm,MVWm2

Complex(dp), Intent(in) :: cplcChaChahhL(2,2,2),cplcChaChahhR(2,2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),      & 
& cplcFdFdhhL(3,3,2),cplcFdFdhhR(3,3,2),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),               & 
& cplcFeFehhL(3,3,2),cplcFeFehhR(3,3,2),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),               & 
& cplcFuFuhhL(3,3,2),cplcFuFuhhR(3,3,2),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),               & 
& cplcgWmgWmhh(2),cplcgWmgWmVZ,cplcgWpCgWpChh(2),cplcgWpCgWpCVZ,cplChiChihhL(4,4,2),     & 
& cplChiChihhR(4,4,2),cplChiChiVZL(4,4),cplChiChiVZR(4,4),cplcHpmVWmVZ(2),               & 
& cplcVWmVWmVZ,cplhhcHpmVWm(2,2),cplhhcVWmVWm(2),cplhhHpmcHpm(2,2,2),cplhhHpmcVWm(2,2),  & 
& cplhhSdcSd(2,6,6),cplhhSecSe(2,6,6),cplhhSucSu(2,6,6),cplhhSvcSv(2,3,3),               & 
& cplHpmcHpmVZ(2,2),cplHpmcVWmVZ(2),cplSdcSdVZ(6,6),cplSecSeVZ(6,6),cplSucSuVZ(6,6),     & 
& cplSvcSvVZ(3,3)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MCha2(i1),MCha2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVZL(i1,i2)
coupR1 = cplcChaChaVZR(i1,i2)
coupL2 = cplcChaChahhL(i2,i1,gO2)
coupR2 = cplcChaChahhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MCha(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MCha(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MChi2(i1),MChi2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MChi2(i1),MChi2(i2)),dp) 
coupL1 = cplChiChiVZL(i1,i2)
coupR1 = cplChiChiVZR(i1,i2)
coupL2 = cplChiChihhL(i1,i2,gO2)
coupR2 = cplChiChihhR(i1,i2,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MChi(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MChi(i2)*B1m2  
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
coupL2 = cplcFdFdhhL(i2,i1,gO2)
coupR2 = cplcFdFdhhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFd(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
coupL2 = cplcFeFehhL(i2,i1,gO2)
coupR2 = cplcFeFehhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFe(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVZL(i1,i2)
coupR1 = cplcFuFuVZR(i1,i2)
coupL2 = cplcFuFuhhL(i2,i1,gO2)
coupR2 = cplcFuFuhhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFu(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWmgWmVZ
coup2 = cplcgWmgWmhh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWpCgWpCVZ
coup2 = cplcgWpCgWpChh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MHpm2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MHpm2(i2),MHpm2(i1)),dp) 
coup1 = cplHpmcHpmVZ(i2,i1)
coup2 = cplhhHpmcHpm(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MHpm2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MHpm2(i2)),dp) 
coup1 = cplHpmcVWmVZ(i2)
coup2 = cplhhcHpmVWm(gO2,i2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSd2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVZ(i2,i1)
coup2 = cplhhSdcSd(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSe2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVZ(i2,i1)
coup2 = cplhhSecSe(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSu2(i2),MSu2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSu2(i2),MSu2(i1)),dp) 
coup1 = cplSucSuVZ(i2,i1)
coup2 = cplhhSucSu(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSv2(i2),MSv2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSv2(i2),MSv2(i1)),dp) 
coup1 = cplSvcSvVZ(i2,i1)
coup2 = cplhhSvcSv(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VWm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MHpm2(i1)),dp) 
coup1 = cplcHpmVWmVZ(i1)
coup2 = cplhhHpmcVWm(gO2,i1)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MVWm2),dp)
B1m2 = Real(SA_B1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcVWmVWmVZ
coup2 = cplhhcVWmVWm(gO2)
    SumI = coup1*coup2*(3._dp/2._dp*B0m2+3._dp*B1m2) 
End do 
res = res +1._dp* SumI  
res = oo16pi2*res 
 
End Subroutine Pi1LoopVZhh 
 
Subroutine DerPi1LoopVZhh(p2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,MFu,             & 
& MFu2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,MVWm,MVWm2,cplcChaChahhL,          & 
& cplcChaChahhR,cplcChaChaVZL,cplcChaChaVZR,cplcFdFdhhL,cplcFdFdhhR,cplcFdFdVZL,         & 
& cplcFdFdVZR,cplcFeFehhL,cplcFeFehhR,cplcFeFeVZL,cplcFeFeVZR,cplcFuFuhhL,               & 
& cplcFuFuhhR,cplcFuFuVZL,cplcFuFuVZR,cplcgWmgWmhh,cplcgWmgWmVZ,cplcgWpCgWpChh,          & 
& cplcgWpCgWpCVZ,cplChiChihhL,cplChiChihhR,cplChiChiVZL,cplChiChiVZR,cplcHpmVWmVZ,       & 
& cplcVWmVWmVZ,cplhhcHpmVWm,cplhhcVWmVWm,cplhhHpmcHpm,cplhhHpmcVWm,cplhhSdcSd,           & 
& cplhhSecSe,cplhhSucSu,cplhhSvcSv,cplHpmcHpmVZ,cplHpmcVWmVZ,cplSdcSdVZ,cplSecSeVZ,      & 
& cplSucSuVZ,cplSvcSvVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),               & 
& MFu2(3),MHpm(2),MHpm2(2),MSd(6),MSd2(6),MSe(6),MSe2(6),MSu(6),MSu2(6),MSv(3),          & 
& MSv2(3),MVWm,MVWm2

Complex(dp), Intent(in) :: cplcChaChahhL(2,2,2),cplcChaChahhR(2,2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),      & 
& cplcFdFdhhL(3,3,2),cplcFdFdhhR(3,3,2),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),               & 
& cplcFeFehhL(3,3,2),cplcFeFehhR(3,3,2),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),               & 
& cplcFuFuhhL(3,3,2),cplcFuFuhhR(3,3,2),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),               & 
& cplcgWmgWmhh(2),cplcgWmgWmVZ,cplcgWpCgWpChh(2),cplcgWpCgWpCVZ,cplChiChihhL(4,4,2),     & 
& cplChiChihhR(4,4,2),cplChiChiVZL(4,4),cplChiChiVZR(4,4),cplcHpmVWmVZ(2),               & 
& cplcVWmVWmVZ,cplhhcHpmVWm(2,2),cplhhcVWmVWm(2),cplhhHpmcHpm(2,2,2),cplhhHpmcVWm(2,2),  & 
& cplhhSdcSd(2,6,6),cplhhSecSe(2,6,6),cplhhSucSu(2,6,6),cplhhSvcSv(2,3,3),               & 
& cplHpmcHpmVZ(2,2),cplHpmcVWmVZ(2),cplSdcSdVZ(6,6),cplSecSeVZ(6,6),cplSucSuVZ(6,6),     & 
& cplSvcSvVZ(3,3)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MCha2(i1),MCha2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVZL(i1,i2)
coupR1 = cplcChaChaVZR(i1,i2)
coupL2 = cplcChaChahhL(i2,i1,gO2)
coupR2 = cplcChaChahhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MCha(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MCha(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MChi2(i1),MChi2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MChi2(i1),MChi2(i2)),dp) 
coupL1 = cplChiChiVZL(i1,i2)
coupR1 = cplChiChiVZR(i1,i2)
coupL2 = cplChiChihhL(i1,i2,gO2)
coupR2 = cplChiChihhR(i1,i2,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MChi(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MChi(i2)*B1m2  
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFd2(i1),MFd2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
coupL2 = cplcFdFdhhL(i2,i1,gO2)
coupR2 = cplcFdFdhhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFd(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFe2(i1),MFe2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
coupL2 = cplcFeFehhL(i2,i1,gO2)
coupR2 = cplcFeFehhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFe(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFu2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVZL(i1,i2)
coupR1 = cplcFuFuVZR(i1,i2)
coupL2 = cplcFuFuhhL(i2,i1,gO2)
coupR2 = cplcFuFuhhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFu(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWmgWmVZ
coup2 = cplcgWmgWmhh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWpCgWpCVZ
coup2 = cplcgWpCgWpChh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MHpm2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MHpm2(i2),MHpm2(i1)),dp) 
coup1 = cplHpmcHpmVZ(i2,i1)
coup2 = cplhhHpmcHpm(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MHpm2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MHpm2(i2)),dp) 
coup1 = cplHpmcVWmVZ(i2)
coup2 = cplhhcHpmVWm(gO2,i2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSd2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVZ(i2,i1)
coup2 = cplhhSdcSd(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSe2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVZ(i2,i1)
coup2 = cplhhSecSe(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSu2(i2),MSu2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSu2(i2),MSu2(i1)),dp) 
coup1 = cplSucSuVZ(i2,i1)
coup2 = cplhhSucSu(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sv], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSv2(i2),MSv2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSv2(i2),MSv2(i1)),dp) 
coup1 = cplSvcSvVZ(i2,i1)
coup2 = cplhhSvcSv(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VWm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MHpm2(i1)),dp) 
coup1 = cplcHpmVWmVZ(i1)
coup2 = cplhhHpmcVWm(gO2,i1)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VWm 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVWm2),dp)
B1m2 = Real(SA_DerB1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcVWmVWmVZ
coup2 = cplhhcVWmVWm(gO2)
    SumI = coup1*coup2*(3._dp/2._dp*B0m2+3._dp*B1m2) 
End do 
res = res +1._dp* SumI  
res = oo16pi2*res 
 
End Subroutine DerPi1LoopVZhh 
 
Subroutine Pi1LoopVZAh(p2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,MFe2,           & 
& MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MVWm,MVWm2,MVZ,MVZ2,           & 
& cplAhAhhh,cplAhcHpmVWm,cplAhhhVZ,cplAhHpmcHpm,cplAhHpmcVWm,cplAhSdcSd,cplAhSecSe,      & 
& cplAhSucSu,cplcChaChaAhL,cplcChaChaAhR,cplcChaChaVZL,cplcChaChaVZR,cplcFdFdAhL,        & 
& cplcFdFdAhR,cplcFdFdVZL,cplcFdFdVZR,cplcFeFeAhL,cplcFeFeAhR,cplcFeFeVZL,               & 
& cplcFeFeVZR,cplcFuFuAhL,cplcFuFuAhR,cplcFuFuVZL,cplcFuFuVZR,cplcgWmgWmAh,              & 
& cplcgWmgWmVZ,cplcgWpCgWpCAh,cplcgWpCgWpCVZ,cplChiChiAhL,cplChiChiAhR,cplChiChiVZL,     & 
& cplChiChiVZR,cplcHpmVWmVZ,cplhhVZVZ,cplHpmcHpmVZ,cplHpmcVWmVZ,cplSdcSdVZ,              & 
& cplSecSeVZ,cplSucSuVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),               & 
& MFe2(3),MFu(3),MFu2(3),Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MSd(6),MSd2(6),MSe(6),          & 
& MSe2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2

Complex(dp), Intent(in) :: cplAhAhhh(2,2,2),cplAhcHpmVWm(2,2),cplAhhhVZ(2,2),cplAhHpmcHpm(2,2,2),cplAhHpmcVWm(2,2),& 
& cplAhSdcSd(2,6,6),cplAhSecSe(2,6,6),cplAhSucSu(2,6,6),cplcChaChaAhL(2,2,2),            & 
& cplcChaChaAhR(2,2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),cplcFdFdAhL(3,3,2),         & 
& cplcFdFdAhR(3,3,2),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFeFeAhL(3,3,2),               & 
& cplcFeFeAhR(3,3,2),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),cplcFuFuAhL(3,3,2),               & 
& cplcFuFuAhR(3,3,2),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcgWmgWmAh(2),cplcgWmgWmVZ,     & 
& cplcgWpCgWpCAh(2),cplcgWpCgWpCVZ,cplChiChiAhL(4,4,2),cplChiChiAhR(4,4,2),              & 
& cplChiChiVZL(4,4),cplChiChiVZR(4,4),cplcHpmVWmVZ(2),cplhhVZVZ(2),cplHpmcHpmVZ(2,2),    & 
& cplHpmcVWmVZ(2),cplSdcSdVZ(6,6),cplSecSeVZ(6,6),cplSucSuVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! hh, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MAh2(i2),Mhh2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MAh2(i2),Mhh2(i1)),dp) 
coup1 = cplAhhhVZ(i2,i1)
coup2 = cplAhAhhh(gO2,i2,i1)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MCha2(i1),MCha2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVZL(i1,i2)
coupR1 = cplcChaChaVZR(i1,i2)
coupL2 = cplcChaChaAhL(i2,i1,gO2)
coupR2 = cplcChaChaAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MCha(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MCha(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MChi2(i1),MChi2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MChi2(i1),MChi2(i2)),dp) 
coupL1 = cplChiChiVZL(i1,i2)
coupR1 = cplChiChiVZR(i1,i2)
coupL2 = cplChiChiAhL(i1,i2,gO2)
coupR2 = cplChiChiAhR(i1,i2,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MChi(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MChi(i2)*B1m2  
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFd2(i1),MFd2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
coupL2 = cplcFdFdAhL(i2,i1,gO2)
coupR2 = cplcFdFdAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFd(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFe2(i1),MFe2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
coupL2 = cplcFeFeAhL(i2,i1,gO2)
coupR2 = cplcFeFeAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFe(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFu2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVZL(i1,i2)
coupR1 = cplcFuFuVZR(i1,i2)
coupL2 = cplcFuFuAhL(i2,i1,gO2)
coupR2 = cplcFuFuAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFu(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWmgWmVZ
coup2 = cplcgWmgWmAh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWpCgWpCVZ
coup2 = cplcgWpCgWpCAh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! VZ, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVZ2,Mhh2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MVZ2,Mhh2(i2)),dp) 
coup1 = cplhhVZVZ(i2)
coup2 = cplAhhhVZ(gO2,i2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +0.5_dp* SumI  
    End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MHpm2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MHpm2(i2),MHpm2(i1)),dp) 
coup1 = cplHpmcHpmVZ(i2,i1)
coup2 = cplAhHpmcHpm(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MHpm2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MHpm2(i2)),dp) 
coup1 = cplHpmcVWmVZ(i2)
coup2 = cplAhcHpmVWm(gO2,i2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSd2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVZ(i2,i1)
coup2 = cplAhSdcSd(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSe2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVZ(i2,i1)
coup2 = cplAhSecSe(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSu2(i2),MSu2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSu2(i2),MSu2(i1)),dp) 
coup1 = cplSucSuVZ(i2,i1)
coup2 = cplAhSucSu(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VWm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MHpm2(i1)),dp) 
coup1 = cplcHpmVWmVZ(i1)
coup2 = cplAhHpmcVWm(gO2,i1)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 res = oo16pi2*res 
 
End Subroutine Pi1LoopVZAh 
 
Subroutine DerPi1LoopVZAh(p2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,             & 
& MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MVWm,MVWm2,               & 
& MVZ,MVZ2,cplAhAhhh,cplAhcHpmVWm,cplAhhhVZ,cplAhHpmcHpm,cplAhHpmcVWm,cplAhSdcSd,        & 
& cplAhSecSe,cplAhSucSu,cplcChaChaAhL,cplcChaChaAhR,cplcChaChaVZL,cplcChaChaVZR,         & 
& cplcFdFdAhL,cplcFdFdAhR,cplcFdFdVZL,cplcFdFdVZR,cplcFeFeAhL,cplcFeFeAhR,               & 
& cplcFeFeVZL,cplcFeFeVZR,cplcFuFuAhL,cplcFuFuAhR,cplcFuFuVZL,cplcFuFuVZR,               & 
& cplcgWmgWmAh,cplcgWmgWmVZ,cplcgWpCgWpCAh,cplcgWpCgWpCVZ,cplChiChiAhL,cplChiChiAhR,     & 
& cplChiChiVZL,cplChiChiVZR,cplcHpmVWmVZ,cplhhVZVZ,cplHpmcHpmVZ,cplHpmcVWmVZ,            & 
& cplSdcSdVZ,cplSecSeVZ,cplSucSuVZ,kont,res)

Implicit None 
Real(dp), Intent(in) :: MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),               & 
& MFe2(3),MFu(3),MFu2(3),Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MSd(6),MSd2(6),MSe(6),          & 
& MSe2(6),MSu(6),MSu2(6),MVWm,MVWm2,MVZ,MVZ2

Complex(dp), Intent(in) :: cplAhAhhh(2,2,2),cplAhcHpmVWm(2,2),cplAhhhVZ(2,2),cplAhHpmcHpm(2,2,2),cplAhHpmcVWm(2,2),& 
& cplAhSdcSd(2,6,6),cplAhSecSe(2,6,6),cplAhSucSu(2,6,6),cplcChaChaAhL(2,2,2),            & 
& cplcChaChaAhR(2,2,2),cplcChaChaVZL(2,2),cplcChaChaVZR(2,2),cplcFdFdAhL(3,3,2),         & 
& cplcFdFdAhR(3,3,2),cplcFdFdVZL(3,3),cplcFdFdVZR(3,3),cplcFeFeAhL(3,3,2),               & 
& cplcFeFeAhR(3,3,2),cplcFeFeVZL(3,3),cplcFeFeVZR(3,3),cplcFuFuAhL(3,3,2),               & 
& cplcFuFuAhR(3,3,2),cplcFuFuVZL(3,3),cplcFuFuVZR(3,3),cplcgWmgWmAh(2),cplcgWmgWmVZ,     & 
& cplcgWpCgWpCAh(2),cplcgWpCgWpCVZ,cplChiChiAhL(4,4,2),cplChiChiAhR(4,4,2),              & 
& cplChiChiVZL(4,4),cplChiChiVZR(4,4),cplcHpmVWmVZ(2),cplhhVZVZ(2),cplHpmcHpmVZ(2,2),    & 
& cplHpmcVWmVZ(2),cplSdcSdVZ(6,6),cplSecSeVZ(6,6),cplSucSuVZ(6,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! hh, Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MAh2(i2),Mhh2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MAh2(i2),Mhh2(i1)),dp) 
coup1 = cplAhhhVZ(i2,i1)
coup2 = cplAhAhhh(gO2,i2,i1)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Cha 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MCha2(i1),MCha2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MCha2(i1),MCha2(i2)),dp) 
coupL1 = cplcChaChaVZL(i1,i2)
coupR1 = cplcChaChaVZR(i1,i2)
coupL2 = cplcChaChaAhL(i2,i1,gO2)
coupR2 = cplcChaChaAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MCha(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MCha(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! Chi, Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 4
       Do i2 = 1, 4
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MChi2(i1),MChi2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MChi2(i1),MChi2(i2)),dp) 
coupL1 = cplChiChiVZL(i1,i2)
coupR1 = cplChiChiVZR(i1,i2)
coupL2 = cplChiChiAhL(i1,i2,gO2)
coupR2 = cplChiChiAhR(i1,i2,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MChi(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MChi(i2)*B1m2  
End do 
res = res +0.5_dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFd2(i1),MFd2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFd2(i1),MFd2(i2)),dp) 
coupL1 = cplcFdFdVZL(i1,i2)
coupR1 = cplcFdFdVZR(i1,i2)
coupL2 = cplcFdFdAhL(i2,i1,gO2)
coupR2 = cplcFdFdAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFd(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fe 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFe2(i1),MFe2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFe2(i1),MFe2(i2)),dp) 
coupL1 = cplcFeFeVZL(i1,i2)
coupR1 = cplcFeFeVZR(i1,i2)
coupL2 = cplcFeFeAhL(i2,i1,gO2)
coupR2 = cplcFeFeAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFe(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fu], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFu2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFu2(i1),MFu2(i2)),dp) 
coupL1 = cplcFuFuVZL(i1,i2)
coupR1 = cplcFuFuVZR(i1,i2)
coupL2 = cplcFuFuAhL(i2,i1,gO2)
coupR2 = cplcFuFuAhR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFu(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gWm], gWm 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWmgWmVZ
coup2 = cplcgWmgWmAh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWmC], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVWm2),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MVWm2),dp) 
coup1 = cplcgWpCgWpCVZ
coup2 = cplcgWpCgWpCAh(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! VZ, hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVZ2,Mhh2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MVZ2,Mhh2(i2)),dp) 
coup1 = cplhhVZVZ(i2)
coup2 = cplAhhhVZ(gO2,i2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +0.5_dp* SumI  
    End Do 
 !------------------------ 
! conj[Hpm], Hpm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MHpm2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MHpm2(i2),MHpm2(i1)),dp) 
coup1 = cplHpmcHpmVZ(i2,i1)
coup2 = cplAhHpmcHpm(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], Hpm 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MHpm2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MHpm2(i2)),dp) 
coup1 = cplHpmcVWmVZ(i2)
coup2 = cplAhcHpmVWm(gO2,i2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Sd 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSd2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSd2(i2),MSd2(i1)),dp) 
coup1 = cplSdcSdVZ(i2,i1)
coup2 = cplAhSdcSd(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Se 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSe2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSe2(i2),MSe2(i1)),dp) 
coup1 = cplSecSeVZ(i2,i1)
coup2 = cplAhSecSe(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Su], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSu2(i2),MSu2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSu2(i2),MSu2(i1)),dp) 
coup1 = cplSucSuVZ(i2,i1)
coup2 = cplAhSucSu(gO2,i1,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VWm 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MHpm2(i1)),dp) 
coup1 = cplcHpmVWmVZ(i1)
coup2 = cplAhHpmcVWm(gO2,i1)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 res = oo16pi2*res 
 
End Subroutine DerPi1LoopVZAh 
 
Subroutine Pi1LoopVWmSe(p2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MSd,MSd2,MSe,MSe2,              & 
& MSu,MSu2,MSv,MSv2,cplcFdFuVWmL,cplcFdFuVWmR,cplcFeFvVWmL,cplcFeFvVWmR,cplcFuFdcSeL,    & 
& cplcFuFdcSeR,cplFvFecSeL,cplFvFecSeR,cplSdcSecSu,cplSecSecSv,cplSucSdVWm,              & 
& cplSvcSeVWm,kont,res)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MSd(6),MSd2(6),MSe(6),MSe2(6),           & 
& MSu(6),MSu2(6),MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplcFdFuVWmL(3,3),cplcFdFuVWmR(3,3),cplcFeFvVWmL(3,3),cplcFeFvVWmR(3,3),              & 
& cplcFuFdcSeL(3,3,6),cplcFuFdcSeR(3,3,6),cplFvFecSeL(3,3,6),cplFvFecSeR(3,3,6),         & 
& cplSdcSecSu(6,6,6),cplSecSecSv(6,6,3),cplSucSdVWm(6,6),cplSvcSeVWm(3,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(6) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! bar[Fd], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,6 
B0m2 = Real(SA_B0(p2,MFd2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFd2(i1),MFu2(i2)),dp) 
coupL1 = cplcFdFuVWmL(i1,i2)
coupR1 = cplcFdFuVWmR(i1,i2)
coupL2 = cplcFuFdcSeL(i2,i1,gO2)
coupR2 = cplcFuFdcSeR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,6 
B0m2 = Real(SA_B0(p2,MFe2(i1),0._dp),dp) 
B1m2 = Real(SA_B1(p2,MFe2(i1),0._dp),dp) 
coupL1 = cplcFeFvVWmL(i1,i2)
coupR1 = cplcFeFvVWmR(i1,i2)
coupL2 = cplFvFecSeL(i2,i1,gO2)
coupR2 = cplFvFecSeR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*0.*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sd], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,6 
B0m2 = Real(SA_B0(p2,MSu2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSu2(i2),MSd2(i1)),dp) 
coup1 = cplSucSdVWm(i2,i1)
coup2 = cplSdcSecSu(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 3
 Do gO2=1,6 
B0m2 = Real(SA_B0(p2,MSv2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSv2(i2),MSe2(i1)),dp) 
coup1 = cplSvcSeVWm(i2,i1)
coup2 = cplSecSecSv(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 res = oo16pi2*res 
 
End Subroutine Pi1LoopVWmSe 
 
Subroutine DerPi1LoopVWmSe(p2,MFd,MFd2,MFe,MFe2,MFu,MFu2,MSd,MSd2,MSe,MSe2,           & 
& MSu,MSu2,MSv,MSv2,cplcFdFuVWmL,cplcFdFuVWmR,cplcFeFvVWmL,cplcFeFvVWmR,cplcFuFdcSeL,    & 
& cplcFuFdcSeR,cplFvFecSeL,cplFvFecSeR,cplSdcSecSu,cplSecSecSv,cplSucSdVWm,              & 
& cplSvcSeVWm,kont,res)

Implicit None 
Real(dp), Intent(in) :: MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MSd(6),MSd2(6),MSe(6),MSe2(6),           & 
& MSu(6),MSu2(6),MSv(3),MSv2(3)

Complex(dp), Intent(in) :: cplcFdFuVWmL(3,3),cplcFdFuVWmR(3,3),cplcFeFvVWmL(3,3),cplcFeFvVWmR(3,3),              & 
& cplcFuFdcSeL(3,3,6),cplcFuFdcSeR(3,3,6),cplFvFecSeL(3,3,6),cplFvFecSeR(3,3,6),         & 
& cplSdcSecSu(6,6,6),cplSecSecSv(6,6,3),cplSucSdVWm(6,6),cplSvcSeVWm(3,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(6) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! bar[Fd], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,6 
B0m2 = Real(SA_DerB0(p2,MFd2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFd2(i1),MFu2(i2)),dp) 
coupL1 = cplcFdFuVWmL(i1,i2)
coupR1 = cplcFdFuVWmR(i1,i2)
coupL2 = cplcFuFdcSeL(i2,i1,gO2)
coupR2 = cplcFuFdcSeR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,6 
B0m2 = Real(SA_DerB0(p2,MFe2(i1),0._dp),dp) 
B1m2 = Real(SA_DerB1(p2,MFe2(i1),0._dp),dp) 
coupL1 = cplcFeFvVWmL(i1,i2)
coupR1 = cplcFeFvVWmR(i1,i2)
coupL2 = cplFvFecSeL(i2,i1,gO2)
coupR2 = cplFvFecSeR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*0.*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Sd], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,6 
B0m2 = Real(SA_DerB0(p2,MSu2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSu2(i2),MSd2(i1)),dp) 
coup1 = cplSucSdVWm(i2,i1)
coup2 = cplSdcSecSu(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 3
 Do gO2=1,6 
B0m2 = Real(SA_DerB0(p2,MSv2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSv2(i2),MSe2(i1)),dp) 
coup1 = cplSvcSeVWm(i2,i1)
coup2 = cplSecSecSv(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 res = oo16pi2*res 
 
End Subroutine DerPi1LoopVWmSe 
 
Subroutine Pi1LoopVWmHpm(p2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,MFe,              & 
& MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,MVWm,            & 
& MVWm2,MVZ,MVZ2,cplAhcHpmVWm,cplAhHpmcHpm,cplcChaChiVWmL,cplcChaChiVWmR,cplcFdFuVWmL,   & 
& cplcFdFuVWmR,cplcFeFvVWmL,cplcFeFvVWmR,cplcFuFdcHpmL,cplcFuFdcHpmR,cplcgAgWpCVWm,      & 
& cplcgWmgZVWm,cplcgWpCgAcHpm,cplcgWpCgZcHpm,cplcgZgWmcHpm,cplcgZgWpCVWm,cplChiChacHpmL, & 
& cplChiChacHpmR,cplcHpmVPVWm,cplcHpmVWmVZ,cplcVWmVPVWm,cplcVWmVWmVZ,cplFvFecHpmL,       & 
& cplFvFecHpmR,cplhhcHpmVWm,cplhhcVWmVWm,cplhhHpmcHpm,cplHpmcHpmVP,cplHpmcHpmVZ,         & 
& cplSdcHpmcSu,cplSecHpmcSv,cplSucSdVWm,cplSvcSeVWm,kont,res)

Implicit None 
Real(dp), Intent(in) :: MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),               & 
& MFe2(3),MFu(3),MFu2(3),Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MSd(6),MSd2(6),MSe(6),          & 
& MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3),MVWm,MVWm2,MVZ,MVZ2

Complex(dp), Intent(in) :: cplAhcHpmVWm(2,2),cplAhHpmcHpm(2,2,2),cplcChaChiVWmL(2,4),cplcChaChiVWmR(2,4),        & 
& cplcFdFuVWmL(3,3),cplcFdFuVWmR(3,3),cplcFeFvVWmL(3,3),cplcFeFvVWmR(3,3),               & 
& cplcFuFdcHpmL(3,3,2),cplcFuFdcHpmR(3,3,2),cplcgAgWpCVWm,cplcgWmgZVWm,cplcgWpCgAcHpm(2),& 
& cplcgWpCgZcHpm(2),cplcgZgWmcHpm(2),cplcgZgWpCVWm,cplChiChacHpmL(4,2,2),cplChiChacHpmR(4,2,2),& 
& cplcHpmVPVWm(2),cplcHpmVWmVZ(2),cplcVWmVPVWm,cplcVWmVWmVZ,cplFvFecHpmL(3,3,2),         & 
& cplFvFecHpmR(3,3,2),cplhhcHpmVWm(2,2),cplhhcVWmVWm(2),cplhhHpmcHpm(2,2,2),             & 
& cplHpmcHpmVP(2,2),cplHpmcHpmVZ(2,2),cplSdcHpmcSu(6,2,6),cplSecHpmcSv(6,2,3),           & 
& cplSucSdVWm(6,6),cplSvcSeVWm(3,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
res = 0._dp 
 
!------------------------ 
! conj[Hpm], Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MAh2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MAh2(i2),MHpm2(i1)),dp) 
coup1 = cplAhcHpmVWm(i2,i1)
coup2 = cplAhHpmcHpm(i2,i1,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 4
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MCha2(i1),MChi2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MCha2(i1),MChi2(i2)),dp) 
coupL1 = cplcChaChiVWmL(i1,i2)
coupR1 = cplcChaChiVWmR(i1,i2)
coupL2 = cplChiChacHpmL(i2,i1,gO2)
coupR2 = cplChiChacHpmR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MCha(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MChi(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFd2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MFd2(i1),MFu2(i2)),dp) 
coupL1 = cplcFdFuVWmL(i1,i2)
coupR1 = cplcFdFuVWmR(i1,i2)
coupL2 = cplcFuFdcHpmL(i2,i1,gO2)
coupR2 = cplcFuFdcHpmR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MFe2(i1),0._dp),dp) 
B1m2 = Real(SA_B1(p2,MFe2(i1),0._dp),dp) 
coupL1 = cplcFeFvVWmL(i1,i2)
coupR1 = cplcFeFvVWmR(i1,i2)
coupL2 = cplFvFecHpmL(i2,i1,gO2)
coupR2 = cplFvFecHpmR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*0.*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gP], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,0._dp),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,0._dp),dp) 
coup1 = cplcgAgWpCVWm
coup2 = cplcgWpCgAcHpm(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gZ], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MVZ2),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,MVZ2),dp) 
coup1 = cplcgZgWpCVWm
coup2 = cplcgWpCgZcHpm(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWm], gZ 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVZ2,MVWm2),dp) 
B1m2 = Real(SA_B1(p2,MVZ2,MVWm2),dp) 
coup1 = cplcgWmgZVWm
coup2 = cplcgZgWmcHpm(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,Mhh2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,Mhh2(i2),MHpm2(i1)),dp) 
coup1 = cplhhcHpmVWm(i2,i1)
coup2 = cplhhHpmcHpm(i2,i1,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,Mhh2(i2)),dp) 
B1m2 = Real(SA_B1(p2,MVWm2,Mhh2(i2)),dp) 
coup1 = cplhhcVWmVWm(i2)
coup2 = cplhhcHpmVWm(i2,gO2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSu2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSu2(i2),MSd2(i1)),dp) 
coup1 = cplSucSdVWm(i2,i1)
coup2 = cplSdcHpmcSu(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MSv2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MSv2(i2),MSe2(i1)),dp) 
coup1 = cplSvcSeVWm(i2,i1)
coup2 = cplSecHpmcSv(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VP 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,0._dp,MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,0._dp,MHpm2(i1)),dp) 
coup1 = cplcHpmVPVWm(i1)
coup2 = cplHpmcHpmVP(i1,gO2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VP 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,0._dp),dp)
B1m2 = Real(SA_B1(p2,MVWm2,0._dp),dp) 
coup1 = cplcVWmVPVWm
coup2 = cplcHpmVPVWm(gO2)
    SumI = coup1*coup2*(3._dp/2._dp*B0m2+3._dp*B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], VZ 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVZ2,MHpm2(i1)),dp) 
B1m2 = Real(SA_B1(p2,MVZ2,MHpm2(i1)),dp) 
coup1 = cplcHpmVWmVZ(i1)
coup2 = cplHpmcHpmVZ(i1,gO2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VZ 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_B0(p2,MVWm2,MVZ2),dp)
B1m2 = Real(SA_B1(p2,MVWm2,MVZ2),dp) 
coup1 = cplcVWmVWmVZ
coup2 = cplcHpmVWmVZ(gO2)
    SumI = coup1*coup2*(3._dp/2._dp*B0m2+3._dp*B1m2) 
End do 
res = res +1._dp* SumI  
res = oo16pi2*res 
 
End Subroutine Pi1LoopVWmHpm 
 
Subroutine DerPi1LoopVWmHpm(p2,MAh,MAh2,MCha,MCha2,MChi,MChi2,MFd,MFd2,               & 
& MFe,MFe2,MFu,MFu2,Mhh,Mhh2,MHpm,MHpm2,MSd,MSd2,MSe,MSe2,MSu,MSu2,MSv,MSv2,             & 
& MVWm,MVWm2,MVZ,MVZ2,cplAhcHpmVWm,cplAhHpmcHpm,cplcChaChiVWmL,cplcChaChiVWmR,           & 
& cplcFdFuVWmL,cplcFdFuVWmR,cplcFeFvVWmL,cplcFeFvVWmR,cplcFuFdcHpmL,cplcFuFdcHpmR,       & 
& cplcgAgWpCVWm,cplcgWmgZVWm,cplcgWpCgAcHpm,cplcgWpCgZcHpm,cplcgZgWmcHpm,cplcgZgWpCVWm,  & 
& cplChiChacHpmL,cplChiChacHpmR,cplcHpmVPVWm,cplcHpmVWmVZ,cplcVWmVPVWm,cplcVWmVWmVZ,     & 
& cplFvFecHpmL,cplFvFecHpmR,cplhhcHpmVWm,cplhhcVWmVWm,cplhhHpmcHpm,cplHpmcHpmVP,         & 
& cplHpmcHpmVZ,cplSdcHpmcSu,cplSecHpmcSv,cplSucSdVWm,cplSvcSeVWm,kont,res)

Implicit None 
Real(dp), Intent(in) :: MAh(2),MAh2(2),MCha(2),MCha2(2),MChi(4),MChi2(4),MFd(3),MFd2(3),MFe(3),               & 
& MFe2(3),MFu(3),MFu2(3),Mhh(2),Mhh2(2),MHpm(2),MHpm2(2),MSd(6),MSd2(6),MSe(6),          & 
& MSe2(6),MSu(6),MSu2(6),MSv(3),MSv2(3),MVWm,MVWm2,MVZ,MVZ2

Complex(dp), Intent(in) :: cplAhcHpmVWm(2,2),cplAhHpmcHpm(2,2,2),cplcChaChiVWmL(2,4),cplcChaChiVWmR(2,4),        & 
& cplcFdFuVWmL(3,3),cplcFdFuVWmR(3,3),cplcFeFvVWmL(3,3),cplcFeFvVWmR(3,3),               & 
& cplcFuFdcHpmL(3,3,2),cplcFuFdcHpmR(3,3,2),cplcgAgWpCVWm,cplcgWmgZVWm,cplcgWpCgAcHpm(2),& 
& cplcgWpCgZcHpm(2),cplcgZgWmcHpm(2),cplcgZgWpCVWm,cplChiChacHpmL(4,2,2),cplChiChacHpmR(4,2,2),& 
& cplcHpmVPVWm(2),cplcHpmVWmVZ(2),cplcVWmVPVWm,cplcVWmVWmVZ,cplFvFecHpmL(3,3,2),         & 
& cplFvFecHpmR(3,3,2),cplhhcHpmVWm(2,2),cplhhcVWmVWm(2),cplhhHpmcHpm(2,2,2),             & 
& cplHpmcHpmVP(2,2),cplHpmcHpmVZ(2,2),cplSdcHpmcSu(6,2,6),cplSecHpmcSv(6,2,3),           & 
& cplSucSdVWm(6,6),cplSvcSeVWm(3,6)

Integer, Intent(inout) :: kont 
Real(dp) :: B0m2, F0m2, G0m2, B1m2, H0m2, B22m2, m1, m2 
Real(dp), Intent(in) :: p2 
Complex(dp) :: A0m2 
Complex(dp), Intent(inout) :: res(2) 
Complex(dp) :: coupL1, coupR1, coupL2,coupR2, coup1,coup2, coup3, temp, sumI 
Integer :: i1,i2,i3,i4, gO1, gO2, ierr 
 
 
Real(dp) ::MVG,MVP,MVG2,MVP2
MVG = Mass_Regulator_PhotonGluon 
MVP = Mass_Regulator_PhotonGluon 
MVG2 = Mass_Regulator_PhotonGluon**2 
MVP2 = Mass_Regulator_PhotonGluon**2 

res = 0._dp 
 
!------------------------ 
! conj[Hpm], Ah 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MAh2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MAh2(i2),MHpm2(i1)),dp) 
coup1 = cplAhcHpmVWm(i2,i1)
coup2 = cplAhHpmcHpm(i2,i1,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Cha], Chi 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 4
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MCha2(i1),MChi2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MCha2(i1),MChi2(i2)),dp) 
coupL1 = cplcChaChiVWmL(i1,i2)
coupR1 = cplcChaChiVWmR(i1,i2)
coupL2 = cplChiChacHpmL(i2,i1,gO2)
coupR2 = cplChiChacHpmR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MCha(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MChi(i2)*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fd], Fu 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFd2(i1),MFu2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MFd2(i1),MFu2(i2)),dp) 
coupL1 = cplcFdFuVWmL(i1,i2)
coupR1 = cplcFdFuVWmR(i1,i2)
coupL2 = cplcFuFdcHpmL(i2,i1,gO2)
coupR2 = cplcFuFdcHpmR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFd(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*MFu(i2)*B1m2  
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[Fe], Fv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 3
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MFe2(i1),0._dp),dp) 
B1m2 = Real(SA_DerB1(p2,MFe2(i1),0._dp),dp) 
coupL1 = cplcFeFvVWmL(i1,i2)
coupR1 = cplcFeFvVWmR(i1,i2)
coupL2 = cplFvFecHpmL(i2,i1,gO2)
coupR2 = cplFvFecHpmR(i2,i1,gO2)
    SumI = (coupL1*coupL2+coupR1*coupR2)*MFe(i1)*(B0m2+B1m2) & 
  & + (coupL1*coupR2+coupR1*coupL2)*0.*B1m2  
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! bar[gP], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVP2),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MVP2),dp) 
coup1 = cplcgAgWpCVWm
coup2 = cplcgWpCgAcHpm(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gZ], gWmC 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVZ2),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,MVZ2),dp) 
coup1 = cplcgZgWpCVWm
coup2 = cplcgWpCgZcHpm(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! bar[gWm], gZ 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVZ2,MVWm2),dp) 
B1m2 = Real(SA_DerB1(p2,MVZ2,MVWm2),dp) 
coup1 = cplcgWmgZVWm
coup2 = cplcgZgWmcHpm(gO2) 
   SumI = -0.5_dp*coup1*coup2*(B0m2+B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], hh 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
       Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,Mhh2(i2),MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,Mhh2(i2),MHpm2(i1)),dp) 
coup1 = cplhhcHpmVWm(i2,i1)
coup2 = cplhhHpmcHpm(i2,i1,gO2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[VWm], hh 
!------------------------ 
sumI = 0._dp 
 
      Do i2 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,Mhh2(i2)),dp) 
B1m2 = Real(SA_DerB1(p2,MVWm2,Mhh2(i2)),dp) 
coup1 = cplhhcVWmVWm(i2)
coup2 = cplhhcHpmVWm(i2,gO2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
    End Do 
 !------------------------ 
! conj[Sd], Su 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 6
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSu2(i2),MSd2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSu2(i2),MSd2(i1)),dp) 
coup1 = cplSucSdVWm(i2,i1)
coup2 = cplSdcHpmcSu(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +3._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Se], Sv 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 6
       Do i2 = 1, 3
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MSv2(i2),MSe2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MSv2(i2),MSe2(i1)),dp) 
coup1 = cplSvcSeVWm(i2,i1)
coup2 = cplSecHpmcSv(i1,gO2,i2)
    SumI = -coup1*coup2*(B1m2+0.5_dp*B0m2)
End do 
res = res +1._dp* SumI  
      End Do 
     End Do 
 !------------------------ 
! conj[Hpm], VP 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVP2,MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MVP2,MHpm2(i1)),dp) 
coup1 = cplcHpmVPVWm(i1)
coup2 = cplHpmcHpmVP(i1,gO2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VP 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVP2),dp)
B1m2 = Real(SA_DerB1(p2,MVWm2,MVP2),dp) 
coup1 = cplcVWmVPVWm
coup2 = cplcHpmVPVWm(gO2)
    SumI = coup1*coup2*(3._dp/2._dp*B0m2+3._dp*B1m2) 
End do 
res = res +1._dp* SumI  
!------------------------ 
! conj[Hpm], VZ 
!------------------------ 
sumI = 0._dp 
 
    Do i1 = 1, 2
 Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVZ2,MHpm2(i1)),dp) 
B1m2 = Real(SA_DerB1(p2,MVZ2,MHpm2(i1)),dp) 
coup1 = cplcHpmVWmVZ(i1)
coup2 = cplHpmcHpmVZ(i1,gO2)
    SumI = coup1*coup2*(B1m2-B0m2) 
End do 
res = res +1._dp* SumI  
      End Do 
 !------------------------ 
! conj[VWm], VZ 
!------------------------ 
sumI = 0._dp 
 
Do gO2=1,2 
B0m2 = Real(SA_DerB0(p2,MVWm2,MVZ2),dp)
B1m2 = Real(SA_DerB1(p2,MVWm2,MVZ2),dp) 
coup1 = cplcVWmVWmVZ
coup2 = cplcHpmVWmVZ(gO2)
    SumI = coup1*coup2*(3._dp/2._dp*B0m2+3._dp*B1m2) 
End do 
res = res +1._dp* SumI  
res = oo16pi2*res 
 
End Subroutine DerPi1LoopVWmHpm 
 
End Module LoopMasses_MSSMTriLnV 

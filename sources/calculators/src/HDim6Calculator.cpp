/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        HDim6Calculator.cpp                                              *
*                                                                              *
* Description Wrapper class for HDim6                                          *
*                                                                              *
* Authors     Bastian Heller    <bastian.heller@rwth-aachen.de>                *
*             Bjoern  Sarrazin  <sarrazin@physik.uni-bonn.de>                  *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include <boost/foreach.hpp>

#include "TMath.h"

#include <LHAPDF/LHAPDF.h>

#include "HDim6/CrossSection_had.h"
#include "HDim6/decay.h"
#include "HDim6/gaugecpl.h"
#include "HDim6/VBF.h"

#include "HDim6Calculator.h"
#include "Quantity.h"
#include "ModelBase.h"
#include "SimplePrediction.h"

Fittino::HDim6Calculator::HDim6Calculator(const ModelBase *model, boost::property_tree::ptree &ptree)
    : CalculatorBase( model, ptree ),
      _calculate_Gamma_hWW     ( ptree.get<bool>( "calculate_Gamma_hWW"      ) ),
      _calculate_Gamma_hZZ     ( ptree.get<bool>( "calculate_Gamma_hZZ"      ) ),
      _calculate_xs_qqh_2flavor( ptree.get<bool>( "calculate_xs_qqh_2flavor" ) ),
      _calculate_xs_qqh_5flavor( ptree.get<bool>( "calculate_xs_qqh_5flavor" ) ),
      _calculate_xs_Wh         ( ptree.get<bool>( "calculate_xs_Wh"          ) ),
      _calculate_xs_Zh         ( ptree.get<bool>( "calculate_xs_Zh"          ) ),
      _effsmvalues ( new effinputs() ),
      _effvalues   ( new effinputs() ),
      _first       ( true ),
      _pdfSet      ( "CT10" ),  // lhapdf-getdata CT10.LHgrid
      _pdfDirectory( "" ),
      _smvalues ( new sminputs() ) {

    if ( _name.empty() )  _name = "HDim6Calculator";
    if ( _tag.empty() ) _tag = "HDim6";

    AddInput( "Mass_h", "Mass_h" );
    AddInput( "f_t", "f_t" );
    AddInput( "f_b", "f_b" );
    AddInput( "f_tau", "f_tau" );
    AddInput( "f_GG", "f_GG" );
    AddInput( "f_W", "f_W" );
    AddInput( "f_B", "f_B" );
    AddInput( "f_BB", "f_BB" );
    AddInput( "f_WW", "f_WW" );
    AddInput( "f_Phi_2", "f_Phi_2" );


    double r = ptree.get<double>( "UnitarityCoefficientR", 1 );
    _effvalues->rbb  = r*r; //ALEX: I now put r*r instead of r so we can use GeV units in the input file while effvalues->r__ is GeV^2 like s 
    _effvalues->rww  = r*r;
    _effvalues->rgg  = r*r;
    _effvalues->rb   = r*r;
    _effvalues->rw   = r*r;
    _effvalues->rbw  = r*r;
    _effvalues->rp1  = r*r;
    _effvalues->rp2  = r*r;
    _effvalues->rp4  = r*r;
    _effvalues->rtop = r*r;
    _effvalues->rbot = r*r;
    _effvalues->rtau = r*r;

    double n = ptree.get<double>( "UnitarityCoefficientN", 0 );
    _effvalues->nbb  = n;
    _effvalues->nww  = n;
    _effvalues->ngg  = n;
    _effvalues->nb   = n;
    _effvalues->nw   = n;
    _effvalues->nbw  = n;
    _effvalues->np1  = n;
    _effvalues->np2  = n;
    _effvalues->np4  = n;
    _effvalues->ntop = n;
    _effvalues->nbot = n;
    _effvalues->ntau = n;

    _effvalues->override_unitarity = ptree.get<bool>( "OverrideUnitarity", true );

    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_g_g",         "",      _normSM_Gamma_hgg     ) );
    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_tau_tau",     "",      _normSM_Gamma_htautau ) );
    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_mu_mu",       "",      _normSM_Gamma_hmumu   ) );
    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_gamma_gamma", "",      _normSM_Gamma_hgaga   ) );
    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_Z_ga",        "",      _normSM_Gamma_hZga    ) );
    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_b_b",         "",      _normSM_Gamma_hbb     ) );
    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_c_c",         "",      _normSM_Gamma_hcc     ) );
    AddQuantity( new SimplePrediction( "NormSM_Gamma_h_s_s",         "",      _normSM_Gamma_hss     ) );
    AddQuantity( new SimplePrediction( "NormSM_xs_ggh",              "",      _normSM_xs_ggh        ) );
    AddQuantity( new SimplePrediction( "NormSM_xs_bbh",              "",      _normSM_xs_bbh        ) );
    AddQuantity( new SimplePrediction( "NormSM_xs_tth",              "",      _normSM_xs_tth        ) );
    AddQuantity( new SimplePrediction( "NormSM_xs_bh",               "",      _normSM_xs_bh         ) );

    if ( _calculate_Gamma_hZZ      ) AddQuantity( new SimplePrediction( "NormSM_Gamma_h_Z_Z",    "", _normSM_Gamma_hZZ      ) );
    if ( _calculate_Gamma_hWW      ) AddQuantity( new SimplePrediction( "NormSM_Gamma_h_W_W",    "", _normSM_Gamma_hWW      ) );
    if ( _calculate_xs_qqh_2flavor ) AddQuantity( new SimplePrediction( "NormSM_xs_qqh_2flavor", "", _normSM_xs_qqh_2flavor ) );
    if ( _calculate_xs_qqh_5flavor ) AddQuantity( new SimplePrediction( "NormSM_xs_qqh_5flavor", "", _normSM_xs_qqh_5flavor ) );
    if ( _calculate_xs_Wh          ) AddQuantity( new SimplePrediction( "NormSM_xs_Wh",          "", _normSM_xs_Wh          ) );
    if ( _calculate_xs_Zh          ) AddQuantity( new SimplePrediction( "NormSM_xs_Zh",          "", _normSM_xs_Zh          ) );

    Messenger::GetInstance()<<Messenger::ALWAYS<<" HDim6Calculator configured with the following unitarity settings: "<<Messenger::Endl;
    Messenger::GetInstance()<<Messenger::ALWAYS<<" OverrideUnitarity: "<<_effvalues->override_unitarity<<Messenger::Endl;
    Messenger::GetInstance()<<Messenger::ALWAYS<<" UnitarityCoefficientR: "<<_effvalues->rbb<<Messenger::Endl;
    Messenger::GetInstance()<<Messenger::ALWAYS<<" UnitarityCoefficientN: "<<_effvalues->nbb<<Messenger::Endl;

}

Fittino::HDim6Calculator::~HDim6Calculator() {

    delete _effvalues;
    delete _effsmvalues;
    delete _smvalues;

}

void Fittino::HDim6Calculator::CalculatePredictions() {

    ConfigureInput();
    CallFunction();

}

void Fittino::HDim6Calculator::Initialize() {

    if ( _pdfDirectory != "" ) {

        LHAPDF::setPDFPath( _pdfDirectory );

    }

    LHAPDF::initPDFSet( _pdfSet, LHAPDF::LHGRID, 0 );

}

void Fittino::HDim6Calculator::CallFunction() {

    bool new_mh = ( _first || _previous_mass_h != _smvalues ->mh );

    bool new_gridParameters = (
                                  _first
                                  || new_mh
                                  || _previous_f_B     != _effvalues->fb
                                  || _previous_f_BB    != _effvalues->fbb
                                  || _previous_f_WW    != _effvalues->fww
                                  || _previous_f_W     != _effvalues->fw
                                  || _previous_f_Phi_2 != _effvalues->fp2
                              );

    if ( _first )  _first = false;

    if ( new_mh )  _previous_mass_h = _smvalues ->mh;

    if ( new_gridParameters ) {

        _previous_f_B     = _effvalues->fb;
        _previous_f_W     = _effvalues->fw;
        _previous_f_BB    = _effvalues->fbb;
        _previous_f_WW    = _effvalues->fww;
        _previous_f_Phi_2 = _effvalues->fp2;

    }

    _Delta_g1_WW       = HDim6::d_g1_ww  ( _smvalues, _effvalues );
    _Delta_g2_WW       = HDim6::d_g2_ww  ( _smvalues, _effvalues );
    _Delta_g1_gaga     = HDim6::d_g1_yy  ( _smvalues, _effvalues );
    _Delta_g2_gaga     = HDim6::d_g2_yy  ( _smvalues, _effvalues );
    _Delta_g1_ZZ       = HDim6::d_g1_zz  ( _smvalues, _effvalues );
    _Delta_g2_ZZ       = HDim6::d_g2_zz  ( _smvalues, _effvalues );
    _Delta_g1_Zga      = HDim6::d_g1_zy  ( _smvalues, _effvalues );
    _Delta_g2_Zga      = HDim6::d_g2_zy  ( _smvalues, _effvalues );
    _Delta_kappa_Gamma = HDim6::d_kappa_y( _smvalues, _effvalues );
    _Delta_kappa_Z     = HDim6::d_kappa_z( _smvalues, _effvalues );
    _Delta_g1_Gamma    = HDim6::d_g1_y   ( _smvalues, _effvalues );
    _Delta_g1_Z        = HDim6::d_g1_z   ( _smvalues, _effvalues );

    double error, chi2;

    if ( new_mh ) {

        hglgl_( _smvalues, _effsmvalues, &_SM_Gamma_hgg,     &error );
        hgaga_( _smvalues, _effsmvalues, &_SM_Gamma_hgaga,   &error );
        hgaz_ ( _smvalues, _effsmvalues, &_SM_Gamma_hZga,    &error );
        hmumu_( _smvalues, _effsmvalues, &_SM_Gamma_hmumu,   &error );
        htata_( _smvalues, _effsmvalues, &_SM_Gamma_htautau, &error );
        hchch_( _smvalues, _effsmvalues, &_SM_Gamma_hcc,     &error );
        hstst_( _smvalues, _effsmvalues, &_SM_Gamma_hss,     &error );
        hbobo_( _smvalues, _effsmvalues, &_SM_Gamma_hbb,     &error );

    }

    hglgl_( _smvalues, _effvalues, &_Gamma_hgg,     &error );
    hgaga_( _smvalues, _effvalues, &_Gamma_hgaga,   &error );
    hgaz_ ( _smvalues, _effvalues, &_Gamma_hZga,    &error );
    hmumu_( _smvalues, _effvalues, &_Gamma_hmumu,   &error );
    htata_( _smvalues, _effvalues, &_Gamma_htautau, &error );
    hchch_( _smvalues, _effvalues, &_Gamma_hcc,     &error );
    hstst_( _smvalues, _effvalues, &_Gamma_hss,     &error );
    hbobo_( _smvalues, _effvalues, &_Gamma_hbb,     &error );

    _normSM_Gamma_hgg     = _Gamma_hgg     / _SM_Gamma_hgg;
    _normSM_Gamma_hgaga   = _Gamma_hgaga   / _SM_Gamma_hgaga;
    _normSM_Gamma_hZga    = _Gamma_hZga    / _SM_Gamma_hZga;
    _normSM_Gamma_hmumu   = _Gamma_hmumu   / _SM_Gamma_hmumu;
    _normSM_Gamma_htautau = _Gamma_htautau / _SM_Gamma_htautau;
    _normSM_Gamma_hcc     = _Gamma_hcc     / _SM_Gamma_hcc;
    _normSM_Gamma_hss     = _Gamma_hss     / _SM_Gamma_hss;
    _normSM_Gamma_hbb     = _Gamma_hbb     / _SM_Gamma_hbb;

    if ( new_gridParameters )  {

        if ( new_mh ) {

            if ( _calculate_Gamma_hZZ ) hzz_( _smvalues, _effsmvalues, &_SM_Gamma_hZZ, &error, &chi2 );
            if ( _calculate_Gamma_hWW ) hww_( _smvalues, _effsmvalues, &_SM_Gamma_hWW, &error, &chi2 );

        }

        if ( _calculate_Gamma_hZZ ) hzz_( _smvalues, _effvalues, &_Gamma_hZZ, &error, &chi2 );
        if ( _calculate_Gamma_hWW ) hww_( _smvalues, _effvalues, &_Gamma_hWW, &error, &chi2 );

        if ( _calculate_Gamma_hZZ ) _normSM_Gamma_hZZ = _Gamma_hZZ / _SM_Gamma_hZZ;
        if ( _calculate_Gamma_hWW ) _normSM_Gamma_hWW = _Gamma_hWW / _SM_Gamma_hWW;

    }

    // todo: calc xs in own class, possible to set E

    k_ggh_      ( _smvalues, _effvalues, &_normSM_xs_ggh, &error, &chi2 );
    ratio_bb_h_ ( _smvalues, _effvalues, &_normSM_xs_bbh, &error, &chi2 );
    ratio_tth_  ( _smvalues, _effvalues, &_normSM_xs_tth, &error, &chi2 );
    ratio_bg_bh_( _smvalues, _effvalues, &_normSM_xs_bh,  &error, &chi2 );

    if ( new_gridParameters ) {

        if ( new_mh ) {

            if ( _calculate_xs_Wh )          HWRadiation_( _smvalues, _effsmvalues, &_SM_xs_Wh,          &error, &chi2 );
            if ( _calculate_xs_Zh )          HZRadiation_( _smvalues, _effsmvalues, &_SM_xs_Zh,          &error, &chi2 );
            if ( _calculate_xs_qqh_2flavor ) ud_jjh_     ( _smvalues, _effsmvalues, &_SM_xs_qqh_2flavor, &error, &chi2 );
            if ( _calculate_xs_qqh_5flavor ) udcsb_jjh_  ( _smvalues, _effsmvalues, &_SM_xs_qqh_5flavor, &error, &chi2 );

        }

        if ( _calculate_xs_Wh )          HWRadiation_( _smvalues, _effvalues, &_xs_Wh,          &error, &chi2 );
        if ( _calculate_xs_Zh )          HZRadiation_( _smvalues, _effvalues, &_xs_Zh,          &error, &chi2 );
        if ( _calculate_xs_qqh_2flavor ) ud_jjh_     ( _smvalues, _effvalues, &_xs_qqh_2flavor, &error, &chi2 );
        if ( _calculate_xs_qqh_5flavor ) udcsb_jjh_  ( _smvalues, _effvalues, &_xs_qqh_5flavor, &error, &chi2 );

        if ( _calculate_xs_Wh )          _normSM_xs_Wh          = _xs_Wh          / _SM_xs_Wh;
        if ( _calculate_xs_Zh )          _normSM_xs_Zh          = _xs_Zh          / _SM_xs_Zh;
        if ( _calculate_xs_qqh_2flavor ) _normSM_xs_qqh_2flavor = _xs_qqh_2flavor / _SM_xs_qqh_2flavor;
        if ( _calculate_xs_qqh_5flavor ) _normSM_xs_qqh_5flavor = _xs_qqh_5flavor / _SM_xs_qqh_5flavor;

    }

}

void Fittino::HDim6Calculator::ConfigureInput() {

    UpdateInput();

    _effvalues->fb   = 1e-6 * GetInput( "f_B" );
    _effvalues->fbb  = 1e-6 * GetInput( "f_BB" );
    _effvalues->fw   = 1e-6 * GetInput( "f_W" );
    _effvalues->fww  = 1e-6 * GetInput( "f_WW" );
    _effvalues->fgg  = 1e-6 * GetInput( "f_GG" );
    _effvalues->fp2  = 1e-6 * GetInput( "f_Phi_2" );
    _effvalues->fboh = 1e-6 * GetInput( "f_b" );
    _effvalues->ftoh = 1e-6 * GetInput( "f_t" );
    _effvalues->ftah = 1e-6 * GetInput( "f_tau" );

    _smvalues ->mh   = GetInput( "Mass_h" );

    _f_g  =  - GetInput( "f_GG" ) * 8 * TMath::Pi() / ( _smvalues->alphas ); // f_g as defined in 1211.4580v4.pdf eq 38 but without factor of vev ( because of units ).

}

/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        HECModel.cpp                                                     *
*                                                                              *
* Description Implementation of the HEC model                                  *
*                                                                              *
* Authors     Sebastian Heer  <s6seheer@uni-bonn.de>                           *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include "HECModel.h"
#include "InputException.h"
#include "HiggsSignalsSLHAModelCalculator.h"
#include "PhysicsParameter.h"
#include "SLHAChi2Contribution.h"
#include "SLHAPrediction.h"

Fittino::HECModel::HECModel() {

    _name = "HEC model";

    //_parameterVector.push_back( new PhysicsParameter( "AlphaEM",           "#alpha_{em}",       128.952, "GeV", "GeV", 0.1,  50., 150., true  ) );
    //_parameterVector.push_back( new PhysicsParameter( "AlphaS",            "#alpha_{S}",         0.1184, "GeV", "GeV", 0.1, -50.,  50., true  ) );
    //_parameterVector.push_back( new PhysicsParameter( "G_F",               "G_{F}",        1.1663787e-5,    "",    "", 0.1,  -5.,   5., true  ) );

    _parameterVector.push_back( new PhysicsParameter( "Massh",             "m_{h}",                126., "GeV", "GeV", 0.02, 125., 128., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_s_hss",       "#Delta_{s}hss",          0.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_p_hss",       "#Delta_{p}hss",         -1.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_s_hcc",       "#Delta_{s}hcc",          0.,    "",    "",  0.1,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_p_hcc",       "#Delta_{p}hcc",         -1.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_s_hbb",       "#Delta_{s}hbb",          0.,    "",    "",  0.1,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_p_hbb",       "#Delta_{p}hbb",         -1.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_s_htt",       "#Delta_{s}htt",          0.,    "",    "",  0.2,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_p_htt",       "#Delta_{p}htt",         -1.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_s_hmumu",     "#Delta_{s}h#mu#mu",      0.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_p_hmumu",     "#Delta_{p}h#mu#mu",     -1.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_s_htautau",   "#Delta_{s}h#tau#tau",    0.,    "",    "",  0.1,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_p_htautau",   "#Delta_{p}h#tau#tau",   -1.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_hWW",         "#Delta hWW",             0.,    "",    "",  0.2,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_hZZ",         "#Delta hZZ",             0.,    "",    "",  0.2,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_hZgamma",     "#Delta hZ#gamma",        0.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_hgammagamma", "#Delta h#gamma#gamma",   0.,    "",    "",  0.1,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_hgg",         "#Delta hgg",             0.,    "",    "",  0.2,  -7.,   5., false ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_hggZ",        "#Delta hggZ",            0.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "Delta_hihjZ",       "#Delta h_{i}h_{j}",     -1.,    "",    "",  0.1,  -7.,   5., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "BR_hjhihi",         "BR_{hjhihi}",            0.,    "",    "",  0.1,   0.,   1., true  ) );
    _parameterVector.push_back( new PhysicsParameter( "BR_hInvisible",     "BR_{hInvisible}",        0.,    "",    "",  0.1,   0.,   1., true  ) );

    //SLHAModelCalculatorBase* slhaModelCalculator = new HiggsSignalsSLHAModelCalculator();
    //_modelCalculatorVector.push_back( slhaModelCalculator );

    //_predictionVector.push_back( new SLHAPrediction( "Gamma_hTotal",                             "#Gamma-hTotal",                                "GeV", "GeV", slhaModelCalculator, "HiggsSignalsAdditionalPredictions", "7", 1 ) );
    ////_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_inclusive_ATL",              "#mu_{h#gamma#gamma}inclusive-ATL",                "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "1", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_1lep_ATL",                   "#mu_{h#gamma#gamma}1lep-ATL",                     "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "1", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_conv_central_highPTt_ATL",   "#mu_{h#gamma#gamma}conv-central-highPTt-AT",      "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "2", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_conv_central_lowPTt_ATL",    "#mu_{h#gamma#gamma}conv-central-lowPTt-ATL",      "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "3", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_conv_rest_highPTt_ATL",      "#mu_{h#gamma#gamma}conv-rest-highPTt-ATL",        "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "4", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_conv_rest_lowPTt_ATL",       "#mu_{h#gamma#gamma}conv-rest-lowPTt-ATL",         "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "5", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_conv_trans_ATL",             "#mu_{h#gamma#gamma}conv-trans-ATL",               "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "6", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_ETmiss_ATL",                 "#mu_{h#gamma#gamma}ETmiss-ATL",                   "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "7", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_highmass2jetloose_ATL",      "#mu_{h#gamma#gamma}highmass2jetloose-ATL",        "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "8", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_highmass2jettight_ATL",      "#mu_{h#gamma#gamma}highmass2jettight-ATL",        "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables",  "9", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_lowmass2jet_ATL",            "#mu_{h#gamma#gamma}lowmass2jet-ATL",              "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "10", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_unconv_central_highPTt_ATL", "#mu_{h#gamma#gamma}unconv-central-highPTt-ATL",   "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "11", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_unconv_central_lowPTt_ATL",  "#mu_{h#gamma#gamma}unconv-central-lowPTt-ATL",    "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "12", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_unconv_rest_highPTt_ATL",    "#mu_{h#gamma#gamma}unconv-rest-highPTt-ATL",      "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "13", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammgamma_unconv_rest_lowPTt_ATL",     "#mu_{h#gamma#gamma}unconv-rest-lowPTt-ATL",       "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "14", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_htautau_ATL",                           "#mu_{h#tau#tau}ATL",                              "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "15", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hWWlnulnu_ATL",                         "#mu_{hWWlnunu}ATL",                               "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "16", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_VBFWWlnulnu_ATL",                       "#mu_{VBFWWlnulnu}ATL",                            "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "17", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hZZ4l_ATL",                             "#mu_{hZZ4l}ATL",                                  "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "18", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_VhVbb_ATL",                             "#mu_{VhVbb}ATL",                                  "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "19", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_CDF",                       "#mu_{h#gamma#gamma}CDF",                          "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "20", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_htautau_CDF",                           "#mu_{h#tau#tau}CDF",                              "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "21", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hWW_CDF",                               "#mu_{hWW}CDF",                                    "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "22", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_tthttbb_CDF",                           "#mu_{tthttb}CDF",                                 "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "23", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_VhVbb_CDF",                             "#mu_{VhVbb}CDF",                                  "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "24", "17", 2 ) );
    ////_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_inclusive_CMS",             "#mu_{h#gamma#gamma}inclusive-CMS",                "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "12", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_2jet_CMS",                  "#mu_{h#gamma#gamma}2jet-CMS",                     "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "25", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_2jetloose_CMS",             "#mu_{h#gamma#gamma}2jetloose-CMS",                "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "26", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_2jettight_CMS",             "#mu_{h#gamma#gamma}2jettight-CMS",                "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "27", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_e_CMS",                     "#mu_{h#gamma#gamma}e-CMS",                        "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "28", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_inclusive_CMS",             "#mu_{h#gamma#gamma}inclusive-CMS",                "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "29", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_ETmiss_CMS",                "#mu_{h#gamma#gamma}ETmiss-CMS",                   "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "30", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_mu_CMS",                    "#mu_{h#gamma#gamma}#mu-CMS",                      "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "31", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_htautau_0_1jet_CMS",                    "#mu_{h#tau#tau}0-1jet-CMS",                       "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "32", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hWW_0_1jet_CMS",                        "#mu_{hWW}0-1jet-CMS",                             "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "33", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hZZ4l_0_1jet_CMS",                      "#mu_{hZZ4l}0-1jet-CMS",                           "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "34", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hZZ4l_2jet_CMS",                        "#mu_{hZZ4l}2jet-CMS",                             "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "35", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_tthttbb_CMS",                           "#mu_{tthttb}CMS",                                 "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "36", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_VBFtautau_CMS",                         "#mu_{VBF#tau#tau}CMS",                            "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "37", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_VBFWW_CMS",                             "#mu_{VBFWW}CMS",                                  "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "38", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_VhVbb_CMS",                             "#mu_{VhVbb}CMS",                                  "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "39", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_VhVtautau_CMS",                         "#mu_{VhV#tau#tau}CMS",                            "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "40", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_WhWWW_CMS",                             "#mu_{WhWWW}CMS",                                  "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "41", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hbb_D0",                                "#mu_{hbb}DO",                                     "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "42", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hgammagamma_D0",                        "#mu_{h#gamma#gamma}DO",                           "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "43", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_htautau_D0",                            "#mu_{h#tau#tau}DO",                               "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "44", "17", 2 ) );
    //_predictionVector.push_back( new SLHAPrediction( "Mu_hWW_D0",                                "#mu_{hWw}DO",                                     "",    "", slhaModelCalculator, "HiggsSignalsPeakObservables", "45", "17", 2 ) );

    //_predictionVector.push_back( new SLHAPrediction( "R_hWW",                                    "R-hWW",                                           "",    "", slhaModelCalculator, "HiggsSignalsAdditionalPredictions", "1", 1 ) );
    //_predictionVector.push_back( new SLHAPrediction( "R_hZZ",                                    "R-hZZ",                                           "",    "", slhaModelCalculator, "HiggsSignalsAdditionalPredictions", "2", 1 ) );
    //_predictionVector.push_back( new SLHAPrediction( "R_hgammagamma",                            "R-h#gamma#gamma",                                 "",    "", slhaModelCalculator, "HiggsSignalsAdditionalPredictions", "3", 1 ) );
    //_predictionVector.push_back( new SLHAPrediction( "R_htautau",                                "R-h#tau#tau",                                     "",    "", slhaModelCalculator, "HiggsSignalsAdditionalPredictions", "4", 1 ) );
    //_predictionVector.push_back( new SLHAPrediction( "R_hbb",                                    "R-hbb",                                           "",    "", slhaModelCalculator, "HiggsSignalsAdditionalPredictions", "5", 1 ) );
    //_predictionVector.push_back( new SLHAPrediction( "R_Vhbb",                                   "R-Vhbb",                                          "",    "", slhaModelCalculator, "HiggsSignalsAdditionalPredictions", "6", 1 ) );
    //
    ////_predictionVector.push_back( new SLHAPrediction( "Masss",   "GeV", slhaModelCalculator, "", "", 0 ) );
    ////_predictionVector.push_back( new SLHAPrediction( "Massc",   "GeV", slhaModelCalculator, "", "", 1 ) );
    ////_predictionVector.push_back( new SLHAPrediction( "Massb",   "GeV", slhaModelCalculator, "", "", 1 ) );
    ////_predictionVector.push_back( new SLHAPrediction( "Masst",   "GeV", slhaModelCalculator, "", "", 1 ) );
    ////_predictionVector.push_back( new SLHAPrediction( "Massmu",  "GeV", slhaModelCalculator, "", "", 1 ) );
    ////_predictionVector.push_back( new SLHAPrediction( "Masstau", "GeV", slhaModelCalculator, "", "", 1 ) );

    //_chi2ContributionVector.push_back( new SLHAChi2Contribution( "HiggsSignals", slhaModelCalculator, "HiggsSignalsResults", "12", 1 ) );

    PhysicsModelBase::Initialize();

}

Fittino::HECModel::~HECModel() {

}

Fittino::HECModel* Fittino::HECModel::Clone() const {

    return new HECModel( *this );

}

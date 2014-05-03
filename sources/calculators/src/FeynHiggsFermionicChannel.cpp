/* $Id$ */

/*******************************************************************************
 *                                                                              *
 * Project     Fittino - A SUSY Parameter Fitting Package                       *
 *                                                                              *
 * File        FeynHiggsFermionicChannel.cpp                                    *
 *                                                                              *
 * Description Wrapper class for FeynHiggs                                      *
 *                                                                              *
 * Authors     Bjoern Sarrazin  <sarrazin@physik.uni-bonn.de>                   *
 *                                                                              *
 * Licence     This program is free software; you can redistribute it and/or    *
 *             modify it under the terms of the GNU General Public License as   *
 *             published by the Free Software Foundation; either version 3 of   *
 *             the License, or (at your option) any later version.              *
 *                                                                              *
 *******************************************************************************/

#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>
#include "boost/format.hpp"
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include "TVector2.h"

#include "CFeynHiggs.h"
#include "CSLHA.h"

#include "CalculatorException.h"
#include "FeynHiggsFermionicChannel.h"
#include "ModelParameter.h"
#include "PhysicsModel.h"
#include "Redirector.h"
#include "SimplePrediction.h"

Fittino::FeynHiggsFermionicChannel::FeynHiggsFermionicChannel( FHRealType* gammas, FHRealType* gammasms, FHComplexType* couplings, FHComplexType* couplingsms, std::string higgsName, std::string channelName, int channelNumber, bool SM )
: FeynHiggsChannel( gammas, gammasms, couplings, couplingsms, higgsName, channelName, channelNumber, SM ) {

    AddQuantity( new SimplePrediction( "Abs_gs_" + higgsName + "_" + channelName, "", _model_gs2 ) );
    AddQuantity( new SimplePrediction( "Abs_gp_" + higgsName + "_" + channelName, "", _model_gp2 ) );
    AddQuantity( new SimplePrediction( "Arg_gs_" + higgsName + "_" + channelName, "", _model_gsPhi ) );
    AddQuantity( new SimplePrediction( "Arg_gp_" + higgsName + "_" + channelName, "", _model_gpPhi ) );

    if ( _doSM ) {

        AddQuantity( new SimplePrediction( "SM_Abs_gs_" + higgsName + "_" + channelName, "", _sm_gs2 ) );
        AddQuantity( new SimplePrediction( "SM_Abs_gp_" + higgsName + "_" + channelName, "", _sm_gp2 ) );
        AddQuantity( new SimplePrediction( "SM_Arg_gs_" + higgsName + "_" + channelName, "", _sm_gsPhi ) );
        AddQuantity( new SimplePrediction( "SM_Arg_gp_" + higgsName + "_" + channelName, "", _sm_gpPhi ) );

        AddQuantity( new SimplePrediction( "NormSM_Abs_gs_" + higgsName + "_" + channelName, "", _normSM_gs2 ) );
        AddQuantity( new SimplePrediction( "NormSM_Abs_gp_" + higgsName + "_" + channelName, "", _normSM_gp2 ) );
        AddQuantity( new SimplePrediction( "DiffSM_Arg_gs_" + higgsName + "_" + channelName, "", _normSM_gsPhi ) );
        AddQuantity( new SimplePrediction( "DiffSM_Arg_gp_" + higgsName + "_" + channelName, "", _normSM_gpPhi ) );

    }

}

Fittino::FeynHiggsFermionicChannel::~FeynHiggsFermionicChannel() {

}

void Fittino::FeynHiggsFermionicChannel::CalculatePredictions() {

    FeynHiggsChannel::CalculatePredictions();

    FHComplexType coup;

    coup = RCoupling( _channel ) + LCoupling( _channel );
    coup = coup / FHComplexType( 2, 0 );

    _model_gs2   = std::abs( coup );
    _model_gsPhi = std::arg( coup );

    coup = RCoupling( _channel) - LCoupling( _channel );
    coup = coup / FHComplexType( 0, 2 );

    _model_gp2   = std::abs( coup );
    _model_gpPhi = std::arg( coup );

    if ( _doSM ) {
        
        coup  = RCouplingSM( _channel ) + LCouplingSM( _channel );
        coup  = coup / FHComplexType( 2, 0 );

        _sm_gs2   = std::abs( coup );
        _sm_gsPhi = std::arg( coup );

        _normSM_gs2   = _model_gs2 / _sm_gs2;
        _normSM_gsPhi = TVector2::Phi_mpi_pi( _model_gsPhi - _sm_gsPhi );

        coup = RCouplingSM( _channel) - LCouplingSM( _channel );
        coup = coup / FHComplexType( 0, 2 );

        _sm_gp2   = std::abs( coup );
        _sm_gpPhi = std::arg( coup );

        _normSM_gp2   = _model_gp2 / _sm_gp2;
        _normSM_gpPhi = TVector2::Phi_mpi_pi( _model_gpPhi - _sm_gpPhi );

    }

}


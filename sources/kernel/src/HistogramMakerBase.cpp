/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        HistogramMakerBase.cpp                                           *
*                                                                              *
* Description Base class for histogram makers                                  *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/


//   _    _       _ _         __  __        _   _     _             
//  | |  | |     | | |       |  \/  |      | | | |   (_)            
//  | |__| | __ _| | | ___   | \  / | __ _ | |_| |__  _  __ _ ___   
//  |  __  |/ _` | | |/ _ \  | |\/| |/ _` || __| '_ \| |/ _` / __|  
//  | |  | | (_| | | | (_) | | |  | | (_| || |_| | | | | (_| \__ \_ 
//  |_|  |_|\__,_|_|_|\___/  |_|  |_|\__,_| \__|_| |_|_|\__,_|___( )
//                                                               |/ 
//           _                  _   _                               _ 
//          (_)                | | | |                             | |
// __      ___ _ __   ___  ___ | | | |_ ___ _ __    _ __ ___   __ _| |
// \ \ /\ / / | '__| / __|/ _ \| | | __/ _ \ '_ \  | '_ ` _ \ / _` | |
//  \ V  V /| | |    \__ \ (_) | | | ||  __/ | | | | | | | | | (_| | |
//   \_/\_/ |_|_|    |___/\___/|_|_|\__\___|_| |_| |_| |_| |_|\__,_|_|
//
//  _   _ _                    _                     _  __         
// (_) (_) |                   | (_)                | |/ / |                   
//  _   _| |__   ___ _ __    __| |_  ___  ___  ___  | ' /| | __ _ ___ ___  ___ 
// | | | | '_ \ / _ \ '__|  / _` | |/ _ \/ __|/ _ \ |  < | |/ _` / __/ __|/ _ \
// | |_| | |_) |  __/ |    | (_| | |  __/\__ \  __/ | . \| | (_| \__ \__ \  __/
//  \__,_|_.__/ \___|_|     \__,_|_|\___||___/\___| |_|\_\_|\__,_|___/___/\___|
//              _               
//              | |             
//  _ __ ___  __| | ___ _ __    
// | '__/ _ \/ _` |/ _ \ '_ \   
// | | |  __/ (_| |  __/ | | |_ 
// |_|  \___|\__,_|\___|_| |_(_)
//
//                                
// Wie kann ich dich erreichen? Wann passt es dir?
// Viele Gruesse, Bjoern                            
//

#include <boost/foreach.hpp>

#include "TGaxis.h"
#include "TH1.h"
#include "TMath.h"

#include "Factory.h"
#include "HistogramMakerBase.h"
#include "ModelBase.h"
#include "ModelParameter.h"
#include "PlotterBase.h"

Fittino::HistogramMakerBase::HistogramMakerBase( ModelBase* model, const boost::property_tree::ptree& ptree )
    : Tool( model, ptree ),
      _iEntry( 0 ),
      _globalAxis( new TGaxis() ) {

    // Determine the number of iterations.
    // If the number of tree iterations is larger than the difference between the upper bound
    // and lower bound of the model parameter P_Iteration, this difference is used.
    // also, in that case, set _iEntry to the lower bound.

    _numberOfPointsToProcess = _model->GetCollectionOfQuantities().At( "TreeIterations" )->GetValue();
    int difference = _model->GetCollectionOfQuantities().At( "P_Iteration" )->GetUpperBound() - _model->GetCollectionOfQuantities().At( "P_Iteration" )->GetLowerBound();

    if ( difference < _numberOfPointsToProcess ) {

        _numberOfPointsToProcess = difference;
        _iEntry = _model->GetCollectionOfQuantities().At( "P_Iteration" )->GetLowerBound();

    }

    _globalAxis->SetMaxDigits( ptree.get<int>( "AxisMaxDigits", 3 ) );

    // Loop over the quantities scheduled in the Fittino input file and fill the histogram
    // information into appropriate lists for later use.

    BOOST_FOREACH( const boost::property_tree::ptree::value_type & node, ptree ) {

        if ( node.first == "Histogram" ) {

            _quantityName.push_back( node.second.get<std::string> ( "Name"               ) );
            _plotName.push_back    ( node.second.get<std::string> ( "PlotName"           ) );
            _numberOfBins.push_back( node.second.get<unsigned int>( "NumberOfBins"       ) );
            _logScale.push_back    ( node.second.get<bool>        ( "LogScale",    false ) );

            // Depending on whether or not LogScale is set, the histogram upper and lower bounds
            // are recomputed.

            if ( node.second.get<bool>( "LogScale", false ) ) {

                if ( node.second.get<double>( "LowerBound" ) <= 0. || node.second.get<double>( "UpperBound" ) <= 0. ) {

                    throw ConfigurationException( "Histogram axis set to log scale but boundaries are not positive." );

                }
                else {

                    _lowerBound.push_back( TMath::Log10( node.second.get<double>( "LowerBound" ) ) );
                    _upperBound.push_back( TMath::Log10( node.second.get<double>( "UpperBound" ) ) );

                }

            }
            else {

                _lowerBound.push_back( node.second.get<double>( "LowerBound" ) );
                _upperBound.push_back( node.second.get<double>( "UpperBound" ) );

            }

        }

    }

    // Detemine the quantity index of each scheduled quantity in the list of model quantities. This
    // is done because lookup via index is usually faster than lookup via name.

    for ( unsigned int iScheduledQuantity = 0; iScheduledQuantity < _quantityName.size(); ++iScheduledQuantity ) {

        for ( unsigned int iQuantity = 0; iQuantity < _model->GetCollectionOfQuantities().GetNumberOfElements(); ++iQuantity ) {

            if ( _quantityName[iScheduledQuantity] == _model->GetCollectionOfQuantities().At( iQuantity )->GetName() ) {

                _quantityIndex.push_back( iQuantity );

            }

        }

    }

    // Set the plotter.

    const Factory factory;

    const boost::property_tree::ptree::value_type& plotterNode = *( ptree.get_child( "Plotter" ).begin() );
    std::string plotterType = plotterNode.first;
    const boost::property_tree::ptree& plotterTree = plotterNode.second;
    _plotter = factory.CreatePlotter( plotterType, _histogramVector, plotterTree );

}

Fittino::HistogramMakerBase::~HistogramMakerBase() {

    delete _globalAxis;
    _globalAxis = 0;

}

void Fittino::HistogramMakerBase::PrintSteeringParameters() const {

    for ( unsigned int iHistogram = 0; iHistogram < _histogramVector.size(); ++iHistogram ) {

        PrintItem( "Histogram", _histogramVector[iHistogram]->GetName() );

    }

}

void Fittino::HistogramMakerBase::Execute() {

    while ( _iterationCounter < _numberOfPointsToProcess ) {

        _iterationCounter++;

        UpdateModel();

        Tool::PrintStatus();

    }

}

void Fittino::HistogramMakerBase::InitializeTool() {

}

void Fittino::HistogramMakerBase::PrintResult() const {

    Messenger& messenger = Messenger::GetInstance();
    messenger << Messenger::ALWAYS << Messenger::_dashedLine << Messenger::Endl;

}

void Fittino::HistogramMakerBase::Terminate() {

    if ( _plotter ) {

        _plotter->MakePlots();

    }

}

void Fittino::HistogramMakerBase::WriteResultToFile() const {

    for ( Int_t i = 0; i < _histogramVector.size(); ++i ) {

        _histogramVector[i]->Write();

    }

}

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        LHCChi2Contribution.cpp                                          *
*                                                                              *
* Description Class for LHC chi2 contributions                                 *
*                                                                              *
* Authors     Matthias Hamer  <mhamer@gwdg.de>                                 *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	          published by the Free Software Foundation; either version 3 of   *
*	          the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include <vector>
#include <sstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "TRandom3.h"

#include "LHCChi2Contribution.h"
#include "LHCModelCalculator.h"
#include "SimpleDataStorage.h"
#include "Collection.h"

Fittino::LHCChi2Contribution::LHCChi2Contribution( std::string name,
        std::string fileName,
        std::string histogramName,
        std::vector<std::string> relevantParameters,
        int nObs,
        double nExpSM,
        double nExpBestFit,
        double systematicErrorBG,
        double systematicErrorSignal,
        LHCModelCalculator* lhcModelCalculator )
    : _fileName( fileName ),
      _histogramName( histogramName ),
      _nObs( nObs ),
      _nExpSM( nExpSM ),
      _nExpBestFit( nExpBestFit ),
      _systematicErrorBG( systematicErrorBG ),
      _systematicErrorSignal( systematicErrorSignal ),
      _lhcModelCalculator( lhcModelCalculator ),
      Chi2ContributionBase( name ) {

      _relevantParameters = relevantParameters;
      std::stringstream nObs_ss;
      nObs_ss << nObs;
      std::string actualHistogramName = _histogramName + "_nObs" + nObs_ss.str();
      _lhcModelCalculator->AddAnalysis( _name, _fileName, actualHistogramName, _relevantParameters );

}

Fittino::LHCChi2Contribution::LHCChi2Contribution( const boost::property_tree::ptree& ptree, 
                                                   Fittino::Collection<ModelCalculatorBase*>* calculators )
    : _fileName( ptree.get<std::string>( "fileName" ) ),
      _histogramName( ptree.get<std::string>( "histogramName" ) ),
      _nObs( ptree.get<int>( "nObs") ),
      _nExpSM( ptree.get<double>( "nExpSM" ) ),
      _nExpBestFit( ptree.get<double>( "nExpBestFit" ) ),
      _systematicErrorBG( ptree.get<double>( "systematicErrorBG" ) ),
      _systematicErrorSignal( ptree.get<double>( "systematicErrorSignal" ) ),
      _lhcModelCalculator( static_cast<LHCModelCalculator*>(calculators->At( ptree.get<std::string>( "calculatorName", "LHCModelCalculator" ) ) ) ),
      Chi2ContributionBase( ptree.get<std::string>( "name" ) ) {

      BOOST_FOREACH( const boost::property_tree::ptree::value_type& node, ptree ) {
        if( node.first == "RelevantParameter" ) {
            _relevantParameters.push_back( node.second.data() );
        }
      }
      std::stringstream nObs_ss;
      nObs_ss << _nObs;
      std::string actualHistogramName = _histogramName + "_nObs" + nObs_ss.str();
      _lhcModelCalculator->AddAnalysis( _name, _fileName, actualHistogramName, _relevantParameters );
       
}

Fittino::LHCChi2Contribution::~LHCChi2Contribution() {

}

void Fittino::LHCChi2Contribution::UpdateValue() {

    _chi2 = _lhcModelCalculator->GetSimpleOutputDataStorage()->GetMap()->at( _name );

}

void Fittino::LHCChi2Contribution::SmearNObs( TRandom3* randomGenerator ) {

    double nObsNew = randomGenerator->Poisson( _nExpSM*randomGenerator->Gaus(1., _systematicErrorBG) + _nExpBestFit*randomGenerator->Gaus(1., _systematicErrorSignal ) );
    std::stringstream nObsNew_ss;
    nObsNew_ss << nObsNew;
    
    std::string actualHistogramName = _histogramName + "_nObs" + nObsNew_ss.str(); 
    _lhcModelCalculator->UpdateAnalysisHistogram( _name, _fileName, actualHistogramName, _relevantParameters );

}

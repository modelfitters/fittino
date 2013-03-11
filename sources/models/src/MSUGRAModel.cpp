/* $Id: MSUGRAModel.cpp 613 2010-05-26 09:42:00Z uhlenbrock $ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        MSUGRAModel.cpp                                                  *
*                                                                              *
* Description Implementation of the MSUGRA model                               *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include "Configuration.h"
#include "MSUGRAModel.h"
#include "ObservableBase.h"
#include "SLHAChi2Contribution.h"
#include "SLHAObservable.h"
#include "SLHAParameter.h"
#include "SPhenoSLHAModelCalculator.h"

Fittino::MSUGRAModel::MSUGRAModel() {

    Configuration* configuration = Configuration::GetInstance();

    _numberOfParameters = 4;
    _name = "MSUGRA model";
    _parameterVector.push_back( new SLHAParameter( "A0"     , configuration->GetSteeringParameter( "A0" ,     0. ), 1., -1.e5, 1.e5, 5 ) );
    _parameterVector.push_back( new SLHAParameter( "M0"     , configuration->GetSteeringParameter( "M0" ,     0. ), 1.,    0., 1.e5, 1 ) );
    _parameterVector.push_back( new SLHAParameter( "M12"    , configuration->GetSteeringParameter( "M12" ,    0. ), 1.,    0., 1.e5, 2 ) );
    _parameterVector.push_back( new SLHAParameter( "TanBeta", configuration->GetSteeringParameter( "TanBeta", 0. ), 1.,    0., 1.e3, 3 ) );

    SLHAModelCalculatorBase* slhaModelCalculator = new SPhenoSLHAModelCalculator();

    _modelCalculatorVector.push_back( slhaModelCalculator );

    _observableVector.push_back( new SLHAObservable( "BR(b -> s gamma)" , 3.55e-04, 0.53e-04, slhaModelCalculator, "SPhenoLowEnergy",  1, 1 ) );
    _observableVector.push_back( new SLHAObservable( "BR(Bs -> mu+ mu-)",  4.5e-09,  0.8e-09, slhaModelCalculator, "SPhenoLowEnergy",  4, 1 ) );
    _observableVector.push_back( new SLHAObservable( "Delta(M_Bs)"      ,   17.719,    4.200, slhaModelCalculator, "SPhenoLowEnergy",  7, 1 ) );
    _observableVector.push_back( new SLHAObservable( "Delta(g-2)_muon"  , 28.7e-10,  8.2e-10, slhaModelCalculator, "SPhenoLowEnergy", 11, 1 ) );

    //_chi2ContributionVector.push_back( new SLHAChi2Contribution( "TestContribution", slhaModelCalculator, "SPhenoLowEnergy", 1, 1 ) );

    ModelBase::InitializeModel();

}

Fittino::MSUGRAModel::~MSUGRAModel() {

}

Fittino::MSUGRAModel* Fittino::MSUGRAModel::Clone() const {

    return new MSUGRAModel( *this );

}

/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        ModelBase.cpp                                                    *
*                                                                              *
* Description Base class for Fittino models                                    *
*                                                                              *
* Authors     Sebastian Heer        <s6seheer@uni-bonn.de>                     *
*             Mathias   Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>            *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include "Chi2ContributionBase.h"
#include "Configuration.h"
#include "ConfigurationException.h"
#include "ModelBase.h"
#include "ModelParameterBase.h"
#include "PredictionBase.h"
#include "SLHAParameter.h"
#include <iostream>

Fittino::ModelBase::ModelBase()
    : _name( "" ) {

    InitializeParameters();

}

Fittino::ModelBase::~ModelBase() {

}

double Fittino::ModelBase::GetChi2() {

    bool evaluate = false;

    for ( unsigned int i = 0; i < GetNumberOfParameters(); i++ ) {

        if ( GetParameterVector()->at( i )->IsUpdated() ) {

            evaluate = true;
            GetParameterVector()->at( i )->SetUpdated( false );

        }

    }

    if ( evaluate ) {

        _chi2 = Evaluate();

    }

    return _chi2;

}

int Fittino::ModelBase::GetNumberOfChi2Contributions() const {

    return _chi2ContributionVector.size();

}

int Fittino::ModelBase::GetNumberOfParameters() const {

    return _parameterVector.size();

}

int Fittino::ModelBase::GetNumberOfPredictions() const {

    return _predictionVector.size();

}

void Fittino::ModelBase::AddParameter( ModelParameterBase* parameter ) {

    _parameterVector.push_back( parameter );

    if ( !_parameterMap.insert( std::make_pair( parameter->GetName(), parameter ) ).second ) {

        std::string message = "Parameter with name " + parameter->GetName() + " has already been added to model " + GetName();

        throw ConfigurationException( message ); //TODO: Dedicated exception class ?

    }

}

void Fittino::ModelBase::AddPrediction( PredictionBase* prediction ) {

    _collectionOfPredictions.AddElement( prediction->GetName(), prediction );
    _predictionVector.push_back( prediction );

}

std::string Fittino::ModelBase::GetName() const {

    return _name;

}

const std::map<std::string, Fittino::ModelParameterBase*>* Fittino::ModelBase::GetParameterMap() const {

    return &_parameterMap;

}

const std::vector<Fittino::Chi2ContributionBase*>* Fittino::ModelBase::GetChi2ContributionVector() const {

    return &_chi2ContributionVector;

}

const std::vector<Fittino::ModelParameterBase*>* Fittino::ModelBase::GetParameterVector() const {

    return &_parameterVector;

}

const std::vector<Fittino::PredictionBase*>* Fittino::ModelBase::GetPredictionVector() const {

    return &_predictionVector;

}

const Fittino::Collection<Fittino::PredictionBase*>& Fittino::ModelBase::GetCollectionOfPredictions() const {

    return _collectionOfPredictions;

}

void Fittino::ModelBase::InitializeParameters() {

    Configuration *configuration = Configuration::GetInstance();
    const boost::property_tree::ptree* propertyTree = configuration->GetPropertyTree();

    BOOST_FOREACH( const boost::property_tree::ptree::value_type & v, propertyTree->get_child( "InputFile" ) ) {
        if( v.first == "Parameter" ) {
            
            std::string name = v.second.get<std::string>( "<xmlattr>.Name" );
            std::string plotName = v.second.get<std::string>( "<xmlattr>.plotName", name );
            std::string unit = v.second.get<std::string>( "<xmlattr>.Unit", "" );
            std::string plotUnit = v.second.get<std::string>( "<xmlattr>.PlotUnit", unit );
            std::string id = v.second.get<std::string>( "<xmlattr>.ID", "" );
            double value = v.second.get<double>( "<xmlattr>.Value" );
            double error = v.second.get<double>( "<xmlattr>.Error" );
            double lowerBound = v.second.get<double>( "<xmlattr>.LowerBound", 0. );
            double upperBound = v.second.get<double>( "<xmlattr>.UpperBound", 0. );
            double plotLowerBound = v.second.get<double>( "<xmlattr>.PlotLowerBound", 0. );
            double plotUpperBound = v.second.get<double>( "<xmlattr>.PlotUpperBound", 0. );
            bool fixed = v.second.get<bool>( "<xmlattr>.Fixed", false );

            if( v.second.get<std::string>( "<xmlattr>.Type" ) == "SLHA" ) {
                
                AddParameter( new SLHAParameter( name, plotName, value, unit, plotUnit, error, lowerBound, upperBound, plotLowerBound, plotUpperBound, id, fixed ) );

            }
        
        }
    
    }

}

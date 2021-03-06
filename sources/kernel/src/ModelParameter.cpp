/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        ModelParameter.cpp                                               *
*                                                                              *
* Description Class for model parameters                                       *
*                                                                              *
* Authors     Sebastian Heer  <s6seheer@uni-bonn.de>                           *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include <iomanip>

#include <boost/property_tree/ptree.hpp>

#include "ConfigurationException.h"
#include "Messenger.h"
#include "ModelParameter.h"

Fittino::ModelParameter::ModelParameter( boost::property_tree::ptree& ptree )
    : Quantity ( ptree ),
      _fixed   ( ptree.get<bool>  ( "Fixed", false ) ),
      _updated ( true ),
      _ptree   ( ptree ),
      _error   ( ptree.get<double>( "Error", 0.1   ) ),
      _minError( ptree.get<double>( "MinError", 0. ) ) {

}

Fittino::ModelParameter::~ModelParameter() {

}

bool Fittino::ModelParameter::IsFixed() const {

    return _fixed;

}

bool Fittino::ModelParameter::IsUpdated() const {

    return _updated;

}

double Fittino::ModelParameter::GetError() const {

    return _error;

}

double Fittino::ModelParameter::GetMinError() const {

    return _minError;

}

void Fittino::ModelParameter::Initialize() const {

    // ConsistencyCheck();

    PrintStatus();

}

void Fittino::ModelParameter::PrintStatus() const {

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::INFO
              << "    "
              << std::left
              << std::setw( 50 )
              << _name
              << std::right
              << std::setw( 9 )
              << std::setprecision( 2 )
              << std::scientific
              << _value;

    if ( _unit != "" && !_fixed ) {

        messenger << std::right
                  << std::setw( 6 )
                  << _unit;

    }
    else if ( _unit != "" && _fixed ) {

        messenger << std::right
                  << std::setw( 6 )
                  << _unit
                  << std::right
                  << std::setw( 10 )
                  << "(fixed)";

    }
    else if ( _unit == "" && _fixed ) {

        messenger << std::right
                  << std::setw( 16 )
                  << "(fixed)";

    }

    messenger << Messenger::Endl;

}

void Fittino::ModelParameter::SetUpdated( bool updated ) {

    _updated = updated;

}

void Fittino::ModelParameter::SetValue( double value ) {

    if ( !_fixed ) {

        _value = value;
        _updated = true;

    }

}

void Fittino::ModelParameter::UpdatePropertyTree() {

    _ptree.put( "Value", GetValue() );

}

void Fittino::ModelParameter::ConsistencyCheck() const {

    // Check if starting values are within bounds.

    if ( _value < _lowerBound || _value > _upperBound ) {

        throw ConfigurationException( "Parameter " + _name + ": Starting value out of bounds." );

    }

}

/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        PredictionBase.cpp                                               *
*                                                                              *
* Description Base class for predictions                                       *
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

#include "Messenger.h"
#include "PredictionBase.h"

Fittino::PredictionBase::PredictionBase( std::string name, 
                                         std::string plotName,
                                         std::string unit,
                                         std::string plotUnit )
         : _predictedValue( 0. ),
           _name( name ),
           _plotName( plotName ),
           _unit( unit ),
           _plotUnit( plotUnit ){

}

Fittino::PredictionBase::~PredictionBase() {

}

std::string Fittino::PredictionBase::GetName() const {

    return _name;

}

std::string Fittino::PredictionBase::GetPlotName() const {

    return _plotName;

}

std::string Fittino::PredictionBase::GetUnit() const {

    return _unit;

}

std::string Fittino::PredictionBase::GetPlotUnit() const {

    return _plotUnit;

}

void Fittino::PredictionBase::PrintStatus() const {

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::INFO
              << "    "
              << std::left
              << std::setw( 43 )
              << _name
              << std::right
              << std::setw( 9 )
              << std::setprecision( 2 )
              << std::scientific
              << _predictedValue;

    if ( _unit != "" ) {

        messenger << std::right
                  << std::setw( 6 )
                  << _unit;

    }

    messenger << Messenger::Endl;

}

double Fittino::PredictionBase::GetPredictedValue() const {

    return _predictedValue;

}

/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        RandomGenerator.cpp                                              *
*                                                                              *
* Description Singleton wrapper class for random number generator              *
*                                                                              *
* Authors     Matthias Hamer  <mhamer@cbpf.br>                                 *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include "TRandom3.h"

#include "RandomGenerator.h"

Fittino::RandomGenerator* Fittino::RandomGenerator::GetInstance() {

    if ( !_instance ) {

        _instance = new RandomGenerator;

    }

    return _instance;

}

double Fittino::RandomGenerator::Gaus( double mean, double sigma ) {

    return _generator->Gaus( mean, sigma );

}

double Fittino::RandomGenerator::Poisson( double lambda ) {

    return _generator->Poisson( lambda );

}

double Fittino::RandomGenerator::Uniform( double x ) {

    return _generator->Uniform( x );

}

double Fittino::RandomGenerator::Uniform( double x1, double x2 ) {

    return _generator->Uniform( x1, x2 );

}

unsigned int Fittino::RandomGenerator::GetSeed() {

    return _generator->GetSeed();

}

unsigned int Fittino::RandomGenerator::Integer( unsigned int imax ) {

    return _generator->Integer( imax );

}

void Fittino::RandomGenerator::SetSeed( unsigned int seed ) {

    _generator->SetSeed( seed );

}

TRandom3* Fittino::RandomGenerator::GetGenerator() {

    return _generator;

}

Fittino::RandomGenerator* Fittino::RandomGenerator::_instance = 0;

Fittino::RandomGenerator::RandomGenerator()
    : _generator( new TRandom3() ) {

}

Fittino::RandomGenerator::~RandomGenerator() {

    delete _generator;

}

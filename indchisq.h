/***************************************************************************
                                  indchisq.h
                             -------------------    
    Main class which does Fittino's calculations.
                             -------------------
    $Id$
    $Name$
                             -------------------
    begin                : September 10, 2003
    authors              : Philip Bechtle, Peter Wienemann
    email                : philip.bechtle@desy.de, peter.wienemann@desy.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef FITTINO_INDCHISQ_H
#define FITTINO_INDCHISQ_H

#include <misc.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <yy.h>

class Input;
class SmearedInput;

extern vector <parameter_t> indchisq_vec;

class IndChisq {

 public:
  
    IndChisq(const Input* input);
    virtual ~IndChisq();
    void CalcIndChisq();

private:
    const Input*  fInput;

};

#endif


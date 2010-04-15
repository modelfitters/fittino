/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        OptimizerException.cpp                                           *
*                                                                              *
* Description Fittino optimizer exception class. It is thrown in case of       *
*             problems with optimizers.                                        *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include "OptimizerException.h"

Fittino::OptimizerException::OptimizerException( const std::string& message )
    : ExceptionBase( message ) {

}

const char* Fittino::OptimizerException::what() const throw() {

    return ( "Fittino::OptimizerException: " + _message ).c_str();

}

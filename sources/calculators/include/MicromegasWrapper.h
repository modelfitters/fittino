/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        MicromegasWrapper.h                                              *
*                                                                              *
* Description Wrapper class for Micromegas.                                    *
*                                                                              *
* Authors     Bjoern Sarrazin  <sarrazin@physik.uni-bonn.de>                   *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_MICROMEGASWRAPPER_H
#define FITTINO_MICROMEGASWRAPPER_H

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  /*!
   *  \ingroup calculators
   *  \brief Wrapper class for Micromegas.
   *
   *  This class introduces an additional "boost free" layer between Micromegas and\n
   *  MicromegasCalculator. It is needed because Micromegas uses a function boost conflicting with\n
   *  the boost namespace used in MicromegasCalculator.
   */
  class MicromegasWrapper {

    public:
      MicromegasWrapper( std::string name );
      void        CalculatePredictions();

    public:
      virtual     ~MicromegasWrapper();

    protected:
      double      _gmin2;
      double      _omegah2;
      std::string _inputFile;
      std::string _mcname;

  };

}

#endif // FITTINO_MICROMEGASWRAPPER_H

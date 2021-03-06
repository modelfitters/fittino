/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        Simple2DHistogramMaker.h                                         *
*                                                                              *
* Description Class for simple 2D histogram maker                              *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_SIMPLE2DHISTOGRAMMAKER_H
#define FITTINO_SIMPLE2DHISTOGRAMMAKER_H

#include "HistogramMakerBase.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  /*!
   *  \ingroup histogrammakers
   *  \brief Class for simple 2D histogram maker.
   */
  class Simple2DHistogramMaker : public HistogramMakerBase {

    public:
      /*!
       *  Standard constructor.
       */
      Simple2DHistogramMaker( ModelBase* model, const boost::property_tree::ptree& ptree );
      /*!
       *  Standard destructor.
       */
      ~Simple2DHistogramMaker();

      /*! \cond UML */
    private:
      void UpdateModel();

      /*! \endcond UML */

  };

}

#endif // FITTINO_SIMPLE2DHISTOGRAMMAKER_H

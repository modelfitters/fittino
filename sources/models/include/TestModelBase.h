/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        TestModelBase.h                                                  *
*                                                                              *
* Description Base class for Fittino test models                               *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_TESTMODELBASE_H
#define FITTINO_TESTMODELBASE_H

#include "ModelBase.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  /*!
   *  \ingroup models
   *  \brief Base class for Fittino test models.
   */
  class TestModelBase : public ModelBase {

    public:
      /*!
       *  Standard constructor.
       */
                     TestModelBase();
      /*!
       *  Standard destructor.
       */
                     ~TestModelBase();
      virtual void   PrintStatus() const;
			virtual void	 SmearObservables( TRandom3* ) {}

    protected:
      virtual void   Initialize() const;

      /*! \cond UML */
    private:
      virtual double Evaluate();
      virtual double TestModelFunction() = 0;

      /*! \endcond UML */

  };

}

#endif // FITTINO_TESTMODELBASE_H

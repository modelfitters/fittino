/* $Id: Chi2ContributionBase.h 851 2011-01-10 18:01:22Z uhlenbrock $ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        Chi2ContributionBase.h                                           *
*                                                                              *
* Description Base class for chi2 contributions                                *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_CHI2CONTRIBUTIONBASE_H
#define FITTINO_CHI2CONTRIBUTIONBASE_H

#include <string>

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  /*!
   *  \ingroup models
   *  \brief Base class for chi2 contributions.
   */
  class Chi2ContributionBase {

    public:
      /*!
       *  Standard constructor.
       */
                   Chi2ContributionBase( std::string name );
      /*!
       *  Standard destructor.
       */
                   ~Chi2ContributionBase();
      virtual void UpdateValue() = 0;
      void         PrintStatus() const;
      double       GetValue() const;

    protected:
      double       _chi2;
      std::string  _name;

  };

}

#endif // FITTINO_CHI2CONTRIBUTIONBASE_H
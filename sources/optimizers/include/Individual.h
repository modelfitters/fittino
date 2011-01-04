/* $Id: Individual.h,v 1.0 2010/11/03 22:00:00 sarrazin $ */

/*******************************************************************************
*                                                                              *
o* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        Individual.h                                                     *
*                                                                              *
* Description Class for Individual (needed by genetic algorithm optimizer)     *
*                                                                              *
* Authors     Bjoern Sarrazin  <bjoern.sarrazin@desy.de>                       *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_INDIVIDUAL_H
#define FITTINO_INDIVIDUAL_H

#include <iostream>
#include <vector>

#include "TRandom.h"

#include "ModelBase.h"

/*!
 *  \brief Fittino namespace
 */
namespace Fittino {

  /*!
   *  \brief Class for individual (needed by genetic algorithm optimizer)
   */
  class Individual {

    public:
                          Individual( ModelBase* model, double mutationRate, int seed );
                          ~Individual();
      bool                UpdatedChi2();
      double              GetGene( int index );
      void                Mutation();
      void                UpdateChi2(); //calls UpdateModel and than sets the chi2 according to that model
      void                SetGene( int index, double newValue );
      void                UpdateModel(); //sets the pointed model according to the genes

    public:
      bool operator       <( const Individual& individual ) const;

    private:
      bool                _updatedChi2;
      double              _chi2;
      double              _mutationRate;
      std::vector<double> _genes;
      TRandom*            _randomGenerator;
      ModelBase*          _model;

  };

}

#endif // FITTINO_INDIVIDUAL_H
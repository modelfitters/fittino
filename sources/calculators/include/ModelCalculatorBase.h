/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        ModelCalculatorBase.h                                            *
*                                                                              *
* Description Base class for model calculators                                 *
*                                                                              *
* Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*	      published by the Free Software Foundation; either version 3 of   *
*	      the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_MODELCALCULATORBASE_H
#define FITTINO_MODELCALCULATORBASE_H

#include <string>
#include <map>

#include "TStopwatch.h"

#include "Collection.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  class PredictionBase;
  class PhysicsModelBase;
  class SimpleDataStorage;

  /*!
   *  \defgroup calculators
   */
  /*!
   *  \ingroup calculators
   *  \brief Base class for model calculators.
   */
  class ModelCalculatorBase {

    public:
      /*!
       *  Standard constructor.
       */
                                    ModelCalculatorBase( const PhysicsModelBase* model );
      /*!
       *  Standard destructor.
       */
                                    ~ModelCalculatorBase();
      std::string                   GetName() const;

    public:  
      const Collection<const PredictionBase*>&   GetCollectionOfQuantities() const;
      const SimpleDataStorage*      GetSimpleOutputDataStorage() const;

    public:
      virtual void                  CalculatePredictions() = 0;
      virtual void                  Initialize() const = 0;

    protected:
      enum                          CallMethod { EXECUTABLE, FUNCTION };
      void                          AddQuantity( const PredictionBase* prediction );

    protected:
      std::string                   _executableName;
      std::string                   _name;
      TStopwatch                    _stopwatch;
      CallMethod                    _callMethod;
      const PhysicsModelBase*       _model;
      SimpleDataStorage*            _simpleOutputDataStorage;
      Collection<const PredictionBase*>        _collectionOfQuantities;

    protected: 
      void                          StopTime( std::string name ); 

    protected:
      virtual void                  CallExecutable() = 0;
      virtual void                  CallFunction() = 0;
      virtual void                  ConfigureInput() = 0;
      
  };

}

#endif // FITTINO_MODELCALCULATORBASE_H

/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        ModelBase.h                                                      *
*                                                                              *
* Description Base class for Fittino models                                    *
*                                                                              *
* Authors     Sebastian Heer        <s6seheer@uni-bonn.de>                     *
*             Mathias   Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>            *
*             Matthias Hamer        <mhamer@gwdg.de>                           *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_MODELBASE_H
#define FITTINO_MODELBASE_H

#include "TRandom3.h"

#include "Chi2ContributionBase.h"
#include "Collection.h"
#include "PredictionBase.h"
#include "VariableBase.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  class ModelParameter;
  class PredictionBase;
  class CalculatorBase;
  class Observable;
  class Chi2ContributionBase;

  /*!
   *  \ingroup kernel
   *  \brief Base class for Fittino models.
   */
  class ModelBase {

    public:
      /*!
       *  Constructor using a property tree.
       */
      ModelBase( const boost::property_tree::ptree& ptree );
      /*!
       *  Standard destructor.
       */
      ~ModelBase();
      /*!
       *  Returns the chi2 of the comparison between the predicted observables of the model and\n
       *  the measured observables. In the case of a test model simply returns the function value.
       */
      double                                     GetChi2();
      /*!
       *  Returns the number of parameters of the model.
       */
      int                                        GetNumberOfParameters() const;
      /*!
       * Update the property tree.
       */
      void                                       UpdatePropertyTree();
      /*!
       *  Returns the name of the model.
       */
      std::string                                GetName() const;
      /*!
       *  Returns the property tree.
       */
      boost::property_tree::ptree                GetPropertyTree();
      const Collection<Chi2ContributionBase*>&   GetCollectionOfChi2Contributions() const;
      /*!
       *  Returns the parameters as a collection.
       */
      const Collection<ModelParameter*>&         GetCollectionOfParameters() const;
      const Collection<PredictionBase*>&         GetCollectionOfPredictions() const;
      const Collection<const Quantity*>&         GetCollectionOfQuantities() const;
      const Collection<const VariableBase<std::string>*>&  GetCollectionOfStringVariables() const;

    public:
      virtual void                               PrintStatus() const = 0;
      /*!
       *  Smears the observables (if existent).
       */
      virtual void                               SmearObservations( TRandom3* ) = 0;
      /*!
       *  Returns a pointer to a copy of the model.
       */
      virtual ModelBase*                         Clone() const = 0;
      virtual const std::vector<Observable*>*    GetObservableVector() const = 0;
      virtual const Collection<CalculatorBase*>& GetCollectionOfCalculators() const = 0;

    protected:
      /*!
       *  Name of the model.
       */
      std::string                                _name;
      boost::property_tree::ptree                _ptree;
      Collection<Chi2ContributionBase*>          _collectionOfChi2Contributions;
      Collection<const VariableBase<std::string>*> _collectionOfStringVariables;
      /*!
       *  Stores the predictions.
       */
      Collection<PredictionBase*>                _collectionOfPredictions;

    protected:
      /*!
       *  Adds a prediction to the model.
       */
      void                                       AddPrediction( PredictionBase* prediction );

    protected:
      /*!
       *  This function is mainly used to provide a dedicated command to print the model\n
       *  configuration to the screen and is usually called at the end of the model constructor\n
       *  after the actual "initialization" of the model is already finished. Therefore, it should\n
       *  not be used to alter the state of the model. However, since third party code (like e.g.\n
       *  model calculators) sometimes does not differ between initialization and printing, this\n
       *  is also the place where third party code is initialized.
       */
      virtual void                               Initialize() const = 0;

      /*! \cond UML */
    private:
      /*!
       *  Value returned by Evaluate().
       */
      double                                     _chi2;
      /*!
       *  Stores the parameters.
       */
      Collection<ModelParameter*>                _collectionOfParameters;
      Collection<const Quantity*>                _collectionOfQuantities;

    private:
      /*!
       *  Returns the number of predictions of the model.
       */
      int                                        GetNumberOfPredictions() const;
      /*!
       *  Adds a parameter to the model.
       */
      void                                       AddParameter( ModelParameter* parameter );
      /*!
       *  Setup all parameters using a ptree.
       */
      void                                       InitializeParameters( const boost::property_tree::ptree& ptree );

    private:
      /*!
       *  Evaluates the chi2 function.
       */
      virtual double                             Evaluate() = 0;

      /*! \endcond UML */

  };

}

#endif // FITTINO_MODELBASE_H

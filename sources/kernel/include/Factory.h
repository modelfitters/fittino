/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        Factory.h                                                        *
*                                                                              *
* Description Factory class for creating file handlers, models, optimizers     *
*             and samplers                                                     *
*                                                                              *
* Authors     Bjoern  Sarrazin    <sarrazin@physik.uni-bonn.de>                *
*             Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_FACTORY_H
#define FITTINO_FACTORY_H

#include "Collection.h"
#include "Configuration.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  class AnalysisTool;
  class Chi2ContributionBase;
  class Configuration;
  class DataStorageBase;
  class ModelBase;
  class ModelCalculatorBase;
  class SLHADataStorageBase;
  class SLHAModelCalculatorBase;
  class LHCModelCalculator;
  class PhysicsModelBase;
  class PlotterBase;
  class SamplerBase;

  /*!
   *  \ingroup kernel
   *  \brief Factory class for creating file handlers, models, optimizers and samplers.
   */
  class Factory {

    public:
      /*!
       *  Standard constructor.
       */
      Factory();
      /*!
       *  Standard destructor.
       */
      ~Factory();

    public:
      ModelCalculatorBase*        CreateCalculator( const std::string& type, const PhysicsModelBase* model, const boost::property_tree::ptree& ptree ) const; 
      /*!
       *  Returns a concrete analysis tool\n
       *  Supported analysis tools are
       *  <ul>
       *    <li> GENETICALGORITHM\n
       *    <li> MINUIT\n
       *    <li> PARTICLESWARM\n
       *    <li> SIMULATEDANNEALING\n
       *  </ul>
       */
      AnalysisTool* const         CreateAnalysisTool( const std::string& type, ModelBase* model, const boost::property_tree::ptree& ptree ) const;
      /*!
       *  Returns a concrete data storage according to the file format passed as an argument.\n
       *  Supported file formats are
       *  <ul>
       *    <li> SLHA\n
       *    <li> XML\n
       *  </ul>
       */
      ModelBase* const            CreateModel( const std::string& type, const boost::property_tree::ptree& ptree ) const;
      /*!
       *  Returns a Prediction according to the type specified in the ptree.
       */
      PredictionBase* const       CreatePrediction( const boost::property_tree::ptree& ptree, const Fittino::ModelCalculatorBase* calculator );
      PredictionBase* const       CreatePrediction( const boost::property_tree::ptree& ptree, Fittino::SLHAModelCalculatorBase* calculator );
      /*!
       *  Returns an Observable with a prediction according to the type specified in the ptree.
       */
      Observable* const           CreateObservable( const boost::property_tree::ptree& ptree, const Fittino::Collection<Fittino::ModelCalculatorBase*>& calculators ) const;
      Observable* const           CreateObservable( const boost::property_tree::ptree& ptree, const Fittino::Collection<Fittino::PredictionBase*>& predictions, const Fittino::Collection<Fittino::ModelCalculatorBase*>& calculators ) const;
      /*!
       *  Returns a Chi2Contribution according to the type passed as an argument.
       */
      Chi2ContributionBase* const CreateChi2Contribution( const std::string& type, const boost::property_tree::ptree& ptree, const Fittino::Collection<Fittino::ModelCalculatorBase*>& calculators ) const;
      Chi2ContributionBase* const CreateChi2Contribution( const boost::property_tree::ptree& ptree, Fittino::LHCModelCalculator* calculator );
      SLHADataStorageBase*        CreateSLHAeaSLHADataStorage();

  };



}

#endif // FITTINO_FACTORY_H

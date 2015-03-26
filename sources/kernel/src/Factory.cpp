/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        Factory.cpp                                                      *
*                                                                              *
* Description Factory class for creating calculators, models, optimizers,      *
*             plotters and samplers                                            *
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

#include "AbsoluteUncertainty.h"
#include "AstroCalculator.h"
#include "AstroUncertainty.h"
#include "AstroFitCalculator.h"
#include "CheckMATECalculator.h"
#include "CheckVacuumCalculator.h"
#include "Chi2Calculator.h"
#include "ContourHistogramMaker.h"
#include "ContourPlotter.h"
#include "EdgeDetectionTool.h"
#include "Factory.h"

#ifdef FEYNHIGGS
#include "FeynHiggsCalculator.h"
#endif

#include "FeynHiggsRescalingCalculator.h"
#include "FormulaCalculator.h"
#include "GeneticAlgorithmOptimizer.h"
#include "HDim6Calculator.h"
#include "HECCalculator.h"
#include "HiggsBoundsHadXSCalculator.h"
#include "HiggsBoundsHECCalculator.h"
#include "HiggsBoundsSLHACalculator.h"
#include "HiggsSignalsHadXSCalculator.h"
#include "HiggsSignalsPartXSCalculator.h"
#include "HiggsSignalsSLHACalculator.h"
#include "HigherOrderMarkovChainSampler.h"
#include "LHCChi2Calculator.h"
#include "LHCLimitCalculator.h"
#include "MadGraphCalculator.h"
#include "MarkovChainSampler.h"
#include "MicromegasCalculator.h"
#include "MinuitOptimizer.h"
#include "NewHiggsSignalsHadXSCalculator.h"
#include "Observable.h"
#include "ParticleSwarmOptimizer.h"
#include "PhysicsModel.h"
#include "ProfileHistogramMaker.h"
#include "ProfilePlotter.h"
#include "RegressionCalculator.h"
#include "RelativeTheoryUncertainty.h"
#include "RosenbrockCalculator.h"
#include "RosenbrockModel.h"
#include "Simple1DHistogramMaker.h"
#include "Simple2DHistogramMaker.h"
#include "Simple3DHistogramMaker.h"
#include "SimpleDoubleCut.h"
#include "SimpleOptimizer.h"
#include "SimplePlotter.h"
#include "SimplePrediction.h"
#include "SimpleSampler.h"
#include "SimpleStringCut.h"
#include "SimulatedAnnealingOptimizer.h"
#include "SLHAeaSLHADataStorage.h"
#include "SPhenoSLHACalculator.h"
#include "SplineCut.h"
#include "SummaryHistogramMaker.h"
#include "SummaryPlotter.h"
#include "SuperIsoCalculator.h"
#include "TreeCalculator.h"
#include "TreeSampler.h"

Fittino::Factory::Factory() {

}

Fittino::Factory::~Factory() {

}

Fittino::CalculatorBase* const Fittino::Factory::CreateCalculator( const std::string& type, const ModelBase* model, boost::property_tree::ptree& ptree ) const {

    if ( type == "AstroCalculator" ) {

        return new AstroCalculator( model, ptree );

    }
    else if ( type == "AstroFitCalculator" ) {

        return new AstroFitCalculator( model, ptree );

    }
    else if ( type == "CheckMATECalculator" ) {

        return new CheckMATECalculator( model, ptree );

    }
    else if ( type == "CheckVacuumCalculator" ) {

        return new CheckVacuumCalculator( model, ptree );

    }
    else if ( type == "Chi2Calculator" ) {

        return new Chi2Calculator( model, ptree );

    }
    else if ( type == "FeynHiggsNativeCalculator" )  {

#if defined FEYNHIGGS

        return new FeynHiggsCalculator( model, ptree, "FeynHiggs" );

#else

        throw ConfigurationException( "Trying to use FeynHiggsCalculator but Fittino was built without FeynHiggs." );

#endif

    }
    else if ( type == "FeynHiggsRescalingCalculator" ) {

#if defined HIGGSSIGNALS_FOUND

        return new FeynHiggsRescalingCalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use FeynHiggsRescalingCalculator but Fittino was built without HiggsBounds." );

#endif

    }
    else if ( type == "FeynHiggsSLHACalculator" ) {

#if defined FEYNHIGGS

        return new FeynHiggsCalculator( model, ptree, "SLHA" );

#else

        throw ConfigurationException( "Trying to use FeynHiggsSLHACalculator but Fittino was built without FeynHiggs." );

#endif

    }
    else if ( type == "FormulaCalculator" ) {

        return new FormulaCalculator( model, ptree );

    }
    else if ( type == "HDim6Calculator" ) {

#if defined LHAPDF_FOUND  && defined GSL

        return new HDim6Calculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use HDim6Calculator but Fittino was built without LHAPDF or GSL." );

#endif

    }
    else if ( type == "HECCalculator" ) {

        return new HECCalculator( model, ptree );

    }
    else if ( type == "HiggsBoundsHadXSCalculator" ) {

#if defined HIGGSBOUNDS_FOUND

        return new HiggsBoundsHadXSCalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use HiggsBoundsHadXSCalculator but Fittino was built without HiggsBounds." );

#endif

    }
    else if ( type == "HiggsBoundsHECCalculator" ) {

#if defined HIGGSBOUNDS_FOUND

        return new HiggsBoundsHECCalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use HiggsBoundsHECCalculator but Fittino was built without HiggsBounds." );

#endif

    }
    else if ( type == "HiggsBoundsSLHACalculator" ) {

#if defined HIGGSBOUNDS_FOUND

        return new HiggsBoundsSLHACalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use HiggsBoundsSLHACalculator but Fittino was built without HiggsBounds." );

#endif

    }
    else if ( type == "HiggsSignalsCalculator" ) {

#if defined HIGGSBOUNDS_FOUND && defined HIGGSSIGNALS_FOUND

        return new HiggsSignalsHadXSCalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use HiggsSignalsHadXSCalculator but Fittino was built without HiggsBounds or HiggsSignals." );

#endif

    }
    else if ( type == "HiggsSignalsPartXSCalculator" ) {

#if defined HIGGSBOUNDS_FOUND && defined HIGGSSIGNALS_FOUND

        return new HiggsSignalsPartXSCalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use HiggsSignalsPartXSCalculator but Fittino was built without HiggsBounds or HiggsSignals." );

#endif

    }
    else if ( type == "HiggsSignalsSLHACalculator" ) {

#if defined HIGGSBOUNDS_FOUND && defined HIGGSSIGNALS_FOUND

        return new HiggsSignalsSLHACalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use HiggsSignalsSLHACalculator but Fittino was built without HiggsBounds or HiggsSignals." );

#endif

    }
    else if ( type == "LHCCalculator" ) {

        return new LHCLimitCalculator( model, ptree );

    }
    else if ( type == "LHCChi2Calculator" ) {

        return new LHCChi2Calculator( model, ptree );

    }
    else if ( type == "MadGraphCalculator" ) {

        return new MadGraphCalculator( model, ptree );

    }
    else if ( type == "MicromegasCalculator" ) {

#ifdef MICROMEGAS

        return new MicromegasCalculator( model, ptree );

#else

        throw ConfigurationException( " Trying to use MicromegasCalculator but Fittino was built without Micromegas." );

#endif

    }
    else if ( type == "NewHiggsSignalsHadXSCalculator" ) {

#if defined HIGGSBOUNDS_FOUND && defined HIGGSSIGNALS_FOUND

        return new NewHiggsSignalsHadXSCalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use NewHiggsSignalsHadXSCalculator but Fittino was built without HiggsBounds or HiggsSignals." );

#endif

    }
    else if ( type == "RegressionCalculator" ) {

        return new RegressionCalculator( model, ptree );

    }
    else if ( type == "RosenbrockCalculator" ) {

        return new RosenbrockCalculator( model, ptree );

    }
    else if ( type == "SPhenoCalculator" ) {

#ifdef SLHAEA

        return new SPhenoSLHACalculator( model, ptree );

#else

        throw ConfigurationException( " Trying to use SPhenoSLHACalculator but Fittino was built without SLHAea." );

#endif

    }
    else if ( type == "SuperIsoCalculator" ) {

#ifdef SUPERISO

        return new SuperIsoCalculator( model, ptree );

#else

        throw ConfigurationException( "Trying to use SuperIsoCalculator but Fittino was built without SuperIso." );

#endif

    }
    else if ( type == "TreeCalculator" ) {

        return new TreeCalculator( model, ptree );

    }
    else {

        throw ConfigurationException( "Calculator type " + type + " not known." );

    }

}

Fittino::CutBase* const Fittino::Factory::CreateCut( const std::string& type, Fittino::ModelBase* model, const boost::property_tree::ptree& ptree ) const {

    if ( type == "SimpleDoubleCut" ) {

        return new SimpleDoubleCut( model, ptree );

    }
    else if ( type == "SimpleStringCut" ) {

        return new SimpleStringCut( model, ptree );

    }
    else if ( type == "SplineCut" ) {

        return new SplineCut( model, ptree );

    }
    else {

        throw ConfigurationException( "Cut type " + type + " not known." );

    }

}

Fittino::ModelBase* const Fittino::Factory::CreateModel( const std::string& type,  boost::property_tree::ptree& ptree ) const {

    if ( type == "ModelBase" ) {

        return new ModelBase( ptree );

    }
    else if ( type == "PhysicsModel" ) {

        return new PhysicsModel( ptree );

    }
    else if ( type == "RosenbrockModel" ) {

        return new RosenbrockModel( ptree );

    }
    else {

        throw ConfigurationException( "Model type " + type + " not known." );

    }

}

Fittino::Observable* const Fittino::Factory::CreateObservable( const boost::property_tree::ptree& ptree, const Fittino::Collection<Fittino::CalculatorBase*>& calculators ) const {

    std::string type = ptree.get<std::string>( "PredictionType" );
    std::string calculatorName = ptree.get<std::string>( "CalculatorName", "NONE" );

    CalculatorBase *calculator = NULL;
    if ( calculatorName != "NONE" ) {

        calculator = calculators.At( ptree.get<std::string>( "CalculatorName", "NONE" ) );

    }

    if ( type == "Simple" ) {

        if ( calculator ) {

            return new Observable( ptree, new SimplePrediction( ptree, calculator ) );

        }
        else {

            return new Observable( ptree, new SimplePrediction( ptree.get<std::string>( "Name" ), "", 0. ) );

        }

    }
    else {

        throw ConfigurationException( "Prediction type" + type + " not known." );

    }

}

Fittino::Observable* const Fittino::Factory::CreateObservable( const boost::property_tree::ptree& ptree, const Fittino::Collection<Fittino::Quantity*>& predictions, const Fittino::Collection<Fittino::CalculatorBase*>& calculators ) const {

    std::string name = ptree.get<std::string>( "Name", "NONE" );
    for ( unsigned int i = 0; i < predictions.GetNumberOfElements(); ++i ) {

        if ( name == predictions.At( i )->GetName() ) {

            return new Observable( ptree, predictions.At( i ) );

        }

    }

    return CreateObservable( ptree, calculators );

}

Fittino::PlotterBase* const Fittino::Factory::CreatePlotter( const std::string& type, std::vector<TH1*>& histogramVector, const boost::property_tree::ptree& ptree ) const {

    if ( type == "ContourPlotter" ) {

        return new ContourPlotter( histogramVector, ptree );

    }
    else if ( type == "ProfilePlotter" ) {

        return new ProfilePlotter( histogramVector, ptree );

    }
    else if ( type == "SimplePlotter" ) {

        return new SimplePlotter( histogramVector, ptree );

    }
    else if ( type == "SummaryPlotter" ) {

        return new SummaryPlotter( histogramVector, ptree );

    }
    else {

        throw ConfigurationException( "Plotter type " + type + " not known." );

    }

}

Fittino::SLHADataStorageBase* const Fittino::Factory::CreateSLHAeaSLHADataStorage() const {

#if defined SLHAEA

    return new SLHAeaSLHADataStorage();

#else

    throw ConfigurationException( "Trying to use SLHAeaSLHADataStorage but Fittino was built without SLHAEA." );

#endif

}

Fittino::Tool* const Fittino::Factory::CreateTool( const std::string& type, ModelBase* model,
                                                   const boost::property_tree::ptree& ptree ) const {

    if ( type == "ContourHistogramMaker" ) {

        return new ContourHistogramMaker( model, ptree );

    }
    else if ( type == "EdgeDetectionTool" ) {

        return new EdgeDetectionTool( model, ptree );

    }
    else if ( type == "GeneticAlgorithmOptimizer" ) {

        return new GeneticAlgorithmOptimizer( model, ptree );

    }
    else if ( type == "HigherOrderMarkovChainSampler" ) {

        return new HigherOrderMarkovChainSampler( model, ptree );

    }
    else if ( type == "MarkovChainSampler" ) {

        return new MarkovChainSampler( model, ptree );

    }
    else if ( type == "MinuitOptimizer" ) {

        return new MinuitOptimizer( model, ptree );

    }
    else if ( type == "ParticleSwarmOptimizer" ) {

        return new ParticleSwarmOptimizer( model, ptree );

    }
    else if ( type == "ProfileHistogramMaker" ) {

        return new ProfileHistogramMaker( model, ptree );

    }
    else if ( type == "Simple1DHistogramMaker" ) {

        return new Simple1DHistogramMaker( model, ptree );

    }
    else if ( type == "Simple2DHistogramMaker" ) {

        return new Simple2DHistogramMaker( model, ptree );

    }
    else if ( type == "Simple3DHistogramMaker" ) {

        return new Simple3DHistogramMaker( model, ptree );

    }
    else if ( type == "SimpleOptimizer" ) {

        return new SimpleOptimizer( model, ptree );

    }
    else if ( type == "SimpleSampler" ) {

        return new SimpleSampler( model, ptree );

    }
    else if ( type == "SimulatedAnnealingOptimizer" ) {

        return new SimulatedAnnealingOptimizer( model, ptree );

    }
    else if ( type == "SummaryHistogramMaker" ) {

        return new SummaryHistogramMaker( model, ptree );

    }
    else if ( type == "TreeSampler" ) {

        return new TreeSampler( model, ptree );

    }
    else {

        throw ConfigurationException( "Tool type " + type + " not known." );

    }

}

Fittino::UncertaintyBase* const Fittino::Factory::CreateUncertainty( const std::string& type, const ModelBase* model, const Measurement* measurement, const boost::property_tree::ptree& ptree ) const {

    if ( type == "AbsoluteUncertainty" ) {

        return new AbsoluteUncertainty( model, measurement, ptree );

    }
    else if ( type == "AstroUncertainty" ) {

        return new AstroUncertainty( model, measurement, ptree );

    }
    else if ( type == "RelativeTheoryUncertainty" ) {

        return new RelativeTheoryUncertainty( model, measurement, ptree );

    }
    else {

        throw ConfigurationException( "Uncertainty of type " + type + " not known." );

    }

}

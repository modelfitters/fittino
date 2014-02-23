/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        AnalysisTool.h                                                   *
*                                                                              *
* Description Abstract class which provides basic functionality for analysis   *
*             tools like optimizers or samplers                                *
*                                                                              *
* Authors     Sebastian Heer        <s6seheer@uni-bonn.de>                     *
*             Mathias   Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>            *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include <iomanip>

#include "TTree.h"

#include "AnalysisTool.h"
#include "Messenger.h"
#include "ModelBase.h"
#include "Quantity.h"
#include "Observable.h"

Fittino::AnalysisTool::AnalysisTool( ModelBase *model, const boost::property_tree::ptree& ptree )
    : _chi2( ptree.get<double>( "Chi2", 1.e99 ) ),
      _iterationCounter( ptree.get<unsigned int>( "IterationCounter", 0 ) ),
      _name( ptree.get<std::string>( "Name", "" ) ),
      _model( model ),
      _outputFile( ( ptree.get<std::string>( "OutputFile", "Fittino.out.root" ) ).c_str(), "RECREATE" ),
      _randomGenerator( ptree.get<int>( "RandomSeed", 0 ) ),
      _metaDataTree( new TTree( ( ptree.get<std::string>( "MetaDataTree", "MetaDataTree" ) ).c_str(), ( ptree.get<std::string>( "MetaDataTree", "MetaDataTree" ) ).c_str() ) ),
      _tree( new TTree( ( ptree.get<std::string>( "OutputTree", "Tree" ) ).c_str(), ( ptree.get<std::string>( "OutputTree", "Tree" ) ).c_str() ) ) {

    _statusParameterVector.push_back( new Quantity( "Chi2",             "#chi^2",           _chi2,             0., 100.  ) );
    _statusParameterVector.push_back( new Quantity( "IterationCounter", "IterationCounter", _iterationCounter, 0., 1.e10 ) );

    _ptree = ptree;

}

Fittino::AnalysisTool::~AnalysisTool() {

}

void Fittino::AnalysisTool::PerformAnalysis() {

    AnalysisTool::InitializeAnalysisTool();
    AnalysisTool::ExecuteAnalysisTool();
    AnalysisTool::TerminateAnalysisTool();

}

boost::property_tree::ptree Fittino::AnalysisTool::GetPropertyTree() {

    UpdatePropertyTree();
    return _ptree;

}

int Fittino::AnalysisTool::GetNumberOfStatusParameters() const {

    return _statusParameterVector.size();

}

void Fittino::AnalysisTool::FillMetaDataTree() {

    _metaDataTree->Fill();

}

void Fittino::AnalysisTool::FillTree() {

    _tree->Fill();

}

void Fittino::AnalysisTool::PrintStatus() const {

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::INFO << Messenger::_dashedLine << Messenger::Endl;
    messenger << Messenger::INFO << Messenger::Endl;
    messenger << Messenger::INFO << "  Iteration step " << _iterationCounter << Messenger::Endl;

    _model->PrintStatus();

    messenger << Messenger::INFO << Messenger::Endl;
    messenger << Messenger::INFO << std::scientific << "    ----------------------------------------------------" << Messenger::Endl;
    messenger << Messenger::INFO << Messenger::Endl;
    messenger << Messenger::INFO << std::scientific << std::setprecision( 2 ) << "    Sum                                         " << _model->GetChi2() << Messenger::Endl;
    messenger << Messenger::INFO << Messenger::Endl;

}

void Fittino::AnalysisTool::UpdatePropertyTree() {

    _ptree.put( "Chi2", _chi2 );
    _ptree.put( "IterationCounter", _iterationCounter );

}

const std::vector<Fittino::Quantity*>* Fittino::AnalysisTool::GetStatusParameterVector() const {

    return &_statusParameterVector;

}

void Fittino::AnalysisTool::ExecuteAnalysisTool() {

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::ALWAYS << Messenger::_dashedLine << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;
    messenger << Messenger::ALWAYS << "  Running the " << _name << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;

    this->Execute();

}

void Fittino::AnalysisTool::InitializeAnalysisTool() {

    InitializeBranches();

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::ALWAYS << Messenger::_dashedLine << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;
    messenger << Messenger::ALWAYS << "  Initializing the " << _name << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;

    AnalysisTool::PrintConfiguration();

}

void Fittino::AnalysisTool::InitializeBranches() {

    // First initialize the branches of the output tree.

    for ( unsigned int i = 0; i < _model->GetCollectionOfQuantities().GetNumberOfElements(); ++i ) {

        _tree->Branch( _model->GetCollectionOfQuantities().At( i )->GetName().c_str(),
                       const_cast<double*>( &_model->GetCollectionOfQuantities().At( i )->GetValue() ) );
    }

    for ( unsigned int i = 0; i < _model->GetCollectionOfChi2Contributions().GetNumberOfElements(); ++i ) {

        _tree->Branch( _model->GetCollectionOfChi2Contributions().At( i )->GetName().c_str(),
                       const_cast<double*>( &_model->GetCollectionOfChi2Contributions().At( i )->GetValue() ) );
    }

    for ( unsigned int i = 0; i < GetNumberOfStatusParameters(); ++i ) {

        _tree->Branch( GetStatusParameterVector()->at( i )->GetName().c_str(),
                       const_cast<double*>( &GetStatusParameterVector()->at( i )->GetValue() ) );
    }

    // Now initialize the branches of the meta data tree.

    if ( _model->GetObservableVector() ) {

        for ( unsigned int i = 0; i < _model->GetObservableVector()->size(); ++i ) {

            std::string measuredValueBranchName = "measuredValue_" + _model->GetObservableVector()->at( i )->GetPrediction()->GetName();
            std::string uncertaintyBranchName = "uncertainty_" + _model->GetObservableVector()->at( i )->GetPrediction()->GetName();

            _metaDataTree->Branch( measuredValueBranchName.c_str(),
                                   const_cast<double*>( &( _model->GetObservableVector()->at( i )->GetMeasuredValue() ) ) );

            _metaDataTree->Branch( uncertaintyBranchName.c_str(),
                                   const_cast<double*>( &( _model->GetObservableVector()->at( i )->GetMeasuredError() ) ) );

        }

    }

}

void Fittino::AnalysisTool::PrintConfiguration() const {

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::ALWAYS << "   Configuration" << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;

    this->PrintSteeringParameters();

    messenger << Messenger::ALWAYS << Messenger::Endl;

}

void Fittino::AnalysisTool::TerminateAnalysisTool() {

    this->Terminate();

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::ALWAYS << Messenger::_dashedLine << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;
    messenger << Messenger::ALWAYS << "  Terminating the " << _name << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;

    this->PrintResult();

    this->WriteResultToFile();

    _outputFile.Close();

}

void Fittino::AnalysisTool::WriteResultToFile() const {

    _tree->Write();
    _metaDataTree->Write();

}

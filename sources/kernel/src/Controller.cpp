/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        Controller.cpp                                                   *
*                                                                              *
* Description Singleton class for controlling the execution flow of Fittino    *
*                                                                              *
* Authors     Philip  Bechtle     <philip.bechtle@desy.de>                     *
*             Klaus   Desch       <desch@physik.uni-bonn.de>                   *
*             Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              *
*             Peter   Wienemann   <wienemann@physik.uni-bonn.de>               *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#include <getopt.h>

#include <cstdlib>

#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Controller.h"
#include "Factory.h"
#include "InputException.h"
#include "ModelBase.h"
#include "Tool.h"
#include "RandomGenerator.h"

Fittino::Controller& Fittino::Controller::GetInstance() {

    static Controller instance;

    return instance;

}

void Fittino::Controller::ExecuteFittino() const {

    try {

        const Factory factory;

        boost::property_tree::ptree::value_type& modelNode = *( _inputPtree->get_child( "InputFile.Model" ).begin() );
        std::string modelType = modelNode.first;
        boost::property_tree::ptree& modelTree = modelNode.second;
        ModelBase* model = factory.CreateModel( modelType, modelTree );

        boost::property_tree::ptree::value_type& toolNode = *( _inputPtree->get_child( "InputFile.Tool" ).begin() );
        std::string toolType = toolNode.first;
        boost::property_tree::ptree& toolTree = toolNode.second;
        Tool* tool = factory.CreateTool( toolType, model, toolTree );

        if ( !_lockFileName.empty() ) {

            boost::property_tree::xml_writer_settings<char> settings( '\t', 1 );
            boost::property_tree::write_xml( _inputFileName, *_inputPtree, std::locale(), settings );

            _scopedLock->unlock();
            
        }
        
        
        tool->PerformTask();

        _outputPtree->put( "InputFile.VerbosityLevel", _inputPtree->get<std::string>( "InputFile.VerbosityLevel" ) );
        _outputPtree->put_child( "InputFile.Model." + modelType, model->GetPropertyTree() );
        _outputPtree->put_child( "InputFile.Tool." + toolType , tool->GetPropertyTree() );

        delete tool;
        delete model;

    }
    catch ( const ConfigurationException& configurationException ) {

        std::cout << "\n" << configurationException.what() << "\n" << std::endl;
        exit( EXIT_FAILURE );

    }

}

void Fittino::Controller::InitializeFittino( int argc, char** argv ) {

    try {

        // If Fittino is called without options or arguments print a help text with further
        // instructions and exit.

        if ( argc == 1 ) {

            Controller::PrintHelp();
            exit( EXIT_SUCCESS );

        }

        // Otherwise handle given command line options.

        Controller::HandleOptions( argc, argv );

        // Print a welcome logo.

        Controller::PrintLogo();

        // A given XML input file is parsed, and the configuration items are stored in the input
        // property tree.

        Controller::CheckInputFileFormat();

        if ( ! _lockFileName.empty() ) {

          _fileLock = new boost::interprocess::file_lock( _lockFileName.c_str() );
          _scopedLock = new boost::interprocess::scoped_lock<boost::interprocess::file_lock>( *_fileLock );

        }

        boost::property_tree::read_xml( _inputFileName,
                                        *_inputPtree,
                                        boost::property_tree::xml_parser::trim_whitespace |
                                        boost::property_tree::xml_parser::no_comments );

        // Set the verbosity level of the output text.

        std::string verbosityLevel = _inputPtree->get<std::string>( "InputFile.VerbosityLevel" );
        Messenger::GetInstance().SetVerbosityLevel( verbosityLevel );

        double randomSeed = _inputPtree->get<double>( "InputFile.RandomSeed", 0 );
        if( randomSeed != 0 ) {
            
            RandomGenerator::GetInstance()->SetSeed( randomSeed );
        
        }

    }
    catch ( const InputException& inputException ) {

        std::cout << "\n" << inputException.what() << "\n" << std::endl;
        exit( EXIT_FAILURE );

    }

}

void Fittino::Controller::TerminateFittino() const {

    boost::property_tree::xml_writer_settings<char> settings( '\t', 1 );
    boost::property_tree::write_xml( "FittinoInterfaceFile.xml", *_outputPtree, std::locale(), settings );

}

Fittino::Controller* Fittino::Controller::_instance = 0;

Fittino::Controller::Controller()
    : _inputFileName( "" ),
      _inputPtree( new boost::property_tree::ptree() ),
      _outputPtree( new boost::property_tree::ptree() ) {

  _fileLock = 0;
  _scopedLock = 0;
  _lockFileName = "";

}

Fittino::Controller::~Controller() {

    delete _inputPtree;
    delete _outputPtree;
    delete _scopedLock;
    delete _fileLock;

}

void Fittino::Controller::CheckInputFileFormat() const {

    try {

        if ( _inputFileName.length() == 0 ) {

            throw InputException( "The input file has to be specified with the option flag -i/--input-file" );

        }
        else if ( _inputFileName.compare( _inputFileName.length() - 4, 4, ".xml" ) ) {

            throw InputException( "Input file suffix must be .xml (XML format)." );

        }

    }
    catch ( const InputException& inputException ) {

        std::cout << "\n" << inputException.what() << "\n" << std::endl;
        exit( EXIT_FAILURE );

    }

}

void Fittino::Controller::HandleOptions( int argc, char** argv ) {

    // Use getopt() to handle given command line options. For more informations on getopt() have a
    // look at the manpage of getopt(3).

    static struct option options[] = {

        {"help",       no_argument,       0, 'h'},
        {"input-file", required_argument, 0, 'i'},
        {"lock-file",  required_argument, 0, 'l'},
        {0,            0,                 0,  0 }

    };

    int optionIndex = 0;
    int optionCode = -1;
    opterr = 0;

    while ( true ) {

        optionCode = getopt_long( argc, argv, ":hi:l:", options, &optionIndex );

        if ( optionCode == -1 ) break;

        switch ( optionCode ) {

            case 'h':
                Controller::PrintHelp();
                exit( EXIT_SUCCESS );

            case 'i':
                _inputFileName = std::string( optarg );
                continue;

            case 'l':  
              _lockFileName = std::string( optarg );
              continue;

            case ':':
                throw InputException( "Missing option parameter." );

            default:
                throw InputException( "Unknown option(s)" );

        }

    }

}

void Fittino::Controller::PrintHelp() const {

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::ALWAYS << Messenger::Endl;
    messenger << Messenger::ALWAYS << "Usage: fittino [OPTION(S)]" << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;
    messenger << Messenger::ALWAYS << "Supported options are:" << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;
    messenger << Messenger::ALWAYS << "  -h, --help" << Messenger::Endl;
    messenger << Messenger::ALWAYS << "      Fittino prints this message." << Messenger::Endl;
    messenger << Messenger::ALWAYS << "  -i, --input-file=FILE" << Messenger::Endl;
    messenger << Messenger::ALWAYS << "      Fittino uses the input file FILE. The input file suffix must" << Messenger::Endl;
    messenger << Messenger::ALWAYS << "      be .xml (XML format)." << Messenger::Endl;
    messenger << Messenger::ALWAYS << "  -l, --lock-file=FILE" << Messenger::Endl;
    messenger << Messenger::ALWAYS << "      Fittino uses the file FILE for inter process locking." << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;

}

void Fittino::Controller::PrintLogo() const {

    Messenger& messenger = Messenger::GetInstance();

    messenger << Messenger::ALWAYS << Messenger::_dashedLine << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;
    messenger << Messenger::ALWAYS << "  Welcome to Fittino" << Messenger::Endl;
    messenger << Messenger::ALWAYS << Messenger::Endl;

}

/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        Controller.h                                                     *
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

/*!
 *  \todo Long-term: Document the code.
 *  \todo Long-term: Write a comprehensive developer's guide.
 *  \todo Long-term: Write a comprehensive user's guide.
 *  \todo Short-term: In xml input files use subtags instead of attributes
 */

#ifndef FITTINO_CONTROLLER_H
#define FITTINO_CONTROLLER_H

#include <string>

#include "PtreeForwardDeclaration.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

  /*!
   *  \defgroup kernel
   */
  /*!
   *  \ingroup kernel
   *  \brief Singleton class for controlling the execution flow of Fittino.
   *
   *  The instance of the Controller class is the first object that is created at the beginning of\n
   *  the execution of Fittino. Depending on the user configuration the controller creates\n
   *  instances of further objects (either directly or with the help of a factory) and advises\n
   *  them to perform the specified tasks, hereby controlling the program's overall execution\n
   *  flow.
   *
   *  The controller adresses its purpose in three distinct phases: At first, the controller\n
   *  initializes Fittino by letting an input file interpreter read in user specified options.\n
   *  After that further objects are created and put into action according to the user's\n
   *  configuration. At the end, the controller provides the controlled termination of Fittino.
   */
  class Controller {

    public:
      /*!
       *  Returns a static pointer to the unique instance of this class.
       */
      static Controller*           GetInstance();

    public:
      /*!
       *  Executes Fittino. For that purpose a model inheriting from ModelBase and the required\n
       *  analysis tools inheriting from AnalysisTool are created and put into action.
       */
      void                         ExecuteFittino() const;
      /*!
       *  Initializes Fittino. Takes as input the command line arguments specified by the user\n
       *  while invocing Fittino and stores them for future reference.
       */
      void                         InitializeFittino( int argc, char** argv );
      /*!
       *  Provides the controlled termination of Fittino.
       */
      void                         TerminateFittino() const;

    private:
      /*!
       *  Pointer to the unique instance of this class.
       */
      static Controller*           _instance;

      /*! \cond UML */
    private:
      /*!
       *  External seed as passed as an argument to Fittino.
       */
      int                          _randomSeed;
      /*!
       *  The name of the input file.
       */
      std::string                  _inputFileName;
      /*!
       *  The name of the data file.
       */
      std::string                  _dataFileName;
      boost::property_tree::ptree* _inputPtree;
      boost::property_tree::ptree* _outputPtree;

    private:
      /*!
       *  Standard constructor.
       */
      Controller();
      /*!
       *  Standard destructor.
       */
      ~Controller();
      /*!
       *  Check the format of the input file. The supported input file format is XML.
       */
      void                         CheckInputFileFormat() const;
      void                         HandleOptions( int argc, char** argv );
      /*!
       *  When Fittino is called without arguments or with the -h/--help option this method\n
       *  prints a help screen with further instructions on how to use Fittino.
       */
      void                         PrintHelp() const;
      /*!
       *  Prints a welcome logo.
       */
      void                         PrintLogo() const;

      /*! \endcond UML */

  };

}

#endif // FITTINO_CONTROLLER_H

  # $Id$ #

################################################################################
#                                                                              #
# Project     Fittino - A SUSY Parameter Fitting Package                       #
#                                                                              #
# File        CMakeLists.txt                                                   #
#                                                                              #
# Description Cmake configuration file for the calculators directory.          #
#                                                                              #
# Authors     Bjoern  Sarrazin    <sarrazin@physik.uni-bonn.de>                #
#             Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              #
#                                                                              #
# Licence     This program is free software; you can redistribute it and/or    #
#             modify it under the terms of the GNU General Public License as   #
#             published by the Free Software Foundation; either version 3 of   #
#             the License, or (at your option) any later version.              #
#                                                                              #
################################################################################

FILE(GLOB includes ./include/*.h)
LIST( APPEND includes "include/AxionCalculator.h" )
LIST( APPEND includes "include/Correlation.h" )
LIST( APPEND includes "include/GM2CalcCalculator.h" )
LIST( APPEND includes "include/InputFile.h" )
LIST( APPEND includes "include/LinearInterpolationCalculator.h" )
LIST( APPEND includes "include/LinearInterpolationCalculatorBase.h" )
LIST( APPEND includes "include/HerwigppCalculator.h" )
LIST( APPEND includes "include/SModelSCalculator.h" )
LIST( APPEND includes "include/HiggsSignalsSMCalculator.h" )
LIST( APPEND includes "include/HiggsSignalsCalculator.h" )
LIST( APPEND includes "include/CheckMATE2Calculator.h" )
LIST( APPEND includes "include/NLLFastCalculator.h" )
LIST( APPEND includes "include/HepMCSplitCalculator.h" )
LIST( APPEND includes "include/MadGraph2Calculator.h" )
LIST( APPEND includes "include/MadGraph5Calculator.h" )
LIST( APPEND includes "include/MadGraph5EWCalculator.h" )
LIST( APPEND includes "include/ClsCalculator_2.h" )
LIST( APPEND includes "include/MinimumCalculator.h" )
LIST( APPEND includes "include/ProspinoCalculator.h" )
LIST( APPEND includes "include/NLLFast2Calculator.h" )
LIST( APPEND includes "include/CheckMATE3Calculator.h" )
LIST( APPEND includes "include/FormulaPropCalculator.h" )
LIST( APPEND includes "include/EventNumber2Calculator.h" )
LIST( APPEND includes "include/EventNumber3Calculator.h" )
LIST( APPEND includes "include/LHCNeuralNet2Calculator.h" )
LIST( APPEND includes "include/CheckMATE4Calculator.h" )
LIST( APPEND includes "include/HepMCSplitCalculator.h" )
LIST( APPEND includes "include/RPVCalculator.h" )


  AUX_SOURCE_DIRECTORY(src sources)

INCLUDE_DIRECTORIES(../../external)
INCLUDE_DIRECTORIES(./include)
INCLUDE_DIRECTORIES(../exceptions/include)
INCLUDE_DIRECTORIES(../datastorages/include)
INCLUDE_DIRECTORIES(../models/include)
INCLUDE_DIRECTORIES(../kernel/include)
INCLUDE_DIRECTORIES(../variables/include)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIRS})


LIST( APPEND sources "src/AbsoluteUncertainty.cpp" )
LIST( APPEND includes "src/AxionCalculator.cpp" )
LIST( APPEND sources "src/Correlation.cpp" )
LIST( APPEND sources "src/InputFile.cpp" )
LIST( APPEND sources "src/FormulaQuantity.cpp" )
LIST( APPEND sources "src/HepMCSplitCalculator.cpp" )
LIST( APPEND sources "src/RelativeTheoryUncertainty.cpp" )
LIST( APPEND sources "src/LinearInterpolationCalculator.cpp" )
LIST( APPEND sources "src/LinearInterpolationCalculatorBase.cpp" )
LIST( APPEND sources "src/AstroUncertainty.cpp" )
LIST( APPEND sources "src/SModelSCalculator.cpp" )
LIST( APPEND sources "src/HiggsSignalsSMCalculator.cpp" )
LIST( APPEND sources "src/HiggsSignalsCalculator.cpp" )
LIST( APPEND sources "src/CheckMATE2Calculator.cpp" )
LIST( APPEND sources "src/HerwigppCalculator.cpp" )
LIST( APPEND sources "src/MinimumCalculator.cpp" )
LIST( APPEND sources "src/NLLFastCalculator.cpp" )
LIST( APPEND sources "src/HepMCSplitCalculator.cpp" ) 
LIST( APPEND sources "src/MadGraph2Calculator.cpp" )
LIST( APPEND sources "src/MadGraph5Calculator.cpp" )
LIST( APPEND sources "src/MadGraph5EWCalculator.cpp" )
LIST( APPEND sources "src/ClsCalculator_2.cpp" )
LIST( APPEND sources "src/ProspinoCalculator.cpp" )
LIST( APPEND sources "src/NLLFast2Calculator.cpp" )
LIST( APPEND sources "src/CheckMATE3Calculator.cpp" )
LIST( APPEND sources "src/FormulaPropCalculator.cpp" )
LIST( APPEND sources "src/EventNumber2Calculator.cpp" )
LIST( APPEND sources "src/EventNumber3Calculator.cpp" )
LIST( APPEND sources "src/LHCNeuralNet2Calculator.cpp" )
LIST( APPEND sources "src/CheckMATE4Calculator.cpp" )
LIST( APPEND sources "src/GM2CalcCalculator.cpp" )
LIST( APPEND sources "src/RPVCalculator.cpp" )

IF(PYTHON2_FOUND OR PYTHON3_FOUND)

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})

ENDIF()

  IF( NOT PYTHON2_FOUND)

    list( APPEND python2files "src/SModelSCalculator.cpp" )
    list( APPEND python2files "src/LHCNeuralNetCalculator.cpp" )

    set_source_files_properties( ${python2files} PROPERTIES HEADER_FILE_ONLY TRUE )


ENDIF()

  IF( NOT PYTHON3_FOUND )

      set_source_files_properties( "src/FlavioCalculator.cpp" PROPERTIES HEADER_FILE_ONLY TRUE )

  ENDIF()

IF(FEYNHIGGS_FOUND)

    INCLUDE_DIRECTORIES(BEFORE ${FEYNHIGGS_INCLUDE_DIRS})

    IF(FEYNHIGGS_QUAD)

        ADD_DEFINITIONS(-DQUAD=1)

    ENDIF()

ELSE()

    LIST(REMOVE_ITEM sources "src/FeynHiggsBosonicChannel.cpp")
    LIST(REMOVE_ITEM sources "src/FeynHiggsCalculator.cpp")
    LIST(REMOVE_ITEM sources "src/FeynHiggsChannel.cpp")
    LIST(REMOVE_ITEM sources "src/FeynHiggsCrossSection.cpp")
    LIST(REMOVE_ITEM sources "src/FeynHiggsFermionicChannel.cpp")
    LIST(REMOVE_ITEM sources "src/FeynHiggsNeutralHiggsCrossSection.cpp")

ENDIF()

  IF( GM2CALC_FOUND AND EIGEN3_FOUND )

      INCLUDE_DIRECTORIES( ${GM2CALC_INCLUDE_DIRS} )
      INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR}   )

  ELSE()

      LIST( REMOVE_ITEM sources "src/GM2CalcCalculator.cpp")

  ENDIF()

IF( HEPMC_FOUND )

    INCLUDE_DIRECTORIES( ${HEPMC_INCLUDE_DIRS} )

ELSE()

    LIST(REMOVE_ITEM sources "src/HepMCSplitCalculator.cpp" )

ENDIF()

IF(NOT HIGGSBOUNDS_FOUND OR NOT HIGGSSIGNALS_FOUND)

    LIST(REMOVE_ITEM sources "src/HiggsSignalsHadXSCalculator.cpp")
    LIST(REMOVE_ITEM sources "src/NewHiggsSignalsHadXSCalculator.cpp")

ENDIF()

IF( NOT HIGGSSIGNALS_FOUND )

    LIST(REMOVE_ITEM sources "src/HiggsSignalsSMCalculator.cpp")
    LIST(REMOVE_ITEM sources "src/HiggsSignalsCalculator.cpp")


ENDIF()

IF(NOT HIGGSBOUNDS_FOUND)

    LIST(REMOVE_ITEM sources "src/HiggsBoundsCalculator.cpp")
    LIST(REMOVE_ITEM sources "src/HiggsBoundsSLHACalculator.cpp")
    LIST(REMOVE_ITEM sources "src/HiggsBoundsHECCalculator.cpp")
    LIST(REMOVE_ITEM sources "src/HiggsBoundsHadXSCalculator.cpp")
    LIST(REMOVE_ITEM sources "src/FeynHiggsRescalingCalculator.cpp")

ENDIF()

IF(NOT MICROMEGAS_FOUND)

    LIST(REMOVE_ITEM sources "src/MicromegasCalculator.cpp")
    LIST(REMOVE_ITEM sources "src/MicromegasWrapper.cpp")
    
ENDIF()

IF(NOT SLHAEA_FOUND)

    LIST(REMOVE_ITEM sources "src/SPhenoSLHACalculator.cpp")

ENDIF()

IF(NOT HIGGSBOUNDS_FOUND OR NOT HIGGSSIGNALS_FOUND)

    LIST(REMOVE_ITEM sources "src/HiggsSignalsSLHACalculator.cpp")

ENDIF()

IF(NOT LHAPDF_FOUND OR NOT GSL_FOUND)

    LIST(REMOVE_ITEM sources "src/HDim6Calculator.cpp")

ENDIF()

#IF(SUPERISO_FOUND)
#
#    INCLUDE_DIRECTORIES(${SUPERISO_INCLUDE_DIRS})
#
#ELSE()
#
#    LIST(REMOVE_ITEM sources "src/SuperIsoCalculator.cpp")
#
#ENDIF()

ADD_LIBRARY(calculators ${includes} ${sources})

IF(LHAPDF_FOUND)

    INCLUDE_DIRECTORIES(${LHAPDF_INCLUDE_DIRS})

ENDIF()

IF(MICROMEGAS_FOUND)

  INCLUDE_DIRECTORIES(${MICROMEGAS_INCLUDE_DIRS})  

ENDIF()

if ( USE_WCxf )

    add_definitions( -DWCxf )

endif()





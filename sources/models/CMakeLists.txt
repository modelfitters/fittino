# $Id$ #

################################################################################
#                                                                              #
# Project     Fittino - A SUSY Parameter Fitting Package                       #
#                                                                              #
# File        CMakeLists.txt                                                   #
#                                                                              #
# Description Cmake configuration file for the models directory.               #
#                                                                              #
# Authors     Mathias Uhlenbrock  <uhlenbrock@physik.uni-bonn.de>              #
#                                                                              #
# Licence     This program is free software; you can redistribute it and/or    #
#             modify it under the terms of the GNU General Public License as   #
#	      published by the Free Software Foundation; either version 3 of   #
#	      the License, or (at your option) any later version.              #
#                                                                              #
################################################################################

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(./include)
INCLUDE_DIRECTORIES(../calculators/include)
INCLUDE_DIRECTORIES(../datastorages/include)
INCLUDE_DIRECTORIES(../exceptions/include)
INCLUDE_DIRECTORIES(../kernel/include)
INCLUDE_DIRECTORIES(../optimizers/include)
INCLUDE_DIRECTORIES(../samplers/include)
INCLUDE_DIRECTORIES(../variables/include)

FILE( GLOB includes ./include/*.h )
AUX_SOURCE_DIRECTORY(src sources)

IF(NOT HIGGSBOUNDS_FOUND AND NOT HIGGSSIGNALS_FOUND)
  LIST(REMOVE_ITEM sources "src/HECModel.cpp")
ENDIF(NOT HIGGSBOUNDS_FOUND AND NOT HIGGSSIGNALS_FOUND)

ADD_LIBRARY(models ${includes} ${sources})



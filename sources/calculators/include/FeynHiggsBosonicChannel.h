/* $Id$ */

/*******************************************************************************
 *                                                                              *
 * Project     Fittino - A SUSY Parameter Fitting Package                       *
 *                                                                              *
 * File        FeynHiggsFermionicChannel                                        *
 *                                                                              *
 * Description Wrapper class for FeynHiggs                                      *
 *                                                                              *
 * Authors     Bjoern Sarrazin  <sarrazin@physik.uni-bonn.de>                   *
 *                                                                              *
 * Licence     This program is free software; you can redistribute it and/or    *
 *             modify it under the terms of the GNU General Public License as   *
 *             published by the Free Software Foundation; either version 3 of   *
 *             the License, or (at your option) any later version.              *
 *                                                                              *
 *******************************************************************************/

#ifndef FITTINO_FEYNHIGGSBOSONICCHANNEL_H
#define FITTINO_FEYNHIGGSBOSONICCHANNEL_H

#include "FeynHiggsChannel.h"
#include "FeynHiggsTypes.h"

/*!
 *  \brief Fittino namespace.
 */

namespace Fittino {

    class SimplePrediction;
    /*!
     *  \ingroup calculators
     *  \brief Wrapper class for FeynHiggs.
     */
    class FeynHiggsBosonicChannel : public FeynHiggsChannel {

    public:
        /*!
         *  Standard constructor.
         */
        FeynHiggsBosonicChannel( FHRealType* gammas, FHRealType* gammasms, FHComplexType* couplings, FHComplexType* couplingsms, std::string higgsName, std::string channelName, int channelNumber, bool SM );
        /*!
         *  Standard destructor.
         */
        virtual ~FeynHiggsBosonicChannel();

        void CalculatePredictions();
        /*! \cond UML */
    private:

        double _model_g_Abs;
        double _model_g_Arg;

        double _sm_g_Abs;
        double _sm_g_Arg;

        double _normSM_g_Abs;
        double _normSM_g_Abs2;

        /*! \endcond UML */
        
    };
    
}

#endif

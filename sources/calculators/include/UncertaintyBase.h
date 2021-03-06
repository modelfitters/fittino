/* $Id$ */

/*******************************************************************************
*                                                                              *
* Project     Fittino - A SUSY Parameter Fitting Package                       *
*                                                                              *
* File        UncertaintyBase.h                                                *
*                                                                              *
* Description                                                                  *
*                                                                              *
* Authors                                                                      *
*                                                                              *
*                                                                              *
* Licence     This program is free software; you can redistribute it and/or    *
*             modify it under the terms of the GNU General Public License as   *
*             published by the Free Software Foundation; either version 3 of   *
*             the License, or (at your option) any later version.              *
*                                                                              *
*******************************************************************************/

#ifndef FITTINO_UNCERTAINTYBASE_H
#define FITTINO_UNCERTAINTYBASE_H

#include "Measurement.h"
#include "Quantity.h"

/*!
 *  \brief Fittino namespace.
 */
namespace Fittino {

    class Measurement;
    class ModelBase;

/*!
 *  \ingroup calculators
 *  \brief .
 */

    class UncertaintyBase : public Quantity {

    public:
        UncertaintyBase(ModelBase const *model, const Measurement *observable, const boost::property_tree::ptree &ptree);
        virtual ~UncertaintyBase();

        virtual void Update() = 0;

    protected:
        const Measurement* _observable;
        double _value;

    };

}

#endif

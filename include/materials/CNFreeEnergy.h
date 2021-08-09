//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeParsedMaterialHelper.h"
#include "ExpressionBuilder.h"

// Forward Declarations

/**
 * Material class that creates regular solution free energy with the expression builder
 * and uses automatic differentiation to get the derivatives
 * \f$ F = (eta-1)^2 (E_f c + k_bT (c\log c + (1 - c)\log(1 - c))) + eta^2 ((c-1)^2)\f$.
 */
class CNFreeEnergy : public DerivativeParsedMaterialHelper, public ExpressionBuilder
{
public:
  static InputParameters validParams();

  CNFreeEnergy(const InputParameters & parameters);

protected:
  /// Coupled variable value for the concentration \f$ c \f$.
  EBTerm _c;
  EBTerm _eta;

  /// Coupled temperature variable \f$ T \f$
  EBTerm _T;

  /// Prefactor
  const Real _Ef;
  
  //const Real _Em;

  //const Real _N;

  //const Real _D0;

  /// Boltzmann constant
  const Real _kB;
};

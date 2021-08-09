//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CNFreeEnergy.h"

registerMooseObject("kookaburraApp", CNFreeEnergy);

InputParameters
CNFreeEnergy::validParams()
{
  InputParameters params = DerivativeParsedMaterialHelper::validParams();
  params.addClassDescription("Material that implements the free energy of a regular solution");
  params.addRequiredCoupledVar("c", "Concentration variable");
  params.addRequiredCoupledVar("eta", "Order parameter variable");
  params.addCoupledVar("T", 300, "Temperature variable");
  params.addParam<Real>("Ef", 1.0, "Formation energy in eV");
  //params.addParam<Real>("Em", 1.0, "Migration energy in eV");
  //params.addParam<Real>("D0", 1.0, "D0 in m^2/s");
  params.addParam<Real>("kB", 8.6173324e-5, "Boltzmann constant");
  params.addParam<Real>(
      "log_tol", "If specified logarithms are evaluated using a Taylor expansion below this value");
  return params;
}

CNFreeEnergy::CNFreeEnergy(const InputParameters & parameters)
  : DerivativeParsedMaterialHelper(parameters),
    _c("c"),
    _eta("eta"),
    _T("T"),
    _Ef(getParam<Real>("Ef")),
    //_Em(getParam<Real>("Em")),
    //_N(1.0),
    //_D0(getParam<Real>("D0")),
    _kB(getParam<Real>("kB"))
    //_M(declareProperty<Real>("M")),
    //_grad_M(declareProperty<RealGradient>("grad_M")),
    //_kappa(declareProperty<Real>("kappa")),
    //_c_eq(declareProperty<Real>("c_eq")),

{
  EBFunction free_energy;
  // Definition of the free energy for the expression builder
  free_energy(_c,_eta,_T) = (_eta-1)*(_eta-1)*(_Ef * _c + _kB * _T * (_c * log(_c) + (1.0 - _c) * log(1.0 - _c)))
                            + _eta*_eta*( (_c-1) * (_c-1) );

  // Use Taylor expanded logarithm?
  if (isParamValid("log_tol"))
    free_energy.substitute(EBLogPlogSubstitution(getParam<Real>("log_tol")));

  // Parse function for automatic differentiation
  functionParse(free_energy);
}

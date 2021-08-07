//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CNInserterAux.h"
#include "CNInserterBase.h"

registerMooseObject("kookaburraApp", CNInserterAux);

InputParameters
CNInserterAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Project the DiscreteNucleationMap state onto an AuxVariable");
  params.addRequiredParam<UserObjectName>("map", "DiscreteNucleationMap user object");
  params.addParam<Real>("no_nucleus_value", 0.0, "Variable value indicating no nucleus is present");
  params.addParam<Real>(
      "nucleus_value", 1.0, "Variable value indicating the presence of a nucleus");
  return params;
}

CNInserterAux::CNInserterAux(const InputParameters & params)
  : AuxKernel(params),
    _map(getUserObject<CNNucleationMap>("map")),
    _v0(getParam<Real>("no_nucleus_value")),
    _v1(getParam<Real>("nucleus_value"))
{
  if (isNodal())
    paramError("variable", "This kernel must operate on an elemental AuxVariable.");
}

void
CNInserterAux::precalculateValue()
{
  // check if a nucleation event list is available for the current element
  _nucleus = &_map.nuclei(_current_elem);
}

Real
CNInserterAux::computeValue()
{
  return (*_nucleus)[_qp] * (_v1 - _v0) + _v0;
}

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IontrackNucleationForce.h"
#include "IontrackNucleationMap.h"
#include "Function.h"

registerMooseObject("kookaburraApp", IontrackNucleationForce);

InputParameters
IontrackNucleationForce::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Term for inserting grain nuclei or phases in non-conserved order parameter fields");
  params.addRequiredParam<UserObjectName>("map", "DiscreteNucleationMap user object");
  params.addParam<MaterialPropertyName>("mask", "Material property defining the mask");
  // params.addCoupledVar("args", "Vector of nonlinear variable arguments this object depends on");
  params.addParam<Real>("no_nucleus_value", 0.0, "Variable value indicating no nucleus is present");
  params.addParam<Real>(
      "nucleus_value", 1.0, "Variable value indicating the presence of a nucleus");
  return params;
}

IontrackNucleationForce::IontrackNucleationForce(const InputParameters & params)
  : Kernel(params),
    _map(getUserObject<IontrackNucleationMap>("map")),
    _mask(getMaterialProperty<Real>("mask")),
    //_dmaskdv(getMaterialPropertyDerivative<Real>("mask", _v_name)),
    //_dmaskdarg(_n_args),
    _v0(getParam<Real>("no_nucleus_value")),
    _v1(getParam<Real>("nucleus_value"))
{
  //  // Get derivatives of mask wrt coupled variables
  //for (unsigned int i = 0; i < _n_args; ++i)
  //  _dmaskdarg[i] = &getMaterialPropertyDerivative<Real>("mask", i);
}


void
IontrackNucleationForce::initialSetup()
{
  //validateNonlinearCoupling<Real>("mask");
}

void
IontrackNucleationForce::precalculateResidual()
{
  // check if a nucleation event list is available for the current element
  _nucleus = &_map.nuclei(_current_elem);
}

Real
IontrackNucleationForce::computeQpResidual()
{
  return -((*_nucleus)[_qp] * (_v1 - _v0) + _v0) * _test[_i][_qp] * _mask[_qp];
}

/*
Real
IontrackNucleationForce::computeQpJacobian()
{
  return _dmaskdv[_qp] * (-((*_nucleus)[_qp] * (_v1 - _v0) + _v0) * _test[_i][_qp]) * _phi[_j][_qp];
}

Real
IontrackNucleationForceForce::computeQpOffDiagJacobian(unsigned int jvar)
{
  const unsigned int cvar = mapJvarToCvar(jvar);
  return (*_dmaskdarg[cvar])[_qp] * (-((*_nucleus)[_qp] * (_v1 - _v0) + _v0) * _test[_i][_qp]) * _phi[_j][_qp];
}
*/

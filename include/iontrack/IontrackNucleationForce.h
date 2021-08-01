//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

//#include "LangevinNoise.h"
//#include "ConservedNoiseBase.h"

// Forward declaration
class IontrackNucleationMap;

/**
 * Free energy penalty contribution to force the nucleation of subresolution particles
 */
//class MaskedBodyForce : public DerivativeMaterialInterface<JvarMapKernelInterface<BodyForce>>
class IontrackNucleationForce : public Kernel
{
public:
  static InputParameters validParams();

  IontrackNucleationForce(const InputParameters & params);
  virtual void initialSetup();
  
  void precalculateResidual() override;
  Real computeQpResidual() override;
  //virtual Real computeQpJacobian();
  //virtual Real computeQpOffDiagJacobian(unsigned int jvar);
protected:
  /// UserObject providing a map of currently active nuclei
  const IontrackNucleationMap & _map;
  const MaterialProperty<Real> & _mask;
  
  /// nucleus data for the current element
  const std::vector<Real> * _nucleus;
  ///@{ Bounds for the returned values
  const Real _v0;
  const Real _v1;
  ///@}

  ///Noise to add intisde the tracks
  //const ConservedNoiseInterface & _noise;
  //const Real &_amplitude;
  
  /*
    /// derivative of the mask wrt the kernel's nonlinear variable
  const MaterialProperty<Real> & _dmaskdv;

  ///  Reaction rate derivatives w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dmaskdarg;
  */
};

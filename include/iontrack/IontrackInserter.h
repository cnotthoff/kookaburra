//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DiscreteNucleationInserterBase.h"

#include "GeneralPostprocessor.h"
// Forward Declarations
class IonFluencePostprocessor;

/**
 * This UserObject manages the insertion and expiration of nuclei in the simulation
 * domain it manages a list of nuclei with their insertion times and their center
 * positions. A DiscreteNucleationMap is needed to enable the DiscreteNucleation
 * material to look up if a nucleus is present at a given element/qp.
 */
class IontrackInserter : public DiscreteNucleationInserterBase
{
public:
  static InputParameters validParams();

  IontrackInserter(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual void finalize();

  const Real & getRate() const { return _nucleation_rate; }
  
  const Real & getFluence() const { return _N_ions_old; }
protected:
  /// Adds a nucleus to the list containing nucleus information
  virtual void addNucleus(unsigned int & qp);
  /// Nucleation rate density (should be a material property implementing nucleation theory)
  const MaterialProperty<Real> & _probability;

  /// Duration of time each nucleus is kept active after insertion
  Real _hold_time;

  /// the local nucleus list of nuclei centered in the domain of the current processor
  NucleusList & _local_nucleus_list;

  /// total nucleation rate
  Real _nucleation_rate;
  
  /// store the local nucleus radius
  const MaterialProperty<Real> & _local_radius;

  /// indicates whether time-dependent statistics are used or not
  const bool _time_dep_stats;
  
  /// track the number of ions during the simulation
  Real _N_ions;
  Real _N_ions_old;
};


//IonFluence
template <>
InputParameters validParams<IonFluencePostprocessor>();

/**
 * Creates a cumulative sum of a post-processor value over a transient.
 *
 * This is useful, for example, for counting the total number of linear or
 * nonlinear iterations during a transient.
 */
class IonFluencePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  IonFluencePostprocessor(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  /// cumulative sum of the post-processor value
  Real _fluence;

  /// UserObject that manages nucleus insertin and deletion
  const IontrackInserter & _inserter;

  /// cumulative sum of the post-processor value from the old time step */
  //const PostprocessorValue & _sum_old;

  /// current post-processor value to be added to the cumulative sum
  //const PostprocessorValue & _pps_value;
};

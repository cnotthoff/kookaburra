//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "BlockRestrictable.h"
#include "MaterialPropertyInterface.h"
#include "Coupleable.h"
#include "MooseVariableDependencyInterface.h"
#include "TransientInterface.h"
#include "RandomInterface.h"

#include "ElementUserObject.h"

#include "libmesh/mesh_tools.h"
/**
 * This UserObject manages the insertion and expiration of nuclei in the simulation
 * domain it manages a list of nuclei with their insertion times and their center
 * positions. A DiscreteNucleationMap is needed to enable the DiscreteNucleation
 * material to look up if a nucleus is present at a given element/qp.
 * This class can take a variable nucleus radius.
 */
class CNInserterBase : public GeneralUserObject,
                       public BlockRestrictable,
//public MaterialPropertyInterface,
//public Coupleable,
//public MooseVariableDependencyInterface,
//public TransientInterface,
                       public RandomInterface
{
public:
  static InputParameters validParams();

  CNInserterBase(const InputParameters & parameters);

  /// A nucleus has an expiration time, a location, and a size.
  // using NucleusLocation = std::tuple<Real, Point, Real>;
  struct NucleusLocation
  {
    Real time;
    Point center;
    Real radius;
  };

  /// Every MPI task should keep a full list of nuclei (in case they cross domains with their finite radii)
  using NucleusList = std::vector<NucleusLocation>;

  // counter pair to track insertions and deletion in the current timestep
  using NucleusChanges = std::pair<unsigned int, unsigned int>;

  virtual bool isMapUpdateRequired() const { return _update_required; }
  virtual const NucleusList & getNucleusList() const { return _global_nucleus_list; }
  virtual const NucleusChanges & getInsertionsAndDeletions() const { return _changes_made; }

  virtual const Real & getRate() const = 0;
  virtual const Point & getVect() const = 0;
protected:
  /// the global list of all tracks over all processors
  NucleusList & _global_nucleus_list;

  /// count the number of nucleus insertions and deletions
  NucleusChanges _changes_made;

  /// is a map update required
  bool _update_required;
};

// Used for Restart
template <>
void dataStore(std::ostream & stream,
               CNInserterBase::NucleusLocation & nl,
               void * context);
template <>
void dataLoad(std::istream & stream,
              CNInserterBase::NucleusLocation & nl,
              void * context);

class CNInserter : public CNInserterBase
{
public:
  static InputParameters validParams();

  CNInserter(const InputParameters & parameters);

  virtual void initialize(){
    //const Real oneoverC=6241509074460762607.776;
    _nucleation_rate=_current*_area*_scale*_scale*1e4;
    _changes_made = {0, 0};
    if(_is_prime){
      // expire entries from the local track list (if the current time step converged)
      if (_fe_problem.converged())
	{
	  unsigned int i = 0;
	  while (i < _global_nucleus_list.size())
	    {
	      if (_global_nucleus_list[i].time <= _fe_problem.time())
		{
		  // remove entry (by replacing with last element and shrinking size by one)
		  _global_nucleus_list[i] = _global_nucleus_list.back();
		  _global_nucleus_list.pop_back();
		  _changes_made.second++;
		}
	      else
		++i;
	    }
	}
    }
  }
  virtual void execute(){
    if(_is_prime){
      //Real nIdt = _nucleation_rate * _fe_problem.dt();
      //if(nIdt > 1){
      //	do{
      //	  addTrack();
      //	}while(nIdt-- >1);
      //      }
      const Real random = getRandomReal();
      // branch the operation for using time-dependent statistics or
      // time-independent probability (e.g., recrystallization fraction)
      // If time-dependent, `rate` refers to a probability rate density
      // If time-independent, `rate` refers to a probability density
      if (!_time_dep_stats)
	{
	  // if using time-independent statistics, this would be more like a nucleation fraction
	  if (random < _nucleation_rate)
	    addTrack();
	}
      else
	{
	  // We check the random number against the inverse of the zero probability.
	  // for performance reasons we do a quick check against the linearized form of
	  // that probability, which is always strictly larger than the actual probability.
	  // The expression below should short circuit and the expensive exponential
	  // should rarely get evaluated
	  if (random < _nucleation_rate * _fe_problem.dt() && random < (1.0 - std::exp(-_nucleation_rate * _fe_problem.dt())))
	    addTrack();
	}
    }
    //_communicator.broadcast(_hold_time);
    
  };
  
  //virtual void threadJoin(const UserObject & y){std::cerr<<"t tid="<<_tid<<"\n";};
  virtual void finalize(){
    /**
     * Pack the _global_nucleus_list into a simple vector of Real.
     * libMesh's allgather does not portably work on the original
     * _global_nucleus_list data structure!
     */
    std::vector<Real> comm_buffer(_global_nucleus_list.size() * 5);
    for (unsigned i = 0; i < _global_nucleus_list.size(); ++i)
      {
	comm_buffer[i * 5 + 0] = _global_nucleus_list[i].time;
	comm_buffer[i * 5 + 1] = _global_nucleus_list[i].center(0);
	comm_buffer[i * 5 + 2] = _global_nucleus_list[i].center(1);
	comm_buffer[i * 5 + 3] = _global_nucleus_list[i].center(2);
	comm_buffer[i * 5 + 4] = _global_nucleus_list[i].radius;
      }
    
    // combine _global_nucleus_lists from all MPI ranks
    _communicator.broadcast(comm_buffer);
    
    // unpack the gathered _global_nucleus_list
    unsigned int n = comm_buffer.size() / 5;
    mooseAssert(comm_buffer.size() % 5 == 0,
		"Communication buffer has an unexpected size (not divisible by 5)");
    if(!_is_prime){
      _global_nucleus_list.resize(n);
      for (unsigned i = 0; i < n; ++i)
	{
	  _global_nucleus_list[i].time = comm_buffer[i * 5 + 0];
	  _global_nucleus_list[i].center(0) = comm_buffer[i * 5 + 1];
	  _global_nucleus_list[i].center(1) = comm_buffer[i * 5 + 2];
	  _global_nucleus_list[i].center(2) = comm_buffer[i * 5 + 3];
	  _global_nucleus_list[i].radius = comm_buffer[i * 5 + 4];
	}
    }
    // get the global number of changes (i.e. changes to _global_nucleus_list)
    _communicator.broadcast(_changes_made.first);
    _communicator.broadcast(_changes_made.second);
    _communicator.broadcast(_N_ions);
    
    // gather the total nucleation rate
    //_communicator.broadcast(_nucleation_rate);
    
    _update_required = _changes_made.first > 0 || _changes_made.second > 0;
    //std::cerr<<"update("<<_communicator.rank()<<")="<<_update_required<<"\n";
    //    std::cerr<<"rad="<<_radius<<"\n";
  };
  
  const Real & getRate() const { return _nucleation_rate; };
  const Point & getVect() const { return _e1; };
protected:
  /// Adds a nucleus to the list containing nucleus information
  virtual void addTrack(){
      NucleusLocation new_nucleus;
      const Real ran1 = getRandomReal();
      const Real ran2 = getRandomReal();
      const Real ran3 = getRandomReal();
      
      new_nucleus.time = _fe_problem.time() + _hold_time;
      new_nucleus.center = Point(ran1*_lwd(0),ran2*_lwd(1),ran3*_lwd(2));
      new_nucleus.center += _box_min;
      new_nucleus.radius = _radius;
      
      _N_ions++;
      _global_nucleus_list.push_back(new_nucleus);
      _changes_made.first++;
  };

  /// Nucleation rate density (should be a material property implementing nucleation theory)
  //const MaterialProperty<Real> & _probability;

  /// Duration of time each nucleus is kept active after insertion
  Real _hold_time;

  /// the local nucleus list of nuclei centered in the domain of the current processor
  //NucleusList & _local_nucleus_list;

  /** Total nucleation rate.
   * For time-dependent statistics, this is probability rate density,
   * for time-independent statistics, it is probability density
   */
  Real _nucleation_rate;
  Real _area;
  Point _lwd;
  Point _box_min;
  Point _e1;

  int  _N_ions;
  /// store the local nucleus radius
  //  const MaterialProperty<Real> & _local_radius;

  /// indicates whether time-dependent statistics are used or not
  const bool _time_dep_stats;

  const bool _is_prime;
  const int _dir;
  const Real _current;
  const Real _scale;
  const Real _radius;
  const Real _angle;
};

class CNNucleationMap : public ElementUserObject
{
public:
  static InputParameters validParams();

  CNNucleationMap(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual void finalize() {};

  virtual void meshChanged();

  const std::vector<Real> & nuclei(const Elem *) const;

  const CNInserterBase & getInserter() const { return _inserter; }

  //Real getWidth() const { return _int_width; }
  //Real getPeriodic() const { return _periodic; }

protected:
  /// Did the mesh change since the last execution of this PP?
  bool _mesh_changed;

  /// Do we need to rebuild the map during this timestep?
  bool _rebuild_map;

  /// Buffer for building the per QP map
  std::vector<Real> _elem_map;

  /// Dummy map for elements without nuclei
  std::vector<Real> _zero_map;

  /// UserObject that manages nucleus insertin and deletion
  const CNInserterBase & _inserter;

  /// variable number to use for minPeriodicDistance calls (i.e. use the periodicity of this variable)
  int _periodic;

  /// Nucleus interface width
  const Real _int_width;

  /// list of nuclei maintained bu the inserter object
  const CNInserterBase::NucleusList & _nucleus_list;

  ///@{ Per element list with 0/1 flags indicating the presence of a nucleus
  using NucleusMap = std::unordered_map<dof_id_type, std::vector<Real>>;
  NucleusMap _nucleus_map;
  ///@}

  const Point _e1;
};

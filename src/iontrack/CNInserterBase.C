//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CNInserterBase.h"
#include "libmesh/parallel_algebra.h"

InputParameters
CNInserterBase::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += BlockRestrictable::validParams();
  //params += MaterialPropertyInterface::validParams();
  params += TransientInterface::validParams();
  params += RandomInterface::validParams();
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
  //params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
  return params;
}

CNInserterBase::CNInserterBase(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    BlockRestrictable(this),
    //TransientInterface(this),
    RandomInterface(parameters, _fe_problem, _tid, false),
    _global_nucleus_list(declareRestartableData("global_nucleus_list", NucleusList(0))),
    _changes_made(0, 0)//,
    //_update_required(_app.isRecovering() || _app.isRestarting())
{
  setRandomResetFrequency(EXEC_TIMESTEP_END);
}

template <>
void
dataStore(std::ostream & stream,
          CNInserterBase::NucleusLocation & nl,
          void * context)
{
  storeHelper(stream, nl.time, context);
  storeHelper(stream, nl.center, context);
  storeHelper(stream, nl.radius, context);
}

template <>
void
dataLoad(std::istream & stream,
         CNInserterBase::NucleusLocation & nl,
         void * context)
{
  loadHelper(stream, nl.time, context);
  loadHelper(stream, nl.center, context);
  loadHelper(stream, nl.radius, context);
}

registerMooseObject("kookaburraApp", CNInserter);
registerMooseObject("kookaburraApp", CNNucleationMap);

InputParameters
CNInserter::validParams()
{
  InputParameters params = CNInserterBase::validParams();
  params.addClassDescription("Manages the list of currently active nucleation sites and adds new "
                             "sites according to a given probability function.");
  //  params.addRequiredParam<MaterialPropertyName>(
  //      "probability", "Probability density for inserting a discrete nucleus");
  params.addRequiredParam<Real>("hold_time", "Time to keep each nucleus active");
  //  params.addParam<MaterialPropertyName>("radius",
  //					"r_crit",
  //					"variable radius material property name, supply a value if "
  //					"radius is constant in the simulation");
  params.addParam<bool>("time_dependent_statistics",
                        true,
                        "flag if time-dependent or time-independent statistics are used");
  params.addParam<int>("dir", 1, "Ion track direction (x=0, y=1, z=2)");
  params.addRequiredParam<Real>("flux", "Ion flux in ions per cm^2 s");
  params.addRequiredParam<Real>("scale", "Lenght scale used in the simulation");
  params.addRequiredParam<Real>("radius", "base radius of the track");
  params.addParam<Real>("angle", 0, "Ion incidence angle");

  return params;
}

CNInserter::CNInserter(const InputParameters & parameters)
  : CNInserterBase(parameters),
    //_probability(getMaterialProperty<Real>("probability")),
    //_hold_time(getParam<Real>("hold_time")),
    //_local_nucleus_list(declareRestartableData("local_nucleus_list", NucleusList(0))),
    //_local_radius(getMaterialProperty<Real>("radius")),
    _time_dep_stats(getParam<bool>("time_dependent_statistics")),
    _is_prime(_communicator.rank() == 0),
    _dir(getParam<int>("dir")),
    _current(getParam<Real>("flux")),
    _scale(getParam<Real>("scale")),
    _radius(getParam<Real>("radius")),
    _angle(getParam<Real>("angle")*M_PI/180.)
{
  MeshBase *mesh = &_fe_problem.mesh().getMesh();
  BoundingBox bbox = MeshTools::create_bounding_box(*mesh);
  const unsigned int dim = _fe_problem.mesh().dimension();
  _box_min = bbox.min();
  _lwd  = bbox.max()-bbox.min();
  int dir = 1<<_dir;
  int not_dir=0;
  
  if(_lwd(0) == 0) not_dir=1;
  else if(_lwd(1) == 0) not_dir=1<<1;
  else if(_lwd(2) == 0) not_dir=1<<2;
  
  if((1<<_dir) & not_dir){
    std::cerr<<"That dose not work!\n";
    exit(0);
  }
  not_dir = not_dir | dir; 
  dir = not_dir | dir;
  dir = ~dir;
  
  _area=( (int)(dir&0b001) * _lwd(0) + (int)((not_dir)&0b001)) * ( (int)((dir&0b010)>>1) * _lwd(1) + (int)( ( ( (not_dir)&0b010 ) )>>1 ) )  * ( (int)((dir&0b100)>>2) * _lwd(2) + (int)( ( ( (not_dir)&0b100 ) )>>2 ));

  Point tmp = Point((1<<_dir)&0b001,((1<<_dir)&0b010)>>1,((1<<_dir)&0b100)>>2);
  _e1(0) = cos(_angle)*tmp(0) - sin(_angle)*tmp(1);
  _e1(1) = cos(_angle)*tmp(1) + sin(_angle)*tmp(0);
  _e1(2) = tmp(2);
}

InputParameters
CNNucleationMap::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription("Generates a spatial smoothed map of all nucleation sites with the "
                             "data of the DiscreteNucleationInserter for use by the "
                             "DiscreteNucleation material.");
  params.addParam<Real>("int_width", 0.0, "Nucleus interface width for smooth nuclei");
  params.addRequiredParam<UserObjectName>("inserter", "CNInserter user object");
  params.addCoupledVar("periodic",
		       "Use the periodicity settings of this variable to populate the grain map");
  // the mapping needs to run at timestep begin, which is after the adaptivity
  // run of the previous timestep.
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;
  return params;
}

CNNucleationMap::CNNucleationMap(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _mesh_changed(false),
    _inserter(getUserObject<CNInserterBase>("inserter")),
    _periodic(isCoupled("periodic") ? coupled("periodic") : -1),
    _int_width(getParam<Real>("int_width")),
    _nucleus_list(_inserter.getNucleusList()),
    _e1(_inserter.getVect())
{
  _zero_map.assign(_fe_problem.getMaxQps(), 0.0);
}

void CNNucleationMap::initialize()
{
  if (_inserter.isMapUpdateRequired() || _mesh_changed)
  {
    _rebuild_map = true;
    _nucleus_map.clear();
  }
  else
    _rebuild_map = false;

  _mesh_changed = false;
}

void CNNucleationMap::execute()
{
  if (_rebuild_map)
  {
    // reserve space for each quadrature point in the element
    _elem_map.assign(_qrule->n_points(), 0);

    // store a random number for each quadrature point
    unsigned int active_nuclei = 0;
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    {
      Real r = std::numeric_limits<Real>::max();

      // find the distance to the closest nucleus
      Real local_radius = 0.0;
      for (unsigned i = 0; i < _nucleus_list.size(); ++i)
      {
	Point tmp1 = _nucleus_list[i].center - _q_point[qp];;
	Real t = -tmp1 * _e1;
	tmp1 = _nucleus_list[i].center + t*_e1;
	
        // use a non-periodic or periodic distance
	r = _periodic < 0
                ? (tmp1 - _q_point[qp]).norm()
                : _mesh.minPeriodicDistance(_periodic, tmp1, _q_point[qp]);
	//r = _periodic < 0
        //        ? (q_tmp - _nucleus_list[i].center).norm()
        //        : _mesh.minPeriodicDistance(_periodic, q_tmp, _nucleus_list[i].center);
        
        // grab the radius of the nucleus that this qp is closest to
        local_radius = _nucleus_list[i].radius;

        // compute intensity value with smooth interface
        Real value = 0.0;
        if (r <= local_radius - _int_width / 2.0) // Inside circle
        {
          active_nuclei++;
          value = 1.0;
        }
        else if (r < local_radius + _int_width / 2.0) // Smooth interface
        {
          Real int_pos = (r - local_radius + _int_width / 2.0) / _int_width;
          active_nuclei++;
          value = (1.0 + std::cos(int_pos * libMesh::pi)) / 2.0;
        }
        if (value > _elem_map[qp])
          _elem_map[qp] = value;
      }
    }

    // if the map is not empty insert it
    if (active_nuclei > 0)
      _nucleus_map.insert(
          std::pair<dof_id_type, std::vector<Real>>(_current_elem->id(), _elem_map));
  }
}

void
CNNucleationMap::threadJoin(const UserObject & y)
{
  // if the map needs to be updated we merge the maps from all threads
  if (_rebuild_map)
  {
    const CNNucleationMap & uo = static_cast<const CNNucleationMap &>(y);
    _nucleus_map.insert(uo._nucleus_map.begin(), uo._nucleus_map.end());
  }
}

void
CNNucleationMap::meshChanged()
{
  _mesh_changed = true;
}

const std::vector<Real> &
CNNucleationMap::nuclei(const Elem * elem) const
{
  NucleusMap::const_iterator i = _nucleus_map.find(elem->id());

  // if no entry in the map was found the element contains no nucleus
  if (i == _nucleus_map.end())
    return _zero_map;

  return i->second;
}

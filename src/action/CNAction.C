#include "CNAction.h"
// MOOSE includes
#include "Conversion.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseObjectAction.h"
#include "MooseMesh.h"
#include "AddVariableAction.h"

#include "libmesh/string_to_enum.h"

registerMooseAction("kookaburraApp", CNAction, "add_variable");

registerMooseAction("kookaburraApp", CNAction, "add_kernel");

InputParameters
CNAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Set up the variable(s) and the kernels needed for a conserved phase field variable."
      " Note that for a direct solve, the element family and order are overwritten with hermite "
      "and third.");
  MooseEnum solves("DIRECT REVERSE_SPLIT FORWARD_SPLIT");
  params.addRequiredParam<MooseEnum>("solve_type", solves, "Split or direct solve?");
  // Get MooseEnums for the possible order/family options for this variable
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  params.addParam<MooseEnum>("family",
                             families,
                             "Specifies the family of FE "
                             "shape functions to use for this variable");
  params.addParam<MooseEnum>("order",
                             orders,
                             "Specifies the order of the FE "
                             "shape function to use for this variable");
  //params.addRequiredParam<std::vector<NonlinearVariableName>>("variables",
  //"The names of the convection and diffusion variables in the simulation");
  params.addParam<Real>("scaling", 1.0, "Specifies a scaling factor to apply to this variable");
  params.addParam<bool>("implicit", true, "Whether kernels are implicit or not");
  params.addParam<bool>("use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  params.addParamNamesToGroup("scaling implicit use_displaced_mesh", "Advanced");
  params.addRequiredParam<MaterialPropertyName>("mobilityc", "The mobility used with the kernel");
  params.addParam<MaterialPropertyName>("mobilityl", "L", "The mobility used with the kernel");
  params.addParam<std::vector<VariableName>>("argsx",
                                             "Vector of variable arguments this kernel depends on");
  params.addRequiredParam<MaterialPropertyName>("free_energy", "Base name of the free energy function F defined in a free energy material");
  params.addRequiredParam<MaterialPropertyName>("kappac", "The kappa used with the kernel");
  params.addParam<MaterialPropertyName>("kappal", "kappa_op", "The kappa used with the kernel");
  params.addParam<bool>("variable_mobility",
  			true,
  			"The mobility is a function of any MOOSE variable (if "
  			"this is set to false, L must be constant over the "
  			"entire domain!)");

  return params;
}

CNAction::CNAction(const InputParameters & params)
  : Action(params),
    _solve_type(getParam<MooseEnum>("solve_type").getEnum<SolveType>()),
    _var_name(name()),
    _scaling(getParam<Real>("scaling")),
    _argsx(getParam<std::vector<VariableName>>("argsx"))
{
  switch (_solve_type)
  {
    case SolveType::DIRECT:
      _fe_type = FEType(Utility::string_to_enum<Order>("THIRD"),
                        Utility::string_to_enum<FEFamily>("HERMITE"));
      if (!parameters().isParamSetByAddParam("order") &&
          !parameters().isParamSetByAddParam("family"))
        mooseWarning("Order and family autoset to third and hermite in ConservedAction");
      break;
    case SolveType::REVERSE_SPLIT:
    case SolveType::FORWARD_SPLIT:
      _fe_type = FEType(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
                        Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family")));
      // Set name of chemical potential variable
      _chempot_name = "chem_pot_" + _var_name;
      break;
    default:
      paramError("solve_type", "Incorrect solve_type in ConservedAction");
  }
}

void
CNAction::act()
{
  //
  // Add variable(s)
  //
  if (_current_task == "add_variable")
  {
    auto type = AddVariableAction::determineType(_fe_type, 1);
    auto var_params = _factory.getValidParams(type);
    var_params.set<MooseEnum>("family") = Moose::stringify(_fe_type.family);
    var_params.set<MooseEnum>("order") = _fe_type.order.get_order();
    var_params.set<std::vector<Real>>("scaling") = {_scaling};

    // Create conserved variable _var_name
    _problem->addVariable(type, _var_name, var_params);

    // Create chemical potential variable for split form
    switch (_solve_type)
    {
      case SolveType::DIRECT:
        break;
      case SolveType::REVERSE_SPLIT:
      case SolveType::FORWARD_SPLIT:
        _problem->addVariable(type, _chempot_name, var_params);
    }
  }

  //
  // Add Kernels
  //
  else if (_current_task == "add_kernel")
  {
    {
      // Add time derivative kernel
      std::string kernel_type_ac = "TimeDerivative";
      
      std::string kernel_name_ac = _var_name + "_AC" + kernel_type_ac;
      InputParameters params1 = _factory.getValidParams(kernel_type_ac);
      params1.set<NonlinearVariableName>("variable") = "eta";
      params1.applyParameters(parameters());
      
      _problem->addKernel(kernel_type_ac, kernel_name_ac, params1);
      
      // Add AllenCahn kernel
      kernel_type_ac = "AllenCahn";
      
      kernel_name_ac = _var_name + "_AC" + kernel_type_ac;
      InputParameters params2 = _factory.getValidParams(kernel_type_ac);
      params2.set<NonlinearVariableName>("variable") = "eta";
      params2.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mobilityl");
      params2.set<MaterialPropertyName>("f_name") = getParam<MaterialPropertyName>("free_energy");
      params2.set<std::vector<VariableName>>("args")=_argsx;
      params2.applyParameters(parameters());
      
      _problem->addKernel(kernel_type_ac, kernel_name_ac, params2);

      // Add ACInterface kernel
      kernel_type_ac = "ACInterface";
      
      kernel_name_ac = _var_name + "_AC" + kernel_type_ac;
      InputParameters params3 = _factory.getValidParams(kernel_type_ac);
      params3.set<NonlinearVariableName>("variable") = "eta";
      params3.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mobilityl");
      params3.set<MaterialPropertyName>("kappa_name") = getParam<MaterialPropertyName>("kappal");
      params3.set<bool>("variable_L") = getParam<bool>("variable_mobility");
      params3.applyParameters(parameters());
      
      _problem->addKernel(kernel_type_ac, kernel_name_ac, params3);
    }
    
    
    switch (_solve_type)
    {
      case SolveType::DIRECT:
        // Add time derivative kernel
        {
          std::string kernel_type = "TimeDerivative";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _var_name;
	  params.set<std::vector<VariableName> >("args") = _argsx;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }

        // Add CahnHilliard kernel
        {
          std::string kernel_type = "CahnHilliard";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _var_name;
          params.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mobilityc");
          params.set<MaterialPropertyName>("f_name") =
              getParam<MaterialPropertyName>("free_energy");
	  params.set<std::vector<VariableName> >("args") = _argsx;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }

        // Add CHInterface kernel
        {
          std::string kernel_type = "CHInterface";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _var_name;
          params.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mobilityc");
          params.set<MaterialPropertyName>("kappa_name") = getParam<MaterialPropertyName>("kappac");
	  params.set<std::vector<VariableName> >("args") = _argsx;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }
        break;

      case SolveType::REVERSE_SPLIT:
        // Add time derivative kernel
        {
          std::string kernel_type = "CoupledTimeDerivative";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _chempot_name;
          params.set<std::vector<VariableName>>("v") = {_var_name};
	  //params.set<std::vector<VariableName> >("args") = _argsx;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }

        // Add SplitCHWRes kernel
        {
          std::string kernel_type = "SplitCHWRes";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _chempot_name;
          params.set<MaterialPropertyName>("mob_name") = getParam<MaterialPropertyName>("mobilityc");
	  params.set<std::vector<VariableName> >("args") = _argsx;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }

        // Add SplitCHParsed kernel
        {
          std::string kernel_type = "SplitCHParsed";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _var_name;
          params.set<std::vector<VariableName>>("w") = {_chempot_name};
          params.set<MaterialPropertyName>("f_name") =
	    getParam<MaterialPropertyName>("free_energy");
          params.set<MaterialPropertyName>("kappa_name") = getParam<MaterialPropertyName>("kappac");
	  params.set<std::vector<VariableName> >("args") = _argsx;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }
        break;

      case SolveType::FORWARD_SPLIT:
        // Add time derivative kernel
        {
          std::string kernel_type = "TimeDerivative";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _var_name;
	  params.set<std::vector<VariableName> >("args") = _argsx;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }

        // Add MatDiffusion kernel for c residual
        {
          std::string kernel_type = "MatDiffusion";

          std::string kernel_name = _var_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _var_name;
          params.set<std::vector<VariableName>>("v") = {_chempot_name};
          params.set<MaterialPropertyName>("diffusivity") =
              getParam<MaterialPropertyName>("mobilityc");
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }
        // Add MatDiffusion kernel for chemical potential residual
        {
          std::string kernel_type = "MatDiffusion";

          std::string kernel_name = _chempot_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _chempot_name;
          params.set<std::vector<VariableName>>("v") = {_var_name};
          params.set<MaterialPropertyName>("diffusivity") = getParam<MaterialPropertyName>("kappac");
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }

        // Add CoupledMaterialDerivative kernel
        {
          std::string kernel_type = "CoupledMaterialDerivative";

          std::string kernel_name = _chempot_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _chempot_name;
          params.set<std::vector<VariableName>>("v") = {_var_name};
          params.set<MaterialPropertyName>("f_name") =
              getParam<MaterialPropertyName>("free_energy");
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }

        // Add CoefReaction kernel
        {
          std::string kernel_type = "CoefReaction";

          std::string kernel_name = _chempot_name + "_" + kernel_type;
          InputParameters params = _factory.getValidParams(kernel_type);
          params.set<NonlinearVariableName>("variable") = _chempot_name;
          params.set<Real>("coefficient") = -1.0;
          params.applyParameters(parameters());

          _problem->addKernel(kernel_type, kernel_name, params);
        }
    }
  }
}

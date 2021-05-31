#include "kookaburraApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
kookaburraApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

kookaburraApp::kookaburraApp(InputParameters parameters) : MooseApp(parameters)
{
  kookaburraApp::registerAll(_factory, _action_factory, _syntax);
}

kookaburraApp::~kookaburraApp() {}

void
kookaburraApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"kookaburraApp"});
  Registry::registerActionsTo(af, {"kookaburraApp"});

  /* register custom execute flags, action syntax, etc. here */
  
  registerSyntax("CNAction", "Modules/CN/CNAC/*");
  
}

void
kookaburraApp::registerApps()
{
  registerApp(kookaburraApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
kookaburraApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  kookaburraApp::registerAll(f, af, s);
}
extern "C" void
kookaburraApp__registerApps()
{
  kookaburraApp::registerApps();
}

#
# s41313-017-0008-y.pdf
# DOI 10.1186/s41313-017-0008-y
# cR = ceq*exp(2*gamma*Omega/(R*k*T))

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmax = 120
  ymax = 120
  #elem_type = QUAD4
[]

[GlobalParams]
#  #polynomial_order = 4
  block = 0
  enable_jit = true
[]
[Variables]
  [./eta]
    order = FIRST
    family = LAGRANGE
  #  #initial_condition = 0.4
  [../]
  #[./c]
  #  order = FIRST
  #  family = LAGRANGE
  #  #initial_condition = 0.4
  #[../]
[]
[Modules]
  [./CN]
    [./CNAC]
      [./c] 
        free_energy = F
        kappac = 0.0
        mobilityc = 0.024
        kappal =    11.622
        mobilityl = 0.0012
        argsx = 'eta c'
        solve_type = REVERSE_SPLIT
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./c_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./c_aux]
    type = IontrackNucleationAux
    map = map
    variable = c_aux
    no_nucleus_value = -1
    nucleus_value = 2
    execute_on = TIMESTEP_END
  [../]
#  [./eta_norm_k]
#    type = NormalizationAux
#    variable = eta_norm
#    source_variable = eta
#    normalization = vol
#    # this coefficient will affect the eigenvalue.
#    normal_factor = 1.0
#    execute_on = 'initial timestep_end'
#  [../]
[]

  
[ICs]
  [./cIC]
    type = RandomIC
    variable = c
    min = 0.00
    #max = 0.00021
    max = 4e-6
  [../]
  [./eatIC]
    type = RandomIC
    variable = eta
    min = 0.0
    max = 0.0001
  [../] 
#  [./cIC]
#    type = SmoothCircleIC
#    x1 = 30.0
#    y1 = 30.0
#    radius = 15.0
#    invalue = 1.0
#    outvalue = 0.019
#    int_width = 1.0
#    variable = c
#  [../]
#  [./eta_IC]
#    type = SmoothCircleIC
#    x1 = 30.0
#    y1 = 30.0
#    radius = 15.0
#    invalue = 1.0
#    outvalue = 0.0
#    int_width = 1.0
#    variable = eta
#  [../]
#  [./eta_IC]
#     type = MultiSmoothCircleIC
#     radius = 1
#     invalue = 1
#     outvalue = 0
#     int_width=0.1
#     numbub = 5
#     bubspac = 3
     #Lx = 60
     #Ly = 60
     #Lz = 0
#     rand_seed = 2000
#     radius_variation=0
#     variable = eta
#  [../]
#  [./c_IC]
#     type = MultiSmoothCircleIC
#     radius = 1
#     invalue = 1
#     outvalue = 0.1
#     int_width=0.1
#     numbub = 5
#     bubspac = 3
     #Lx = 60
     #Ly = 60
     #Lz = 0
#     rand_seed = 2000
#     radius_variation=0
#     variable = c
#  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  # dummy energy
  [./Fdummy]
      # simple double well free energy
    type = DerivativeParsedMaterial
    f_name = F
    args = 'c eta'
    constant_names       = 'A   B'
    constant_expressions = '46.5 465'
    function = 'A*eta^2*(1-eta)^2 + B*(c-(3*eta^2-2*eta^3))^2'
    derivative_order = 2
    outputs = exodus
  [../]  
#  [./precipitate_indicator]  # Returns 1/vol if precipitate
#    type = ParsedMaterial
#    f_name = prec_indic
#    args = 'eta_norm eta'
#    #material_property_names = 'vol'
#    function = if(eta>0.8,eta_norm,0)
#  [../]
  # forcing function mask
  [./mask]
    type = ParsedMaterial
    f_name = irrad_mask
    function = 'if(eta > 0.8, 0, 1)'
    #material_property_names = 'vol'
    args = 'eta'
  [../]
  [./const_pfm]
    type = GenericConstantMaterial
    prop_names = 'radius'
    prop_values = ' 2'
  [../]
[]

[Kernels]
  #active = 'PhaseConc ChemPotSolute CHBulk ACBulkF ACBulkC ACInterface dcdt detadt ckernel'
  [./c_source]
    type = IontrackNucleationForce
    variable = chem_pot_c
    map = map
    mask = irrad_mask
    no_nucleus_value = 0
    nucleus_value = 0.0005
  [../]
#  [./ffn]
#    type = MaskedBodyForce # BodyForce
#    variable = chem_pot_c
#    function = forcing_fn
#    mask = irrad_mask
#  [../]
  [./conserved_langevin]
    type = ConservedLangevinNoise
    amplitude = 1e-5
    variable = eta
    noise = normal_masked_noise
  []
  #  [./noisesource1]
#    type = SourceNoiseKernel
#    amplitude = 1e-2
#    variable = chem_pot_c
#    noise = source_masked_noise
#  [../]
#  [./noisesource2]
#    type = SourceNoiseKernel
#    amplitude = 1e-2
#    variable = eta
#    noise = source_masked_noise
#  [../]

[]

[UserObjects]
  [./inserter]
    # The inserter runs at the end of each time step to add nucleation events
    # that happend during the timestep (if it converged) to the list of nuclei
    type = IontrackInserter
    hold_time = 0.1
    probability = 0.0001
  [../]
  [./map]
    # The map UO runs at the beginning of a timestep and generates a per-element/qp
    # map of nucleus locations. The map is only regenerated if the mesh changed or
    # the list of nuclei was modified.
    # The map converts the nucleation points into finite area objects with a given radius.
    type = IontrackNucleationMap
    radius = 2
    int_width = 2
    angle = 0
    periodic = c
    inserter = inserter
  [../]

  [./normal_masked_noise]
    type = ConservedMaskedNormalNoise
    mask = irrad_mask
  [../]
#  [./source_masked_noise]
#    type = SourceMaskedNoise
#    Pcasc = 0.001
#    mask = irrad_mask
#  [../]
[]

[Functions]
#  [./forcing_fn]
#    type = ParsedFunction
#    # dudt = 3*t^2*(x^2 + y^2)
#    # value = 3*t*t*((x*x)+(y*y))-(4*t*t*t)
#    value = 'if(t<0,10.0e-4,0)'
#  [../]
[]

[Postprocessors]
#  [./void_fraction]      # Fraction of surface devoted to precipitates
#      type = ElementIntegralMaterialProperty
#      mat_prop = prec_indic
#  [../]
#  [./vol]
#      type = VolumePostprocessor
#      outputs = 'none'
#      execute_on = 'initial timestep_end'
#  [../]
  [./fluence]
    type = IonFluencePostprocessor
    inserter = inserter
  []
  [./area_eta_v]
    type = LevelSetVolume
    threshold = 0.5
    variable = eta
    location = outside
    execute_on = 'initial timestep_end'
  [../]
  [./area_eta_m]
    type = LevelSetVolume
    threshold = 0.5
    variable = eta
    location = inside
    execute_on = 'initial timestep_end'
  [../]

[]

[Preconditioning]
  [./cw_coupling]
    type = SMP
    full = true
  [../]
[]


[Adaptivity]
  marker = errorfrac # this specifies which marker from 'Markers' subsection to use
  #steps = 2 # run adaptivity 2 times, recomputing solution, indicators, and markers each time

  # Use an indicator to compute an error-estimate for each element:
  [./Indicators]
    # create an indicator computing an error metric for the convected variable
    [./error] # arbitrary, use-chosen name
      type = GradientJumpIndicator
      variable = c
      outputs = none
    [../]
  [../]
  initial_steps = 4
  initial_marker = imarker
  # Create a marker that determines which elements to refine/coarsen based on error estimates
  # from an indicator:
  [./Markers]
    [./imarker]
      type = UniformMarker
      mark = REFINE
      outputs = none
    [../]
    [./errorfrac] # arbitrary, use-chosen name (must match 'marker=...' name above
      type = ErrorFractionMarker
      indicator = error # use the 'error' indicator specified above
      refine = 0.5 # split/refine elements in the upper half of the indicator error range
      coarsen = 0.1 # dont do any coarsening
      #outputs = none
    [../]
  [../]
  cycles_per_step = 3
  recompute_markers_during_cycles = true
  max_h_level = 5
  start_time = 500
[]
    
[Executioner]

  type = Transient
  solve_type = 'PJFNK'
  scheme = bdf2
  
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'

  l_max_its = 50
  nl_max_its = 50

  end_time = 2000.0

  start_time=0
  dtmax = 5
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 8
    iteration_window = 4
  [../]

[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = false
[]

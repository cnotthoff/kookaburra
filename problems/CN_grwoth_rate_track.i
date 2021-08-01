#
# s41313-017-0008-y.pdf
# DOI 10.1186/s41313-017-0008-y
# cR = ceq*exp(2*gamma*Omega/(R*k*T))

#description of the energy model:
#https://mooseframework.inl.gov/modules/phase_field/Quantitative.html

#energies:
#https://repository.kaust.edu.sa/bitstream/handle/10754/315753/Antisites+and+anisotropic+diffusion+in+GaAs+and+GaSb.pdf?sequence=1

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 200
  ny = 100
  xmax = 200
  ymax = 100
  #elem_type = QUAD4
[]

[GlobalParams]
  block = 0
  polynomial_order = 4
  enable_jit = true
[]
  
[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]
[]
  
[ICs]
  #inactive = 'cIC'
  [./cIC]
    type = MultiSmoothCircleIC
    variable = c
    invalue = 1.0
    outvalue = 0.02
    bubspac = 20.0 # This spacing is from bubble center to bubble center
    numbub = 30 # for Ayu: please play with this parameter
    radius = 5
    int_width = 6.0
    rand_seed = 2000
    radius_variation = 0.01 #This is the standard deviation
    radius_variation_type = normal
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./GaSb]
    type = PFParamsPolyFreeEnergy
    c = c
    T = 400 # K
    int_width = 3.0
    length_scale = 1.0e-9
    time_scale = 1.0
    #D0 = 1.0e-13 # m^2/s,#https://www-nature-com.virtual.anu.edu.au/articles/35040526.pdf
    D0 = 1.0e-2 # m^2/s, # too large but looks like we need this
    Em = 1.2 # in eV, #https://aip-scitation-org.virtual.anu.edu.au/doi/pdf/10.1063/1.3010300 
    Ef = 2.0 # in eV, 
    surface_energy = 0.7 # Total guess
  [../]
  [./free_energy]
    type = PolynomialFreeEnergy
    c = c
    derivative_order = 2
  [../]
  
  [./mask]
    type = ParsedMaterial
    f_name = irrad_mask
    function = 'if(c > 0.8, 0, 1)'
    args = 'c'
  [../]
  [./const_pfm]
    type = GenericConstantMaterial
    prop_names = 'radius'
    prop_values = ' 10'
  [../]
[]

[Kernels]
  [./c_source]
    type = IontrackNucleationForce
    variable = w
    map = map
    mask = irrad_mask
    no_nucleus_value = 0
    nucleus_value = 0.2 # for Ayu: please play with this parameter to see if you can reproduce the experimental grwoth rate of about 0.02 per ion
# plot area_c_v/(area_c_v+area_c_m) as a function of ion impacts and compare this to rho = (1-exp(-0.02*N_ion + off1) + off2
# where off1 and off2 are arbritray fit parameters to adjust for best fit.
  [../]
  [./c_res]
    type = SplitCHParsed
    variable = c
    kappa_name = kappa
    w = w
    f_name = F
  [../]
  [./w_res]
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
  [./time]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
[]

[UserObjects]
  [./inserter]
    # The inserter runs at the end of each time step to add nucleation events
    # that happend during the timestep (if it converged) to the list of nuclei
    type = IontrackInserter
    hold_time = 0.01
    probability = 0.0000015
  [../]
  [./map]
    # The map UO runs at the beginning of a timestep and generates a per-element/qp
    # map of nucleus locations. The map is only regenerated if the mesh changed or
    # the list of nuclei was modified.
    # The map converts the nucleation points into finite area objects with a given radius.
    type = IontrackNucleationMap
    radius = 10
    int_width = 5
    angle = 0
    periodic = c
    inserter = inserter
  [../]
[]

[Postprocessors]
  [./fluence]
    type = IonFluencePostprocessor
    inserter = inserter
  []
  [./area_c_v]
    type = LevelSetVolume
    threshold = 0.5
    variable = c
  location = outside
    execute_on = 'initial timestep_end'
  [../]
  [./area_c_m]
    type = LevelSetVolume
    threshold = 0.5
    variable = c
    location = inside
    execute_on = 'initial timestep_end'
  [../]
  [./dtnuc]
    type = IontrackNucleationTimeStep
    inserter = inserter
    p2nucleus = 0.02
    dt_max = 0.001
    outputs = none
  [../]
[]

[Preconditioning]
  [./cw_coupling]
    type = SMP
    full = true
  [../]
[]
    
[Executioner]

  type = Transient
  solve_type = 'PJFNK'
  scheme = bdf2
  
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_factor_mat_solver_package'
  petsc_options_value = 'asm      lu          nonzero mumps'

  l_max_its = 10
  nl_max_its = 10

  nl_rel_tol  = 1e-10
  nl_abs_tol  = 1e-10
  
  end_time = 2000.0

  start_time=0
  dtmax = 5
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 10.
    timestep_limiting_postprocessor = dtnuc
  [../]
[]

[Outputs]
  checkpoint = true

  #exodus = false
  exodus = true
  
  interval = 1
  file_base = track_growth_rate_1_out #for Ayu: here you can define the output file name
  [console]
    type = Console
    interval = 1
  []
  [./csv]
    type = CSV
    append_restart=true
    interval = 1
  [../]
  perf_graph = false
[]

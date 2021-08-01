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
    outvalue = 0.032
    bubspac = 20.0 # This spacing is from bubble center to bubble center
    numbub = 30 # for Ayu: please play with this parameter
    radius = 5
    int_width = 8.0
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
    int_width = 4.0 # I had to increase the width to stabelise the simulation, I guess we have to change to c, eta at some point.
    length_scale = 1.0e-9
    time_scale = 1.0
    #D0 = 1.0e-13 # m^2/s,#https://www-nature-com.virtual.anu.edu.au/articles/35040526.pdf
    D0 = 1.0e-2 # m^2/s, # for Ayu: if you have some time may be play with this parameter too.
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
    type = MaskedBodyForce
    variable = w
    value = 0.0001
    mask = irrad_mask
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

[Postprocessors]
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
    dt = 0.1
    #growth_factor = 1.5
    #cutback_factor = 0.5
    #optimal_iterations = 8
    #iteration_window = 4
  [../]
[]

[Outputs]
  checkpoint = true

  #exodus = false
  exodus = true
  
  interval = 10
  file_base = cont_growth_rate_1_out #for Ayu: here you can define the output file name
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

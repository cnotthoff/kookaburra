[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 100
  nz = 0
  xmin = -50
  xmax = 50
  ymin = -200
  ymax = 0
  zmin = 0
  zmax = 0
  elem_type = QUAD8
[]

[GlobalParams]
  block = 0
  enable_jit = true
[]
[Variables]
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]
  [./c]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Functions]
  [./forcing_fn]
    type = ParsedFunction
    value = 'if(y<-10,0,if(y>10,0,1))'
  [../]
[]
  
[ICs]
[]

[BCs]
  #active = 'Periodic '
  [./Periodic]
    [./all]
      auto_direction = 'x'
    [../]
  [../]
  [./dlc]
     type = NeumannBC
     variable = c
     boundary = 'bottom'
     value = 0.0
  [../]
  [./drc]
     type = DirichletBC
     variable = c
     boundary = 'top'
     value = 0.0
  [../]
  [./dle]
     type = NeumannBC
     variable = eta
     boundary = 'bottom'
     value = 0
  [../]
  [./dre]
     type = DirichletBC
     variable = eta
     boundary = 'top'
     value = 0.0
  [../]
[]

[Materials]
  # free energy
  [./Fdummy]
    type = CNFreeEnergy # (Ef*c + kbT( c*ln(c)+ (1-c)*ln(1-c))*g(eta) + h(eta)*(c-1)^2
    Ef = 0.7 #interesting to vary, maybe from 1 eV down 
    T = 600  #might be also interesting to vary
    c = c
    eta = eta
    log_tol = 1e-2
    f_name = F
  [../]  
  # forcing function mask
  [./mask]
    type = ParsedMaterial
    f_name = irrad_mask
    function = 'if(eta > 0.8, 0, 1)'
    args = 'eta'
  [../]

  # sink mask to model free surface at the top
  [./s_mask]
    type = ParsedMaterial
    f_name = s_mask
    function = 'c'
    args = 'c'
  [../]

  [./tmp_pfm]
    type = GenericConstantMaterial
    # interesting to vary D
    prop_names  = '   D        N     kappa_c    kappa_eta    l_scale t_scale'
    prop_values = ' 5.5e-18  1.6e28   0.1        0.1          1e-9      1'
			     
    #prop_names  = ' Em     D0        N     kappa_c    kappa_eta  T    l_scale t_scale'
    #prop_values = ' 0.5  1.10e-11  1.6e28    1             1    400      1e-9      1'
    #prop_values = ' 1.0  1.10e-11  1.6e28  0.05e2     0.1e2    600      1e-9      1'
  []

  [./mob_c]
    type = DerivativeParsedMaterial
    f_name = M
    args = 'c eta'
    #constant_names       = 'kb'
    #constant_expressions = '8.6173324e-5'
    #material_property_names = 'Em T D0 l_scale t_scale d2F:=D[F(c,eta),c,c,]'
    #function = 'D0*exp(-Em/(kb*T))/(l_scale*l_scale/t_scale)/d2F'
    material_property_names = 'D l_scale t_scale d2F:=D[F(c,eta),c,c,]'
    function = 'D/(l_scale*l_scale/t_scale)/d2F'
    derivative_order = 1
    #outputs = exodus
  [../]  
  
  [./mob_eta]
    type = DerivativeParsedMaterial
    f_name = L
    args = 'c eta'
    material_property_names = 'M l_scale t_scale N'
    function = 'M*(l_scale*l_scale*l_scale)*N*t_scale'
    derivative_order = 1
    #outputs = exodus
  [../]  
[]

[Kernels]
  #sink term to model free surface
  [./sink]
    type = MaskedBodyForce # BodyForce
    variable = w
    function = forcing_fn
    mask = s_mask
    value = -1e3 #might need adjustment to keep concentration zero at the boundary
  [../]

  [./source]
    type = MaskedBodyForce # BodyForce
    variable = w
    mask = irrad_mask
    value = 0.00025 #interesting
  [../]
  
  [./conserved_langevin_eta]
    type = ConservedLangevinNoise
    amplitude = 1e-4
    variable = eta
    noise = normal_masked_noise
  []
  # one run with conserved_langevin_c and one without
  #inactive = 'conserved_langevin_c'
  [./conserved_langevin_c]
    type = ConservedLangevinNoise
    amplitude = 1e-4
    variable = c
    noise = normal_masked_noise
  []
			     
  #Cahn Hilliard
  [./c_res]
    type = SplitCHParsed
    variable = c
    kappa_name = kappa_c
    w = w
    f_name = F
    args = eta
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
  
  #Allen Cahn
  [./detadt]
    type = TimeDerivative
    variable = eta
  [../]
  [./ACBulk]
    type = AllenCahn
    variable = eta
    f_name = F
    mob_name = L
    args = c
  [../]
  [./ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = kappa_eta
    variable_L = true
    mob_name = L
  [../]
[]
			     
[UserObjects]
  [./normal_masked_noise]
    type = ConservedMaskedNormalNoise
    mask = irrad_mask
  [../]
[]

[Postprocessors]
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
  [./flood_count_pp]
    type = FeatureFloodCount
    variable = eta
    threshold = 0.5
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

  l_max_its = 20
  nl_max_its = 20

  end_time = 200000.0

  start_time=0
  dtmax = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.2
    cutback_factor = 0.5
    optimal_iterations = 10
    iteration_window = 3
  [../]

[]

[Outputs]
  checkpoint = true

  #exodus = false
  exodus = true
  
  interval = 10
  #file_base = nucleation_test_1_out #for Ayu: here you can define the output file name
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

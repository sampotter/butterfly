examples = [
  'cheb_cov',
  'lbo_cov',
  'save_FEM_matrices'
]

example_deps = [argtable3_dep, butterfly_dep, m_dep]

foreach example : examples
  executable(
    example,
    [example + '.c'],
    dependencies : example_deps
  )
endforeach

scripts = [
  'lbo_cov_test.sh',
  'run_bf_lbo_test.sh',
  'run_lbo_cov.sh',
  'plot_test_output.jl',
  'covariance_sphere.jl',
  'sphere_cov.jl',
  'Project.toml'
]

foreach script : scripts
  custom_target(
    script,
    input : script, output : script,
    command : ['cp', '@INPUT@', '@OUTPUT@'],
    install : false,
    build_by_default : true
  )
endforeach
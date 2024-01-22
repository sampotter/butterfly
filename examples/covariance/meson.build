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
  # 'lbo_cov_test.sh',
  'run_lbo_cov.sh',
  'exact_cov_test.sh',
  'fast_cov_test.sh',
  'plot_sample.jl',
  'plot_exact_cov_test.jl',
  'plot_fast_cov_test.jl',
  'utils.jl',
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
examples = [
  'test_hankel_evaluator',
]

example_deps = [argtable3_dep, butterfly_dep, m_dep]

foreach example : examples
  executable(
    example,
    example + '.c',
    dependencies : example_deps
  )
endforeach

scripts = [
  'make_hankel_evaluator_plots.py',
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

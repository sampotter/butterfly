examples = [
  'bf_lbo'
]

example_deps = [argtable3_dep, butterfly_dep, m_dep]

foreach example : examples
  executable(example, example + '.c', dependencies : example_deps)
endforeach

scripts = [
  'bf_lbo_plot_blocks.py',
  'collect_lbo_data.py'
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

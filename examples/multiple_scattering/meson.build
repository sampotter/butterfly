example_deps = [argtable3_dep, butterfly_dep, m_dep]

executable('multiple_scattering', ['multiple_scattering.c',
  'multiple_scattering_opts.c', 'multiple_scattering_context.c'],
  dependencies : example_deps + [fmm2d_dep])

scripts = [
  'collect_multiple_scattering_data.py',
  'convert_raw_multiple_scattering_data_to_csv.py',
  'make_multiple_scattering_plots.py',
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

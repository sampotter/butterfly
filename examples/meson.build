examples = [
  'bf_one_block',
  'bf_all_blocks',
  'bf_lbo',
  'helm2_bie',
  'test_hankel_evaluator',
  'test_helm2_adapt_fac'
]

foreach example : examples
  example_deps = [argtable3_dep, butterfly_dep, m_dep]
  executable(
    example,
    example + '.c',
    dependencies : example_deps
  )
endforeach

executable('multiple_scattering', ['multiple_scattering.c',
  'multiple_scattering_opts.c', 'multiple_scattering_context.c'],
  dependencies : example_deps + [fmm2d_dep])

scripts = [
  'bf_lbo_plot_blocks.py',
  'collect_multiple_scattering_data.py',
  'convert_raw_multiple_scattering_data_to_csv.py',
  'ellipse.py',
  'make_circle_test_data.py',
  'make_ellipse_test_data.py',
  'make_hankel_evaluator_plots.py',
  'make_multiple_scattering_plots.py',
  'plot_bf_all_blocks.py',
  'plot_bf_one_block.py',
  'plot_blocks.py',
  'plot_mesh_and_octree.py',
  'plot_points_and_quadtree.py',
  'run_bf_all_blocks_test.sh',
  'run_bf_lbo_test.sh',
  'run_bf_one_block_test.sh',
  'run_helm2_bie_test.sh',
  'scratch.py',
  'test.py',
  'util.py',
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

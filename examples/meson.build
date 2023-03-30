examples = [
  'bf_one_block',
  'bf_all_blocks',
  'bf_lbo',
  'helm2_bie',
  'multiple_scattering',
  'test_hankel_evaluator',
  'test_helm2_adapt_fac'
]

foreach example : examples
  executable(
    example,
    example + '.c',
    dependencies : [argtable3_dep, butterfly_dep, m_dep]
  )
endforeach

scripts = [
  'bf_lbo_plot_blocks.py',
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

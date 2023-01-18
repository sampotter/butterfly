examples = [
  'bf_one_block',
  'bf_all_blocks',
  'bf_lbo',
  'helm2_bie'
]

foreach example : examples
  executable(
    example,
    example + '.c',
    dependencies : [butterfly_dep, m_dep]
  )
endforeach

scripts = [
  'ellipse.py',
  'make_circle_test_data.py',
  'make_ellipse_test_data.py',
  'plot_bf_all_blocks.py',
  'plot_bf_one_block.py',
  'plot_blocks.py',
  'plot_points_and_quadtree.py',
  'run_bf_all_blocks_test.sh',
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

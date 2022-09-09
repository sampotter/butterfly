examples = [
  'bf_one_block',
  'bf_all_blocks',
  'helm2_bie'
]

foreach example : examples
  executable(
    example,
    example + '.c',
    dependencies : butterfly_dep,
  )
endforeach

python_scripts = [
  'ellipse.py',
  'make_test_data.py',
  'plot_bf_one_block.py',
  'plot_blocks.py',
  'plot_points_and_quadtree.py',
  'scratch.py',
  'test.py',
  'util.py',
]

foreach python_script : python_scripts
  custom_target(
    python_script,
    input : python_script, output : python_script,
    command : ['cp', '@INPUT@', '@OUTPUT@'],
    install : false,
    build_by_default : true
  )
endforeach

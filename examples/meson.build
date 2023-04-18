subdir('lbo')
subdir('multiple_scattering')
subdir('simple')
subdir('tree_evaluator')

# Extra files that should be copied over, such as common files,
# resources, etc.
extra_files = [
  'util.py',
]

foreach path : extra_files
  custom_target(
    path,
    input : path, output : path,
    command : ['cp', '@INPUT@', '@OUTPUT@'],
    install : false,
    build_by_default : true
  )
endforeach

subdir('fast_direct_solver')
subdir('fiedler_tree')
subdir('lbo')
if get_option('mfem').enabled()
  subdir('lbo_MFEM')
endif
subdir('multiple_scattering')
subdir('radiosity')
subdir('simple')
subdir('tree_evaluator')
subdir('covariance')

# Extra files that should be copied over, such as common files,
# resources, etc.
extra_files = [
  'util.py',
]

foreach extra_file : extra_files
  fs.copyfile(extra_file, extra_file)
endforeach

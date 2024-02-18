assert(get_option('python').enabled(),
       'The meson.build file for the Python wrapper should only be run if ' +
       'the "python" feature is enabled')

python3 = import('python').find_installation('python3')
python3_dep = python3.dependency(embed : true)

python3.extension_module(
  'butterfly',
  'butterfly.pyx',
  dependencies : [python3_dep, butterfly_dep],
  install : true
)

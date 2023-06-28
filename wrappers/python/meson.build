python = import('python').find_installation()

python_dep = python.dependency()

python.extension_module(
  'butterfly',
  'butterfly.pyx',
  dependencies : [python_dep, butterfly_dep]
)

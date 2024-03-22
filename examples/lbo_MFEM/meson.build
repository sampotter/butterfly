hypre_root = '/home/sfp/Build/hypre/src/hypre'
metis_root = '/home/sfp/Build/metis-4.0'
mfem_root = '/home/sfp/Build/mfem-4.5.2'

cpp = meson.get_compiler('cpp')

hypre_dep = cpp.find_library('HYPRE', dirs : [hypre_root + '/lib'])
metis_dep = cpp.find_library('metis', dirs : [metis_root])
mfem_dep = cpp.find_library('mfem', dirs : [mfem_root])
mpi_dep = dependency('mpi', language : 'cpp')

executable(
  'lbo_MFEM',
  'lbo_MFEM.cpp',
  dependencies : [hypre_dep, metis_dep, mfem_dep, butterfly_dep, mpi_dep],
  include_directories : [hypre_root + '/include', mfem_root]
)

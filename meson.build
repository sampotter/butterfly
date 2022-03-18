project(
  'butterfly',
  'c',
  default_options : [
    'c_std=c11'
  ]
)

add_global_arguments('-DBF_DEBUG', language : 'c')
add_global_arguments('-DBF_DOUBLE', language : 'c')

host_system = host_machine.system()
if host_system == 'linux'
  add_global_arguments('-DBF_LINUX', language : 'c')
elif host_system == 'darwin'
  add_global_arguments('-DBF_DARWIN', language : 'c')
else
  error('unsupported platform:', host_system)
endif

cc = meson.get_compiler('c')

cblas_dep = cc.find_library('openblas', required : true)
m_dep = cc.find_library('m', required : false)

butterfly_lib_src = [
  'bessel.c',
  'cheb.c',
  'dtype.c',
  'error.c',
  'fac.c',
  'geom.c',
  'helm2.c',
  'mat.c',
  'mat_block.c',
  'mat_block_coo.c',
  'mat_block_dense.c',
  'mat_block_diag.c',
  'mat_dense_complex.c',
  'mat_diag_real.c',
  'mat_product.c',
  'ptr_array.c',
  'quadtree.c',
  'rand.c',
  'splitmix64.c',
  'util.c',
  'vec.c',
  'xoshiro256plus.c',
]

butterfly_lib = library(
  'butterfly',
  butterfly_lib_src,
  dependencies : [cblas_dep, m_dep]
)

executable(
  'bf_one_block',
  'bf_one_block.c',
  dependencies : cblas_dep,
  link_with : butterfly_lib
)

executable(
  'bf_all_blocks',
  'bf_all_blocks.c',
  dependencies : cblas_dep,
  link_with : butterfly_lib
)

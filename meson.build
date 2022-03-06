project(
  'butterfly',
  'c',
  default_options : [
    'c_std=c11'
  ]
)

add_global_arguments('-DBF_DEBUG', language : 'c')
add_global_arguments('-DBF_DOUBLE', language : 'c')

cc = meson.get_compiler('c')

cblas_dep = dependency('openblas')
gsl_dep = dependency('gsl')
m_dep = cc.find_library('m', required : false)

butterfly_lib_src = [
  'dtype.c',
  'error.c',
  'fac.c',
  'geom.c',
  'helm2.c',
  'mat.c',
  'mat_block.c',
  'mat_block_coo.c',
  'mat_block_diag.c',
  'mat_dense_complex.c',
  'mat_diag_real.c',
  'mat_product.c',
  'ptr_array.c',
  'quadtree.c',
  'rand.c',
  'splitmix64.c',
  'util.c',
  'util_math.c',
  'vec.c',
  'xoshiro256plus.c',
]

butterfly_lib = library(
  'butterfly',
  butterfly_lib_src,
  dependencies : [cblas_dep, gsl_dep, m_dep]
)

executable('bf_one_block', 'bf_one_block.c', link_with : butterfly_lib)
executable('bf_all_blocks', 'bf_all_blocks.c', link_with : butterfly_lib)

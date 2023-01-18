project('butterfly', 'c', default_options : ['c_std=c11'])

cmake = import('cmake')

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
cholmod_dep = cc.find_library('cholmod', required : true)
m_dep = meson.get_compiler('c').find_library('m', required : false)
umfpack_dep = cc.find_library('umfpack', required : true)

arpack_ng_dep = dependency('arpack', required : true)
gsl_dep = dependency('gsl', required : true)

butterfly_lib_src = [
  'src/bbox.c',
  'src/bessel.c',
  'src/cheb.c',
  'src/circle.c',
  'src/cholmod.c',
  'src/dtype.c',
  'src/error.c',
  'src/fac.c',
  'src/fac_streamer.c',
# 'src/fiedler_tree.c',
  'src/helm2.c',
  'src/interval_tree.c',
  'src/interval_tree_node.c',
  'src/lbo.c',
  'src/mat.c',
  'src/mat_block.c',
  'src/mat_block_coo.c',
  'src/mat_block_dense.c',
  'src/mat_block_diag.c',
  'src/mat_coo_complex.c',
  'src/mat_coo_real.c',
  'src/mat_csr_real.c',
  'src/mat_csr_real.cholmod.c',
  'src/mat_csr_real.umfpack.c',
  'src/mat_dense.c',
  'src/mat_dense_complex.c',
  'src/mat_dense_real.c',
  'src/mat_diag_real.c',
  'src/mat_eig.arpack.c',
  'src/mat_givens.c',
  'src/mat_product.c',
  'src/mat_solve.c',
  'src/mat_sum.c',
  'src/mat_util.c',
  'src/mat_zero.c',
  'src/octree.c',
  'src/octree_node.c',
  'src/perm.c',
  'src/points.c',
  'src/ptr_array.c',
  'src/quadrature.c',
  'src/quadtree.c',
  'src/quadtree_node.c',
  'src/rand.c',
  'src/splitmix64.c',
  'src/tree.c',
  'src/tree_iter.c',
  'src/tree_level_iter.c',
  'src/tree_node.c',
  'src/trimesh.c',
  'src/types.c',
  'src/util.c',
  'src/vec.c',
  'src/vec_complex.c',
  'src/vec_real.c',
  'src/vec_zero.c',
  'src/vectors.c',
  'src/xoshiro256plus.c',
]

butterfly_inc = include_directories('include')

butterfly_lib = library(
  'butterfly',
  butterfly_lib_src,
  dependencies : [arpack_ng_dep, cblas_dep, cholmod_dep, gsl_dep, m_dep],
  include_directories : [butterfly_inc]
)

butterfly_dep = declare_dependency(
  link_with : butterfly_lib,
  include_directories : butterfly_inc
)

subdir('examples')

project('butterfly', 'c', default_options : ['c_std=c11'])

pkg = import('pkgconfig')

if get_option('buildtype') == 'debug'
  add_project_arguments('-DBF_DEBUG', language : 'c')
endif

add_project_arguments('-DBF_DOUBLE', language : 'c')

host_system = host_machine.system()
if host_system == 'linux'
  add_project_arguments('-DBF_LINUX', language : 'c')
elif host_system == 'darwin'
  add_project_arguments('-DBF_DARWIN', language : 'c')
else
  error('unsupported platform:', host_system)
endif

cc = meson.get_compiler('c')

openblas_dep = cc.find_library('openblas', required : true)
# cholmod_dep = cc.find_library('cholmod', required : true)
m_dep = meson.get_compiler('c').find_library('m', required : false)
umfpack_dep = cc.find_library('umfpack', required : true)

fmm2d_dep = dependency('fmm2d', required : true)
argtable3_dep = dependency('argtable3', required : true)
arpack_ng_dep = dependency('arpack', required : true)
gsl_dep = dependency('gsl', required : true)
# openblas_dep = dependency('openblas', required : true)
# umfpack_dep = dependency('umfpack', required : true)

butterfly_lib_src = [
  'src/array.c',
  'src/bbox.c',
  'src/bessel.c',
  'src/cheb.c',
  'src/circle.c',
#  'src/cholmod.c',
  'src/dtype.c',
  'src/ellipse.c',
  'src/error.c',
  'src/eval_tree.c',
  'src/fac.c',
  'src/fac_helm2.c',
  'src/fac_span.c',
  'src/fac_streamer.c',
  'src/fiedler_tree.c',
  'src/fiedler_tree_node.c',
  'src/helm2.c',
  'src/indexed_mat.c',
  'src/interval_tree.c',
  'src/interval_tree_node.c',
  'src/lbo.c',
  'src/linalg.c',
  'src/logging.c',
  'src/lu.c',
  'src/lu_csr_real.c',
  'src/lu_dense_complex.c',
  'src/mat.c',
  'src/mat_block.c',
  'src/mat_block_coo.c',
  'src/mat_block_dense.c',
  'src/mat_block_diag.c',
  'src/mat_coo_complex.c',
  'src/mat_coo_real.c',
  'src/mat_csr_real.c',
#  'src/mat_csr_real.cholmod.c',
  'src/mat_csr_real.umfpack.c',
  'src/mat_dense.c',
  'src/mat_dense_complex.c',
  'src/mat_dense_real.c',
  'src/mat_diag_real.c',
  'src/mat_func.c',
  'src/mat_givens.c',
  'src/mat_identity.c',
  'src/mat_perm.c',
  'src/mat_product.c',
  'src/mat_sum.c',
  'src/mat_zero.c',
  'src/mem.c',
  'src/node_array.c',
  'src/octree.c',
  'src/octree_node.c',
  'src/perm.c',
  'src/points.c',
  'src/poisson_disk_sampling.c',
  'src/ptr_array.c',
  'src/quadrature.c',
  'src/quadtree.c',
  'src/quadtree_node.c',
  'src/rand.c',
  'src/real_array.c',
  'src/size_array.c',
  'src/splitmix64.c',
  'src/tree.c',
  'src/tree_iter.c',
  'src/tree_iter_post_order.c',
  'src/tree_level_iter.c',
  'src/tree_node.c',
  'src/trimesh.c',
  'src/trimesh.get_level_set_submesh.c',
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

butterfly_dependencies = [
  arpack_ng_dep,
#  cholmod_dep,
  gsl_dep,
  m_dep,
  openblas_dep,
  umfpack_dep
]

butterfly_include_directories = [
  butterfly_inc
]

butterfly_lib = library(
  'butterfly',
  butterfly_lib_src,
  dependencies : butterfly_dependencies,
  include_directories : butterfly_include_directories,
  install : true
)

install_subdir('include/bf', install_dir : 'include', install_tag : 'devel')

butterfly_dep = declare_dependency(
  link_with : butterfly_lib,
  include_directories : butterfly_inc
)

pkg.generate(butterfly_lib)

subdir('examples')

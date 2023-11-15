cmocka_dep = dependency('cmocka')

test_files = [
  'sphere.obj',
  'sphere_Lam.txt',
  'sphere_Phi.txt'
]

foreach test_file : test_files
  fs.copyfile(test_file, test_file)
endforeach

# Simple tests with no arguments:

test_names = [
  'test_disjoint_interval_list',
  'test_interval',
]

test_deps = [butterfly_dep, cmocka_dep, m_dep]

foreach test_name : test_names
  test_exe = executable(test_name, test_name + '.c', dependencies : test_deps)
  test(test_name, test_exe)
endforeach

# Tests with arguments require more involved setup:

test_linalg = executable('test_linalg', 'test_linalg.c', dependencies : test_deps)
test('test_linalg', test_linalg)

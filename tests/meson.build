cmocka_dep = dependency('cmocka')

test_names = [
  'test_interval',
  'test_disjoint_interval_list'
]

test_deps = [butterfly_dep, cmocka_dep, m_dep]

foreach test_name : test_names
  test_exe = executable(test_name, test_name + '.c', dependencies : test_deps)
  test(test_name, test_exe)
endforeach

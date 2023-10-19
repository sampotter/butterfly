cmocka_dep = dependency('cmocka')

test_interval = executable('test_interval', 'test_interval.c', dependencies : [butterfly_dep, cmocka_dep, m_dep])

test('interval', test_interval)

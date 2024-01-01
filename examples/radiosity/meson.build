examples = [
  'radiosity'
]

example_deps = [butterfly_dep, m_dep]

foreach example : examples
  executable(example, example + '.c', dependencies : example_deps)
endforeach

fs.copyfile('radiosity.py', 'radiosity.py')
fs.copyfile('67p.obj', '67p.obj')

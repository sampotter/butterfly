executable('fast_direct_solver', ['fast_direct_solver.c'],
           dependencies : [butterfly_dep])

scripts = [
  'fast_direct_solver.py'
]

foreach script : scripts
  custom_target(
    script,
    input : script, output : script,
    command : ['cp', '@INPUT@', '@OUTPUT@'],
    install : false,
    build_by_default : true
  )
endforeach

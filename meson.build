project('gmpy2', 'c',
        version: run_command('python', '-m', 'setuptools_scm',
                             check: true).stdout().strip())
add_global_arguments('-DSHARED=1', language : 'c')
py = import('python').find_installation(pure: true)
libgmp = dependency('gmp', version: '>= 6.2.0')
libmpfr = dependency('mpfr', version: '>= 4.1.0')
libmpc = dependency('mpc', version: '>= 1.2.1')
py.install_sources(['src/gmpy2.pxd', 'src/gmpy2.h'])
py.extension_module('gmpy2', ['src/gmpy2.c'], install: true,
                    dependencies: [libgmp, libmpfr, libmpc])

#!/bin/bash

set -e -x

wheel=$WHEELNAME
dest_dir=$WHEELHOUSE

delvewheel repair ${wheel} -w ${dest_dir} --add-path .local/bin --no-mangle-all

cp .local/bin/{gmp,mpfr,mpc}.lib ${dest_dir}
(cd ${dest_dir}; wheel unpack --dest . gmpy2-*.whl; rm gmpy2-*.whl; mv *.lib gmpy2-*/gmpy2.libs; wheel pack gmpy2-*)

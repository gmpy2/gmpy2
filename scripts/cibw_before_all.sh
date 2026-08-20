#!/bin/bash

set -e -x

GMP_VERSION=6.3.0
MPFR_VERSION=4.2.2
MPC_VERSION=1.4.1

PREFIX="$(pwd)/.local/"

CURL_OPTS="--fail --location --retry 4 --connect-timeout 32"

download () {
  sleep_ivl=16
  until curl ${CURL_OPTS} --remote-name $1
  do
    sleep ${sleep_ivl}
    sleep_ivl=$((${sleep_ivl}*2))
  done
}
genlib () {
  cd .local/bin/
  dll_file=$1
  lib_name=$(basename -s .dll ${dll_file})
  name=$(echo ${lib_name}|sed 's/^lib//;s/-[0-9]\+//')

  gendef ${dll_file}
  dlltool -d ${lib_name}.def -l ${name}.lib

  cp ${name}.lib ../lib/
  cp ${dll_file} ../lib/
  cd ../../
}

# -- build GMP --
download https://ftp.gnu.org/gnu/gmp/gmp-${GMP_VERSION}.tar.xz
tar -xf gmp-${GMP_VERSION}.tar.xz
cd gmp-${GMP_VERSION}
# Patch the mp_bitcnt_t to "unsigned long long int" on WINDOWS AMD64:
#patch -N -Z -p0 < ../scripts/mp_bitcnt_t.diff

patch -N -Z -p0 < ../scripts/fat_build_fix.diff
patch -N -Z -p0 < ../scripts/dll-importexport.diff
patch -N -Z -p1 < ../scripts/gcc15.diff

CONFIG_ARGS="--enable-shared --disable-static --with-pic --prefix=$PREFIX"
if [ "$OSTYPE" = "msys" ] || [ "$OSTYPE" = "cygwin" ]
then
  if [ "${RUNNER_ARCH}" = "ARM64" ]
  then
    autoreconf -fi
    CONFIG_ARGS="${CONFIG_ARGS} --disable-assembly"
  else
    CONFIG_ARGS="${CONFIG_ARGS} --enable-fat"
  fi
else
  CONFIG_ARGS="${CONFIG_ARGS} --enable-fat"
fi
# config.guess uses microarchitecture and configfsf.guess doesn't
# We replace config.guess with configfsf.guess to avoid microarchitecture
# specific code in common code.
rm config.guess && mv configfsf.guess config.guess && chmod +x config.guess
./configure ${CONFIG_ARGS}
make -j6
make install
cd ../

# -- build MPFR --
download https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz
tar -xf mpfr-${MPFR_VERSION}.tar.gz
cd mpfr-${MPFR_VERSION}
./configure --enable-shared \
            --disable-static \
            --with-pic \
            --with-gmp=$PREFIX \
            --prefix=$PREFIX
make -j6
make install
cd ../
# -- build MPC --
download https://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.xz
tar -xf mpc-${MPC_VERSION}.tar.xz
cd mpc-${MPC_VERSION}

./configure --enable-shared \
            --disable-static \
            --with-pic \
            --with-gmp=$PREFIX \
            --with-mpfr=$PREFIX \
            --prefix=$PREFIX
make -j6
make install

cd ../

# -- copy headers --
cp $PREFIX/include/{gmp,mpfr,mpc}.h gmpy2/

# -- generate *.lib files from *.dll on M$ Windows --
if [ "$OSTYPE" = "cygwin" ]
then
  genlib libgmp-10.dll
  genlib libmpfr-6.dll
  genlib libmpc-3.dll
fi

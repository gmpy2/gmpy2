:orphan:

Building for Windows
====================

A new build approach is being tested for the final gmpy2 v2.2.x release.

The GMP, MPFR, and MPC DLLs are created using the MinGW64 compiler suite and
are linked to UCRT.  The final gmpy2 binaries are created using Visual Studio.


MSYS2 Steps
-----------

1.  Install MSYS2

    1. Upgrade the initial installation::

        pacman -Syuu

    2. Reboot and launch the "MSYS2 UCRT" terminal window.  Upgrade the user
       system::

        pacman -Syu

    3. Install the MinGW64 compiler suite that links with ucrt::

        pacman -Sy mingw-w64-ucrt-x86_64-gcc

    4. Install rquired tools::

        pacman -Sy patch m4 lzip wget tar make diffutils git

2.  Create a temporary working directory for the build::

      mkdir ~/temp

3.  Clone (or download) the full gmpy2 source code::

      cd ~/temp
      git clone https://github.com/gmpy2/gmpy2.git
      cd gmpy

4.  Download and compile GMP, MPFR, and MPC:

      bash scripts/cibw_before_all.sh


Visual Studio Steps
-------------------

1.  Launch the "Visual Studio x64 Native Tools Command Prompt" shell.

2.  Change to location of the build directory::

      cd c:\msys64\home\<username>\temp\gmpy

3.  Run::

      scripts\cibw_before_all_windows.bat

4.  Compile and install Windows binary wheels, e.g. for CPython 3.12::

      py -3.12 -m pip install --upgrade build pytest hypothesis cython mpmath
      set CIBUILDWHEEL=1
      py -3.12 -m build --wheel --no-isolation
      py -3.12 -m pip install dist\gmpy2-*312-*.whl --force-reinstall

5.  Run the test suite::

      py -3.12 -m pytest test/

6.  Compile the C-API demo program::

      cd demo
      py -3.12 setup.py build_ext
      cd build\lib.win-amd64-cpython-312
      py -3.12
      >>> import gmpy2
      >>> import gmpy2_demo
      >>> dir(gmpy2_demo)
      ['__doc__', '__file__', '__loader__', '__name__', '__package__', '__spec__', 'factor']
      >>> gmpy2_demo.factor(123456789)
      [(mpz(3), mpz(2)), (mpz(3607), mpz(1)), (mpz(3803), mpz(1))]

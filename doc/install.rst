============
Installation
============

Easy install
============

* Required packages in python: `numpy`, `scipy`, `h5py`

  * we suggest using Anaconda_ distribute, which includes most packages, and 
    provides you a user specific environment, i.e., all your 
    python packages go to a single folder. Thus, you don't need the root to 
    install packages.

  .. _Anaconda: http://continuum.io/downloads

* You can install `Vireo` simply via pypi in terminal (suggested), or upgrade 
  by add ``--upgrade`` as follows:

  ::

    pip install vireoSNP

    pip install --upgrade --no-deps vireoSNP


Source code
===========

* Alternatively, you also could download the source code via GitHub with the 
  latest version) and run python setup in terminal:

* GitHub: https://github.com/huangyh09/vireo

  ::

    python setup.py install

* In any case, if had the permission error for installation as you are not 
  root, add ``--user``.


Test
====

In order to test the installation, you could type ``vireo``. If successful, you
will see the following output.

.. code-block:: html

  Welcome to vireoSNP v0.1.1!

  use -h or --help for help on argument.

If installation is sucessful, but can't run it, then check whether the directory 
which contains the executable binary file is added to PATH environment. 

.. code-block:: html

  vireo: command not found

Usually the directory is ``~/.local/bin`` if you don't use Anaconda. You could 
add the path into PATH environment variable, by write the following line into 
``.profile`` or ``.bashrc`` file.

:: 
  
  export PATH="~/.local/bin:$PATH"


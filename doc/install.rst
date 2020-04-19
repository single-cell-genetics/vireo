============
Installation
============

**Required packages in python:** ``numpy>=1.9.0``, ``scipy>=1.0``, ``matplotlib``

**Environment:** we only tested Vireo in ``Python 3`` environment, so if it 
fails in Python 2, please try it in Python 3 before reporting the issue. 

We recommend using Anaconda_ distribute to set up the environment. It not only
includes all dependent packages, but also provides a user controlled 
environment, namely, you will have the root permission for this distribution, 
including installation of any package.

.. _Anaconda: http://continuum.io/downloads

Easy install from PyPI
======================

You can install `Vireo` simply via PyPI_ in terminal (**suggested**), or upgrade 
by adding ``--upgrade`` as follows:

::
  
  pip install vireoSNP

  pip install --upgrade --no-deps vireoSNP

.. _PyPI: https://pypi.org/project/vireoSNP


Install from source code
========================

Alternatively, you can download the source code from GitHub_ (for the 
latest version) and run python setup in terminal:

.. _GitHub: https://github.com/huangyh09/vireo

::
  
  wget https://github.com/huangyh09/vireo/archive/master.zip
  unzip master.zip
  cd vireo-master

  python setup.py install

You can also use the following shortcut

.. code-block:: bash

  pip install -U git+https://github.com/single-cell-genetics/vireo

In any case, if had the permission error for installation as you are not root, 
add ``--user``.


Quick check
===========

In order to test the installation, you could type ``vireo`` in terminal. If 
successfully installed, you will see the following output.

.. code-block:: html

  Welcome to vireoSNP v0.1.1!

  use -h or --help for help on argument.

If installation is sucessful, but can't run it (e.g., message below), then 
check whether the directory which contains the executable binary file is added 
to PATH environment. 

.. code-block:: html

  vireo: command not found

If using Anaconda, the executable ``vireo`` is located in 
``$anaconda3/bin/vireo``. 
If not using Anaconda, it is usually located in directory ``~/.local/bin``. You 
could add the path into PATH environment variable, by write the following line 
into ``.profile`` or ``.bashrc`` file.

:: 
  
  export PATH="~/.local/bin:$PATH"


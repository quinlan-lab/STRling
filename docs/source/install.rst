Install
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


We recommending downloading the static binary.

Download the `strling` binary from the latest release `here <https://github.com/quinlan-lab/STRling/releases/latest>`_.

Make it executable:
`chmod +x strling`

Install with conda
------------------

We suggest conda installation if you want to perform outlier analysis, which requires a python script.

.. code-block:: bash
  conda config --add channels bioconda
  conda install -y strling

The following should now be in your path:  
`strling`  
`strling-outliers.py`

Note: It should be safe to ignore "Solving environment: failed ..." messages as long as the message "Solving environment: done" appears following the failed message and the install proceeds.

Install from source
-------------------

Install nim:
`curl https://nim-lang.org/choosenim/init.sh -sSf > init.sh && sh init.sh`

Install STRling:

.. code-block:: bash

  git clone <URL>
  cd STRling
  nimble install

Compile options for development:

Compile in fast mode (danger):
`nim c -d:danger -d:release src/strling.nim`

Install
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


We recommending downloading the static binary.

Download the `strling` binary from the latest release `here <https://github.com/quinlan-lab/STRling/releases/latest>`_.

Make it executable:
`chmod +x strling`

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

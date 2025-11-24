Contribute
==========

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Install Python dependencies
--------------------------

Using conda:
.. code-block:: bash

    conda env create -f environment.yml
    conda activate strling

Creating a STRling release
--------------------------

- Rebase dev on master as needed
- Increment version in two places (using the appropriate semantic version https://semver.org/):

`strling.nimble:version       = "0.4.1"`

`src/strpkg/version.nim:const strlingVersion* = "0.4.1"`

- PR/merge dev into master
- Make a release draft
- Add static binaries (see below)

Creating a static binary
------------------------

See https://github.com/brentp/hts-nim#static-binary-with-singularity for more details.

On branch master:

.. code-block:: bash

    git pull
    module load singularity

    #debug binary
    singularity run --bind $(pwd):/load --bind /scratch --bind /uufs 'docker://brentp/musl-hts-nim:latest' /usr/local/bin/nsb -n strling.nimble -s src/strling.nim -- -d:danger -d:release -d:debug
    cp strling strling_debug

    #regular binary
    singularity run --bind $(pwd):/load --bind /scratch --bind /uufs 'docker://brentp/musl-hts-nim:latest' /usr/local/bin/nsb -n strling.nimble -s src/strling.nim -- -d:danger -d:release


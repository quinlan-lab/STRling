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


Continuous Integration
----------------------

STRling now uses GitHub Actions for continuous integration instead of Travis CI.
The workflow file lives at ``.github/workflows/ci.yml`` and runs:

* Nim build and unit tests (``nimble test``)
* Optional end-to-end extraction/call tests on pushes to ``master``
* htslib build from source (matching previous Travis configuration)

To view the current build status, see the CI badge in the README or visit:
``https://github.com/quinlan-lab/STRling/actions``

If you add new dependencies or tests, update ``strling.nimble`` and (if needed)
the workflow steps. Keep heavy end-to-end tests in the integration job so that
pull requests remain fast.


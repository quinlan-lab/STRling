# Package
template thisModuleFile: string = instantiationInfo(fullPaths = true).filename
#

#when fileExists(thisModuleFile.parentDir / "src/strling.nim"):
#  # In the git repository the Nimble sources are in a ``src`` directory.
#  import src/strpkg/version as _
#else:
#  # When the package is installed, the ``src`` directory disappears.
#  import strpkg/version as _

version       = "0.5.2"
author        = "Harriet and Brent"
description   = "Novel STR caller"
license       = "MIT"


# Dependencies

requires "nim >= 0.18.0", "kmer >= 0.2.2", "hts", "itertools", "argparse >= 0.7.0 & < 1.0", "msgpack4nim >= 0.4.2", "lapper"
bin = @["strling"]

srcDir = "src"

skipDirs = @["tests"]

task test, "run the tests":
  exec "nim c -x:on --lineDir:on --debuginfo -r tests/all"



# Package

version       = "0.0.1"
author        = "Harriet and Brent"
description   = "work in progress STR caller"
license       = "MIT"


# Dependencies

requires "nim >= 0.18.0", "kmer", "hts", "itertools", "argparse"
bin = @["str"]

srcDir = "src"

skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"



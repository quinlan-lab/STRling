# Package

version       = "0.0.1"
author        = "Harriet and Brent"
description   = "work in progress STR caller"
license       = "MIT"


# Dependencies

requires "nim >= 0.18.0", "kmer >= 0.2.2", "hts", "itertools", "argparse", "msgpack4nim"
bin = @["str"]

srcDir = "src"

skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"



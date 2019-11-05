# Package

version       = "0.0.1"
author        = "Harriet and Brent"
description   = "Novel STR caller"
license       = "MIT"


# Dependencies

requires "nim >= 0.18.0", "kmer >= 0.2.2", "hts", "itertools", "argparse", "msgpack4nim", "lapper"
bin = @["strling"]

srcDir = "src"

skipDirs = @["tests"]

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"



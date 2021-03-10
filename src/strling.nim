import tables
import os
import strformat

import ./strpkg/extract
import ./strpkg/merge
import ./strpkg/call
import ./strpkg/version
import ./strpkg/extract_region
import ./strpkg/genome_strs

proc main*() =

  type pair = object
    f: proc()
    description: string

  var dispatcher = {
    "extract": pair(f:extract_main, description:"extract informative STR reads from a BAM/CRAM. This is a required first step."),
    "merge": pair(f:merge_main, description:"merge putitive STR loci from multiple samples. Only required for joint calling."),
    "call": pair(f:call_main, description:"call STRs"),
    "index": pair(f:index_main, description:"identify large STRs in the reference genome, to produce ref.fasta.str."),
    "pull_region": pair(f:extract_region_main, description:"for debugging; pull all reads (and mates) for a given regions"),
    }.toOrderedTable
  var args = commandLineParams()

  stderr.write_line &"\nstrling version: {strlingVersion}"
  when defined(debug):
    stderr.write_line &"compiled in debug mode"
  if len(args) == 0 or not (args[0] in dispatcher):
    stderr.write_line "\nCommands: "
    for k, v in dispatcher:
      echo &"  {k:<13}:   {v.description}"
    if len(args) > 0 and (args[0] notin dispatcher):
      if args[0] in @["-h", "--help"]:
        quit(0)
      else:
        echo &"unknown program '{args[0]}'"
    quit &"ERROR: please enter a valid command"

  dispatcher[args[0]].f()

when isMainModule:
  main()


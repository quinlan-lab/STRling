import ./strpkg/extract
import ./strpkg/call
import tables
import os
import strformat

proc main*() =

  type pair = object
    f: proc()
    description: string

  var dispatcher = {
    "extract": pair(f:extract_main, description:"extract informative STR reads from a BAM/CRAM. This is a required first step."),
    "call": pair(f:call_main, description:"call STRs"),
    }.toOrderedTable
  var args = commandLineParams()

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


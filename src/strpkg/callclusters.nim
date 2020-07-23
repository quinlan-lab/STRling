import algorithm
import hts/bam
import strformat

import ./cluster
import ./collect
import ./utils
import ./genotyper
import ./genome_strs
import ./extract

# returns reads assign to locus and removes them from table.
# Updates Bounds with number of supporting reads
proc assign_reads_locus*(locus: var Bounds, treads_by_tid_rep: TableRef[tid_rep, seq[tread_id]]): seq[tread] =
    var key: tid_rep = (locus.tid, locus.repeat.as_array)
    var trs = treads_by_tid_rep.getOrDefault(key, @[])

    var left_most = if locus.left_most == 0: 0'u32 else: locus.left_most - 1'u32
    var li = lowerBound(trs, tread_id(tr: tread(position: left_most)), proc(a, b:tread_id):int =
      return cmp(a.tr.position, b.tr.position)
    )
    var ri = upperBound(trs, tread_id(tr: tread(position: locus.right_most)), proc(a, b:tread_id):int =
      return cmp(a.tr.position, b.tr.position)
    )
    when defined(debug):
      if ri - li > 0:
        stderr.write_line &"[strling] got {ri - li} treads for locus: {locus} with indexes {li}..{ri}"

    if trs.len > 0:
      # now we have a the subset of treads that support the given locus
      for tr_id in trs[li..<ri]:
        result.add(tr_id.tr)

      # remove these from the table
      treads_by_tid_rep[key] = trs[0..<li]
      if ri < trs.high:
        treads_by_tid_rep[key].add(trs[ri + 1..trs.high])

    # Force reporting of loci from input files
    # even if they have insufficient supporting reads
    locus.force_report = true 
    # Update read counts in Bounds
    locus.n_total = 0
    locus.n_left = 0
    locus.n_right = 0
    for r in result:
      locus.n_total.inc
      if r.split == Soft.right:
        locus.n_right.inc
      elif r.split == Soft.left:
        locus.n_left.inc

proc bounds*(c: Cluster, min_clip: uint16, min_clip_total: uint16): (Bounds, bool) =
  if c.reads.len >= uint16.high.int:
    stderr.write_line "More than " & &"{uint16.high.int}" & " reads in cluster with first read:" & $c.reads[0] & " skipping"
    return
  var b = c.bounds
  if b.right - b.left > 1000'u32:
    stderr.write_line "large bounds:" & $b & " skipping"
    return
  # require left and right support unless force_report is set
  # i.e. when Bounds is from an input file rarather than a new cluster
  if not b.force_report:
    if b.n_left < min_clip: return
    if b.n_right < min_clip: return
    if (b.n_right + b.n_left) < min_clip_total: return
  return (b, true)

#TODO: put common genotyping code in a function
proc genotype_bounds*() =
  discard

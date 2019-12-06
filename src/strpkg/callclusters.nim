import algorithm
import hts/bam
import strformat

import ./cluster
import ./collect
import ./utils
import ./genotyper
import ./genome_strs
import ./extract

# returns reads assign to locus and removes them from table
proc assign_reads_locus*(locus: Bounds, treads_by_tid_rep: TableRef[tid_rep, seq[tread]]): seq[tread] =
    var key: tid_rep = (locus.tid, locus.repeat.as_array)
    var trs = treads_by_tid_rep.getOrDefault(key, @[])
    
    var li = lowerBound(trs, tread(position: locus.left_most - 1), proc(a, b:tread):int =
      return cmp(a.position, b.position)
    )
    var ri = upperBound(trs, tread(position: locus.right_most), proc(a, b:tread):int =
      return cmp(a.position, b.position)
    )
    stderr.write_line "[strling] got {ri - li} treads for locus: {locus} with indexes {li}..{ri}"
    
    if trs.len > 0:
      # now we have a the subset of treads that support the given locus
      result = trs[li..ri]

      #XXX TODO: Count up the number of supporting reads and record in Bounds

      # remove these from the table
      treads_by_tid_rep[key] = trs[0..<li]
      if ri < trs.high:
        treads_by_tid_rep[key].add(trs[ri + 1..trs.high])

proc check_cluster*(c: Cluster, min_clip: uint16, min_clip_total: uint16): (Bounds, bool) =
  if c.reads.len >= uint16.high.int:
    stderr.write_line "More than " & &"{uint16.high.int}" & " reads in cluster with first read:" & $c.reads[0] & " skipping"
    return
  var b = c.bounds
  if b.right - b.left > 1000'u32:
    stderr.write_line "large bounds:" & $b & " skipping"
    return
  # require left and right support
  if not b.force_report:
    if b.n_left < min_clip: return
    if b.n_right < min_clip: return
    if (b.n_right + b.n_left) < min_clip_total: return
  return (b, true)

#TODO: put common genotyping code in a function
proc genotype_bounds*() =
  discard

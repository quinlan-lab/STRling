
type Event = enum
   Insertion
   Deletion

type Genotype = object
   event: Event
   copy_number: int
   GQ: float

type Call = object
   chrom: string
   start: int
   stop: int
   # and confidence intervals around position and size
   genotype: Genotype
   quality: float
   # ...

proc genotype*(b:Bounds, tandems: seq[tread]): Call =
   # TBD
   discard

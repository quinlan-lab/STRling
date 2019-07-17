proc complement*(s:char): char {.inline.} =
    if s == 'C':
        return 'G'
    elif s == 'G':
        return 'C'
    elif s == 'A':
        return 'T'
    elif s == 'T':
        return 'A'
    else:
        return s
proc reverse_complement*(xs: string): string =
  result = newString(xs.len)
  for i, x in xs:
    # high == len - 1
    result[xs.high-i] = complement(x)


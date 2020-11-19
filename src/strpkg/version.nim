const strlingVersion* = "0.4.0"

# bin file format version
const thisFmtVersion* = 0'i16

proc asArray9*(s:string): array[9, char] =
  for i, c in s:
    result[i] = c

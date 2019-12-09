const strlingVersion* = "0.0.2"

proc asArray9*(s:string): array[9, char] =
  for i, c in s:
    result[i] = c

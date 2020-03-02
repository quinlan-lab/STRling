import ./cluster
import math
import hts/bam

type cumulative_dist* = array[4096, float32]

proc cumulative*(frag_dist: array[4096, uint32]): cumulative_dist =
  # this returns a cumulative distribution that indicates the proportion of
  # fragments with length less than the given index.
  for i in 0..<frag_dist.len:
    const window = 11
    for j in max(0, i - window)..min(i+window, frag_dist.high):
      result[i] += frag_dist[j].float32

  result.cumsum
  var fmax = result[result.high].float32
  for i, v in result:
    result[i] = v.float32 / fmax

proc expected_spanning_probability*(cd: cumulative_dist, read:Record, event_start: int, event_stop:int=event_start+1, min_spanning_bases:int=20): float =

  var dist: int
  if read.start.int < event_stop - min_spanning_bases:
    if read.flag.reverse: return
    dist = event_start - read.start.int
    if dist < 0: return
    # below we check dist + event size is >= min_spanning_bases so that left can align:
    # event                         |    |
    #                  <------{dist}---->
    # read       AAAAAA>
    if dist + (event_stop - event_start).int < min_spanning_bases:
      return

  else:
    if not read.flag.reverse: return
    dist = read.stop.int - event_stop
    if dist < 0: return
    if dist + (event_stop - event_start).int < min_spanning_bases:
      return

  # require e.g. 20 bases on RHS of event to allow a mapping chunk from bwa mem
  dist += min_spanning_bases
  dist += (event_stop - event_start)
  if dist < 0 or dist > cd.high:
    return

  # cd dist of 0.9 means 90% of fragments less than distance
  # so we do 1 - to get the 10% that would span the dist.
  return 1 - cd[dist]

when isMainModule:
  let frag_sizes: array[4096, uint32] = [0'u32,0,2,1,2,1,2,1,4,0,2,1,4,4,2,1,5,0,5,3,6,6,5,4,6,5,5,0,3,4,6,7,5,6,3,4,9,5,8,7,4,3,4,2,9,1,5,10,8,6,5,12,4,9,11,10,10,17,13,19,11,13,14,22,12,21,12,17,18,6,15,16,24,18,11,34,20,22,31,19,33,30,24,41,34,23,39,34,33,45,48,40,40,57,36,41,52,60,48,42,51,61,57,69,58,73,73,80,79,82,84,106,103,93,95,111,107,121,113,116,140,151,156,152,184,173,176,155,192,177,171,173,202,226,264,205,262,247,225,251,255,286,271,265,297,295,255,290,325,312,301,308,319,310,332,314,319,355,325,354,365,372,388,385,409,349,389,368,385,401,442,427,412,456,410,413,419,453,430,430,456,517,475,436,493,468,453,490,487,467,435,489,473,479,469,530,501,529,516,493,485,511,545,560,532,503,531,512,537,500,519,569,565,529,554,547,550,557,531,556,516,607,588,601,517,568,582,589,564,608,618,620,632,630,641,644,602,607,598,610,636,618,619,670,617,630,659,638,659,643,664,651,712,679,686,679,660,685,671,676,741,690,722,778,698,735,743,767,755,744,796,809,771,786,777,775,756,828,771,833,819,877,804,871,785,862,864,870,880,913,853,922,989,935,978,896,969,978,1003,964,941,995,962,1070,960,1053,1082,1003,1015,1058,1055,1050,1062,1092,1071,1156,1167,1110,1119,1086,1139,1167,1159,1167,1166,1225,1246,1206,1132,1270,1285,1238,1274,1357,1246,1277,1257,1312,1362,1324,1388,1387,1460,1376,1507,1455,1449,1428,1502,1454,1413,1515,1505,1480,1519,1563,1558,1598,1587,1652,1687,1681,1626,1697,1629,1772,1707,1786,1765,1814,1729,1798,1843,1816,1838,1959,1913,1846,1877,1908,1878,1956,2009,2059,1972,1934,2093,2055,2138,2100,2108,2115,2152,2152,2194,2236,2253,2306,2288,2277,2280,2417,2345,2323,2423,2467,2368,2417,2558,2467,2536,2544,2511,2571,2673,2595,2624,2749,2714,2693,2617,2724,2823,2812,2848,2831,2880,2790,2970,2918,2960,2917,3050,2973,3082,3042,3115,3162,3204,3188,3274,3277,3180,3317,3244,3230,3388,3380,3476,3410,3464,3480,3389,3578,3483,3561,3584,3483,3655,3623,3658,3657,3706,3890,3862,3737,3828,3823,3935,4012,3909,4035,4015,3986,4116,4135,4104,4071,4150,4198,4144,4417,4353,4366,4393,4275,4412,4472,4475,4428,4439,4462,4518,4446,4499,4562,4661,4710,4770,4717,4708,4773,4850,4774,4896,5001,4855,4759,4947,5136,5038,5144,5114,5065,5207,5104,5107,5095,5128,5199,5329,5134,5314,5353,5480,5400,5439,5594,5475,5454,5589,5452,5456,5526,5568,5548,5592,5585,5519,5568,5661,5696,5698,5798,5591,5719,5688,5749,5931,5660,5850,5770,5695,5837,5876,5885,5730,5870,5910,5902,5963,5820,5786,5959,5940,5879,5948,5784,5822,5958,5879,6016,6012,6042,5861,5830,5946,5877,5965,5917,5970,5940,5830,5836,6021,5925,5904,5775,6006,5944,6014,5905,5855,5937,5743,5828,5869,5723,5688,5876,5933,5897,5697,5722,5575,5656,5615,5748,5598,5689,5650,5554,5557,5463,5707,5472,5530,5440,5281,5525,5363,5292,5413,5497,5315,5372,5239,5289,5400,5292,5253,5048,5203,5060,5073,5186,5037,5031,5061,5027,4944,4931,4829,4725,4824,4842,4881,4791,4738,4673,4679,4666,4651,4534,4598,4606,4433,4541,4492,4453,4590,4551,4585,4370,4418,4468,4272,4224,4246,4074,4181,4165,4141,4002,4101,4023,3993,4028,3864,4061,3852,3983,3868,3799,3771,3833,3776,3689,3705,3704,3612,3570,3643,3650,3482,3640,3468,3494,3618,3501,3353,3451,3386,3331,3281,3236,3142,3233,3232,3219,3069,3208,3219,2994,3064,3025,2993,2943,2965,2946,2871,2871,2797,2911,2803,2782,2721,2720,2712,2784,2704,2762,2684,2686,2525,2606,2583,2518,2599,2526,2372,2473,2409,2397,2316,2395,2382,2300,2361,2287,2184,2169,2199,2154,2222,2243,2166,2108,2024,2115,2008,2002,2006,2005,1959,1930,1895,1962,1891,1887,1903,1790,1809,1822,1808,1791,1738,1801,1816,1722,1634,1675,1755,1667,1639,1658,1574,1649,1545,1550,1536,1468,1546,1453,1493,1433,1520,1464,1426,1435,1422,1447,1328,1394,1430,1389,1293,1268,1312,1261,1265,1214,1288,1263,1188,1213,1201,1147,1150,1161,1179,1164,1148,1143,1145,1072,1038,1002,1064,1060,1077,990,1037,986,1057,988,1006,967,990,1031,940,968,924,914,929,881,921,851,920,922,870,887,790,839,835,799,770,819,792,828,762,792,825,749,723,754,741,750,746,718,685,690,696,656,671,674,686,707,600,642,643,605,641,591,607,635,588,579,592,580,595,587,550,586,561,537,569,559,532,537,530,534,546,502,477,505,488,490,514,463,467,449,448,458,464,454,413,440,474,416,410,424,405,407,418,389,397,412,378,400,368,399,358,405,372,388,360,383,358,371,367,345,324,351,325,307,330,331,360,317,317,313,348,301,312,287,277,333,311,263,283,292,272,282,271,283,252,259,250,239,251,267,271,257,295,272,226,244,250,259,222,228,243,221,226,231,225,233,205,254,221,235,194,212,214,219,200,205,189,195,188,202,189,172,174,172,170,165,166,158,164,186,157,154,170,182,171,180,161,131,153,154,143,147,158,160,136,134,142,133,132,130,142,118,120,130,122,119,115,133,115,141,132,129,118,110,114,99,107,92,122,108,126,96,101,104,92,88,102,91,115,100,97,114,88,127,84,83,106,92,72,85,88,90,86,86,78,91,83,72,68,96,86,89,85,73,70,85,73,65,73,57,67,69,74,69,69,62,71,72,85,76,74,66,55,67,52,63,59,59,58,56,55,60,50,68,54,48,59,53,60,65,59,42,49,61,51,55,46,46,39,49,38,48,38,42,37,52,42,52,41,43,47,42,37,42,36,30,41,40,29,39,33,34,37,40,27,33,32,27,37,29,34,28,33,24,25,33,39,32,23,21,28,23,27,33,24,25,24,24,36,35,27,
25,25,31,20,30,26,27,34,14,26,17,23,24,31,22,20,23,21,15,31,13,17,25,21,23,23,11,20,11,27,19,17,19,20,16,14,15,15,5,10,15,9,10,11,13,10,6,5,12,4,4,6,18,2,5,2,3,3,6,2,4,0,2,1,2,1,2,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

  var cd = frag_sizes.cumulative


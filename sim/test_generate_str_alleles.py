import sys
#sys.path.append("..")
from generate_str_alleles import *
import pytest

@pytest.mark.parametrize("ref_sequence, repeatunit, delta, expected", [
    ('ATATAT', 'AT', 1, 'ATATATAT'),
    ('ATATAT', 'AT', -2, 'AT'),
    ('CGATATATCG', 'AT', 1, 'CGATATATATCG'),
    ('CGATATATCG', 'AT', -1, 'CGATATCG'),
])
def test_mutate_str(ref_sequence, repeatunit, delta, expected):
    assert mutate_str(ref_sequence, repeatunit, delta, random=False) == expected

@pytest.mark.parametrize("ref, variant, start, stop, expected", [
    ('CGATATATCG', 'AT', 2, 8, 'CGATCG'),
    ('012345678', 'AT', 3, 6, '012AT678'),
    ('012345678', 'AT', 2, None, '01AT2345678'),
])
def test_replace_variant(ref, variant, start, stop, expected):
    assert replace_variant(ref, variant, start, stop) == expected

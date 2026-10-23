from pandasgwas import get_traits, bind, intersect, set_diff, set_xor, union, set_equal


def test_bind():
    one = get_traits(efo_id='EFO_0009713')
    two = get_traits(efo_id='OBA_VT0001321')
    three = get_traits(efo_id='EFO_0920093')
    one_and_two = bind(one, two)
    two_and_three = bind(two, three)
    assert len(one_and_two) == 2
    assert len(two_and_three) == 2


def test_intersect():
    one = get_traits(efo_id='EFO_0009713')
    two = get_traits(efo_id='OBA_VT0001321')
    three = get_traits(efo_id='EFO_0920093')
    one_and_two = bind(one, two)
    two_and_three = bind(two, three)
    result = intersect(one_and_two, two_and_three)
    assert len(result) == 1


def test_set_diff():
    one = get_traits(efo_id='EFO_0009713')
    two = get_traits(efo_id='OBA_VT0001321')
    three = get_traits(efo_id='EFO_0920093')
    one_and_two = bind(one, two)
    two_and_three = bind(two, three)
    result = set_diff(one_and_two, two_and_three)
    assert len(result) == 1


def test_set_xor():
    one = get_traits(efo_id='EFO_0009713')
    two = get_traits(efo_id='OBA_VT0001321')
    three = get_traits(efo_id='EFO_0920093')
    one_and_two = bind(one, two)
    two_and_three = bind(two, three)
    result = set_xor(one_and_two, two_and_three)
    assert len(result) == 2


def test_union():
    one = get_traits(efo_id='EFO_0009713')
    two = get_traits(efo_id='OBA_VT0001321')
    three = get_traits(efo_id='EFO_0920093')
    one_and_two = bind(one, two)
    two_and_three = bind(two, three)
    result = union(one_and_two, two_and_three)
    assert len(result) == 3


def test_set_equal():
    df1 = get_traits(efo_id='EFO_0009713')
    df2 = get_traits(efo_id='EFO_0009713')
    assert set_equal(df1, df2)

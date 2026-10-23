from pandasgwas import is_API_available

def test_is_api_available():
    is_available = is_API_available()
    assert is_available is True

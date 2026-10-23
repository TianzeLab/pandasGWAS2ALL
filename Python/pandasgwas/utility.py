from typing import List

from pandasgwas import client
from requests.adapters import HTTPAdapter
import requests

def is_API_available() -> bool:
    """
    Test if the API is available. If the return value is not True, it means that the API is not available.

    ```Python
    from pandasgwas import is_API_available
    is_available = is_API_available()
    ```

    Returns:
        True or False
    """
    try:
        with requests.Session() as s:
            s.mount('https://', HTTPAdapter(max_retries=5))
            r = s.get("https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits/")
            if r.status_code == 200:
                return True
            return False
    except:
        return False

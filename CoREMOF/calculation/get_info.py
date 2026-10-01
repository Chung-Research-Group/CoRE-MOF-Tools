"""Get publication information from DOI.
"""

from datetime import date as calendar_date
import math
import time
from urllib.parse import quote


def _publication_date(data):
    """Prefer a valid online date, then print, without inventing precision."""
    if not isinstance(data, dict):
        return "unknown"
    for name in ("published-online", "published-print"):
        value = data.get(name)
        if not isinstance(value, dict):
            continue
        parts = value.get("date-parts")
        if not isinstance(parts, list) or len(parts) != 1:
            continue
        parts = parts[0]
        if (not isinstance(parts, list) or not 1 <= len(parts) <= 3
                or any(type(part) is not int for part in parts)):
            continue
        try:
            calendar_date(parts[0], parts[1] if len(parts) > 1 else 1, parts[2] if len(parts) > 2 else 1)
        except ValueError:
            continue
        return "-".join([str(parts[0])] + [f"{part:02d}" for part in parts[1:]])
    # Licence start, metadata creation, deposit and indexing dates are not
    # publication dates. Never use them to assign a release naming year.
    return "unknown"


def get_publication_date(doi, max_retries=3, delay=5):
    """Get Crossref's online publication date, otherwise its print date.

    Args:
        doi (str): DOI identifier, without a resolver URL.
        max_retries (int): maximum request attempts, including the first.
        delay (float): seconds between retryable transport failures.
       
    Returns:
        str:
            -   YYYY, YYYY-MM, YYYY-MM-DD, or ``unknown``. Missing date
                components are not filled. A licence start date is never
                substituted for publication. This lookup does not modify
                existing metadata or release IDs.
    """
    if type(max_retries) is not int or max_retries < 1:
        raise ValueError("max_retries must be a positive integer")
    if isinstance(delay, bool) or not isinstance(delay, (int, float)) or not math.isfinite(delay) or delay < 0:
        raise ValueError("delay must be a finite nonnegative number")
    if not isinstance(doi, str) or not doi.strip():
        return "unknown"
    import requests
    url = "https://api.crossref.org/works/" + quote(doi.strip(), safe="/")
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                payload = response.json()
                return _publication_date(payload.get("message") if isinstance(payload, dict) else None)
            else:
                return "unknown"
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            if attempt + 1 < max_retries:
                time.sleep(delay)
        except Exception as e:
            print(f"[Error] DOI: {doi} failed due to {e}")
            return "unknown"
    return "unknown"


def extract_publication(doi):
    """Get publisher.

    Args:
        doi (str): DOI.
       
    Returns:
        str:
            -   Publisher
    """
    try:
        doi_part1 = doi.split("/")[0]
        part_1 = ["10.1021","10.1039","10.1002","10.1038","10.1126","10.1016","10.1007","10.3390"]
        part_name = ["ACS","RSC","WILEY","Nature","Science","SciDirect","Springer","MDPI"]
        try:
            index = part_1.index(doi_part1)
            return part_name[index]
        except:
            if doi_part1 == "unknown":
                return "unknown"
            else:
                return "other"
    except:
        return "unknown"

"""Get publication information from DOI.
"""

import requests, time


def get_publication_date(doi, max_retries=3, delay=5):
    """Get publicated date.

    Args:
        doi (str): DOI.
       
    Returns:
        str:
            -   time of DOI.
    """
        
    url = f"https://api.crossref.org/works/{doi}"
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                data = response.json().get("message", {})
                if 'published-online' in data:
                    date = data['published-online']['date-parts'][0]
                elif 'published-print' in data:
                    date = data['published-print']['date-parts'][0]
                elif 'license' in data:
                    date = data['license'][0]['start']['date-parts'][0]
                else:
                    return "unknown"
                if len(date) == 3:
                    return f"{date[0]}-{date[1]:02d}-{date[2]:02d}"
                elif len(date) == 2:
                    return f"{date[0]}-{date[1]:02d}"
                else:
                    return str(date[0])
            else:
                return "unknown"
        except requests.exceptions.ReadTimeout:
            print(f"[Timeout] DOI: {doi}, retrying ({attempt + 1}/{max_retries})...")
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
import requests

def get_publication_date(doi):
    url = f"https://api.crossref.org/works/{doi}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'license' in data['message']:
            pub_date = data['message']['license'][0]['start']['date-parts'][0]
        elif 'published-online' in data['message']:
            pub_date = data['message']['published-online']['date-parts'][0]
        elif 'published-print' in data['message']:
            pub_date = data['message']['published-print']['date-parts'][0]
        else:
            pub_date = 'Unknown'
    else:
        pub_date = "Unknown"
    time = pub_date
    if len(time)==3:
        return str(time[0])+"-"+str(time[1])+"-"+str(time[2])
    elif len(time)==2:
        return str(time[0])+"-"+str(time[1])
    else:
        return time
    
def extract_publication(doi):
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
            return "Other"
import requests
import os
import asf_search as asf
import functionality

#check user credentials 
def check_login(username, password):
    #create temporary session for credential check
    try:
        session = asf.ASFSession().auth_with_creds(username, password)
        response = session.get(
            "https://urs.earthdata.nasa.gov/profile",
            timeout=30
        )

        if response.status_code == 200:
            return True

        print("Login failed with status:", response.status_code)
        return False
    
    except Exception as e:
        print("Login failed:", e)   
        return False

#build URL to be send 
def search_asf(params):
    #base URL where we can add parameters
    base_URL = "https://api.daac.asf.alaska.edu/services/search/param"

    #handeling network errors
    try:
        return requests.get(base_URL, params=params, timeout=30)
    except requests.exceptions.RequestException as e:
        print("Request failed:", e)
        return None



#using the urls from the first API call to download the files
def download_asf(information, username, password):
    #create a session with credentials for download
    session = asf.ASFSession().auth_with_creds(username, password)

    #list error, same as in controllor for compatability check 
    if isinstance(information, dict):
        information = [information]


    #path to download folder
    download_folder = os.path.join(
        os.path.dirname(__file__), "..", "PHASE_Preprocessing", "slaves"
        )

    #create folder if it not exist
    os.makedirs(download_folder, exist_ok=True)


    #delting old files, but cheking if we are trying do download the same picture again to then not delte it
    #finding incoming pictures
    expected_files = {
        item["sceneName"] + ".zip"
        for item in information
    }

    # deletion pending
    datanot = "pending"
    functionality.write_json("deletion.json", datanot)

    #removing old files
    for filename in os.listdir(download_folder):
        if filename.endswith(".zip") and filename not in expected_files:
            os.remove(os.path.join(download_folder, filename))
    
    #deletion complete
    dataok = "complete"
    functionality.write_json("deletion.json", dataok)

    #retriving the remaning download URls
    urls = []
    for item in information:
        filename = item["sceneName"] + ".zip"
        filepath = os.path.join(download_folder, filename)

        if not os.path.exists(filepath):
            url = item.get("url")
            if url:
                urls.append(url)


    #handel no url error
    if not urls:
        print("no URLs")
        return
    
    #download using the found url, credentials from session and a path for files to be stored
    asf.download_urls(
        urls=urls,
        path=download_folder,
        session=session

    )

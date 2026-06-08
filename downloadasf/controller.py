import API
import sys
import functionality

#check login for user
def run_login():
    #get username and password
    userInfo = functionality.read_json("login_request.json")

    username = userInfo.get("username")
    password = userInfo.get("password")

    #check login and save result for matlab
    login_result = {}

    if API.check_login(username, password):
        login_result["status"] = "success"
        print("login ok")
    else:
        login_result["status"] = "failure"
        print("login faild, password or username wrong")
        
    functionality.write_json("login_result.json", login_result)


#run the first API call
def run_search():
    #read data from matlab, and build API call
    request_data = functionality.read_json("search_request.json")
    params = functionality.build_search_parameters(request_data)
    response = API.search_asf(params)

    #handel errors
    if response is None:
        print("request faild")
        return
    
    data = response.json()
    print("Number of raw ASF features:", len(data["features"]))

    #handle sampling if user want
    sampling = request_data.get("sampling")
    if sampling:
        data["features"] = functionality.sampling_rate(
            data["features"],
            sampling.get("rate"),
            sampling.get("unit")
        )

    #get all relevant info that is needed for matlab
    info = functionality.relevant_info(data)
    
    #retive the best footprint, to get best path, frame and orbit direction
    aoi_lon, aoi_lat = functionality.find_aoi_center(request_data["aoi"])
    best_product = functionality.best_footprint(info["information"], aoi_lon, aoi_lat)

    #save all data needed for download
    download_info = {
        "information": [
            {
                "sceneName": product["sceneName"],
                "url": product["url"]
            }
            for product in info["information"]
        ]
    }

    functionality.write_json("download_data.json", download_info)

    #write summary for matlab, containing size, footprint, path and frame 
    products = []

    for product in info["information"]:
        products.append({
            "sceneName": product["sceneName"],
            "size": round(product["size"], 2),
            "size_bytes": product["size_bytes"],
            "startTime": product["startTime"],
            "pathNumber": product["pathNumber"],
            "frameNumber": product["frameNumber"],
            "flightDirection": product["flightDirection"],
            "footprint": product["footprint"]
        })

    summary = {
        "product_count": info["product_count"],
        "total_size_gb": round(info["total_size_gb"], 2),
        "total_size_bytes": info["total_size_bytes"],
        "products": products
    }

    #only add best result of there exist best results
    if best_product:
        summary.update({
            "best_path": best_product["pathNumber"],
            "best_frame": best_product["frameNumber"],
            "best_direction": best_product["flightDirection"]
        })
    else:
        summary.update({
            "best_path": None,
            "best_frame": None,
            "best_direction": None,
            "warning" : "no products footptint contains the AOI center"
        })


    functionality.write_json("search_summary.json", summary)
    
#if download confirmed, start downloading data
def run_download():
    print("run_download started")
    #get username and password
    userInfo = functionality.read_json("login_request.json")
    username = userInfo.get("username")
    password = userInfo.get("password")

    #get data from json file
    info = functionality.read_json("download_data.json")

    summary = functionality.read_json("search_summary.json")

    products = summary.get("products", [])

    #handle not a list error if only one product is selected
    if isinstance(products, dict):
        products = [products]


    #check if there are any products for download
    if not products:
        functionality.write_json("search_summary.json", {
            "status": "invalid",
            "message": "No products available for download."
        })
        return
    
    #check if download request is valid, if not send error back to matlab
    compatibility = functionality.check_compatibility(products)
    if not compatibility["valid"]:
        summary = {
            "status": "invalid",
            "message": "Products are not compatible for PHASE.",
            "compatibility": compatibility
        }

        functionality.write_json("search_summary.json", summary)

        print("Products are not compatible for PHASE.")
        return

    #if all checks passed, start download
    API.download_asf(info["information"], username, password)


if __name__ == "__main__":
    #make shure mode has a value
    if len(sys.argv) < 2:
        print("Missing mode. Use: search, download, or login.")
        sys.exit(1)
    
    mode = sys.argv[1]


    #sepreate if we want to run search, download or username/password check 
    if mode == "search":
        run_search()

    elif mode == "download":
        run_download()

    elif mode == "login":
        run_login()

    #incase wrong mode
    else:
        print(f"unknown mode{mode}")
        sys.exit(1)
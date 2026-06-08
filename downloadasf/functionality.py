import json
import os
from datetime import datetime
from collections import defaultdict
from shapely.geometry import Point, Polygon

#read and write to json files
def read_json(filename):
    path = os.path.join(os.path.dirname(__file__), filename)
    
    with open(path, "r") as f:
        return json.load(f)
    
def write_json(filename, data):
    path = os.path.join(os.path.dirname(__file__), filename)

    with open(path, "w") as f:
        json.dump(data, f, indent=4)
    

#setting up AOI for URL
def build_aoi(aoi):
    aoi_type = aoi["type"]
    
    #seperate based on aoi type
    if aoi_type == "point":
        lon, lat = aoi["coordinates"]
        return f"POINT({lon} {lat})"
    
    elif aoi_type == "line":
        coords = aoi["coordinates"]
        line_coords = ", ".join(f"{lon} {lat}" for lon, lat in coords)
        return f"LINESTRING({line_coords})"
    
    elif aoi_type in ["polygon", "rectangle"]:
        #close the polygon
        coords = aoi.get("coordinates", aoi.get("corners"))
        polygon_coords = coords + [coords[0]]
        poly_string = ", ".join(f"{lon} {lat}" for lon, lat in polygon_coords)
        return f"POLYGON(({poly_string}))"
    
    else:
        raise ValueError(f"Unsupported AOI type: {aoi_type}")
    
#handling parameters that can have multiple values
def handle_multi_variable(params, asf_name, value):
    if value:
        params[asf_name] = value

#build search parameter for url search
def build_search_parameters(data):
    #import parameters and harcode sentinel-1
    aoi = data["aoi"]
    start_date = data.get("startDate")
    end_date = data.get("endDate") 
    pathStart = data.get("pathStart")
    pathEnd = data.get("pathEnd")
    frameStart = data.get("frameStart")
    frameEnd = data.get("frameEnd")

    #parameter that always has to be in the search
    params = {
        "dataset": "SENTINEL-1",
        "intersectsWith": build_aoi(aoi),
        "output": "geojson"
    }

    #retrive parameter if they are a part of the search
    if start_date:
        params["start"] = start_date
    if end_date:
        params["end"] = end_date + "T23:59:59Z"
    handle_multi_variable(params, "polarization", data.get("polarization"))
    handle_multi_variable(params, "processingLevel", data.get("processingLevel"))
    handle_multi_variable(params, "beamMode", data.get("beamMode")) 
    handle_multi_variable(params, "flightDirection", data.get("flightDirection"))
    handle_multi_variable(params, "platform", data.get("subtype"))
    handle_multi_variable(params, "groupID", data.get("groupID"))  
     
    if pathStart and pathEnd:
        params["relativeOrbit"] = f"{pathStart}-{pathEnd}"
    elif pathStart:
        params["relativeOrbit"] = pathStart
    elif pathEnd:
        params["relativeOrbit"] = pathEnd
    
    if frameStart and frameEnd:
        params["frame"] = f"{frameStart}-{frameEnd}"
    elif frameStart:
        params["frame"] = frameStart
    elif frameEnd:
        params["frame"] = frameEnd
    return params

#retriving relevant data to show user after first asf search
def relevant_info(data):
    features = data["features"]
    results = []
    total_size = 0
    product_count = 0

    #relevant info for each picture and total size off all pictures
    for feature in features:
        props = feature["properties"]

        url = props.get("url", "")
        if not url.endswith(".zip"):
            continue
        
        size_bytes = props.get("bytes", 0)
        total_size += size_bytes
        startTime = props.get("startTime", "")
        date = startTime[:10]

        results.append({
        "sceneName": props["sceneName"],
        "url": url,
        "size": size_bytes / (1024**3),
        "size_bytes": size_bytes,
        "startTime": date,
        "pathNumber": props["pathNumber"],
        "frameNumber": props["frameNumber"],
        "footprint": feature["geometry"],
        "flightDirection": props["flightDirection"]
    })

    total_size_gb = total_size / (1024**3)
    product_count = len(results)
    return {
        "information": results,
        "total_size_gb": total_size_gb,
        "total_size_bytes": total_size,
        "product_count" : product_count
    }

#Sampling rate if user do not want all of the pictures
def sampling_rate(information, rate, period):
    #chechk if rate is a valid number
    if not rate or rate <= 0 or not period:
        return information
    
    period = str(period).lower()
    
    #group by week/month/year
    groups = defaultdict(list)
    for product in information:
        props = product["properties"]
        
        url = props.get("url", "")
        if not url.endswith(".zip"):
            continue

        date = datetime.fromisoformat(
            props["startTime"].replace("Z", "+00:00")
        )

        #grouping key for different periods
        if period == "week":
            year, week, _ = date.isocalendar()
            key = f"{year}-W{week}"
        elif period == "month":
            key = date.strftime("%Y-%m")
        elif period == "year":
            key = date.strftime("%Y") 
        else:
            raise ValueError(f"Unsuported sampling period: {period}")  

        #add product to the correct group
        groups[key].append(product)

    sampled = []
    #sample evenly inside of each group
    for key in sorted(groups.keys()):
        #sort based on time, oldest to newest
        products = sorted(
            groups[key],
            key=lambda product: product["properties"]["startTime"]
        )

        #if fewer products then requestet, return all products
        if len(products) <= rate:
            sampled.extend(products)
            continue

        #calculate even space between products
        step = len(products) / rate

        for i in range(rate):
            index = int(i * step)
            sampled.append(products[index])
    return sampled


#cheking if paths and orbits are the same
def check_compatibility(information):
    #retrieving information from search summary 
    paths = {product["pathNumber"] for product in information}
    frames = {product["frameNumber"] for product in information}
    directions = {product["flightDirection"] for product in information}

    #check if they are the same
    results = {
        "same_path": len(paths) == 1,
        "same_frame": len(frames) == 1,
        "same_direction": len(directions) == 1,
        "paths": list(paths),
        "frames": list(frames),
        "directions": list(directions)
    }

    results["valid"] = (
        results["same_path"]
        and results["same_frame"]
        and results["same_direction"]
    )

    return results

#finding the best producrt for user
def best_footprint(products, aoi_lon, aoi_lat):
    #use center of aoi to find best product
    aoi_center = Point(aoi_lon, aoi_lat)

    best_product = None
    best_score = -1

    #chechk every product
    for product in products:
        coords = product["footprint"]["coordinates"][0]
        polygon = Polygon(coords)

        if not polygon.contains(aoi_center):
            continue
        
        #finding distance, and updating new best product if it is better
        score = aoi_center.distance(polygon.boundary)

        if score > best_score:
            best_score = score
            best_product = product
    
    return best_product

#help function for finding the center of a aoi and returning coordinates
def find_aoi_center(aoi):
    if aoi["type"] == "point":
        lon, lat = aoi["coordinates"]
        return lon, lat

    elif aoi["type"] == "line":
        coords = aoi["coordinates"]

    elif aoi["type"] in ["polygon", "rectangle"]:
        coords = aoi.get("coordinates", aoi.get("corners"))

    polygon = Polygon(coords)
    center = polygon.centroid

    return center.x, center.y

import folium
import requests as r
import json
import os
import api_key


def geocoder_conn_str(input_addr, asset_pos, asset_name, facility_name, facility_addr):

    QUERY = '%20'.join(input_addr.split(' '))
    BASE = 'https://us1.locationiq.com/v1/search.php?key='

    conn_str = BASE + os.getenv('API_KEY') + '&q=' + QUERY + '&format=json'

    data = r.get(conn_str).json()

    geo_corr = [float(data[0]['lat']), float(data[0]['lon'])]

    m = folium.Map(
        location=geo_corr,
        zoom_start=16,
        tiles='OpenStreetMap')

    folium.Marker(asset_pos, popup=facility_name, tooltip=asset_name).add_to(m)

    return m.save('./templates/map.html')

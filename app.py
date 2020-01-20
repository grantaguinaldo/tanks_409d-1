import math
import pandas as pd
import numpy as np
from tanks_helper import *
import os
from flask import Flask, render_template, request
from mapr import geocoder_conn_str

app = Flask(__name__)


@app.route('/api/vfrtk', methods=['POST'])
def vfrtk():

    df_chem = pd.read_csv('chemical_db.csv')
    df_met = pd.read_csv('met_db.csv')
    df_shade = pd.read_csv('table-7-1-6-solarabs.csv')

    POST_DATA = request.get_json()

    INPUT_CITY = POST_DATA['input_city']
    INPUT_TANK = POST_DATA['input_tank']
    INPUT_CONTENTS = POST_DATA['input_contents']
    CHEM_LIST = POST_DATA['input_chem']
    ANNUAL_QUANTITY = POST_DATA['input_qty']
    DEFAULT_LIST = POST_DATA['input_default']
    CONDITION_LIST = POST_DATA['input_condition']
    INPUT_TANK_TYPE = POST_DATA['input_tank_type']
    INPUT_TANK_NAME = POST_DATA['input_tank_name']
    INPUT_FACILITY_NAME = POST_DATA['input_facility_name']
    INPUT_FACILITY_ADDR = POST_DATA['input_facility_address']
    INPUT_TANK_GEO = POST_DATA['input_tank_geo']

    geocoder_conn_str(input_addr=INPUT_FACILITY_ADDR,
                      asset_pos=INPUT_TANK_GEO,
                      asset_name=INPUT_TANK_NAME,
                      facility_name=INPUT_FACILITY_NAME,
                      facility_addr=INPUT_FACILITY_ADDR)

    MET_LIST = filterMetList(df=df_met, input_city=INPUT_CITY)

    solarabs = solarabsLookUp(df=df_shade,
                              color=CONDITION_LIST[0],
                              shade=CONDITION_LIST[1],
                              condition=CONDITION_LIST[2])

    tank = FixedRoofTank(tkshellht=INPUT_TANK[0],                 # From User Data
                         skliqht=INPUT_TANK[1],                   # From User Data
                         tkrfslope=DEFAULT_LIST[0],               # Default
                         diameter=INPUT_TANK[2],                  # From User Data
                         ins=MET_LIST[0][1],                      # From Met Table
                         solarabs=solarabs,
                         tax=MET_LIST[0][5],                      # From Met Table
                         tan=MET_LIST[0][6],                      # From Met Table
                         atmplocal=MET_LIST[0][0],                # From Met Table
                         throughput=INPUT_CONTENTS[0],            # From User Data
                         productfactor=INPUT_CONTENTS[1],         # From User Data
                         hlx=INPUT_CONTENTS[2],                   # From User Data
                         hln=INPUT_CONTENTS[3],                   # From User Data
                         ventsetting=DEFAULT_LIST[1],
                         tanktype=INPUT_TANK_TYPE)                # Default

    calculation(df=df_chem,
                chem_list=CHEM_LIST,
                annual_qty=ANNUAL_QUANTITY,
                tank=tank,
                file_name='vert_fixed_roof_tk.html',
                tank_name=INPUT_TANK_NAME,
                tank_type=INPUT_TANK_TYPE,
                facility_name=INPUT_FACILITY_NAME)

    return render_template('vert_fixed_roof_tk.html')


@app.route('/is_alive')
def is_alive():
    return 'This is Tanks 4.09_d on the Web'


@app.route('/emissiondetail')
def emissiondetail():
    return render_template('vert_fixed_roof_tk.html')


@app.route('/summary')
def summary():
    return render_template('loss_summary.html')


@app.route('/map')
def map():
    return render_template('map.html')


@app.route('/')
def index():
    return render_template('index.html')


if __name__ == "__main__":
    app.run(debug=True)

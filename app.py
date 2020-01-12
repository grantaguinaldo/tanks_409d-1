import math
import pandas as pd
import numpy as np
from tanks_helper import *
import os
from flask import Flask, render_template, request

app = Flask(__name__)

# TODO: Need to know how to pass the payload from a post request into this route.
# https://stackoverflow.com/questions/20001229/how-to-get-posted-json-in-flask


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

    ######################################################################
    # Needs to come into the app as a POST request payload.
    #
    # INPUT_CITY = 'Denver, Colorado'                     # From User Data
    # INPUT_TANK = [12, 8, 6]                             # tkshellht, skliqht, diameter
    # INPUT_CONTENTS = [8450, 'other stocks', 11.5, 4.5]  # throughput, productfactor, hlx, hln
    # CHEM_LIST = ['Cyclohexane', 'Benzene', 'Toluene']   # From User Data
    # ANNUAL_QUANTITY = [101, 2812, 258]                  # From User Data
    # DEFAULT_LIST = [0.0625, 1491, 1]                    # tkrfslope, ins, ventsetting
    # CONDITION_LIST = ['White', 'None', 'Aged']          # color, shade, condition
    ######################################################################

    MET_LIST = filterMetList(df=df_met, input_city=INPUT_CITY)

    solarabs = solarabsLookUp(df=df_shade,
                              color=CONDITION_LIST[0],
                              shade=CONDITION_LIST[1],
                              condition=CONDITION_LIST[2])

    tank = VerticalFixedRoofTank(tkshellht=INPUT_TANK[0],                 # From User Data
                                 skliqht=INPUT_TANK[1],                   # From User Data
                                 tkrfslope=DEFAULT_LIST[0],               # Default
                                 diameter=INPUT_TANK[2],                  # From User Data
                                 ins=DEFAULT_LIST[1],                     # Default
                                 solarabs=solarabs,
                                 tax=MET_LIST[0][5],                      # From Met Table
                                 tan=MET_LIST[0][6],                      # From Met Table
                                 atmplocal=MET_LIST[0][0],                # From Met Table
                                 throughput=INPUT_CONTENTS[0],            # From User Data
                                 productfactor=INPUT_CONTENTS[1],         # From User Data
                                 hlx=INPUT_CONTENTS[2],                   # From User Data
                                 hln=INPUT_CONTENTS[3],                   # From User Data
                                 ventsetting=DEFAULT_LIST[2])             # Default

    calculation(df=df_chem,
                chem_list=CHEM_LIST,
                annual_qty=ANNUAL_QUANTITY,
                tank=tank,
                file_name='vert_fixed_roof_tk.html')

    print('At the end of the script!')

    return render_template('vert_fixed_roof_tk.html')


@app.route('/is_alive')
def alive():
    return 'This App Is Alive!'


@app.route('/')
def index():
    return 'This is Tanks 4.09_d on the Web'


@app.route('/results')
def results():
    return render_template('vert_fixed_roof_tk.html')


if __name__ == "__main__":
    app.run(debug=True)

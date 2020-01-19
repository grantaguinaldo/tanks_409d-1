import requests as r

POST_JSON = {
    'input_city': 'Denver, Colorado',                     # From User Data
    'input_tank': [12, 8, 6],                             # tkshellht, skliqht, diameter
    'input_contents': [8450, 'other stocks', 11.5, 4.5],  # throughput, productfactor, hlx, hln
    'input_chem': ['Cyclohexane', 'Benzene', 'Toluene'],  # From User Data
    'input_qty': [101, 2812, 258],                        # From User Data
    'input_default': [0.0625, 1491, 1],                   # tkrfslope, ins, ventsetting
    'input_condition': ['White', 'None', 'Average'],      # color, shade, condition
    'input_tank_type': 'Vertical'                         # From User Data
}

r.post('http://127.0.0.1:5000/api/vfrtk', json=POST_JSON)

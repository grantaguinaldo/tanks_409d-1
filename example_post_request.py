import requests as r

POST_JSON = {
    'input_city': 'Cape Hatteras, North Carolina',
    'input_tank': [12, 8, 6],
    'input_contents': [8450, 'other stocks', 11.5, 4.5],
    'input_chem': ['Cyclohexane', 'Benzene', 'Toluene'],
    'input_qty': [101, 2812, 258],
    'input_default': [0.0625, 1491, 1],
    'input_condition': ['White', 'None', 'Average'],
    'input_tank_type': 'Horizontal',
    'input_tank_name': 'Test Tank Name',
    'input_facility_name': 'Test Facility Name'
}

r.post('http://127.0.0.1:5000/api/vfrtk', json=POST_JSON)

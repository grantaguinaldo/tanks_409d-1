from flask import Flask, render_template, request, json

app = Flask(__name__)


@route('/api', methods=['POST'])
def check():
    #data = request.get_json()
    data = request.json

    return jsonify(data)


if __name__ == '__main__':
    app.run(debug=True)

'''
import requests

requests.post('http://127.0.0.1:5000/api', json={
    'input_city_': 'Denver, Colorado',
    'input_tank_': [12, 8, 6],
    'input_contents': [8450, 'other stocks', 11.5, 4.5],
    'input_chem': ['Cyclohexane', 'Benzene', 'Toluene'],
    'input_qty': [101, 2812, 258],
    'input_default': [0.0625, 1491, 1],
    'input_condition': ['White', 'None', 'Aged']})
'''

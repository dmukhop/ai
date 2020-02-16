from flask import Blueprint, jsonify, current_app

import json
import os

bp = Blueprint("drugs", __name__)

drug_response = []

APP_ROOT = os.path.dirname(os.path.abspath(__file__))

with open(APP_ROOT + '/data/drug_response.json') as json_file:
    drug_response = json.load(json_file)



@bp.route("/drugs/<int:id>/responses")
def get_drug_responses(id):

    results = filter(lambda p: p['DRUG_ID'] == id,drug_response)

    return jsonify(list(results))


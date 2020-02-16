from flask import Blueprint, jsonify, current_app, request

import json
import os
from pyhive import presto  # or import hive

bp = Blueprint("patients", __name__)

patients = []
patient_results = []

APP_ROOT = os.path.dirname(os.path.abspath(__file__))

with open(APP_ROOT + '/data/patients.json') as json_file:
    patients = json.load(json_file)

with open(APP_ROOT + '/data/patient_results.json') as json_file:
    patient_results = json.load(json_file)


@bp.route("/patients")
def index():
    return jsonify(patients)



@bp.route("/patients/<string:id>/condition/<string:condition>/results")
def get_patient_results(id, condition):
    print("#############",id,condition)
    cursor = presto.connect(host = 'presto', port='8080').cursor()
    sql_txt = ' SELECT BINARY_RESPONSE, CONDITION, DRUG_ID, DRUG_NAME, LN_IC50, MODEL, PATHWAY, PATIENT_ID, THRESHOLD FROM default.recommendations where PATIENT_ID ='+"'"+id+"' and CONDITION = "+ "'"+ condition + "'"
    print(sql_txt)
    cursor.execute(sql_txt)
    patient_results = []
    while True:
        a = cursor.fetchone()
        if a == None :
            break
        s3_1 = []
        #print("S3_1:",s3_1)
        for c, b in zip(a, [ 'BINARY_RESPONSE', 'CONDITION', 'DRUG_ID', 'DRUG_NAME', 'LN_IC50', 'MODEL', 'PATHWAY', 'PATIENT_ID', 'THRESHOLD']):
            s3_1.append((b,c))
        patient_results.append( dict((s3_1)) )
    cursor.close()
    if len(patient_results) == 0:
        return jsonify([]), 404
    return(jsonify(patient_results) )

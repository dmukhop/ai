#!/usr/bin/env python
from flask import Blueprint, jsonify, current_app, request
import json
import os
from pyhive import presto  # or import hive

bp = Blueprint("query", __name__)



APP_ROOT = os.path.dirname(os.path.abspath(__file__))


@bp.route("/patients/<string:id>/results")
def predict_patient_condition(id, condition):
    process_patient_string = {'event_type': 'Process_patient', "COSMIC_ID": id, "Condition" : condition }
    cursor = presto.connect(host = 'presto', port='8080').cursor()
    sql_txt = ' SELECT BINARY_RESPONSE, CONDITION, DRUG_ID, DRUG_NAME, LN_IC50, MODEL, PATHWAY, PATIENT_ID, THRESHOLD FROM default.recommendations where PATIENT_ID ='+"'"+id+"'"
    #print(sql_txt)
    cursor.execute(sql_txt)
    patient_results = []
    while True:
        a = cursor.fetchone()
        if a == None :
            break
        #s += '{'
        s3_1 = []
        #print("S3_1:",s3_1)
        for c, b in zip(a, [ 'BINARY_RESPONSE', 'CONDITION', 'DRUG_ID', 'DRUG_NAME', 'LN_IC50', 'MODEL', 'PATHWAY', 'PATIENT_ID', 'THRESHOLD']):
            #s += '"'+b+'": "'+ c + '",'
            #s2.append('"'+b+'": "'+ c + '"')
            s3_1.append((b,c))
        #s += '}'
        patient_results.append( dict((s3_1)) )
        #s[len(s)-1]= '}'
    cursor.close()
    if len(patient_results) == 0:
        return jsonify([]), 404
    return(jsonify(patient_results) )

    #return jsonify("Patient :"+id+" Processed")



@bp.route("/query/<string:id>/condition/<string:condition>/pathway/<string:pathway>")
def predict_patient_condition_pathway(id, condition, pathway):

    print("ID - {0}, CONDITION - {1}, PATHWAY - {2}".format(id, condition, pathway))

    return jsonify("OK")





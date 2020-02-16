#!/usr/bin/env python
from flask import Blueprint, jsonify, current_app, request
from kafka import KafkaProducer
import json
import os
from pyhive import presto  # or import hive
import time

bp = Blueprint("predict", __name__)
producer = KafkaProducer(bootstrap_servers='kafka:29092')

patients = []
patient_results = []
cursor = presto.connect(host = 'presto', port='8080').cursor()


APP_ROOT = os.path.dirname(os.path.abspath(__file__))

def log_to_kafka(topic, event):
    event.update(request.headers)
    producer.send(topic, json.dumps(event).encode())


@bp.route("/predict/<string:id>/condition/<string:condition>")
def predict_patient_condition(id, condition):
    process_patient_string = {'event_type': 'Process_patient', "COSMIC_ID": id, "Condition" : condition }
    log_to_kafka('process_patient', process_patient_string)

    #processing = True
    #while processing:
    #    cursor.execute('SELECT * FROM default.recommendations where COSMIC_ID = "'+ id + '"')
    #    if cursor.fetchone() != None :
    #        processing = False
    #    time.sleep(5)
    #cursor.fetchall()

    print("ID - {0}, CONDITION - {1}".format(id, condition))

    return jsonify("Patient :"+id+" Processed")



@bp.route("/predict/<string:id>/condition/<string:condition>/pathway/<string:pathway>")
def predict_patient_condition_pathway(id, condition, pathway):
    process_patient_string = {'event_type': 'Process_patient', "COSMIC_ID": id, "Condition" : condition }
    log_to_kafka('process_patient', process_patient_string)

    #processing = True
    #while processing:
    #    cursor.execute('SELECT * FROM default.recommendations where COSMIC_ID = "'+ id + '"')
    #    if cursor.fetchone() != None :
    #        processing = False
    #    time.sleep(5)
    #cursor.fetchall()

    print("ID - {0}, CONDITION - {1}".format(id, condition))

    return jsonify("Patient :"+id+" Processed")


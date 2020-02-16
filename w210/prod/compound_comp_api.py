#!/usr/bin/env python
import json,jsonify
from kafka import KafkaProducer
from flask import Flask, request
from flask_cors import CORS
import time


app = Flask(__name__)
CORS(app)
producer = KafkaProducer(bootstrap_servers='kafka:29092')


def log_to_kafka(topic, event):
    event.update(request.headers)
    producer.send(topic, json.dumps(event).encode())


@app.route("/")
def default_response():
    default_event = {'event_type': 'default'}
    log_to_kafka('default', default_event)
    return "Welcome to Compound Companion API Gateway. \nFollwing APIs are available: \n\t1. Process Patient \n\t2. Query Patient\n"


@app.route("/Process_patient")
def process_patient():
    payload = request.args.get("patient")
    process_patient_string = {'event_type': 'Process_patient', "payload": payload}
    log_to_kafka('process_patient', process_patient_string)
    return "Patient Processed" 

@app.route("/Query_Patient")
def query_patient():
    payload = request.args.get("patient")
    process_query_string = {'event_type': 'Query_patient', "payload": payload}
    query_string = {'event_type': 'Query Patient'}
    log_to_kafka('query_patient', query_string)
    time.sleep(10)
    ## get the query response from kafka here ##
    return jsonify(data = "Here the drug recommendation!")

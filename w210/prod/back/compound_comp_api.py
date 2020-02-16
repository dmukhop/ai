#!/usr/bin/env python
import json,jsonify
from kafka import KafkaProducer
from flask import Flask, request
from flask_cors import CORS

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


@app.route("/Process_WES")
def process_WES():
    payload = request.args.get("WES")
    process_WES_string = {'event_type': 'Process_WES', "payload": payload}
    log_to_kafka('WES', process_WES_string)
    return "WES Processed" 

@app.route("/Process_CNA")
def process_CNA():
    process_CNA_string = {'event_type': 'Process_CNA'}
    log_to_kafka('CNA', process_CNA_string)
    return jsonify(data = "CNA data process successfully!") 

@app.route("/Query_Patient")
def query_patient():
    query_string = {'event_type': 'Query Patient'}
    log_to_kafka('query', query_string)
    return jsonify(data = "Here the drug recommendation!")

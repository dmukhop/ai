#!/usr/bin/env python
from flask import Blueprint, jsonify, current_app
from kafka import KafkaProducer
import json
import os

from pyhive import presto  # or import hive
cursor = presto.connect(host = 'presto', port='8080').cursor()
sql_txt = 'SELECT BINARY_RESPONSE, CONDITION, DRUG_ID, DRUG_NAME, LN_IC50, MODEL, PATHWAY, PATIENT_ID, THRESHOLD FROM default.recommendations where PATIENT_ID ='+"'910937' and CONDITION="+ "'SCLC'"
print(sql_txt)
s = "["
s2 = []
s3 = []
#cursor.execute("select * from recommendations where cosmic_id = '910937'")
cursor.execute(sql_txt)
while True:
    a = cursor.fetchone()
    if a == None :
        break
    s += '{'
    s3_1 = []
    for c, b in zip(a, ['BINARY_RESPONSE', 'CONDITION', 'DRUG_ID', 'DRUG_NAME', 'LN_IC50', 'MODEL', 'PATHWAY', 'PATIENT_ID', 'THRESHOLD' ]):
        s += '"'+b+'": "'+ c + '",'
        s2.append('"'+b+'": "'+ c + '"')
        s3_1.append((b,c))
    s += '}'
    s3.append( dict((s3_1)) )
    #s[len(s)-1]= '}'
cursor.close()
s += ']'
print(s)
print(json.dumps(s2))
print(s3)
#print(jsonify(data = json.dumps(s2)) )
#print(type(a),cursor.columns, a)
#print json.dumps(cursor.fetchone())

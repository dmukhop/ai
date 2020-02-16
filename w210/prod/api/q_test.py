#!/usr/bin/env python
import json
import os
from pyhive import presto  # or import hive






def predict_patient_condition(id, condition):
    cursor = presto.connect(host = 'presto', port='8080').cursor()
    process_patient_string = {'event_type': 'Process_patient', "COSMIC_ID": id, "Condition" : condition }
    sql_txt = 'SELECT cosmic_id, drug_id, ic50_threshold, predicted_ic50, predicted_resp  FROM default.recommendations where COSMIC_ID ='+"'"+id+"'" 
    #sql_txt += ' where COSMIC_ID ='+id
    print(sql_txt)
    cursor.execute(sql_txt)
    print(1)
    patient_results = []
    while True:
        print(2)
        a = cursor.fetchone()
        print(3)
        if a == None :
            break
        #s += '{'
        s3_1 = []
        #print("S3_1:",s3_1)
        for c, b in zip(a, ["cosmic_id", "drug_id", "ic50_threshold", "predicted_ic50", "predicted_resp" ]):
            #s += '"'+b+'": "'+ c + '",'
            #s2.append('"'+b+'": "'+ c + '"')
            s3_1.append((b,c))
        #s += '}'
        patient_results.append( dict((s3_1)) )
        #s[len(s)-1]= '}'
    cursor.close()
    #return(jsonify(patient_results) )

predict_patient_condition('1299052','Bla')




#!/bin/bash

rundir=/w210/prod


# Start up docker containers - This step will never be called from this script. The docker-compose file will execute this script as the containers are created
docker-compose up -d

# To check Zookeeper logs 
# docker-compose logs zookeeper | grep -i binding

# Quick look at Kafla logs
# docker-compose logs kafka | grep -i started

# Let's now create Kafka topics
docker-compose exec kafka kafka-topics --create --topic WES --partitions 1 --replication-factor 1 --if-not-exists --zookeeper zookeeper:32181

# docker-compose exec kafka kafka-topics --create --topic CNA --partitions 1 --replication-factor 1 --if-not-exists --zookeeper zookeeper:32181

# Sample publisher code using bash
# docker-compose exec mids bash -c “cat /w205/github-example-large.json | jq ‘.’ -c | kafkacat -P -b kafka:29092 -t foo && echo ‘Produced 100 messages.’”


# If we want to see the console using Kafkacat then use following in a new terminal
# docker-compose exec mids kafkacat -C -b kafka:29092 -t WES -o beginning

# Sample Apache Bench stress test command
# docker-compose exec mids ab -n 10 -H "Host: user1.comcast.com" http://localhost:5000/

# Now we invoke WES_feature_write.py. It will run continuously and as we specified in the code below, 
# every 10 seconds it will read kafka, process our data, and write it to the console.

docker-compose exec spark spark-submit $rundir/write_wes_stream.py

# Create Hadoop Directories

docker-compose exec cloudera hdfs dfs -mkdir /user/root/raw
docker-compose exec cloudera hdfs dfs -mkdir /user/root/processed


# Create hive table
# create external table if not exists default.raw_wes (Accept string, Host string, User_Agent string, event_type string, payload string, timestamp string) stored as parquet location '/tmp/wes_process'  tblproperties ("parquet.compress"="SNAPPY");

#Install jsonify
docker-compose exec mids pip install jsonify

#Install flask-cors
docker-compose exec mids pip install -U flask-cors

# Start the API Server - ****** WE NEED TO EXPOSE A PORT *****


docker-compose exec mids env FLASK_APP=$rundir/compound_comp_api.py flask run --host 0.0.0.0


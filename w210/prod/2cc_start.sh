rundir=/w210/prod

# Create Hadoop Directories

#docker-compose exec cloudera hdfs dfs -mkdir /user/root/raw
#docker-compose exec cloudera hdfs dfs -mkdir /user/root/processed


# Create hive table
# create external table if not exists default.raw_wes (Accept string, Host string, User_Agent string, event_type string, payload string, timestamp string) stored as parquet location '/tmp/wes_process'  tblproperties ("parquet.compress"="SNAPPY");
,/hadoop_cp.sh
#Install jsonify
docker-compose exec mids pip install jsonify
#Hive and presto

docker-compose exec mids pip install pyhive
docker-compose exec mids pip install pyhive[presto]
docker-compose exec mids pip install pyspark

#Install flask-cors
docker-compose exec mids pip install -U flask-cors

# Start the API Server - ****** WE NEED TO EXPOSE A PORT *****


docker-compose exec mids env FLASK_APP=$rundir/api flask run --host 0.0.0.0


docker-compose exec cloudera hdfs dfs -rm /user/root/processed/*
docker-compose exec cloudera hdfs dfs -copyFromLocal /w210/ucb-w210-cap-drug-response/data/complete_score.csv /user/root/processed
docker-compose exec cloudera hdfs dfs -copyFromLocal /w210/ucb-w210-cap-drug-response/data/complete_score2.csv /user/root/processed
docker-compose exec cloudera hive -f /w210/prod/create_tables.sql

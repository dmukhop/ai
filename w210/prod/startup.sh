docker-compose exec mids env FLASK_APP=~/w210/prod/compound_comp_api.py flask run --host 0.0.0.0
docker-compose exec spark spark-submit ~/w210/prod/filtered_writes.py

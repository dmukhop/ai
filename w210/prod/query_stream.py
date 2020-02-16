#
## Streaming

##

#!/usr/bin/env python
"""Extract request from kafka and writes back
"""
import json
from pyspark.sql import SparkSession
from pyspark.sql.functions import udf, from_json
from pyspark.sql.types import StructType, StructField, StringType


def query_patient_event_schema():
    """
    root
    |-- Accept: string (nullable = true)
    |-- Host: string (nullable = true)
    |-- User-Agent: string (nullable = true)
    |-- event_type: string (nullable = true)
    |-- payload : string (nullable = true)
    """
    return StructType([
        StructField("Accept", StringType(), True),
        StructField("Host", StringType(), True),
        StructField("User-Agent", StringType(), True),
        StructField("event_type", StringType(), True),
        StructField("payload", StringType(), True),
    ])


@udf('boolean')
def is_query_process(event_as_json):
    """udf for filtering events
    """
    event = json.loads(event_as_json)
    if event['event_type'] == 'Query_patient':
        return True
    return False


def main():
    """main
    """
    spark = SparkSession \
        .builder \
        .appName("QueryEventsJob") \
        .getOrCreate()

    raw_events = spark \
        .readStream \
        .format("kafka") \
        .option("kafka.bootstrap.servers", "kafka:29092") \
        .option("subscribe", "query_patient") \
        .load()

    patient_raw_data = raw_events \
        .filter(is_wes_process(raw_events.value.cast('string'))) \
        .select(raw_events.value.cast('string').alias('raw_event'),
                raw_events.timestamp.cast('string'),
                from_json(raw_events.value.cast('string'),
                          process_patient_event_schema()).alias('json')) \
        .select('raw_event', 'timestamp', 'json.*')

    #convert json to CSV here 

    # Spark SQL here

    # write back to kafka
    sink = patient_raw_data \
        .writeStream \
        .format("csv") \
        .option("checkpointLocation", "/tmp/checkpoints_for_wes_process") \
        .option("path", "/user/root/raw") \
        .trigger(processingTime="5 seconds") \
        .start()

    sink.awaitTermination()


if __name__ == "__main__":
    main()



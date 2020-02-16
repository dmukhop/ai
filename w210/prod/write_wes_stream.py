#
## Streaming

##

#!/usr/bin/env python
"""Extract events from kafka and write them to hdfs
"""
import json
from pyspark.sql import SparkSession
from pyspark.sql.functions import udf, from_json
from pyspark.sql.types import StructType, StructField, StringType


def process_patient_event_schema():
    """
    root
    |-- Accept: string (nullable = true)
    |-- Host: string (nullable = true)
    |-- User-Agent: string (nullable = true)
    |-- event_type: string (nullable = true)
    |-- COSMIC_ID : string (nullable = true)
    |-- condition : string (nullable = true)
    """
    return StructType([
        StructField("Accept", StringType(), True),
        StructField("Host", StringType(), True),
        StructField("User-Agent", StringType(), True),
        StructField("event_type", StringType(), True),
        StructField("COSMIC_ID", StringType(), True),
        StructField("Condition", StringType(), True)
    ])


@udf('boolean')
def is_patient_process(event_as_json):
    """udf for filtering events
    """
    event = json.loads(event_as_json)
    if event['event_type'] == 'Process_patient':
        return True
    return False


def main():
    """main
    """
    print('********************* STARTING SPAK SESSION***********')
    spark = SparkSession \
        .builder \
        .appName("ExtractEventsJob") \
        .getOrCreate()

    raw_events = spark \
        .readStream \
        .format("kafka") \
        .option("kafka.bootstrap.servers", "kafka:29092") \
        .option("subscribe", "process_patient") \
        .load()
    print('********************* RAW EVENTS')

    patient_raw_data = raw_events \
        .filter(is_patient_process(raw_events.value.cast('string'))) \
        .select(raw_events.value.cast('string').alias('raw_event'),
                raw_events.timestamp.cast('string'),
                from_json(raw_events.value.cast('string'),
                          process_patient_event_schema()).alias('json')) \
        .select('raw_event', 'timestamp', 'json.*')

    #convert json to CSV here 
    print('*********************', patient_raw_data)
    sink = patient_raw_data \
        .writeStream \
        .format("json") \
        .option("checkpointLocation", "/tmp/checkpoints_for_patient_process") \
        .option("path", "/user/root/raw") \
        .trigger(processingTime="60 seconds") \
        .start()
    print('********************* WROTE RAW DATA IN SPARK')

    sink.awaitTermination()


if __name__ == "__main__":
    main()



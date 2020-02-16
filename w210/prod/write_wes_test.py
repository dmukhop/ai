#!/usr/bin/env python
"""Extract events from kafka and write them to hdfs
"""
import json
from pyspark.sql import SparkSession, Row
from pyspark.sql.functions import udf


@udf('boolean')
def is_wes_process(event_as_json):
    event = json.loads(event_as_json)
    if event['event_type'] == 'Process_WES':
        return True
    return False


def main():
    """main
    """
    spark = SparkSession \
        .builder \
        .appName("ExtractEventsJob") \
        .getOrCreate()

    raw_events = spark \
        .readStream \
        .format("kafka") \
        .option("kafka.bootstrap.servers", "kafka:29092") \
        .option("subscribe", "WES") \
        .load()

    wes_process = raw_events \
        .select(raw_events.value.cast('string').alias('raw'),
                raw_events.timestamp.cast('string')) \
        .filter(is_wes_process('raw'))


    #extracted_wes_process = wes_process \
    #    .rdd \
    #    .map(lambda r: Row(timestamp=r.timestamp, **json.loads(r.raw))) \
    #    .toDF()
    wes_process.printSchema()
    wes_process.show()

    sink = wes_process \
        .writeStream \
        .format("parquet") \
        .option("checkpointLocation", "/tmp/checkpoints_for_wes_process") \
        .option("path", "/tmp/wes_process") \
        .trigger(processingTime="120 seconds") \
        .start()

    sink.awaitTermination()


if __name__ == "__main__":
    main()

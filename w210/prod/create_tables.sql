drop table raw_input;
create external table raw_input ( value string ) location '/user/root/raw/';

drop view clean_input;
create view clean_input as select get_json_object(value,'$.COSMIC_ID') as cosmic_id, get_json_object(value,'$.Condition') as condition, get_json_object(value,'$.timestamp') as timestamp, current_timestamp as curr_time from raw_input;

drop table recommendations;
create external table if not exists default.recommendations 
(	COSMIC_DRUG_ID string
	,LN_IC50 string
	,DRUG_ID string
	,PATIENT_ID string
	,DRUG_NAME string
	,THRESHOLD string
	,BINARY_RESPONSE string
	,CONDITION string
	,PATHWAY string
	,MODEL string)
ROW FORMAT DELIMITED FIELDS TERMINATED BY ',' STORED AS TEXTFILE LOCATION '/user/root/processed';

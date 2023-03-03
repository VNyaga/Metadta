options nonumber nodate ls=90;

%include "\\sciensano.be\fs\1120_Cancer_Employee\BMH\Victoria\SAS\Metadta\Scripts\Metadas_v1.3.sas";

*=====================================================================;
*Example One;

proc import out=telomerase replace
	datafile='\\sciensano.be\fs\1120_Cancer_Employee\BMH\Victoria\SAS\Metadta\Data\telomerase.csv'
	dbms=csv;
	getnames = yes;
 run;

data telomerase;
set telomerase;
rename  study = study_id;
run;

%metadas(import=n, dsname=telomerase, tp=TP, fp=FP, fn=FN, tn=TN,
method=b, keepds=none);

*=====================================================================;
*Example Two;
proc import out=dementia replace
	datafile="\\sciensano.be\fs\1120_Cancer_Employee\BMH\Victoria\SAS\Metadta\Data\dementia.csv"
	dbms=csv;
	getnames = yes;
 run;

data dementia;
set dementia;
rename  study = study_id;
run;

 %metadas(import=n, dsname=dementia, tp=TP, fp=FP, fn=FN, tn=TN,
method=b, keepds=none);

*=====================================================================;
*Example Three;
proc import out=ascus replace
	datafile="\\sciensano.be\fs\1120_Cancer_Employee\BMH\Victoria\SAS\Metadta\Data\ascus.csv"
	dbms=csv;
	getnames = yes;
 run;

 data ascus;
set ascus;
rename  studyid = study_id;
run;

 %metadas(import=n, dsname=ascus, tp=TP, fp=FP, fn=FN, tn=TN,
method=b, covariate=test, cvref="RepC", cvtype=cat,  keepds=none);

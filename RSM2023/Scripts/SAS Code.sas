options nonumber nodate ls=90;

%include "\Metadta\Scripts\Metadas_v1.3.sas";
*NOTE: Change the directory in the -option datafile- while reading in the data;
*=====================================================================;
*Section 5.1 ;

proc import out=telomerase replace
	datafile='\Metadta\Data\telomerase.csv'
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
*Section 5.2;
proc import out=dementia replace
	datafile="\Metadta\Data\dementia.csv"
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
*Section 5.3;
proc import out=ascus replace
	datafile="\Metadta\Data\ascus.csv"
	dbms=csv;
	getnames = yes;
 run;

 data ascus;
set ascus;
rename  studyid = study_id;
run;

 %metadas(import=n, dsname=ascus, tp=TP, fp=FP, fn=FN, tn=TN,
method=b, covariate=test, cvref="RepC", cvtype=cat,  keepds=none);

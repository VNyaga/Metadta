cap program drop metadta_examples
program metadta_examples
	version 14.0
	`1'
end

program define example_one
	preserve
	di " "
	use "http://fmwww.bc.edu/repec/bocode/t/telomerase.dta", clear
	di  " "
	di `". metadta tp fp fn tn,  	///"'  
	di `"{phang}studyid(study) sumtable(all) smooth	///{p_end}"' 
	di `"{phang}soptions(xtitle("False positive rate")) 	///{p_end}"'
	di `"{phang}foptions(texts(2) xlabel(0, 0.5, 1)	///{p_end}"'
	di `"{pmore} ysize(10) xsize(20) astext(60)) {p_end}"' 
	
	set more off

	
	#delimit ;
	metadta tp fp fn tn, 
		studyid(study) sumtable(all) smooth 
		soptions(xtitle("False positive rate")) 
		foptions(texts(2) xlabel(0, 0.5, 1) ysize(10) xsize(20) astext(60)) 
	;
	#delimit cr
	restore
end


program define example_two_one
	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/a/ascus.dta", clear
	di _n
	
	di `". metadta tp fp fn tn test, 			///{p_end}"'
	di `"{phang}studyid(studyid)  comparative sumtable(none) smooth gof			///{p_end}"'
	di `"{phang}soptions(xtitle("False positive rate") col(red blue)) 			///{p_end}"'	
	di `"{phang}foptions(outplot(abs) texts(2) xlabel(0, 0.5, 1) ysize(10) xsize(20) astext(60)) {p_end}"'
	
	
	set more off
	
	#delimit ;
	metadta tp fp fn tn test, 
		studyid(studyid)  comparative sumtable(none) smooth gof
		soptions(xtitle("False positive rate") col(red blue)) 	
		foptions(texts(2) xlabel(0, 0.5, 1) ysize(10) xsize(20) astext(60)) 
	;
	#delimit cr
	restore
end

program define example_two_two
	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/a/ascus.dta", clear
	di _n

	di `". metadta tp fp fn tn test, 	///{p_end}"' 
	di `"{phangstudyid(studyid) comparative sumtable(all)	smooth gof	///{p_end}"' 
	di `"{phangfoptions(logscale outplot(rr) texts(2) xlabel(0.5, 1, 2)  ysize(10) xsize(20) astext(60)) 	{p_end}"'

	set more off
	
	#delimit ;
	metadta tp fp fn tn test, 
		studyid(studyid) comparative sumtable(all)	smooth gof
		foptions(logscale outplot(rr) texts(2) xlabel(0.5, 1, 2)  ysize(10) xsize(20) astext(60)) 
	;
	#delimit cr
	restore
end

program define example_two_three
	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/a/ascus.dta", clear
	di _n

	di `".metadta tp fp fn tn test, 	///{p_end}"' 
	di `"{phang}studyid(studyid) comparative model(random, laplace) cov(unstructured, independent) sumtable(all) smooth gof 	///{p_end}"' 
	di `"{phang}foptions(logscale outplot(rr) texts(2) xlabel(0.5, 1, 2)  ysize(10) xsize(20) astext(60)) 	{p_end}"'

	set more off
	
	#delimit ;
	metadta tp fp fn tn test, 
		studyid(studyid) comparative model(random, laplace) cov(unstructured, independent) sumtable(all) smooth gof
		foptions(logscale outplot(rr) texts(2) xlabel(0.5, 1, 2)  ysize(10) xsize(20) astext(60)) 
	;
	#delimit cr
	restore
end

program define example_three_one
	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/c/clinself.dta", clear
	di _n

	di `". metadta tp fp fn tn sample Setting, nomc gof 	///{p_end}"'
	di `"{phang}studyid(study) interaction(sesp) 	///{p_end}"'
	di `"{phang}summaryonly comparative sumtable(rr) noitable 	///{p_end}"'
	di `"{phang}foptions(outplot(rr) texts(2) xlabel(0.6, 1) astext(60)) 	{p_end}"'
		
	set more off
		
	#delimit ;
	metadta tp fp fn tn sample Setting, nomc gof 
		studyid(study) interaction(sesp) 
		summaryonly comparative sumtable(rr) noitable 
		foptions(outplot(rr) texts(2) xlabel(0.6, 1) astext(60)) 	
	;
	#delimit cr
	restore
end

program define example_three_two
	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/c/clinself.dta", clear
	di _n
	
	di `". metadta tp fp fn tn sample Setting TA, nomc smooth gof	///{p_end}"'
	di `"{phang}studyid(study) interaction(sesp) 	///{p_end}"'
	di `"{phang}model(random, laplace) comparative noitable sumtable(rr)	///{p_end}"'
	di `"{phang}foptions(logscale outplot(rr) grid  texts(1.5) xlabel(0.5, 1, 2) ysize(10) xsize(20) astext(60)) 	{p_end}"'
	
	set more off
	#delimit ;
	metadta tp fp fn tn sample Setting TA, nomc smooth gof
	studyid(study) interaction(sesp) 
		model(random, laplace) comparative noitable sumtable(rr)
		foptions(logscale outplot(rr) grid  texts(1.5) xlabel(0.5, 1, 2) ysize(10) xsize(20) astext(60)) 
	;
	#delimit cr
	restore
end

program define example_three_three
	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/c/clinself.dta", clear
	di _n
	
	di `". metadta tp fp fn tn sample, nomc smooth gof progress		///{p_end}"'
	di `"{phang}studyid(study) cov(unstructured, independent) 		///{p_end}"'
	di `"{phang}comparative noitable sumtable(rr)		///{p_end}"'
	di `"{phang}foptions(logscale outplot(rr) grid  texts(1.5) xlabel(0.5, 1, 2) ysize(10) xsize(20) astext(60)) 		{p_end}"'
	
	set more off
	#delimit ;
	metadta tp fp fn tn sample, nomc smooth gof progress
	studyid(study) cov(unstructured, independent) 
		comparative noitable sumtable(rr)
		foptions(logscale outplot(rr) grid  texts(1.5) xlabel(0.5, 1, 2) ysize(10) xsize(20) astext(60)) 
	;
	#delimit cr
	restore
end

program define example_four
	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/p/pairedta.dta", clear
	di _n


	di `". metadta tp1 fp1 fn1 tn1 tp2 fp2 fn2 tn2 hpv1 hpv2, smooth gof  progress ///{p_end}"' 
	di `"{phang}studyid(study) model(random) cov(, independent) ///{p_end}"' 
	di `"{phang}cbnetwork sumtable(rr)  ///{p_end}"' 
	di `"{phang}foptions(outplot(rr) grid  texts(1.85) ///{p_end}"' 
	di `"{pmore}xlabel(0.75, 0.90, 1, 1.11, 1.33) logscale lcols(hpv2 setting)  astext(70) ysize(10) xsize(20))   {p_end}"'
 
	// metadta tp1---fn2 index comparator, cbnetwork
	set more off
	
	 metadta tp1 fp1 fn1 tn1 tp2 fp2 fn2 tn2 hpv1 hpv2, smooth gof  progress ///
		studyid(study) model(random) cov(, independent) ///
		cbnetwork sumtable(rr)  ///
		foptions(outplot(rr) grid  texts(1.85) ///
			xlabel(0.75, 0.90, 1, 1.11, 1.33) logscale lcols(hpv2 setting)  astext(70) ysize(10) xsize(20))    
			
	restore
end

program define example_five

	preserve
	di _n
	use "http://fmwww.bc.edu/repec/bocode/n/network.dta", clear
	di _n
	
	di `". metadta  tp fp fn tn test ,  smooth gof  ///{p_end}"'
	di `"{phang}studyid(study)  ///{p_end}"'
	di `"{phang}abnetwork ref(HC2) sumtable(all) ///{p_end}"'
	di `"{phang}foptions(outplot(rr) texts(1.75) ///{p_end}"'
	di `"{pmore}xlabel(0.80, 0.90, 1, 1.11, 2) logscale astext(70) ysize(10) xsize(20)) {p_end}"'


//ab network meta-analysis
	set more off
	metadta  tp fp fn tn test ,  smooth gof  /// 
		studyid(study)  ///
		abnetwork ref(HC2) sumtable(all) ///
		foptions(outplot(rr) texts(1.75) ///
		xlabel(0.80, 0.90, 1, 1.11, 2) logscale astext(70) ysize(10) xsize(20))
	restore
end


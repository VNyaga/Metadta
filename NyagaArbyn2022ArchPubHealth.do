cd "C:\DATA\WIV\Projects\Stata\Metadta\Graphs\" //Change this

set more off


use "http://fmwww.bc.edu/repec/bocode/t/telomerase.dta", clear 
list, noobs clean
*EXAMPLE 1
*1 - Figure 1
set more off

metadta tp fp fn tn,     							///
  studyid(study) model(random) dp(2) sumtable(all) 		///
  soptions(xtitle("False positive rate") ciopt(lpattern(dash_dot))				/// 
	xlabel(0(0.2)1) xscale(range(0 1)) 				///
     	ytitle("Sensitivity") yscale(range(0 1))			/// 
	ylabel(0(0.2)1, nogrid) 						///
     	graphregion(color(white)) plotregion(margin(medium)) 	/// 			
		xsize(15) ysize(15) 							///
     	legend(order(1 "Summary" 5 "Observed data" 2 "SROC" 3 ///				
		"Confidence region" 4 "Prediction region") 		///
            cols(1) ring(0) bplacement(6)) graphsave("s1.gph")) 			///
    foptions(xline(1) graphregion(color(white)) texts(2.5) 		///
		xlabel(0, 0.5, 1) diamopt(color(red))  		/// 				
		olineopt(color(red) lpattern(dash)))
		
estimates restore metadta_modest
estimates store unstructured		
 
*2
set more off
metadta tp fp fn tn,     							///
  studyid(study) model(random) cov(independe) dp(2) sumtable(all) 		///
  soptions(xtitle("False positive rate") ciopt(lpattern(dash_dot)) 				/// 
	xlabel(0(0.2)1) xscale(range(0 1)) 				///
     	ytitle("Sensitivity") yscale(range(0 1))			/// 
	ylabel(0(0.2)1, nogrid) 						///
     	graphregion(color(white)) plotregion(margin(medium)) 	/// 			
		xsize(15) ysize(15) 							///
     	legend(order(1 "Summary" 5 "Observed data" 2 "SROC" 3 ///				
		"Confidence region" 4 "Prediction region") 		///
            cols(1) ring(0) bplacement(6)) graphsave("s2.gph")) 			///
    foptions(graphregion(color(white)) texts(2.5) 		///
		xlabel(0, 0.5, 1) diamopt(color(red))  		/// 				
		olineopt(color(red) lpattern(dash)))
		
estimates restore metadta_modest
estimates store reduced	

lrtest unstructured reduced, stats

//Figure 2
gr combine s1.gph s2.gph, cols(2) graphregion(color(white)) xsize(8) ysize(4) iscale(1)
	
*Example 2	 
use "https://github.com/VNyaga/Metadta/blob/master/clinselfdemo.dta?raw=true", clear
list in 1/10, noobs clean


*1. 
gsort sample

set more off

metadta tp fp fn tn sample if ta=="SA" & setting=="follow-up",                   		///
    studyid(study) sortby(year study) nooverall		///
    noitable sumtable(abs rr) nosroc 		///
    foptions( subline pointopt(msize(1.5)) diamopt(color(red)) olineopt(color(red) 		///
		lpattern(dash)) outplot(abs) fysize(35)						///
		graphregion(color(white)) texts(2) tit("SA, follow-up, by sample, CIN2+") ///
		graphsave("f1.gph")		///
		xlabel(0, .2, .4, .6, .8, 1) arrowopt(msize(1)))
	


*2. 
metadta tp fp fn tn sample if ta=="SA" & setting=="screening" ,                   		///
    studyid(study) sortby(year study) nooverall		///
    noitable sumtable(abs rr) nosroc nomc		///
    foptions(subline pointopt(msize(1.5)) diamopt(color(red)) olineopt(color(red) 		///
		lpattern(dash)) outplot(abs) 				///
		graphregion(color(white)) texts(2) tit("SA, screening, by sample, CIN2+") graphsave("f2.gph")		///
		xlabel(0, .2, .4, .6, .8, 1) arrowopt(msize(1)))
		
gr combine f1.gph f2.gph, graphregion(color(white)) 

*graph save sa.gph, replace
*3. 
gsort sample
metadta tp fp fn tn sample if ta=="TA" & setting=="follow-up",                  		///
    studyid(study) sortby(year study) nooverall	nomc	///
    noitable sumtable(abs rr) nosroc 		///
    foptions(subline pointopt(msize(1.5)) diamopt(color(red)) olineopt(color(red) 		///
		lpattern(dash)) outplot(abs) 						///
		graphregion(color(white)) texts(2) tit("TA, follow-up, by sample, CIN2+") graphsave("f3.gph")			///
		xlabel(0, .2, .4, .6, .8, 1) arrowopt(msize(1)))

		
*4. 
metadta tp fp fn tn sample if ta=="TA" & setting=="screening",                   		///
    studyid(study) sortby(year study) nooverall nomc		///
    noitable sumtable(abs rr) nosroc 		///
    foptions(subline pointopt(msize(1.5)) diamopt(color(red)) olineopt(color(red) 		///
		lpattern(dash)) outplot(abs) fysize(35)						///
		graphregion(color(white)) texts(2) tit("TA, screening, by sample, CIN2+") graphsave("f4.gph")			///
		xlabel(0, .2, .4, .6, .8, 1) arrowopt(msize(1)))

gr combine f1.gph f3.gph, graphregion(color(white)) cols(1)  imargin(0 0 0 0)
graph save 13.gph, replace

gr combine f2.gph f4.gph, graphregion(color(white)) cols(1)  imargin(0 0 0 0)	
graph save 24.gph, replace

//Figure 3
gr combine 13.gph 24.gph, cols(2) graphregion(color(white)) xsize(8) ysize(4) iscale(1)

*RELATIVE ACCURACY SELF VS CLIN BY SETTING WITH SA 
set more off

gsort ta -setting study sample
metadta tp fp fn tn sample setting if ta=="SA",                   		///
    studyid(study)  nomc			///
    comparative noitable sumtable(rr) 						///
    foptions(pointopt(msize(1.5)) diamopt(color(red)) olineopt(color(red) 			///
		lpattern(dash)) outplot(rr) tit("SA, self- vs clin-samples, by setting, CIN2+")	graphsave("f5.gph")	///
		graphregion(color(white)) texts(1.75) 				///
		xlabel(0.3, 0.5, 1, 2, 3) logscale arrowopt(msize(1)))
		
*RELATIVE ACCURACY SELF VS CLIN BY SETTING WITH TA 
sort ta setting study sample
metadta tp fp fn tn sample setting if ta=="TA",                   		///
    studyid(study)  nomc			///
    comparative noitable sumtable(rr) 						///
    foptions(pointopt(msize(1.5)) diamopt(color(red)) olineopt(color(red) 			///
		lpattern(dash)) outplot(rr) tit("TA, self- vs clin-samples, by setting, CIN2+") graphsave("f6.gph")		///
		graphregion(color(white)) texts(1.75) 				///
		xlabel(0.3, 0.5, 1, 2, 3) logscale arrowopt(msize(1)))

//Figure 4
gr combine f5.gph f6.gph, col(1) graphregion(color(white))

*RELATIVE ACCURACY SELF VS CLIN interaction(sesp)
set more off
set trace off 
sort ta setting study sample
metadta tp fp fn tn sample ta setting,                   		///
    studyid(study) interaction(sesp) nomc 			///
    comparative noitable sumtable(rr) nooverall						///
    foptions(subline pointopt(msize(1.3)) diamopt(color(red)) olineopt(color(red) 			///
		lpattern(dash)) outplot(rr) tit("self- vs clin-samples - interaction(sesp)")				///
		graphregion(color(white)) texts(1.7) lcol(setting) graphsave("f7.gph")				///
		xlabel(0.3, 0.5, 1, 2, 3) logscale arrowopt(msize(1))) 
		
estimates restore metadta_modest
estimates store full


*RELATIVE ACCURACY SELF VS CLIN interaction(sesp)		
set more off
sort ta setting study sample
metadta tp fp fn tn sample ta setting,                   		///
    studyid(study) interaction(se) 			///
    comparative noitable sumtable(rr) nooverall 						///
    foptions(subline pointopt(msize(1.3)) diamopt(color(red)) olineopt(color(red) 			///
		lpattern(dash)) outplot(rr) tit("self- vs clin-samples - interaction(se)")					///
		graphregion(color(white)) texts(1.7) lcol(setting)	graphsave("f8.gph")				///
		xlabel(0.3, 0.5, 1, 2, 3) logscale arrowopt(msize(1))) 
//Figure 5
gr combine f7.gph f8.gph, row(1) graphregion(color(white)) xsize(8)
		
estimates restore metadta_modest
estimates store reduced

lrtest full reduced, stats

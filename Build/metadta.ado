/*
CREATED:	8 Sep 2017
AUTHOR:		Victoria N Nyaga
PURPOSE: 	To fit a bivariate random-effects model to diagnostic data and 
			produce a series of graphs(sroc and forestsplots).
VERSION: 	3.0.0
NOTES
1. Variable names should not contain underscore(_)
2. Data should be sorted and no duplicates
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UPDATES
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATE:						DETAILS:
24.08.2020
							grid: Grid lines between studies
							noveral in help file changed to noOVerall
							Print full matrix of rr when data is not repeated
							Print # of studies
							Correct the tick & axis position for sp
							Correct computation of the I2
							graphsave(filename) option included
03.09.2020					Change paired to comparative
14.09.2020					Correct way of counting the distinct groups in the meta-analysis
							Check to ensure variable names do no contain underscore.
15.02.2021					paired data: tp1 fp1 .. fn2 tn2 comparator index covariates, by(byvar)
							comparator, index, byvar need to be string
							Need to test more with more covariates!!!
							vline > xline
29.03.2021					order tp fp fn tn
07.05.2021					Reduce the blank lines in the output
11.05.2021					fixed dp = 4 for p-value, fixed dp=0 for df
14.06.2021					Graph save issues
17.06.2021					Subline fix
15.10.2021					Subgroup analysis with superimposed graphs; stratify option
09.06.2022					Network meta-analysis; 
26.07.2022					REF(label, top|bottom) ; default is bottom
01.08.2022					Renamed network to abnetwork; and paired to cbnetwork
10.10.2022					introduce variance only on se or sp
19.10.2022					Cutoff optimization
15.02.2023					rename independent to general
27.11.2023					Include outplot(OR)	
							gof - option to display goodness of fit
							Simulate posterior distributions
							smooth:Option to generate smooth estimates
12.01.2023					Correction on computation of the sroc predicticn region	
19.01.2023					option optimize - report conditional/exact if simulated summary is not sensible						
							
FUTURE 						Work on the absolutes for the paired analysis; 
							if version 16 or later; use melogit instead of meqrlogit

*/



/*++++++++++++++++++++++	METADTA +++++++++++++++++++++++++++++++++++++++++++
						WRAPPER FUNCTION
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop metadta
program define metadta, eclass sortpreserve byable(recall)
version 14.0

	#delimit ;
	syntax varlist(min=4) [if] [in],  /*tp fp fn tn tp2 fp2 fn2 tn2  */
	STudyid(varname) /*Study Idenfier*/
	[
	STRatify  /*Stratified analysis, requires byvar()*/
	LAbel(string asis) /*namevar=namevar, yearvar=yearvar*/
	DP(integer 2) /*Decimal place*/
	POWer(integer 0) /*Exponentiating power*/  
	MODel(string asis) /*fixed (1 < n < 3) | random */ 
	COV(string asis) /*UNstructured(default)|INdependent | IDentity | EXchangeable, also includes second covariance for the network, comparative analysis*/
	SORTby(varlist) /*order data by varlist. How data appears on the table and forest plot*/
	INteraction(string) /*sesp(default)|se|sp*/ 
	CVeffect(string) /*sesp(default)|se|sp*/ 
	Level(integer 95) /*Significance level*/
	COMParative /*Comparative data or not*/
	SUMtable(string) /*Which summary tables to present:abs|logodds|rr|all*/
	CImethod(string) /*ci method for the study proportions*/
	noFPlot /*No forest plot*/
	noITable /*No study specific summary table*/
	noHTable /*No heterogeneity table*/
	noMC /*No Model comparison - Saves time*/
	PROGress /*See the model fitting*/
	noSRoc /*No SROC*/
	noWT  /*Suppress weights in the Itable & fplot*/
	noOVerall /*Dont report the overall in the Itable & fplot*/ 
	noSUBgroup /*Dont report the subgroup in the Itable & fplot*/ 
	SUMMaryonly /*Present only summary in the Itable, SROC & fplot*/
	DOWNload(string) /*Keep a copy of data used in the plotting*/
	Alphasort /*Sort the categorical variable alphabetically*/
	FOptions(string asis) /*Options specific to the forest plot*/
	SOptions(string asis) /*Options specific to the sroc plot*/
	by(varname)  /*the grouping variable*/
	PAIRed  /*Paired; now called CBnetwork*/
	CBnetwork /*CB network*/
	ABnetwork /*AB network*/
	REF(string asis) /*Reference level in network, comparative analysis */
	TRADEoff	 
	SMooth nsims(integer 800) //max dim in stata ic
	GOF //Goodness of fit
	OPTIMize //Replace population-averaged estimates with Conditional/exact estimates if model has issues e.g complete seperation etc
	] ;
	#delimit cr
	
	preserve
	cap ereturn clear
	marksample touse, strok 
	qui drop if !`touse'

	tempvar rid se sp event total invtotal use id neolabel es lci uci ///
			modeles modellci modeluci grptotal uniq rowid obsid ///
			newobs predevent iw pairid
			
	tempname logodds absout popabsout logoddsi absouti popabsouti rrout rrouti sptestnlrr setestnlrr sptestnlrri setestnlrri ///
		orout orouti sptestnlor setestnlor sptestnlori setestnlori selogodds absoutse popabsoutse selogoddsi ///
		absoutsei popabsoutsei serrout serrouti seorout seorouti popserrout popserrouti popseorout popseorouti ///
		splogodds absoutsp popabsoutsp splogoddsi absoutspi popabsoutspi poprrouti poprrout poporouti poporout ///
		sprrout sprrouti sporout sporouti popsprrout popsprrouti popsporout popsporouti  ///
		coefmat coefvar BVar BVari WVar WVari Esigma omat isq2 isq2i Isq Isq2 Isq2i absexact  absexactse absexactsp ///
		bghet refe bgheti refei lrtestp V Vi dftestnl ptestnl semc semci spmc spmci samtrix semcrow spmcrow ///
		serow sprow  pserow psprow postserow postsprow serrow sprrow pserrow psprrow nltestrr cutmat nltestor matgof matgofi ///
		exactabsoutsei exactabsoutspi exactabsoutse exactabsoutsp exactabsout
	
	if _by() {
		global by_index_ = _byindex()
		if ("`fplot'" == "" | "`sroc'" == "") & "$by_index_" == "1" {
			cap graph drop _all
			global fplotname = 0
			global srocname = 0
		}
	}
	else {
		global by_index_ 
	}
	
	/*Check if variables exist*/
	foreach var of local varlist {
		cap confirm var `var'
		if _rc!=0  {
			di in re "Variable `var' not in the dataset"
			exit _rc
		}
	}
	
	/*Check for se and sp; they are reserved*/
	qui ds
	local vlist = r(varlist)
	foreach v of local vlist {
		if ("`v'" == "se") | ("`v'" == "sp") {
			di as error "se/sp is a reserved variables name; drop or rename se/sp"
			exit _rc
		}
	}

	//define the design of analysis
	if "`paired'" != "" {
		di as res "Use of the option -paired- is deprecated and replaced with -cbnetwork-"
		local cbnetwork "cbnetwork"
	}
	if "`abnetwork'`comparative'`cbnetwork'" != "" {
		cap assert ("`cbnetwork'" != "") + ("`comparative'" != "") + ("`abnetwork'" != "") == 1
		if _rc!=0  {
			di as error "Define 1 option from: `cbnetwork' `comparative' `abnetwork'"
			exit _rc
		}
	}
	
	if "`cbnetwork'" != "" {
		local design = "cbnetwork"
	}
	else if "`comparative'" != "" {
		local design = "comparative"
	}
	else if "`abnetwork'" != "" {
		local design = "abnetwork"
	}
	else if "`abnetwork'`comparative'`cbnetwork'" == "" { 
		local design = "general"  //default design
	}
	
	//General housekeeping
	if 	"`ref'" != "" {
		tokenize "`ref'", parse(",")
		local ref `1'
		local refpos "`3'"
	}
	if "`refpos'" != "" {
		if (strpos("`refpos'", "top") == 1) {
			local refpos "top"
		}
		else if (strpos("`refpos'", "bot") == 0) {
			local refpos "bottom"
		}
		else {
			di as error "Option `refpos' not allowed in ref(`ref', `refpos')"
			exit
		}
	}
	
	//default
	if "`refpos'" == "" {
		local refpos "bottom"
	}
	if 	"`model'" == "" {
		local model random
	}
	else {
		tokenize "`model'", parse(",")
		local model `1'
		local modelopts "`3'"
	}
	if strpos("`model'", "f") == 1 {
		local model "fixed"
	}
	else if strpos("`model'", "r") == 1 {
		local model "random"
	}
	else {
		di as error "Option `model' not allowed in [`model', `modelopts']"
		di as error "Specify either -fixed- or -random-"
		exit
	}
	if "`model'" == "fixed" & strpos("`modelopts'", "ml") != 0 {
		di as error "Option ml not allowed in [`model', `modelopts']"
		exit
	}
	if "`model'" == "fixed" & strpos("`modelopts'", "irls") != 0 {
		di as error "Option irls not allowed in [`model', `modelopts']"
		exit
	}
	qui count
	if `=r(N)' < 2 {
		di as err "Insufficient data to perform meta-analysis"
		exit 
	}
	if `=r(N)' < 3 & "`model'" == "random"  {
		local model fixed //If less than 3 studies, use fixed model
		di as res _n  "Note: Fixed-effects model imposed whenever number of studies is less than 3."
		if "`modelopts'" != "" {
			local modelopts
			di as res _n  "Warning: Model options ignored."
			di as res _n  "Warning: Consider specifying options for the fixed-effects model."
		}
	}
	if "`model'" == "random" {
		if "`cov'" != "" {
			tokenize "`cov'", parse(",")
			if "`1'" != "," {
				local bcov "`1'"
				local wcov "`3'"
			}
			else{
				local bcov 
				local wcov "`2'"
			}
			
			if "`bcov'" != "" {
				if strpos("`bcov'", "un")== 1 {
					local bcov = "unstructured"
				}	
				else if ustrregexm("`bcov'", "ind", 1){
					local bcov = "independent"
				}
				else if ustrregexm("`bcov'", "id", 1){
					local bcov = "identity"
				}
				else if strpos("`bcov'", "ex") == 1 {
					local bcov = "exchangeable"
				}
				else if strpos("`bcov'", "se") == 1 {
					local bcov = "se"
				}
				else if strpos("`bcov'", "sp") == 1 {
					local bcov = "sp"
				}
				else {
					di as error "Allowed covariance structures: se, sp, unstructured, independent, identity, or exchangeable"
					exit
				}
			}
			
			if /*"`abnetwork'`cbnetwork'`comparative'" != "" &*/ "`wcov'" != "" {
				if ustrregexm("`wcov'", "ind", 1){
					local wcov = "independent"
				}
				else if ustrregexm("`wcov'", "id", 1){
					local wcov = "identity"
				}
				else if ustrregexm("`wcov'", "ze", 1){
					local wcov "zero"
				}
				else {
					di as error "Allowed second covariance structures: independent, identity or zero"
					exit
				}
			}
		}
		if "`bcov'" == ""  {
			local bcov = "unstructured"	
		}
		*if "`abnetwork'`cbnetwork'`comparative'" != "" {
			if "`wcov'" == "" & ("`abnetwork'`cbnetwork'" != "") {
				local wcov = "independent"
			}
			if "`wcov'" == "zero" | ("`comparative'" != "" & "`wcov'" == "") {
				local wcov 
			}
		*}
	}
	else {
		local bcov
		local wcov
	}		
	if `level' < 1 {
			local level `level'*100
	}
	if `level'>99 | `level'<10 {
		local level 95
	}

	/*By default the regressor variabels apply to both sensitivity and specificity*/
	if "`cveffect'" == "" {
		local cveffect "sesp"
	}
	else {
		local rc_ = ("`cveffect'"=="sesp") + ("`cveffect'"=="se") + ("`cveffect'"=="sp")
		if `rc_' != 1 {
			di as err "Options cveffect(`cveffect') incorrectly specified"
			di as err "Allowed options: sesp, se sp"
			exit
		}
	}
	if "`interaction'" != "" {
		if ("`interaction'" != "`cveffect'") & ("`cveffect'" != "sesp"){
			di as err "Conflict in cveffect(`cveffect') & interaction(`interaction')"
			exit
		}
	}

	tokenize `varlist'
	if "`design'" != "cbnetwork" {		
		local depvars "`1' `2' `3' `4'" //	tp fp fn tn 
	
	macro shift 4
	}
	else {
		tempvar index byvar assignment idpair
		cap assert "`10'" != ""
		if _rc != 0 {
			di as err "cbnetwork data requires atleast 10 variable"
			exit _rc
		}
		local depvars "`1' `2' `3' `4' `5' `6' `7' `8'" //	tp1 fp1 fn1 tn1 tp2 fp2 fn2 tn2
		local tp = "`1'"
		local fp = "`2'"
		local fn = "`3'"
		local tn = "`4'"
		local Comparator = "`10'"
		local Index = "`9'"
		
		forvalues num = 1/8 {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		cap confirm string variable `9'
		if _rc != 0 {
			di as error "The index variable in cbnetwork analysis should be a string"
			exit, _rc
		}
		cap confirm string variable `10'
		if _rc != 0 {
			di as error "The comparator variable in cbnetwork analysis should be a string"
			exit, _rc
		}
		macro shift 10
		
	}
	local regressors "`*'" 
	gettoken varx confounders : regressors
	local p: word count `regressors'
	*local idpair: word 1 of `regressors'
	
	/*
	if ("`design'" == "general") & ("`stratify'" != "") & ("`by'" != "") & (`p' > 0) {
		di as err "Re-frame your analysis. The options _stratify & by()_ in meta-regression is confusing. Consider using the prefix by: instead"
		exit
	}
	*/
	
	//check no underscore in the variable names
	if strpos("`regressors'", "_") != 0  {
		di as error "Underscore is a reserved character and covariate(s) containing underscore(s) are not allowed"
		di as error "Rename the covariate(s) and remove the underscore(s) character"
		exit	
	}
	
	if `p' < 2 & "`interaction'" !="" {
		di as error "Interactions allowed with atleast 2 covariates"
		exit
	}
	if ("`design'" == "comparative") | ("`design'" == "abnetwork") | ("`wcov'" != "" & "`design'" != "cbnetwork") {
		gettoken first confounders : regressors
		cap assert `p' > 0
		if _rc != 0 {
			di as error "`design' analysis requires at least 1 covariate to be specified"
			exit _rc
		}
		*gettoken varx confounders : regressors
		if "`first'" != "" {
			cap confirm string variable `first'
			if _rc != 0 {
				di as error "The first covariate in `design' analysis should be a string"
				exit, _rc
			}
		}
		local typevarx = "i"
	}
	
	//=======================================================================================================================
	//=======================================================================================================================
	tempfile master
	qui save "`master'"
	
	fplotcheck,`design' `foptions' first(`first') by(`by') //Forest plot advance housekeeping
	local outplot = r(outplot)
	local foptions = r(foptions)
	local lcols = r(lcols)
	if "`lcols'" == " " { //if empty
		local lcols
	}
	
	*declare study labels for display
	if "`label'"!="" {
		tokenize "`label'", parse("=,")
		while "`1'"!="" {
			cap confirm var `3'
			if _rc!=0  {
				di as err "Variable `3' not defined"
				exit
			}
			local `1' "`3'"
			mac shift 4
		}
	}	
	qui {
		*put name/year variables into appropriate macros
		if "`namevar'"!="" {
			local lbnvl : value label `namevar'
			if "`lbnvl'"!=""  {
				quietly decode `namevar', gen(`neolabel')
			}
			else {
				gen str10 `neolabel'=""
				cap confirm string variable `namevar'
				if _rc==0 {
					replace `neolabel'=`namevar'
				}
				else if _rc==7 {
					replace `neolabel'=string(`namevar')
				}
			}
		}
		if "`namevar'"==""  {
			cap confirm numeric variable `studyid'
			if _rc != 0 {
				gen `neolabel' = `studyid'
			}
			if _rc == 0{
				gen `neolabel' = string(`studyid')
			}
		}
		if "`yearvar'"!="" {
			local yearvar "`yearvar'"
			cap confirm string variable `yearvar'
			if _rc==7 {
				local str "string"
			}
			if "`namevar'"=="" {
				replace `neolabel'=`str'(`yearvar')
			}
			else {
				replace `neolabel'=`neolabel'+" ("+`str'(`yearvar')+")"
			}
		}
	}
		
	//Long format
	longsetup `varlist', rid(`rid') se(`se') event(`event') total(`total') `design' rowid(`rowid')  idpair(`idpair') assignment(`assignment') first(`first')
	
	//Gen variables to store weight and model-based estimates
	qui {		
		cap gen _ESAMPLE = 0
		cap drop _WT
		gen _WT = .
		cap gen `modeles' = .
		cap gen `modellci' = .
		cap gen `modeluci' = .
	}
	
	//Index
	if "`design'" == "cbnetwork" {
		tempvar ipair
		qui gen `ipair' = "Yes"
		qui replace `ipair' = "No" if `idpair'
	}

	//byvar
	if "`by'" != "" {
		cap confirm string variable `by'
		if _rc != 0 {
			di as error "The by() variable should be a string"
			exit, _rc
		}
		local found 0
		foreach v of local varlist {
			if "`v'" == "`by'" {
				local found = 1
				continue, break
			}
		}
		if !`found' {
		/*if strpos(`"`varlist'"', "`by'") == 0*/ 
			tempvar byvar
			my_ncod `byvar', oldvar(`by')
			drop `by'
			rename `byvar' `by'
		}
	}
	if ("`design'" == "comparative" | "`design'" == "abnetwork" ) & ("`outplot'" == "abs") {
		cap assert ("`first'" != "`by'") 
		if _rc != 0 { 
				di as error "Remove the option by(`by') or specify a different by-variable"
				exit _rc
		}
	}
		
	buildregexpr `varlist', cveffect(`cveffect') interaction(`interaction') se(`se') sp(`sp') `alphasort' `design'  ipair(`ipair') baselevel(`ref')
	
	local regexpression = r(regexpression)
	local seregexpression = r(seregexpression)
	local spregexpression = r(spregexpression)
	local catreg = r(catreg)
	local contreg = r(contreg)
	local basecode = r(basecode)
	
	if "`interaction'" != "" { 
		local varx = r(varx)
		local typevarx = r(typevarx)		
	}
	if "`design'" == "comparative" {
		*local varx : word 1 of `regressors'
		gettoken varx catreg : catreg
		local typevarx = "i"
		local baselab:label `varx' `basecode'
		if `basecode' == 1 {
			local indexcode "2"
		}
		else {
			local indexcode "1"
		}
		local indexlab:label `varx' `indexcode'
		
		if "`outplot'" != "abs" {
			local varxlabs "`varx' `indexlab' `baselab'"
		}
	}
	local pcat: word count `catreg'
	
	/*if "`cbnetwork'" != "" {
		local varx = "`index'"
		local typevarx  = "i"
	*/
	
	
	/*if ("`cbnetwork'" == "") & ("`comparative'" == "") & ("`interaction'" == "")  {
		local varx 
		local typevarx  
	}*/
	if "`design'" == "cbnetwork" { 
		local varx = "`ipair'"
		local typevarx = "i"		
	}
	
	
	//Ensure varx in comparative is bi-categorical
	if "`comparative'" != "" {
		qui label list `varx'
		local ncat = r(max)
		cap assert `ncat' == 2
		if _rc != 0 {
			di as error "Comparative analysis requires that `varx' only has 2 categories but found `ncat'"
			exit _rc
		} 
	}
	
	local pcont: word count `contreg'
	if "`typevarx'" != "" & "`typevarx'" == "c" {
		local ++pcont
	}

	if `pcont' > 0 {
		local continuous = "continuous"
	}
	
	/*Model presenations*/
	local regressorss "`regressors'"
	/*Overall*/
	local lmuse = "mu_lse"
	local lmusp = "mu_lsp"
	
	if "`design'" == "general" | "`design'" == "comparative" {	
		local nuse = "mu_lse"
		local nusp = "mu_lsp"
	}
	else if "`design'" == "cbnetwork" {
		if "`interaction'" != "" { 
			local nuse = "mu_lse + Ipair*`Comparator' + `Index'"
			local nusp = "mu_lsp + Ipair*`Comparator' + `Index'"
		}
		else {
			local nuse = "mu_lse + Ipair + `Index'"
			local nusp = "mu_lsp + Ipair + `Index'"
		}
	}
	else {
		*abnetwork
		local nuse = "mu.`first'_lse"
		local nusp = "mu.`first'_lsp"
		tokenize `regressors'
		mac shift
		local regressorss "`*'"
	}
	
	//Build the rest of the equation
	*local VarX: word 1 of `regressors'
	local q: word count `regressors'
	forvalues i=1/`q' {
		local c:word `i' of `regressorss'
		
		//se
		if "`cveffect'" != "sp" {
			local nuse = "`nuse' + `c'"
			if ("`interaction'" == "sesp" | "`interaction'" == "se") & `i' > 1 {
				local nuse = "`nuse' + `c'*`varx'"
			}
		}
		
		//sp
		if "`cveffect'" != "se" {
			local nusp = "`nusp' + `c'"
			if ("`interaction'" == "sesp" | "`interaction'" == "sp") & `i' > 1 {
				local nusp = "`nusp' + `c'*`varx'"
			}
		}
	}
	if ("`catreg'" != " " | "`typevarx'" =="i" | ("`design'" == "comparative" | "`design'" == "cbnetwork"))  {

		if "`design'" == "cbnetwork" {
			local catregs = "`catreg' `Index'"
		}

		if "`design'" == "comparative" {
			local catregs = "`catreg' `varx'" 
		}
		if "`design'" == "abnetwork" {
			tokenize `catreg'
			macro shift
			local catregs "`*'"
		}
		if "`design'" == "general" {
			local catregs "`catreg'"
		}
	}

	if "`subgroup'" == "" & ("`catreg'" != "" | "`typevarx'" =="i" ) {
		if "`outplot'" == "abs" {
			if "`typevarx'" =="i" {
				local groupvar = "`varx'"
			}
			else {
				local groupvar : word 1 of `catreg'
			}
		}
		if ("`outplot'" == "rr" | "`outplot'" == "or")& "`varx'" !="" {
			local groupvar : word 1 of `catreg'
		}
	}
	if "`by'" != "" {
		local groupvar "`by'"
		local byvar "`by'"
		*How many times to loop
		qui label list `by'
		local nlevels = r(max)
	}
	if "`design'" == "abnetwork" {
		local groupvar "`first'"
		local overall "nooverall"
		if "`outplot'" == "rr" | "`outplot'" == "or" {
			local itable "noitable"
		}
	}
	if "`by'" == "" & "`cbnetwork'" != "" {
		local groupvar  "`Index'"
		local byvar "`Index'"
	} 
	
	*Stratify not allow in cbnetwork or abnetwork analysis
	if "`stratify'" != "" {
		if ("`design'" == "cbnetwork") | ("`design'" == "abnetwork")  {
			di as error"The option stratify is not allowed in `design' analysis"
			exit
		}
	}
	
	*Check by is active & that there are more than 1 levels
	if "`stratify'" != "" {
		if "`by'" == "" {
			di as error "The by() variable needs to be specified in stratified analysis"
			exit			
		}
		else {
			if `nlevels' < 2 {
				di as error "The by() variable should have atleast 2 categories in stratified analysis"	
				exit
			}
		}
		//nomc 
		local mc
	}
	
	if "`groupvar'" == "" {
		local subgroup nosubgroup
	}
	//nobox if summaruonly
	if "`summaryonly'" != "" {
		local box "nobox"
	}
	
	qui gen `sp' = 1 - `se'
	*qui gen `use' = .
	
	//fit the model
	if "`progress'" != "" {
		local echo noi
	}
	else {
		local echo qui
	}
	
	*Loop should begin here
	if "`stratify'" == "" {
		local nlevels = 0
	}

	local i = 1
	local byrownames 
	local bybirownames
	
	if "`design'" == "abnetwork"  {
		local hetdim 7
	}
	else {
		if (`p' == 0) & ("`model'" == "random") &  ("`design'" != "cbnetwork" )  {
			local hetdim 5
		}
		else {
			local hetdim 4
		}
	}
	
	if "`outplot'" == "abs" {
		local sumstatse "Sensitivity"
		local sumstatsp "Specificity"
	}
	else {
		if "`outplot'" == "or"  {
			local sumstatse "Sensitivity OR"
			local sumstatsp "Specificity OR"
		}
		else {
			local sumstatse "Relative Sensitivity"
			local sumstatsp "Relative Specificity"
		}
	}
	
	//Should run atleast once
	while `i' < `=`nlevels' + 2' {
		local modeli = "`model'"
		local modeloptsi = "`modelopts'"
		local smoothi = "`smooth'"
		local getmodel
		local optimizedi = 0
	
		//don't run last loop if stratify
		if (`i' > `nlevels') & ("`stratify'" != "") & ("`design'" == "comparative") {
			local overall "nooverall"
			continue, break
		}
	
		*Stratify except the last loop for the overall
		if (`i' < `=`nlevels' + 1') & ("`stratify'" != "") {
			local strataif `"if `by' == `i'"'
			local ilab:label `by' `i'
			local stratalab `":`by' = `ilab'"'
			local ilab = ustrregexra("`ilab'", " ", "_")
			local byrownames = "`byrownames' `by':`ilab'"
			local byrowname = "`by'|`ilab'"
			if "`design'" == "comparative" & "`stratify'" != "" {
				local bybirownames = "`bybirownames' `ilab':`baselab' `ilab':`indexlab' `ilab':Overall"
			}
			
			*Check if there is enough data in each strata
			//Number of obs in the analysis
			qui egen `obsid' = group(`rid') if `by' == `i'
			qui summ `obsid'
			local Nobs= r(max)
			drop `obsid'	

			//Number of studies in the analysis
			qui egen `uniq' = group(`studyid') if `by' == `i'
			qui summ `uniq'
			local Nuniq = r(max)
			drop `uniq'	
		}
		else {
			//Skip if overall not needed
			if ("`overall'" != "") & (`i' > `=`nlevels'+1')  & ("`stratify'" != "" | "`design'" == "comparative")  {
				continue, break
			}
			
			//Don't smoothen
			//Don't smoothen after last loop if stratify
			if (`i' > `nlevels') & ("`stratify'" != "") {
				local smoothi
			}
			
			//Nullify
			local strataif 					
			if "`stratify'" != "" {
				local stratalab ": all studies"
				local byrownames = "`byrownames' Overall"	
				local byrowname = "All_studies"
			}
			
			//Number of obs in the analysis
			qui count
			local Nobs= r(N)
			if "`cbnetwork'" != "" {
				local Nobs = `Nobs'*0.25
			}
			else {
				local Nobs = `Nobs'*0.5
			}

			qui egen `uniq' = group(`studyid')
			qui summ `uniq'
			local Nuniq = r(max)
			drop `uniq'
		}
	
		if "`design'" == "comparative" {
			cap assert mod(`Nobs', 2) == 0 
			if _rc != 0 {
				di as error "Comparative analysis requires 2 observations per study"
				exit _rc
			}
		}
		if "`design'" == "abnetwork" {
			cap assert `Nobs'/`Nuniq' >= 2 
			if _rc != 0 {
				di as error "abnetwork design requires atleast 2 observations per study"
				exit _rc
			}
		}		
		if `Nuniq' < 3 & "`modeli'" == "random"  {
			local modeli fixed //If less than 3 studies, use fixed model
			if "`modeloptsi'" != "" {
				local modeloptsi
				di as res _n  "Warning: Random-effects model options ignored."
				di as res _n  "Warning: Fixed-effects model fitted instead."
			}
		}		
	
		*Run model if more than 1 study
		if `Nobs' > 1 {			
			//Requested model
			`echo' fitmodel `event' `total' `se' `sp' `strataif', bcov(`bcov') wcov(`wcov') modelopts(`modeloptsi') model(`modeli') ///
			regexpression(`regexpression') sid(`studyid') `design' ipair(`ipair') level(`level') nested(`first')
			
			//Returned model	
			local getmodel = r(model)
			
			//Returned estimates			
			estimates store metadta_modest

			cap drop _ESAMPLE
			qui gen _ESAMPLE = e(sample)
			
			mat `coefmat' = e(b)
			mat `coefvar' = e(V)
			
			qui estat ic
			mat `matgofi' = r(S)
			local BIC =  `matgofi'[1, 6]
			mat `matgofi' = `matgofi'[1..., 5..6]
			mat rownames `matgofi' = Value	
		}
		else {
			mat `coefmat' = (0 , 0)
			mat `coefvar' = e(V)
			local getmodel = "none"
		}
		
		//Obtain the prediction		
		//if random, needs atleast 7 studies to run predict command
		qui {
			count
			local nobs = r(N)
			if "`getmodel'" == "random" {
				if ((`nobs' < 7) & ("`getmodel'" == "random")) {
					local multipler = int(ceil(7/`nobs'))
					qui expand `multipler', gen(`newobs')
				}
			}
			
			cap drop `predevent'
			predict `predevent', mu
			//Revert to original data if filler data was generated
			 if (("`getmodel'" == "random") & (`nobs' < 7))  {
				keep if !`newobs'
			}
			
			//compute the weight
			cap drop `iw'
			if "`getmodel'" == "random" {
				gen `iw' = `total'*(`predevent'/`total')*(1 - `predevent'/`total')
			}
			else {
				gen `iw' = `total'
			}
			
			//compute the relative weight
			sum `iw' if (_ESAMPLE == 1) & `se'
			local seW = r(sum)
			
			sum `iw' if (_ESAMPLE == 1) & !`se'
			local spW = r(sum)
		
			//compute the weights
			replace _WT = (`iw'/`seW')*100 if (_ESAMPLE == 1) & (_WT == .) & `se'
			replace _WT = (`iw'/`spW')*100 if (_ESAMPLE == 1) & (_WT == .) & !`se'
		}
		
		estcovar, matrix(`coefmat') model(`getmodel') bcov(`bcov') wcov(`wcov') `design'
		local kcov = r(k) //#covariance parameters
		mat `BVari' = r(BVar)  //Between var-cov
		mat `WVari' = r(WVar)  //Within var-cov
		mat colnames `BVari' = logitse logitsp
		mat rownames `BVari' = logitse logitsp
		
		mat colnames `WVari' = logitse logitsp
		mat rownames `WVari' = logitse logitsp
		
		local tausesp 	= `BVari'[1, 2]
		local rho 		= `BVari'[1, 2]/sqrt(`BVari'[1, 1]*`BVari'[2, 2])
		local tau2se 	= `BVari'[1, 1]
		local tau2sp	= `BVari'[2, 2]
		local tau2g		= (1 - (`BVari'[1, 2]/sqrt(`BVari'[1, 1]*`BVari'[2, 2]))^2)*`BVari'[1, 1]*`BVari'[2, 2]
			
		if("`sumtable'" != "") {
			local loddslabel = "Log_odds"
			local abslabel = "Proportion"
			local rrlabel = "Rel_Ratio"
		}
		if `Nobs' > 1 {
			local S_1 = e(N) -  e(k) //df
		}
		else {
			local S_1 = .
		}
		
		local S_2 = . //between study heterogeneity chi2
		local S_3 = . // between study heterogeneity pvalues
		local i2g = . //Isq
		local i2se = . //Isqse
		local i2sp = . //Isqsp
		local S_81 = . //Full vs Null chi2 -- se
		local S_91 = . //Full vs Null  pvalue -- se
		local S_891 = . //Full vs Null  df -- se
		local S_82 = . //Full vs Null chi2 -- sp
		local S_92 = . //Full vs Null  pvalue -- sp
		local S_892 = . //Full vs Null  df -- sp

		//Consider a reduced model	
		if "`getmodel'" == "random" {
			qui estimates restore metadta_modest
			local S_2 = e(chi2_c)
			local S_3 = e(p_c)
		}
	
		if `p' == 0 & "`design'" == "general" {
			/*Compute I2*/
			mat `Esigma' = J(2, 2, 0) /*Expected within study variance*/
			
			if "`strataif'" != "" {
				qui gen `invtotal' = 1/`total'
				
				qui summ `invtotal' if `se' & `by' == `i'
				local invtotalse = r(sum)
				
				qui summ `invtotal' if `sp' & `by' == `i'
				local invtotalsp = r(sum)
				drop `invtotal'
			}
			else {
				qui gen `invtotal' = 1/`total'
				qui summ `invtotal' if `se'
				local invtotalse = r(sum)
				
				qui summ `invtotal' if `sp' 
				local invtotalsp = r(sum)
			}

			mat `Esigma'[1, 1] = (exp(`BVari'[1, 1]*0.5 + `coefmat'[1, 1]) + exp(`BVari'[1, 1]*0.5 - `coefmat'[1, 1]) + 2)*(1/(`Nuniq'))*`invtotalse'
			mat `Esigma'[2, 2] = (exp(`BVari'[2, 2]*0.5 + `coefmat'[1, 2]) + exp(`BVari'[2, 2]*0.5 - `coefmat'[1, 2]) + 2)*(1/(`Nuniq'))*`invtotalsp'
			
			local detEsigma = `Esigma'[1, 1]*`Esigma'[2, 2]
			
			local detSigma = (1 - (`BVari'[2, 1]/sqrt(`BVari'[1, 1]*`BVari'[2, 2]))^2)*`BVari'[1, 1]*`BVari'[2, 2]
			
			local IsqE = sqrt(`detSigma')/(sqrt(`detEsigma') + sqrt(`detSigma'))
			
			local i2g = `IsqE'
			local i2se = (`BVari'[1, 1]/(`Esigma'[1, 1] + `BVari'[1, 1]))  //se
			local i2sp = (`BVari'[2, 2]/(`Esigma'[2, 2] + `BVari'[2, 2])) //sp
		}
		
		local nmc = 0
		if (`p' > 0  & "`mc'" == "") {
			forvalues j=1/2 {
				local S_9`j' = .
				local S_8`j' = .
				local S_89`j' = .
			}
			
			if "`interaction'" !="" {
				local confariates "`confounders'"
			}
			if "`interaction'" ==""  {
				if ("`design'" == "abnetwork" | "`design'" == "comparative") {
					tokenize `regressors'
					macro shift
					local confariates "`*'"
				}
				else {
					local confariates "`regressors'"
				}
				
			}
			local initialse 1
			local initialsp 1
			local rownamesmcse
			local rownamesmcsp
			
			/*if "`confariates'" != "" {
				di "*********************************** ************* ***************************************"
				di as txt "Just a moment - Fitting reduced models for comparisons"
			}*/
			foreach c of local confariates {
				if "`cveffect'" != "sp" {
					if "`interaction'" =="sesp" | "`interaction'" =="se" {
						local xterm = "`c'#`typevarx'.`varx'"
						local xnu = "`c'*`varx'"
					}
					else {
						local xterm = "`c'"
						local xnu = "`c'"
					}
					//Sensivitivity terms
					local nullse		
					foreach term of local seregexpression {
						if ("`term'" != "i.`xterm'#c.`se'")&("`term'" != "c.`xterm'#c.`se'")&("`term'" != "`xterm'#c.`se'") {
							local nullse "`nullse' `term'"
						} 
					}
					
					local nullnuse = subinstr("`nuse'", "+ `xnu'", "", 1)
					
					local textsemc1`nmc' "Ommitted : `xnu' in logit(se)"
					local textsemc2`nmc' "logit(se) = `nullnuse'"
					local textsemc3`nmc' "logit(sp) = `nusp'"
										
					/*
					di as res _n "Ommitted : `xnu' in logit(se)"
					di as res "{phang} logit(se) = `nullnuse'{p_end}"
					di as res "{phang} logit(sp) = `nusp'{p_end}"
					*/
					
					local nullse = "`nullse' `spregexpression'"
					`echo' fitmodel `event' `total' `se' `sp' `strataif',  bcov(`bcov') wcov(`wcov') modelopts(`modeloptsi') model(`getmodel') ///
					regexpression(`nullse') sid(`studyid') `design' ipair(`ipair') level(`level') nested(`first')
					
					estimates store metadta_Nullse
					
					//LR test the model
					qui lrtest metadta_modest metadta_Nullse
					local selrp :di %10.`dp'f chi2tail(r(df), r(chi2))
					local selrchi2 = r(chi2)
					local selrdf = r(df)
					estimates drop metadta_Nullse
					
					if `initialse'  {
						mat `semci' = [`selrchi2', `selrdf', `selrp']
						local initialse 0
					}
					else {
						mat `semci' = [`selrchi2', `selrdf', `selrp'] \ `semci'
					}
					local rownamesmcse "`rownamesmcse' `xnu'"
				}
				if "`cveffect'" != "se" {
					if "`interaction'" =="sesp" | "`interaction'" =="sp" {
						local xterm = "`c'#`typevarx'.`varx'"
						local xnu = "`c'*`varx'"
					}
					else {
						local xterm = "`c'"
						local xnu = "`c'"
					}
					//Specificity terms
					local nullsp		
					foreach term of local spregexpression {
						if ("`term'" != "i.`xterm'#c.`sp'")&("`term'" != "c.`xterm'#c.`sp'")&("`term'" != "`xterm'#c.`sp'") {
							local nullsp "`nullsp' `term'"
						} 
					}
					
					local nullnusp = subinstr("`nusp'", "+ `xnu'", "", 1)
					
					local textspmc1`nmc' "Ommitted : `xnu' in logit(sp)"
					local textspmc2`nmc' "logit(se) = `nuse'"
					local textspmc3`nmc' "logit(sp) = `nullnusp'"
					
					/*
					di as res _n "Ommitted : `xnu' in logit(sp)"
					di as res "{phang} logit(se) = `nuse'{p_end}"
					di as res "{phang} logit(sp) = `nullnusp'{p_end}"
					*/
					local nullsp = "`seregexpression' `nullsp'" 
					`echo' fitmodel `event' `total' `se' `sp' `strataif', bcov(`bcov') wcov(`wcov') modelopts(`modeloptsi') model(`getmodel') ///
					regexpression(`nullsp') sid(`studyid') `design' ipair(`ipair') level(`level') nested(`first')
					estimates store metadta_Nullsp
					
					//LR test the model
					qui lrtest metadta_modest metadta_Nullsp
					local splrp :di %10.`dp'f chi2tail(r(df), r(chi2))
					local splrchi2 = r(chi2)
					local splrdf = r(df)
					estimates drop metadta_Nullsp
					
					if `initialsp' {
						mat `spmci' = [`splrchi2', `splrdf', `splrp']
						local initialsp 0
					}
					else {
						mat `spmci' = [`splrchi2', `splrdf', `splrp'] \ `spmci'
					}
					local rownamesmcsp "`rownamesmcsp' `xnu'"
				}
				local ++nmc
			}
			
			//Ultimate null model if more than one term
			if (`p' > 0) & (`nmc' > 1) {
				if "`cveffect'" != "sp" {
					local nullse `se'		
					local nullse = "`nullse' `spregexpression'"
					`echo' fitmodel `event' `total' `se' `sp' `strataif',  bcov(`bcov') wcov(`wcov') modelopts(`modeloptsi') model(`getmodel') regexpression(`nullse') ///
					sid(`studyid') `design' ipair(`ipair') level(`level') nested(`first')
					estimates store metadta_Nullse
					
					qui lrtest metadta_modest metadta_Nullse
					local selrp :di %10.`dp'f chi2tail(r(df), r(chi2))
					local selrchi2 = r(chi2)
					local selrdf = r(df)
					estimates drop metadta_Nullse
					
					if `initialse'  {
						mat `semci' = [`selrchi2', `selrdf', `selrp']
						local initialse 0
					}
					else {
						mat `semci' = `semci' \ [`selrchi2', `selrdf', `selrp']
					}
					local rownamesmcse "`rownamesmcse' All"
				}
				if "`cveffect'" != "se" {
					local nullsp `sp'
					local nullsp = "`seregexpression' `nullsp'"
					`echo' fitmodel `event' `total' `se' `sp' `strataif',  bcov(`bcov') wcov(`wcov') modelopts(`modeloptsi') model(`getmodel') regexpression(`nullsp') ///
					sid(`studyid') `design' ipair(`ipair') level(`level') nested(`first')
					estimates store metadta_Nullsp
					
					qui lrtest metadta_modest metadta_Nullsp
					local splrp :di %10.`dp'f chi2tail(r(df), r(chi2))
					local splrchi2 = r(chi2)
					local splrdf = r(df)
					estimates drop metadta_Nullsp
					
					if `initialsp' {
						mat `spmci' = [`splrchi2', `splrdf', `splrp']
						local initialsp 0
					}
					else {
						mat `spmci' = `spmci' \ [`splrchi2', `splrdf', `splrp']
					}
					local rownamesmcsp "`rownamesmcsp' All"
				}	
			}
			
			if "`cveffect'" != "sp" & `nmc' > 0 {
				mat roweq `semci' = Sensitivity
				mat rownames `semci' = `rownamesmcse'
				mat colnames `semci' =  chi2 df pval
			}

			if "`cveffect'" != "se" & `nmc' > 0 {
				mat roweq `spmci' = Specificity
				mat rownames `spmci' = `rownamesmcsp'
				mat colnames `spmci' = chi2 df pval
			}
		}
		if `nmc' == 0 {
			local mc "nomc"
		}
		
		mat `BVari' = (`tau2se', `i2se', `tau2sp', `i2sp', `tau2g', `i2g', `tausesp', `rho')
		mat colnames `BVari' = sensitivity:tausq sensitivity:isq specificity:tausq specificity:isq Generalized:tausq Generalized:isq covar rho 
		mat rownames `BVari' = Overall
		
		*mat `isq2i' = (`S_71', `S_7' \ `S_7', `S_72') //Isq
		*mat `Isq2i' = (`S_7' , `S_71', `S_72')
		*mat colnames `Isq2i' = Generalized logitse logitse
		*mat rownames `Isq2i' = Overall
		
		//model comparison
		*mat `mci' = (`S_81', `S_91',  `S_891', `S_82', `S_92' ,  `S_892' ) // Full vs Null 
		*mat colnames `mci' = sensitivity:Chi2 sensitivity:pval sensitivity:df specificity:Chi2 specificity:pval specificity:df
		*mat rownames `bgheti' = overall
		
		mat `refei' = (`S_2', `kcov', `S_3') // chisq re vs fe, df, pv re vs fe
		mat colnames `refei' = Chi2 df pval
		mat rownames `refei' = Overall
		
		//if no df, get exact estimates	
		if `Nobs' > 1 & `p' > 0 & "`getmodel'" == "fixed" {
			estimates restore metadta_modest
			
			local datapoints = e(N) 
			
			local df = e(df) //number of df
				
			if (`=`datapoints'*0.5 - `df'' < 1) | (`df' < 1) {
				 qui {
					
					//Exact inference					
					sum `event' if e(sample) & `se'
					local n1 = r(sum)
					sum `total' if e(sample) & `se'
					local N1 = r(sum)
					
					metadta_absexactci `N1' `n1',  level(`level') //exact ci
					mat `absexactse' = r(absexact)
					
					sum `event' if e(sample) & `sp'
					local n2 = r(sum)
					sum `total' if e(sample) & `sp'
					local N2 = r(sum)
					
					metadta_absexactci `N2' `n2',  level(`level') //exact ci
					mat `absexactsp' = r(absexact)
					
					mat `absexact' = `absexactse' \ `absexactsp'					
					mat rownames `absexact' = Sensitivity Specificity
					
				}
			}
			local optimizedi = 1
		}
		
	
		if `Nobs' > 1 {		
			
			//LOG ODDS			
			estp `strataif', estimates(metadta_modest) sumstat(`loddslabel') depname(Effect) interaction(`interaction') cveffect(`cveffect') ///
				catreg(`catreg') contreg(`contreg') se(`se') sp(`sp') level(`level') dp(`dp') varx(`varx') typevarx(`typevarx') ///
				by(`by') `tradeoff' regexpression(`regexpression') `design'  `stratify'

			mat `Vi' = r(Vmatrix) //var-cov for catreg & overall 

			mat `logoddsi' = r(outmatrix)
			mat `selogoddsi' = r(outmatrixse)
			mat `splogoddsi' = r(outmatrixsp)
			
			if "`tradeoff'" != "" {
				mat `cutmat' = r(cutmat)
			}					
   
			if "`interaction'" == "" {
				//names of the V matrix
				local vnames
				local rnames :rownames `Vi'
				local nrowsv = rowsof(`Vi')
				forvalues r = 1(1)`nrowsv' {
				//Labels
					local rname`r':word `r' of `rnames'
					tokenize `rname`r'', parse("#")					
					
					local left = "`1'"
					local right = "`3'"
					
					tokenize `left', parse(.)
					local parm = substr("`1'", 1, 1)
					if `parm' == 0 {
						local eqlab "sp"
					}
					else {
						local eqlab "se"
					}
					
					if "`right'" == "" {
						local lab = "Overall"
					}
					else {
						tokenize `right', parse(.)
						local rightv = "`3'"
						local rightlabel = substr("`1'", 1, 1)
					
						local rlab:label `rightv' `rightlabel'
						local rlab = ustrregexra("`rlab'", " ", "-")
						local lab = "`rightv'_`rlab'"
					}
					
					local vnames = "`vnames' `eqlab':`lab'"	
				}
				mat rownames `Vi' = `vnames'
				mat colnames `Vi' = `vnames'
			}
				
			//ABS
			estp `strataif', estimates(metadta_modest) sumstat(`abslabel') ///
				depname(Effect) interaction(`interaction') cveffect(`cveffect') ///
				catreg(`catreg') contreg(`contreg')  se(`se') level(`level') ///
				expit power(`power') dp(`dp') varx(`varx') ///
				typevarx(`typevarx') by(`byvar') regexpression(`regexpression') `design' `stratify' 
		
			
			mat `absouti' = r(outmatrix)
			mat `absoutsei' = r(outmatrixse)
			mat `absoutspi' = r(outmatrixsp)
			
			
			//Simulation
			postsim `strataif', orderid(`rid') studyid(`studyid') todo(p) ///
				estimates(metadta_modest) absoutse(`absoutsei') absoutsp(`absoutspi')  ///
				level(`level')  model(`getmodel')  by(`byvar')  se(`se') sp(`sp') ///
				`design'  varx(`varx') bcov(`bcov') wcov(`wcov') nsims(`nsims') link(`link') regressors(`regressors')
					
			mat `popabsoutsei' = r(outmatrixse)
			mat `popabsoutspi' = r(outmatrixsp)
			mat `popabsouti' = r(outmatrix)
			
			
			//Create a matrix with exact abs
			local nrowsepop = rowsof(`popabsoutsei')
			mat `exactabsoutsei' = J(`nrowsepop', 5, .)
			local eserownames : rownames `popabsoutsei'
			mat rownames `exactabsoutsei' = `eserownames'
			
			//Replace overall with exact estimate if necessary
			if `optimizedi' {
				mat `exactabsoutsei'[`nrowsepop', 1] = `absexactse'[1, 1]
				mat `exactabsoutsei'[`nrowsepop', 4] = `absexactse'[1, 5]
				mat `exactabsoutsei'[`nrowsepop', 5] = `absexactse'[1, 6]
			}
						
			local nrowsppop = rowsof(`popabsoutspi')
			mat `exactabsoutspi' = J(`nrowsppop', 5, .)
			local esprownames : rownames `popabsoutspi'
			mat rownames `exactabsoutspi' = `esprownames'
			
			//Replace overall with exact estimate if necessary
			if `optimizedi' {
				mat `exactabsoutspi'[`nrowsppop', 1] = `absexactsp'[1, 1]
				mat `exactabsoutspi'[`nrowsppop', 4] = `absexactsp'[1, 5]
				mat `exactabsoutspi'[`nrowsppop', 5] = `absexactsp'[1, 6]
			}
			
			//RR
			if `pcat' > 0 | "`typevarx'" == "i" {
				estr `strataif', estimates(metadta_modest) sumstat(`rrlabel')  cveffect(`cveffect') ///
				catreg(`catreg') se(`se') level(`level') power(`power') dp(`dp') varx(`varx') ///
				typevarx(`typevarx') by(`byvar') `stratify' regexpression(`regexpression') `design' ///
				baselevel(`basecode') refpos(`refpos')  comparator(`Comparator') 
			
				mat `rrouti' = r(RRoutmatrix)
				mat `serrouti' = r(RRoutmatrixse)
				mat `sprrouti' = r(RRoutmatrixsp)
				
				mat `orouti' = r(ORoutmatrix)
				mat `seorouti' = r(ORoutmatrixse)
				mat `sporouti' = r(ORoutmatrixsp)
				
				local inltest = r(inltest)

				if "`inltest'" == "yes" {
					mat `setestnlrri' = r(setestnlRR) //Equality of RR
					mat `sptestnlrri' = r(sptestnlRR) 
					mat `setestnlori' = r(setestnlOR) //Equality of RR
					mat `sptestnlori' = r(sptestnlOR) 
				}
				
				//Simulation
				postsim `strataif', orderid(`rid') studyid(`studyid') todo(r)  cveffect(`cveffect')  ///
					estimates(metadta_modest) serrout(`serrouti') sprrout(`sprrouti')    ///
					level(`level')  model(`getmodel')  by(`byvar') se(`se') sp(`sp') ///
					`design'  varx(`varx') bcov(`bcov') wcov(`wcov') nsims(`nsims') link(`link') regressors(`regressors')
						
				mat `popserrouti' = r(rroutmatrixse)
				mat `popsprrouti' = r(rroutmatrixsp)
				mat `poprrouti' = r(rroutmatrix)
				
				mat `popseorouti' = r(oroutmatrixse)
				mat `popsporouti' = r(oroutmatrixsp)
				mat `poporouti' = r(oroutmatrix)
			}
			else {
				local rr "norr"
			}
			
			//Smooth estimates
			//simulations
			if "`smoothi'" != "" {
				postsim `strataif', orderid(`rid') studyid(`studyid') todo(smooth) estimates(metadta_modest)  ///
						level(`level')  model(`getmodel')  by(`byvar') se(`se') sp(`sp')    ///
						modeles(`modeles') modellci(`modellci') modeluci(`modeluci') outplot(`outplot')	///
						`design'  varx(`varx')  bcov(`bcov') wcov(`wcov') nsims(`nsims') link(`link') regressors(`regressors')
			}			
		}
		*if one study
		else { 
			mat `absouti' = J(2, 6, 0)
			mat `absoutsei' = J(1, 6, 0)
			mat `absoutspi' = J(1, 6, 0)
			
			mat `popabsoutsei' = J(1, 5, .)
			mat `popabsoutspi' = J(1, 5, .)
			
			mat `exactabsoutsei' = J(1, 5, .)
			mat `exactabsoutspi' = J(1, 5, .)
			
			mat `rrouti' = J(2, 6, 1)
			mat `serrouti' = J(1, 6, 1)
			mat `sprrouti' = J(1, 6, 1)
			
			mat `popserrouti' = J(1, 5, 1)
			mat `popsprrouti' = J(1, 5, 1)
			
			mat `setestnlrri' = J(1, 3, .)
			mat `sptestnlrri' = J(1, 3, .)	

			mat `orouti' = J(2, 6, 1)
			mat `seorouti' = J(1, 6, 1)
			mat `sporouti' = J(1, 6, 1)
			
			mat `popseorouti' = J(1, 5, 1)
			mat `popsporouti' = J(1, 5, 1)
			
			mat `setestnlori' = J(1, 3, .)
			mat `sptestnlori' = J(1, 3, .)				
			
			mat `Vi' = J(2, 2, .) //var-cov for catreg & overall 
			mat `logoddsi' = J(2, 6, 0)
			mat `selogoddsi' = J(1, 6, 0)
			mat `splogoddsi'	= J(1, 6, 0)
			mat	`BVari' = J(1, 8, .)
			mat	`WVari' = J(1, 2, 0)
			if ((`p' > 0 & "`design'" != "abnetwork") | (`p' > 1 & "`design'" == "abnetwork")) & "`mc'" == "" {
				if "`cveffect'" != "sp" {
					mat `semci'	= J(1, 3, .)
				}
				if "`cveffect'" != "se" {
					mat `spmci'	= J(1, 3, .)
				}
			}
			mat `matgofi' = J(1, 2, .)
		}
		/*
		if "`stratify'" != "" {
			if ((`p' > 0 & "`design'" != "abnetwork") | (`p' > 1 & "`design'" == "abnetwork")) & "`mc'" == "" {
				if "`cveffect'" != "sp" {
					mat roweq `semci' = Sensitivity
					mat coleq `semci' = `by'
				}
				if "`cveffect'" != "se" {
					mat roweq `spmci' = Specificity
					mat coleq `spmci' = `by'
				}
			}
		}
		*/
		*Stack the matrices
		//Remove the first empty row
		if ("`stratify'" != "") {
			if (`p' > 0) {
				mat `logoddsi' 	= `logoddsi'[2..., 1...]
				mat `selogoddsi' = `selogoddsi'[2..., 1...]
				mat `splogoddsi' = `splogoddsi'[2..., 1...]
				
				mat `absouti' 	= `absouti'[2..., 1...]
				mat `absoutsei' = `absoutsei'[2..., 1...]
				mat `absoutspi' = `absoutspi'[2..., 1...]
				
				mat `popabsouti' 	= `popabsouti'[2..., 1...]
				mat `popabsoutsei' 	= `popabsoutsei'[2..., 1...]
				mat `popabsoutspi' 	= `popabsoutspi'[2..., 1...]
				
				mat `exactabsoutsei' = `exactabsoutsei'[2..., 1...]
				mat `exactabsoutspi' = `exactabsoutspi'[2..., 1...]
			}
			//Add row for byvar
			/*
			tempname byvarrow
			mat `byvarrow' = J(1, 6, .)
			mat rownames `byvarrow' = `byrowname'
			mat `absoutsei' = `byvarrow' \ `absoutsei'
			mat `absoutspi' = `byvarrow' \ `absoutspi'
			*/
			mat roweq `absoutsei' = `byrowname'
			mat roweq `absoutspi' = `byrowname'
			
			mat roweq `popabsoutsei' = `byrowname'
			mat roweq `popabsoutspi' = `byrowname'
			
			mat roweq `exactabsoutsei' = `byrowname'
			mat roweq `exactabsoutspi' = `byrowname'
			
			if "`rr'" == "" {
				mat roweq `serrouti' = `byrowname'
				mat roweq `sprrouti' = `byrowname'
				
				mat roweq `seorouti' = `byrowname'
				mat roweq `sporouti' = `byrowname'
				
				mat roweq `popserrouti' = `byrowname'
				mat roweq `popsprrouti' = `byrowname'
				
				mat roweq `popseorouti' = `byrowname'
				mat roweq `popsporouti' = `byrowname'
			}
		}
		
		if "`stratify'" != "" {
			if ((`p' > 0 & "`design'" != "abnetwork") | (`p' > 1 & "`design'" == "abnetwork")) & "`mc'" == "" {
				if "`cveffect'" != "sp" {
					mat roweq `semci' = `byrowname'
				}
				if "`cveffect'" != "se" {
					mat roweq `spmci' = `byrowname'
				}
			}
		}
		
		if `i' == 1 {		
			mat `absout' =	`absouti'	
			mat `absoutse' = `absoutsei'	
			mat `absoutsp' = `absoutspi'
			
			mat `popabsoutse' =  `popabsoutsei'	
			mat `popabsoutsp' = `popabsoutspi'
			mat `popabsout' = `popabsouti' 
			
			mat `exactabsoutse' =  `exactabsoutsei'	
			mat `exactabsoutsp' =  `exactabsoutspi'
			
			if "`rr'" == "" {		
				mat `rrout' =	`rrouti'
				mat `serrout' = `serrouti'	
				mat `sprrout' = `sprrouti'
				
				mat `orout' =	`orouti'
				mat `seorout' = `seorouti'	
				mat `sporout' = `sporouti'
				
				mat `popserrout' =  `popserrouti'	
				mat `popsprrout' =  `popsprrouti'
				mat `poprrout' =  `poprrouti'
				
				mat `popseorout' = `popseorouti'	
				mat `popsporout' = `popsporouti'
				mat `poporout' =  `poporouti'
				
				if "`inltest'" == "yes" {
					mat `setestnlrr' = `setestnlrri' //Equality of RR
					mat `sptestnlrr' = `sptestnlrri' 
					
					mat `setestnlor' = `setestnlori' //Equality of RR
					mat `sptestnlor' = `sptestnlori' 
				}
			}
			mat `BVar' = `BVari'
			mat `WVar' = `WVari'
			mat `V' = `Vi' //var-cov for catreg & overall 
			mat `logodds' = `logoddsi'
			mat `selogodds' = `selogoddsi'
			mat `splogodds' = `splogoddsi'	
			*mat `bghet' = `bgheti'
			mat `refe' = `refei'
			*mat `Isq2' =  `Isq2i'
			if ((`p' > 0 & "`design'" != "abnetwork") | (`p' > 1 & "`design'" == "abnetwork")) & "`mc'" == "" {
				if "`cveffect'" != "sp" {
					mat `semc'	= `semci'
				}
				if "`cveffect'" != "se" {
					mat `spmc'	= `spmci'
				}
			}
			mat `matgof' = `matgofi'
		}
		else {
			mat `absout' = `absout' \ `absouti'
			mat `absoutse' = `absoutse' \ `absoutsei'
			mat `absoutsp' =  `absoutsp' \ `absoutspi'
			
			mat `popabsoutse' = `popabsoutse' \ `popabsoutsei'
			mat `popabsoutsp' =  `popabsoutsp' \ `popabsoutspi'
			mat `popabsout' =  `popabsout' \ `popabsouti' 
			
			mat `exactabsoutse' = `exactabsoutse' \ `exactabsoutsei'
			mat `exactabsoutsp' = `exactabsoutsp' \ `exactabsoutspi'
			
			if "`rr'" == "" {
				mat `rrout' = `rrout' \ `rrouti'
				mat `serrout' = `serrout' \ `serrouti'
				mat `sprrout' =  `sprrout' \ `sprrouti'
				
				mat `orout' = `orout' \ `orouti'
				mat `seorout' = `seorout' \ `seorouti'
				mat `sporout' =  `sporout' \ `sporouti'
				
				mat `popserrout' = `popserrout' \ `popserrouti'
				mat `popsprrout' =  `popsprrout' \ `popsprrouti'
				mat `poprrout' =  `poprrout' \ `poprrouti' 
				
				mat `popseorout' = `popseorout' \ `popseorouti'
				mat `popsporout' =  `popsporout' \ `popsporouti'
				mat `poporout' =  `poporout' \ `poporouti' 
				
				if "`inltest'" == "yes" {
					mat `setestnlrr' = `setestnlrr' \ `setestnlrri' //Equality of RR
					mat `sptestnlrr' = `sptestnlrr' \ `sptestnlrri' 
					
					mat `setestnlor' = `setestnlor' \ `setestnlori' //Equality of RR
					mat `sptestnlor' = `sptestnlor' \ `sptestnlori' 
				}
			}
			mat `BVar' = `BVar' \ `BVari'
			mat `WVar' = `WVar' \ `WVari'
			mat `V' = `V' \ `Vi' 
			mat `logodds' = `logodds' \ `logoddsi'
			mat `selogodds' = `selogodds' \ `selogoddsi'
			mat `splogodds' = `splogodds' \ `splogoddsi'
			if ((`p' > 0 & "`design'" != "abnetwork") | (`p' > 1 & "`design'" == "abnetwork")) & "`mc'" == "" {
				if "`cveffect'" != "sp" {
					mat `semc'	= `semc' \ `semci'
				}
				if "`cveffect'" != "se" {
					mat `spmc'	= `spmc' \ `semci'
				}
			}
			*mat `bghet' =`bghet' \ `bgheti' 
			mat `refe' = `refe' \ `refei'
			*mat `Isq2' = `Isq2' \ `Isq2i'
			
			mat `matgof' = `matgof' \ `matgofi'
		}
		if "`stratify'" != "" {
			di as res _n "*********************************** Model for `stratalab' ***************************************" 
		}
		else {
			di as res _n "**************************************************************************" 
		}
		
		tokenize `varlist'
		di "{phang} `1' ~ binomial(se, `1' + `3'){p_end}"
		di "{phang} `4' ~ binomial(sp, `4' + `2'){p_end}"
		if "`getmodel'" == "random" {
			di "{phang} logit(se) = `nuse' + `studyid'_lse{p_end}"
			di "{phang} logit(sp) = `nusp' + `studyid'_lsp{p_end}"
		}
		if "`getmodel'" == "fixed" {
			di "{phang} logit(se) = `nuse'{p_end}"
			di "{phang} logit(sp) = `nusp'{p_end}"
		}
		if "`getmodel'" == "random" {
			if "`bcov'" == "se" {
				di "{phang}`studyid'_lse ~ normal(0, tau2){p_end}"
			}
			else if "`bcov'" == "sp"  {
				di "{phang}`studyid'_lsp ~ normal(0, tau2){p_end}"
			}
			else {
				di "{phang}`studyid'_lse, `studyid'_lsp ~ biv.normal(0, tau2){p_end}"
			}
		}
		if "`design'" == "cbnetwork"  {
			di "{phang} Ipair = 0 if 1st set{p_end}"
			di "{phang} Ipair = 1 if 2nd set{p_end}"
		}
		*if "`design'" == "abnetwork" | "`design'" == "comparative" {
		if "`getmodel'" == "random"  {
			if "`design'" != "cbnetwork"   {
				if "`wcov'" == "identity" {
					di "{phang}`first'_lse, `first'_lsp  ~ N(0, `first'.sigma2){p_end}"
				}
				else if "`wcov'" == "independent"  {
					di "{phang}`first'_lse ~ N(0, `first'_lse.sigma2){p_end}"
					di "{phang}`first'_lsp ~ N(0, `first'_lsp.sigma2){p_end}"
				}
				qui label list `first'
				local nfirst = r(max)
			}
			if "`design'" == "cbnetwork"  {
				if "`wcov'" == "identity" {
					di "{phang}Ipair_lse, Ipair_lsp  ~ N(0, Ipair.sigma2){p_end}"
				}
				else if "`wcov'" == "independent"  {
					di "{phang}Ipair_lse ~ N(0, Ipair_lse.sigma2){p_end}"
					di "{phang}Ipair_lsp ~ N(0, Ipair_lsp.sigma2){p_end}"
				}
			}
		}
		if ("`catreg'" != " " | "`typevarx'" =="i" | ("`design'" == "comparative" | "`design'" == "cbnetwork"))  {
			di _n "{phang}Base levels{p_end}"
			di _n as txt "{pmore} Variable  -- Base Level{p_end}"
		}
		foreach fv of local catregs  {			
			local lab:label `fv' 1
			di "{pmore} `fv'  -- `lab'{p_end}"	
		}
		if "`design'" == "abnetwork" {
			local lab:label `first' `basecode'
			di "{pmore} `first'  -- `lab'{p_end}"
		}
		
		di as txt "{phang}Number of observations = " as res "`Nobs'{p_end}"
		di as txt "{phang}Number of studies = " as res "`Nuniq'{p_end}"
		if "`design'" == "abnetwork" {
			di as txt "{phang}Number of `first's = " as res "`nfirst'{p_end}"
		}
		
		if "`stratify'" != "" & `Nobs' > 1 { 
			di  _n	
			qui estimates restore metadta_modest
			local ilab = ustrregexra("`ilab'", " ", "_")
			local ilab = ustrregexra("`ilab'", "-", "_")
			local ilab = ustrregexra("`ilab'", "+", "")
			local ilab = ustrregexra("`ilab'", "/", "_")
			local ilab = "`ilab'" + "$by_index_"
			qui estimates store metadta_`ilab'
			display `"{stata "estimates replay metadta_`ilab'":Click to show the raw estimates}"'
		}
		else {
			if "$by_index_" !="" {
				di  _n	
				qui estimates restore metadta_modest
				local ilab = "$by_index_"
				qui estimates store metadta_`ilab'
				display `"{stata "estimates replay metadta_`ilab'":Click to show the raw estimates}"'
			}
			else {
				di  _n				
				display `"{stata "estimates replay metadta_modest":Click to show the raw estimates}"'
			}
		}
		
		if "`mc'" =="" {
			if "`confariates'" != "" {
				di as res _n "**************************************************************************" 
				di as txt "Fitted reduced models for comparisons"	
			}
			
			forvalues t=0(1)`=`nmc'-1'  {
				forvalues l=1(1)3 {
					di as res "`textsemc`l'`t''"
				}
				 
				di _n
				forvalues l=1(1)3 {
					di as res "`textsemc`l'`t''"
				}
			}
		}
				
		local ++i
	}
	*Loop should end here
	//Replace model with fitted model in ran once
	if `i' == 2 {
		local model "`getmodel'"
	}
	

	//rownames for the matrix
	if "`stratify'" != "" & `i' > 1 {
		*mat rownames `hetout' = `byrownames'
		mat `serow' = J(1, 6, .)
		mat `sprow' = J(1, 6, .)
		
		mat `pserow' = J(1, 5, .)
		mat `psprow' = J(1, 5, .)
		
		mat `semcrow' = J(1, 3, .)
		mat `spmcrow' = J(1, 3, .)
		
		mat rownames `serow' = "*--Sensitivity--*"
		mat rownames `sprow' = "*--Specificity--*" //19 characters
		mat rownames `pserow' = "*--Sensitivity--*"
		mat rownames `psprow' = "*--Specificity--*" //19 characters
		mat rownames `semcrow' = "*--Sensitivity--*"
		mat rownames `spmcrow' = "*--Specificity--*"
		
		/*
		if "`design'" != "comparative" {
			mat rownames `absoutse' = `byrownames'
			mat rownames `absoutsp' = `byrownames'
			mat rownames `popabsoutse' = `byrownames'
			mat rownames `popabsoutsp' = `byrownames'
			mat rownames `selogodds' = `byrownames'
			mat rownames `splogodds' = `byrownames'
		}
		else {
			mat rownames `absoutse' = `bybirownames'
			mat rownames `absoutsp' = `bybirownames'
			mat rownames `popabsoutse' = `bybirownames'
			mat rownames `popabsoutsp' = `bybirownames'
			mat rownames `selogodds' = `bybirownames'
			mat rownames `splogodds' = `bybirownames'
		}
		*/
		if "`mc'" == "" {
			mat `semc' = `semcrow' \ `semc'
			mat `spmc' = `spmcrow' \ `spmc'	
		}		
		
		mat `selogodds' = `serow' \  `selogodds'
		mat `splogodds' = `sprow' \  `splogodds'
		
		mat `absoutse' = `serow' \  `absoutse'
		mat `absoutsp' = `sprow' \  `absoutsp'
		
		mat `popabsoutse' = `pserow' \  `popabsoutse'
		mat `popabsoutsp' = `psprow' \  `popabsoutsp'
		
		mat `exactabsoutse' = `pserow' \  `exactabsoutse'
		mat `exactabsoutsp' = `psprow' \  `exactabsoutsp'

		mat `absout' = `absoutse' \ `absoutsp'
		mat `popabsout' = `popabsoutse' \ `popabsoutsp'
		
		mat colnames `absout' = Estimate SE(logit) z(logit) P>|z| Lower Upper
		mat colnames `absoutse' = Estimate SE(logit) z(logit) P>|z| Lower Upper
		mat colnames `absoutsp' = Estimate SE(logit) z(logit) P>|z| Lower Upper
		
		mat colnames `popabsout' = Mean SE Median Lower Upper
		mat colnames `popabsoutse' = Mean SE Median Lower Upper
		mat colnames `popabsoutsp' = Mean SE Median Lower Upper
		mat colnames `exactabsoutse' = Mean SE Median Lower Upper
		mat colnames `exactabsoutsp' = Mean SE Median Lower Upper
					
		if "`rr'" == "" {
			mat `serrow' = J(1, 6, .)
			mat `sprrow' = J(1, 6, .)
			
			mat `pserrow' = J(1, 5, .)
			mat `psprrow' = J(1, 5, .)
			
			mat rownames `serrow' = "Relative Sensitivity"
			mat rownames `sprrow' = "Relative Specificity"  //20 chars
			mat rownames `pserrow' = "Relative Sensitivity"
			mat rownames `psprrow' = "Relative Specificity"  //20 chars
			/*
			mat rownames `serrout' = `byrownames'
			mat rownames `sprrout' = `byrownames'
			mat rownames `seorout' = `byrownames'
			mat rownames `sporout' = `byrownames'
			
			mat rownames `popserrout' = `byrownames'
			mat rownames `popsprrout' = `byrownames'
			mat rownames `popseorout' = `byrownames'
			mat rownames `popsporout' = `byrownames'
			*/
			mat `serrout' = `serrow' \  `serrout'
			mat `sprrout' = `sprrow' \  `sprrout'
			mat `rrout' = `serrout' \ `sprrout'
			
			mat `seorout' = `serrow' \  `seorout'
			mat `sporout' = `sprrow' \  `sporout'
			mat `orout' = `seorout' \ `sporout'
			
			mat `popserrout' = `pserrow' \  `popserrout'
			mat `popsprrout' = `psprrow' \  `popsprrout'
			mat `poprrout' = `popserrout' \ `popsprrout'
			
			mat `popseorout' = `pserrow' \  `popseorout'
			mat `popsporout' = `psprrow' \  `popsporout'
			mat `poporout' = `popseorout' \ `popsporout'
			
			mat colnames `rrout' = Estimate SE(log) z(log) P>|z| Lower Upper
			mat colnames `serrout' = Estimate SE(log) z(log) P>|z| Lower Upper
			mat colnames `sprrout' = Estimate SE(log) z(log) P>|z| Lower Upper
			
			mat colnames `poprrout' = Mean SE Median Lower Upper
			mat colnames `popserrout' = Mean SE Median Lower Upper
			mat colnames `popsprrout' = Mean SE Median Lower Upper
			
			mat colnames `orout' = Estimate SE(log) z(log) P>|z| Lower Upper
			mat colnames `seorout' = Estimate SE(log) z(log) P>|z| Lower Upper
			mat colnames `sporout' = Estimate SE(log) z(log) P>|z| Lower Upper
			
			mat colnames `poporout' = Mean SE Median Lower Upper
			mat colnames `popseorout' = Mean SE Median Lower Upper
			mat colnames `popsporout' = Mean SE Median Lower Upper
			
		}
		
		mat rownames `BVar' = `byrownames'
		mat rownames `refe' = `byrownames'
		*mat coleq `spmc' = `byrownames'
		*mat coleq `semc' = `byrownames'
		mat rownames `matgof' = `byrownames'
		mat colnames `matgof' = AIC BIC
	}
	
	if "`rr'" == "" {
		if "`inltest'" == "yes" {
			mat `nltestrr' = `setestnlrr' \ `sptestnlrr'
			mat `nltestor' = `setestnlor' \ `sptestnlor'
		}
	}

	//CI
	if "`outplot'" != "abs" {
		if "`design'" == "comparative" {
			drop `sp'
			gettoken idpair confounders : regressors
			/*tokenize `regressors'
			macro shift
			local confounders `*'*/
			qui count
			*local Nobs = `=r(N)'*0.25
			local Nobs = r(N)
			cap assert mod(`Nobs', 4) == 0 /*Check if the number of studies is half*/
			if _rc != 0 {
				di as error "Some studies cannot be compared properly"
				exit _rc, STATA
			}
			
			
			sort `se' `regressors' `rid'
			egen `id' = seq(), f(1) t(`=`Nobs'*0.25') b(1) 
			sort `se' `id' `rid'
			by `se' `id': egen `pairid' = seq()
			widesetup `event' `total' `confounders' , varx(`varx') idpair(`pairid') se(`se') sid(`id') `design'
			
			gen `sp' = 1 - `se'
			local vlist = r(vlist)
			local cc0 = r(cc0)
			local cc1 = r(cc1)
				
			if "`refpos'" == "bottom" {
				if "`outplot'" == "rr" {
					koopmanci `event'`=`indexcode'-1' `total'`=`indexcode'-1' `event'`=`basecode'-1' `total'`=`basecode'-1', r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01')
				}
				if "`outplot'" == "or" {
					orccci `event'`=`indexcode'-1' `total'`=`indexcode'-1' `event'`=`basecode'-1' `total'`=`basecode'-1', r(`es') upperci(`uci') lowerci(`lci') level(`level') cimethod(`cimethod')
				}
			}
			else {
				if "`outplot'" == "rr" {
					koopmanci `event'`=`basecode'-1' `total'`=`basecode'-1' `event'`=`indexcode'-1' `total'`=`indexcode'-1', r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01')
				}
				if "`outplot'" == "or" {
					orccci `event'`=`basecode'-1' `total'`=`basecode'-1' `event'`=`indexcode'-1' `total'`=`indexcode'-1', r(`es') upperci(`uci') lowerci(`lci') level(`level') cimethod(`cimethod')			
				}
			}
			
			//Rename the varying columns
			foreach v of local vlist {
				rename `v'0 `v'_`cc0'
				label var `v'_`cc0' "`v'_`cc0'"
				rename `v'1 `v'_`cc1'
				label var `v'_`cc1' "`v'_`cc1'"
			}
			
			//Add the weights
			qui {
				cap confirm variable _WT_`cc0'
				if _rc == 0 {
					gen _WT = _WT_`cc0' + _WT_`cc1'
					drop _WT_`cc0'  _WT_`cc1'
				}
				else {
					replace _WT = 2*_WT //Same weights
				}
			}
			
			//Remove unnecessary columns
			if "`smooth'" !="" {
				qui drop `modeles'_`cc0'
				rename `modeles'_`cc1' `modeles'
				
				qui drop `modellci'_`cc0'
				rename `modellci'_`cc1' `modellci'
				
				qui drop `modeluci'_`cc0'
				rename `modeluci'_`cc1' `modeluci'
			}
			
			//make new lcols, rcols
			foreach v of local lcols {
				if strpos("`vlist'", "`v'") != 0 {
					local lcols_rr "`lcols_rr' `v'_`cc0' `v'_`cc1'"
				}
				else {
					local lcols_rr "`lcols_rr' `v'"
				}
			}
			local lcols "`lcols_rr'"
			
			//make new depvars
			local depvars_rr 
			
			foreach v of local depvars {
				if strpos("`vlist'", "`v'") != 0 {
					local depvars_rr "`depvars_rr' `v'_`cc0' `v'_`cc1'"
				}
				else {
					local depvars_rr "`depvars_rr' `v'"
				}
			}
			local depvars "`depvars_rr'"
			
			//make new indvars
			local indvars_rr 
			
			foreach v of local indvars {
				if strpos("`vlist'", "`v'") != 0 {
					local indvars_rr "`indvars_rr' `v'_`cc0' `v'_`cc1'"
				}
				else {
					local indvars_rr "`indvars_rr' `v'"
				}
			}
			local regressors "`indvars_rr'"
		}
		if "`design'" == "cbnetwork" {
			sort `rowid' `se' `rid'
			drop `sp' `rid' `predevent' `iw'
			
			qui reshape wide `event' `total' `ipair' `assignment' _WT `modeles' `modellci' `modeluci', i(`rowid' `se') j(`idpair')	
			
			//Add the weights
			qui {
				cap drop _WT
				gen _WT = _WT0 + _WT1
				drop _WT0  _WT1
			}
			
			//Remove unnecessary columns
			if "`smooth'" !="" {
				qui drop `modeles'0
				rename `modeles'1 `modeles'
				
				qui drop `modellci'0
				rename `modellci'1 `modellci'
				
				qui drop `modeluci'0
				rename `modeluci'1 `modeluci'
			}
			
			if "`outplot'" == "rr" {
				koopmanci `event'1 `total'1 `event'0 `total'0, r(`es') upperci(`uci') lowerci(`lci') alpha(`=1 - `level'*0.01')
			}
			if "`outplot'" == "or" {
				orccci`event'1 `total'1 `event'0 `total'0, r(`es') upperci(`uci') lowerci(`lci') level(`level') cimethod(`cimethod')
			}
			gen `id' = `rowid'
		}
		if "`design'" == "abnetwork" {
			qui {
				gen `id' = _n
				gen `es' = .
				gen `lci' = .
				gen `uci' = .
			}
		}
	}
	else {
		metadta_propci `total' `event', p(`es') lowerci(`lci') upperci(`uci') cimethod(`cimethod') level(`level')
		sort `rid' `se'
		gen `id' = _n
	}
	forvalues l = 1(1)6 {
		local S_`l'1 = .
		local S_`l'2 = .
	}

	//===================================================================================
	//Prepare data for display
	gen `use' = 1  //Individual studies
	
	prep4show `id' `se' `use' `neolabel' `es' `lci' `uci' `modeles' `modellci' `modeluci', ///
		sortby(`sortby') groupvar(`groupvar') grptotal(`grptotal') 	cveffect(`cveffect') ///
		outplot(`outplot') serrout(`serrout') popserrout(`popserrout') popseorout(`popseorout') ///
		absoutse(`absoutse') absoutsp(`absoutsp') popabsoutse(`popabsoutse') popabsoutsp(`popabsoutsp') exactabsoutse(`exactabsoutse') exactabsoutsp(`exactabsoutsp')    ///
		popsprrout(`popsprrout') sprrout(`sprrout') popsporout(`popsporout') `subgroup' `summaryonly' `stratify' ///
		`overall' download(`download') indvars(`regressors') depvars(`depvars') `design' `optimize'
		

	//Extra tables
	//het
	if "`model'" == "random" & "`htable'" == "" {			
		printmat, matrixout(`BVar') type(bhet) dp(`dp') `design'  p(`p')
	}
	//re vs fe
	if "`model'" == "random" & "`htable'" == "" {			
		printmat, matrixout(`refe') type(refe) dp(`dp') `design' 
	}
	
	//Display GOF
	if ("`gof'" != "") {
		printmat, matrixout(`matgof') type(gof) dp(`dp') 
	}
	//logodds
	/*if  (("`sumtable'" == "all") |(strpos("`sumtable'", "logit") != 0)) {
		printmat, matrixout(`logodds') type(logit) dp(`dp') power(`power') `continuous' cveffect(`cveffect')
	}*/
	
	//abs
	if  (("`sumtable'" == "all") |(strpos("`sumtable'", "abs") != 0)) {
		printmat, matrixout(`absout') type(abs) dp(`dp') power(`power') `continuous' cveffect(`cveffect') model(`model') 

		printmat, matrixout(`popabsout') type(popabs) dp(`dp') power(`power')  cveffect(`cveffect') nsims(`nsims') level(`level')
	}
	
	//rr
	if (("`sumtable'" == "all") |(strpos("`sumtable'", "rr") != 0)) & (`pcat' > 0 | "`typevarx'" == "i") {
	*if (("`sumtable'" == "all" & `pcat' > 0) | (strpos("`sumtable'", "rr") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))   {
		//rr
		printmat, matrixout(`rrout') type(rr) dp(`dp') power(`power') cveffect(`cveffect') model(`model') 
		
		//rr equal
		if "`inltest'" == "yes" {
			printmat, matrixout(`nltestrr') type(rre) dp(`dp') cveffect(`cveffect')
		}
		
		//pop rr
		printmat, matrixout(`poprrout') type(poprr) dp(`dp') power(`power') cveffect(`cveffect') nsims(`nsims') level(`level')
	}
	//or
	if (("`sumtable'" == "all") |(strpos("`sumtable'", "or") != 0)) & (`pcat' > 0 | "`typevarx'" == "i") {
	*if (("`sumtable'" == "all" & `pcat' > 0) | (strpos("`sumtable'", "rr") != 0)) & (("`catreg'" != " ") | ("`typevarx'" == "i"))   {
		//or
		printmat, matrixout(`orout') type(or) dp(`dp') power(`power') cveffect(`cveffect') model(`model') 
		
		//rr equal
		if "`inltest'" == "yes" {
			printmat, matrixout(`nltestor') type(rre) dp(`dp') cveffect(`cveffect')
		}

		//pop or
		printmat, matrixout(`poporout') type(popor) dp(`dp') power(`power') cveffect(`cveffect') nsims(`nsims') level(`level')		
	}	
	//model comparison
	if ((`p' > 0 & "`design'" != "abnetwork") | (`p' > 1 & "`design'" == "abnetwork")) & ("`mc'" =="") {
		printmat, matrixoutse(`semc') matrixoutsp(`spmc') type(mc) dp(`dp') cveffect(`cveffect')
	}
	/*
	//Display heterogeneity
	if "`model'" == "random" & "`htable'" == "" {
		*disphetab, `htable' dp(`dp') /*isq2(`Isq2')*/ bshet(`bshet')  bvar(`BVari') wvar(`WVari') p(`p')  
	}
	*/
	//Display the studies
	if "`itable'" == "" {
		disptab `id' `se' `use' `neolabel' `es' `lci' `uci' `grptotal' `modeles' `modellci' `modeluci', `itable' dp(`dp') power(`power') ///
			`subgroup' `overall' `summaryonly' sumstatse(`sumstatse') sumstatsp(`sumstatsp')  	///
			isq2(`isq2') bghet(`bghet') bshet(`bshet') model(`model') bvar(`BVar') 	///
			catreg(`catreg') outplot(`outplot') interaction(`interaction') `smooth' `wt' ///
			se_lrtest(`se_lrtest') sp_lrtest(`sp_lrtest') p(`p') `mc' `design' level(`level') cveffect(`cveffect') 
	}
	
	//Draw the forestplot
	if "`fplot'" == "" {
		fplot `es' `lci' `uci' `use' `neolabel' `grptotal' `id' `se' `modeles' `modellci' `modeluci', ///	
			studyid(`studyid') power(`power') dp(`dp') level(`level') cveffect(`cveffect')  ///
			groupvar(`groupvar') `smooth' `wt' `box'  `summaryonly' ///
			outplot(`outplot') lcols(`lcols') `foptions' `design' varxlabs(`varxlabs')
	}
	
	//Draw the SROC curve
	if "`outplot'" != "abs"  /*| ("`stratify'" != "" & `p' > 0)*/ {
		local sroc "nosroc"
	}
	if "`sroc'" == "" {	
		if "`groupvar'" == "" & `p' > 0  {
			di as res "NOTE: SROC presented for the overall mean."
		}
		use "`master'", clear
		sroc `varlist',  model(`model') selogodds(`selogodds') splogodds(`splogodds') cveffect(`cveffect') ///
			popabsoutsp(`popabsoutsp') popabsoutse(`popabsoutse') v(`V') bvar(`BVar') ///
			groupvar(`groupvar') cimethod(`cimethod') level(`level') p(`p') `soptions' `stratify' `summaryonly'
	}
	
	if "`tradeoff'" != "" {		
		use "`master'", clear
		tradeoff `varlist', sevar(`se') spvar(`sp') cutmat(`cutmat') coefmat(`coefmat') ///
			groupvar(`groupvar') p(`p') `soptions' `stratify'
	}
																																
	cap ereturn clear
	if ((`p' > 0 & "`design'" != "abnetwork") |(`p' > 1 & "`design'" == "abnetwork")) & "`mc'" == "" {
		cap confirm matrix `semc'
		if _rc == 0 {
			ereturn matrix mctestse = `semc' //model comparison se
		}
		cap confirm matrix `pemc'
		if _rc == 0 {
			ereturn matrix mctestsp = `spmc' //model comparison sp
		}
		*ereturn matrix bghet = `bghet' //Full vs Null model
	}
	/*if `p' == 0 & "`cbnetwork'" == "" {
		ereturn matrix isq2 = `isq2' //isq2
	}*/
	/*if "`inltest'" == "yes" {
		ereturn matrix setestnl = `setestnl' //Equality of RR - se
		ereturn matrix sptestnl = `sptestnl'
	}*/	
	cap confirm matrix `logodds'
	if _rc == 0 {
		ereturn matrix logodds = `logodds' //logodds se and sp
		ereturn matrix Vlogodds = `V' //var-cov for catreg & overall 
	}
	cap confirm matrix `absout'
	if _rc == 0 {
		ereturn matrix absout = `absout'
		ereturn matrix absoutse = `absoutse'
		ereturn matrix absoutsp = `absoutsp'
	}
	cap confirm matrix `popabsout'
	if _rc == 0 {
		ereturn matrix popabsout = `popabsout'
		ereturn matrix popabsoutse = `popabsoutse'
		ereturn matrix popabsoutsp = `popabsoutsp'
	}
	cap confirm matrix `rrout'
	if _rc == 0 {
		ereturn matrix rrout = `rrout'
		ereturn matrix serrout = `serrout'
		ereturn matrix sprrout = `sprrout'
	}
	cap confirm matrix `poprrout'
	if _rc == 0 {
		ereturn matrix poprrout = `poprrout'
		ereturn matrix popserrout = `popserrout'
		ereturn matrix popsprrout = `popsprrout'
	}
	cap confirm matrix `orout'
	if _rc == 0 {
		ereturn matrix orout = `orout'
		ereturn matrix seorout = `seorout'
		ereturn matrix sporout = `sporout'
	}
	cap confirm matrix `nltestrr'
	if _rc == 0 {
		ereturn matrix rrtest = `nltestrr'
		ereturn matrix ortest = `nltestor'
	}	
	cap confirm matrix `poporout'
	if _rc == 0 {
		ereturn matrix poporout = `poporout'
		ereturn matrix popseorout = `popseorout'
		ereturn matrix popsporout = `popsporout'
	}
	ereturn matrix refe = `refe' //Re vs FE 
	ereturn matrix vcovar = `BVar' //var-covar between logit se and logit sp
	cap confirm matrix `WVar'
	if _rc == 0 {
		ereturn matrix wvar = `WVar' //2nd variance
	}
	
	restore 
end

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: INDEX +++++++++++++++++++++++++
							Find index of word in a string
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

cap program drop index
program define index, rclass
version 14.0

	syntax, source(string asis) word(string asis)
	local nwords: word count `source'
	local found 0
	local index 1

	while (!`found') & (`index' <= `nwords'){
		local iword:word `index' of `source'
		if "`iword'" == `word' {
			local found 1
		}
		local index = `index' + 1
	}
	
	if `found' {
		local index = `index' - 1
	}
	else{
		local index = 0
	}
	return local index `index'
end

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: myncod +++++++++++++++++++++++++
								Decode by order of data
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/	
cap program drop my_ncod
program define my_ncod
version 14.1

	syntax newvarname(gen), oldvar(varname)
	
	qui {
		cap confirm numeric var `oldvar'
		tempvar by_num 
		
		if _rc == 0 {				
			decode `oldvar', gen(`by_num')
			drop `oldvar'
			rename `by_num' `oldvar'
		}

		* The _by variable is generated according to the original
		* sort order of the data, and not done alpha-numerically

		qui count
		local N = r(N)
		cap drop `varlist'
		gen `varlist' = 1 in 1
		local lab = `oldvar'[1]
		cap label drop `oldvar'
		if "`lab'" != ""{
			label define `oldvar' 1 "`lab'"
		}
		local found1 "`lab'"
		local max = 1
		forvalues i = 2/`N'{
			local thisval = `oldvar'[`i']
			local already = 0
			forvalues j = 1/`max'{
				if "`thisval'" == "`found`j''"{
					local already = `j'
				}
			}
			if `already' > 0{
				replace `varlist' = `already' in `i'
			}
			else{
				local max = `max' + 1
				replace `varlist' = `max' in `i'
				local lab = `oldvar'[`i']
				if "`lab'" != ""{
					label define `oldvar' `max' "`lab'", modify
				}
				local found`max' "`lab'"
			}
		}

		label values `varlist' `oldvar'
		label copy `oldvar' `varlist', replace
		
	}
end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: FITMODEL +++++++++++++++++++
							Fit the logistic model
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop fitmodel
program define fitmodel, rclass
version 14.0

	syntax varlist [if], [ bcov(string) wcov(string) model(string) modelopts(string asis) regexpression(string) sid(varname) p(string) ///
		nested(varname) comparative ipair(varname) level(integer 95) abnetwork cbnetwork general comparative]
		tokenize `varlist'	
		
		marksample touse, strok
		
		tempname varcov coefs initmat
		
		local exactest = 0
		
		if "`wcov'" != "" {
			local nested = `"|| (`nested'`ipair': `3' `4', noc cov(`wcov'))"'
		}
		else {
			local nested
		}
		
		if ("`bcov'" == "se") | ("`bcov'" == "sp"){
			if ("`bcov'" == "se") {
				local re = "`3'"
			}
			else {
				local re = "`4'"
			}
			local cov
		}
		else {
			local re = "`3' `4'"
			local cov = "cov(`bcov')"
		}
		local try = 0
	
		if ("`model'" == "fixed") {
			local ++try
			capture noisily binreg `1' `regexpression' if `touse', noconstant n(`2') ml `modelopts' l(`level')
			local success = _rc
		}
				
		if ("`model'" == "random") {
			//Specify the engine
			if _caller() >= 16 {
				local fitcommand "melogit"
			}
			else {
				local fitcommand "meqrlogit"
			}
		
			if strpos(`"`modelopts'"', "(iterate") == 0  {
				if "`fitcommand'" == "meqrlogit" {
					local iterate = `"iterate(30)"'
				}
				else {
					local iterate = `"iterate(100)"'
				}
			}
			if strpos(`"`modelopts'"', "intpoi") == 0  {
				qui count if `touse'
				if `=r(N)' < 7 {
					local intopts = `"intpoints(`=r(N)')"'
				}
			}
			
			//First trial
			local ++try
			#delim ;
			capture noisily `fitcommand' (`1' `regexpression' if `touse', noc )|| 
			  (`sid': `re', noc `cov') `nested',
			  binomial(`2') `modelopts' `intopts' `iterate'  l(`level');
			#delimit cr 
			
			local success = _rc
			local converged = e(converged)
			
			if (strpos(`"`modelopts'"', "from") == 0) & ((`converged' == 0) | (`success' != 0))  {
				//First fit fixed effects model to get better starting values
				noi di _n"*********************************** ************* ***************************************" 
				noi di as txt _n "Just a moment - Obtaining better initial values "
				noi di   "*********************************** ************* ***************************************" 
				local ++try
				capture noisily binreg `1' `regexpression' if `touse', noconstant n(`2') ml  
				local success = _rc
				if _rc == 0 {				
					qui estimates table
					mat `coefs' = r(coef)
					local initialized = 0
					local ninits = rowsof(`coefs')
					forvalues e = 1(1)`ninits' {
						local init = `coefs'[`e', 1]
						if `init' != .b {
							local value "`init'"
							if `initialized' == 0 {
								local inits = `"`value'"'
								local initialized = 1
							}
							else {
								local inits = `"`inits', `value'"'
							}
						}
					}
				}
				
				local b = 0
				local w = 0	
				if "`wcov'" != "" {
					if strpos("`wcov'", "ind") != 0 {
						local w = 2
					}
					else if strpos("`wcov'", "id") != 0 {
						local w = 1
					}
				}
				if strpos("`bcov'", "uns") != 0 {
					local b = 3
				}
				else if (strpos("`bcov'", "ind") != 0) | (strpos("`bcov'", "exc") != 0) {
					local b = 2
				}
				else if (strpos("`bcov'", "id") != 0) | (strpos("`bcov'", "sp") != 0) | (strpos("`bcov'", "se") != 0)   {
					local b = 1
				}
				
				mat `varcov' = J(1, `b', 0)					
				mat `initmat' = (`inits', `varcov')
				
				local inits = `"from(`initmat', copy)"'
				
				//Laplace first without nested
				local ++try
				#delim ;
				capture noisily `fitcommand' (`1' `regexpression' if `touse', noc )|| 
				  (`sid': `re', noc `cov'),
				  binomial(`2') `inits' `intopts' `iterate'  l(`level') matlog laplace;
				#delimit cr 
				
				local success = _rc
				local converged = e(converged)
				if _rc == 0 {
					//Get the new estimates
					qui estimates table
					mat `coefs' = r(coef)
					local initialized = 0
					local ninits = rowsof(`coefs')
					forvalues e = 1(1)`ninits' {
						local init = `coefs'[`e', 1]
						if `init' != .b {
							local value "`init'"
							if `initialized' == 0 {
								local inits = `"`value'"'
								local initialized = 1
							}
							else {
								local inits = `"`inits', `value'"'
							}
						}
					}
				}
				
				if `w' > 0 {
					mat `varcov' = J(1, `w', 0)					
					mat `initmat' = (`inits', `varcov')
					
					local inits = `"from(`initmat', copy)"'
					//fit once more with nested
					local ++try
					#delim ;
					capture noisily `fitcommand' (`1' `regexpression' if `touse', noc )|| 
					  (`sid': `re', noc `cov') `nested',
					  binomial(`2') `inits' `intopts' `iterate'  l(`level') matlog laplace;
					#delimit cr 
					
					local success = _rc
					local converged = e(converged)
					
					if _rc == 0 {
						//Get the new estimates
						qui estimates table
						mat `coefs' = r(coef)
						local initialized = 0
						local ninits = rowsof(`coefs')
						forvalues e = 1(1)`ninits' {
							local init = `coefs'[`e', 1]
							if `init' != .b {
								local value "`init'"
								if `initialized' == 0 {
									local inits = `"`value'"'
									local initialized = 1
								}
								else {
									local inits = `"`inits', `value'"'
								}
							}
						}
					}
					
				}
				
				mat `initmat' = (`inits')	
				local inits = `"from(`initmat', copy)"'

				//Final fit
				local ++try
				#delim ;
				capture noisily `fitcommand' (`1' `regexpression' if `touse', noc )|| 
				  (`sid': `re', noc `cov') `nested',
				  binomial(`2') `inits' `modelopts' `intopts' `iterate'  l(`level');
				#delimit cr 
				
				local success = _rc
				local converged = e(converged)
			}
			*==========================================================================
			if "`fitcommand'" == "meqrlogit" {				
				//Try to refineopts 2 times
				if strpos(`"`modelopts'"', "refineopts") == 0 {				
					if (`try' < 7) & ((`converged' == 0) | (`success' != 0)) {
						local ++try 
						local refine "refineopts(iterate(`=10 * `try''))"
						#delim ;					
						capture noisily `fitcommand' (`1' `regexpression' if `touse', noc )|| 	
								(`sid': `re', noc `cov') `nested',						
								binomial(`2') `modelopts' l(`level') `refine' `intopts' `iterate' ;
						#delimit cr 
						
						local success = _rc
						
						local converged = e(converged) 
					}
				}
			}
			*==========================================================================
			/*if "`fitcommand'" == "meqrlogit" {
				*Try matlog if still difficult
				if (strpos(`"`modelopts'"', "matlog") == 0) {
					local matlog "matlog"
				}
			}
			else {
				*Try difficult if still difficult
				if (strpos(`"`modelopts'"', "difficult") == 0) {
					local difficult "difficult"
				}
				if strpos(`"`modelopts'"', "(iterate") == 0 {
					local iterate = `"iterate(200)"'
				}
			}
			if ((`converged' == 0) | (`success' != 0)) {
				local ++try 
				#delim ;
				capture noisily `fitcommand' (`1' `regexpression' if `touse', noc )|| 
					(`sid': `re', noc `cov') `nested',
					binomial(`2') `modelopts' l(`level') `refine' `matlog' `difficult' `intopts' `iterate' ;
				#delimit cr
				
				local success = _rc				
				local converged = e(converged)
			}*/
			
		}
		//Revert to FE if ME fails
		if (`success' != 0) & ("`model'" == "random") {
			capture noisily binreg `1' `regexpression' if `touse', noconstant n(`2') ml  l(`level')
			local success = _rc
			local model "fixed"
		}		
		
		if `success' != 0 {
			*display as error "Unexpected error performing regression"
			di as err "Model could not converge after `try' attempts with different options. Try fitting a simpler model"
            exit `success'
		}
		
		return local model =  "`model'"
end
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: metadta_PROPCI +++++++++++++++++++++++++
								CI for proportions
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop metadta_propci
	program define metadta_propci
	version 14.1

		syntax varlist [if] [in], p(name) lowerci(name) upperci(name) [cimethod(string) level(real 95)]
		
		qui {	
			tokenize `varlist'
			gen `p' = .
			gen `lowerci' = .
			gen `upperci' = .
			
			count `if' `in'
			forvalues i = 1/`r(N)' {
				local N = `1'[`i']
				local n = `2'[`i']

				cii proportions `N' `n', `cimethod' level(`level')
				
				replace `p' = r(proportion) in `i'
				replace `lowerci' = r(lb) in `i'
				replace `upperci' = r(ub) in `i'
			}
		}
	end
	
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: absexactci +++++++++++++++++++++++++
								CI for proportions
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop metadta_absexactci
	program define metadta_absexactci, rclass
	version 14.1
		syntax anything(name=data id="data"), [level(real 95)]
		
		tempname absexact
		mat `absexact' = J(1, 6, .)
		
		local len: word count `data'
		if `len' != 2 {
			di as error "Specify full data: N n"
			exit
		}
		
		foreach num of local data {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		tokenize `data'
		cap assert (`2' <= `1') 
		if _rc != 0 {
			di as err "Order should be N n"
			exit _rc
		}
		
		cii proportions `1' `2', `cimethod' level(`level')
		
		mat `absexact'[1, 1] = r(proportion) 
		mat `absexact'[1, 2] = r(se)
		mat `absexact'[1, 5] = r(lb) 
		mat `absexact'[1, 6] = r(ub)
		
		local zvalue = (`absexact'[1, 1] - 0.5)/sqrt(0.25/`1')
		mat `absexact'[1, 3] = `zvalue'
		
		local pvalue = normprob(-abs(`zvalue'))*2
		mat `absexact'[1, 4] = `pvalue'
		
		return matrix absexact = `absexact'
	end
	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: WIDESETUP +++++++++++++++++++++++++
							Transform data to wide format
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop widesetup
	program define widesetup, rclass
	version 14.1

	syntax varlist, [sid(varlist) varx(varname) idpair(varname) se(varname) comparative sortby(varlist) cbnetwork rowid(varname) index(varname) assignment(varname) comparator(varname) ]

		qui{
			tokenize `varlist'
			local event = "`1'"
			local total = "`2'"

			tempvar jvar modey diffy
		
			if "`cbnetwork'" == "" {
				gen `jvar' = `idpair' - 1
				
				/*Check for varying variable and store them*/
				ds
				local vnames = r(varlist)
				local vlist
				foreach v of local vnames {	
					cap drop `modey' `diffy'
					bysort `sid': egen `modey' = mode(`v'), minmode
					egen `diffy' = diff(`v' `modey')
					sum `diffy'
					local sumy = r(sum)
					if (strpos(`"`varlist'"', "`v'") == 0) & (`sumy' > 0) & "`v'" != "`jvar'" & "`v'" != "`idpair'" & "`v'" != "`varx'"{
						if "`se'" != "" & "`v'" == "`se'"{
							local v
						}
						local vlist "`vlist' `v'"
					}
				}
				cap drop `modey' `diffy'
				
				sort `sid' `jvar' `sortby'
				
				/*2 variables per study : n N*/			
				reshape wide `event' `total'  `idpair' `varx' `vlist', i(`sid' `se') j(`jvar')
			}
			if "`cbnetwork'" != "" { 
				reshape wide `varlist', i(`sid') j(`se')
				
				drop `sid'
				
				reshape wide `event'0 `total'0 `event'1 `total'1 `index' `assignment', i(`rowid') j(`comparator')
			
			}
			if "`varx'" != "" {			
				tempvar holder
				decode `varx'0, gen(`holder')
				local cc0 = `holder'[1]
				drop `holder'
				
				decode `varx'1, gen(`holder')
				local cc1 = `holder'[1]

			}

			return local vlist = "`vlist'"
			return local cc0 = "`cc0'"
			return local cc1 = "`cc1'"
		}
	end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: PREP4SHOW +++++++++++++++++++++++++
							Prepare data for display table and graph
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop prep4show
program define prep4show
version 14.0

	#delimit ;
	syntax varlist, [popserrout(name) popsprrout(name) serrout(name) sprrout(name) popseorout(name) popsporout(name) 
		absoutse(name) absoutsp(name) popabsoutse(name) popabsoutsp(name) exactabsoutse(name) exactabsoutsp(name) sortby(varlist) 
		groupvar(varname) summaryonly nooverall nosubgroup outplot(string) grptotal(name) download(string asis) 
		indvars(varlist) depvars(varlist) comparative abnetwork general cbnetwork stratify cveffect(string) optimize ] 
	;
	#delimit cr
	tempvar sp  expand 
	tokenize `varlist'
	 
	local id = "`1'"
	local se = "`2'" 
	local use = "`3'"
	local label = "`4'"
	local es = "`5'"
	local lci = "`6'"
	local uci = "`7'"
	local modeles = "`8'"
	local modellci = "`9'"
	local modeluci = "`10'"
	qui {
		gen `sp' = 1 - `se'

		gen `expand' = 1

		//Groups
		if "`groupvar'" != "" {
			bys `groupvar' `se' : egen `grptotal' = count(`id') //# studies in each group
			gsort `groupvar' `se' `sortby' `id'
			bys `groupvar' `se' : replace `expand' = 1 + 1*(_n==1) + 2*(_n==_N)
			expand `expand'
			gsort `groupvar' `se' `sortby' `id' `expand'
			bys `groupvar' `se' : replace `use' = -2 if _n==1  //group label
			bys `groupvar' `se' : replace `use' = 2 if _n==_N-1  //subgroup
			bys `groupvar' `se' : replace `use' = 0 if _n==_N //blank
			replace `id' = `id' + 1 if `use' == 1
			replace `id' = `id' + 2 if `use' == 2  //subgroup
			replace `id' = `id' + 3 if `use' == 0 //blank
			replace `label' = "Group Mean" if `use' == 2
			replace _WT = . if `use' == 2
			replace `es' = . if `use' == 2
			replace `lci' = . if `use' == 2
			replace `uci' = . if `use' == 2
			replace `modeles' = . if `use' == 2
			replace `modellci' = . if `use' == 2
			replace `modeluci' = . if `use' == 2
			
			qui label list `groupvar'
			local nlevels = r(max)
			local c = 0
			local m = 1
			
			forvalues l = 1/`nlevels' {
				if "`outplot'" == "abs" {
					if  "`general'`comparative'" != "" {
							local c = 1
						}
					if ("`stratify'" !="") {
						local c = 1
						if ("`comparative'" !="") {
							local m = 3
						}
					}
					//ES
					local S_112 = `popabsoutse'[`=`l'*`m' + `c'', 1]
					local S_122 = `popabsoutsp'[`=`l'*`m' + `c'', 1]
					
					//lci
					local S_312 = `popabsoutse'[`=`l'*`m' + `c'', 4]
					local S_322 = `popabsoutsp'[`=`l'*`m' + `c'', 4]
					
					//uci
					local S_412 = `popabsoutse'[`=`l'*`m' + `c'', 5]
					local S_422 = `popabsoutsp'[`=`l'*`m' + `c'', 5]
					
					//Get exact CI
					local E_312 = `exactabsoutse'[`=`l'*`m' + `c'', 4]
					local E_322 = `exactabsoutsp'[`=`l'*`m' + `c'', 4]
					
					local E_412 = `exactabsoutse'[`=`l'*`m' + `c'', 5]
					local E_422 = `exactabsoutsp'[`=`l'*`m' + `c'', 5]
					
					//Get Conditional CI
					//lci
					local C_312 = `absoutse'[`=`l'*`m' + `c'', 5]
					local C_322 = `absoutsp'[`=`l'*`m' + `c'', 5]
					
					//uci
					local C_412 = `absoutse'[`=`l'*`m' + `c'', 6]
					local C_422 = `absoutsp'[`=`l'*`m' + `c'', 6]
					
					if "`optimize'" != "" {
						//if simulated FE variance 5 times larger replace with exact
						//se
						if (`=(`S_412' - `S_312')/(`E_412' - `E_312')' > 5) & (`E_412' != .) & (`E_312' != .)  {
							local S_312 = `E_312'
							local S_412 = `E_412'
						}
						
						//sp
						if (`=(`S_422' - `S_322')/(`E_422' - `E_322')' > 5) & (`E_422' != .) & (`E_322' != .)  {
							local S_322 = `E_322'
							local S_422 = `E_422'
						}
						//if simulated RE more than 5 times larger, replace with conditional stats 
						//se
						if (`=(`S_412' - `S_312')/(`C_412' - `C_312')' > 5) & (`E_412' == .) & (`E_312' == .)  {
							local S_312 = `C_312'
							local S_412 = `C_412'
						}
						
						//sp
						if (`=(`S_422' - `S_322')/(`C_422' - `C_322')' > 5) & (`E_422' == .) & (`E_322' == .)  {
							local S_322 = `C_322'
							local S_422 = `C_422'
						}
					}	
				}
				else {
					if ("`abnetwork'" !="") | ("`stratify'" !="")  {
						local c = 1
					}
					if "`outplot'" == "or" {
						local S_112 = `popseorout'[`=`l' + `c'', 1]
						local S_122 = `popsporout'[`=`l' + `c'', 1]
						
						local S_312 = `popseorout'[`=`l' + `c'', 4]
						local S_322 = `popsporout'[`=`l' + `c'', 4]
						
						local S_412 = `popseorout'[`=`l' + `c'', 5]
						local S_422 = `popsporout'[`=`l' + `c'', 5]
					}
					else {
						//ES
						local S_112 = `popserrout'[`=`l' + `c'', 1]
						local S_122 = `popsprrout'[`=`l' + `c'', 1]
						
						//lci
						local S_312 = `popserrout'[`=`l' + `c'', 4]
						local S_322 = `popsprrout'[`=`l' + `c'', 4]
						
						//uci
						local S_412 = `popserrout'[`=`l' + `c'', 5]
						local S_422 = `popsprrout'[`=`l' + `c'', 5]
						
						//Conditional stats
						//lci
						local C_312 = `serrout'[`=`l' + `c'', 5]
						local C_322 = `sprrout'[`=`l' + `c'', 5]
						
						//uci
						local C_412 = `serrout'[`=`l' + `c'', 6]
						local C_422 = `sprrout'[`=`l' + `c'', 6]
						
						if "`optimize'" != "" {
							//if simulated more than 5 times larger, replace with conditional stats 
							//se
							if `=(`S_412' - `S_312')/(`C_412' - `C_312')' > 5 {
								local S_312 = `C_312'
								local S_412 = `C_412'
							}
							
							//sp
							if `=(`S_422' - `S_322')/(`C_422' - `C_322')' > 5 {
								local S_322 = `C_322'
								local S_422 = `C_422'
							}
						}
					}
				}
				local lab:label `groupvar' `l'
				replace `label'  = "`lab'" if `use' == -2 & `groupvar' == `l'	
				replace `label'  = "`groupvar' = `lab'" if `use' == -2 & `groupvar' == `l' & "`abnetwork'" == ""
				replace `label'  = "`lab'" if `use' == 2 & `groupvar' == `l'	& "`outplot'" == "rr" & "`abnetwork'" != ""
				if "`cveffect'" == "sesp" {
					replace `es'  = `S_112'*`se' + `S_122'*`sp' if `use' == 2 & `groupvar' == `l'	
					replace `lci' = `S_312'*`se' + `S_322'*`sp' if `use' == 2 & `groupvar' == `l'	
					replace `uci' = `S_412'*`se' + `S_422'*`sp' if `use' == 2 & `groupvar' == `l'
				}
				else if "`cveffect'" == "se" {
					replace `es'  = `S_112' if `use' == 2 & `groupvar' == `l' & `se'
					replace `lci' = `S_312' if `use' == 2 & `groupvar' == `l' & `se'
					replace `uci' = `S_412' if `use' == 2 & `groupvar' == `l' & `se'
				}
				else if "`cveffect'" == "sp" {
					replace `es'  =  `S_122' if `use' == 2 & `groupvar' == `l'	& `sp'
					replace `lci' =  `S_322' if `use' == 2 & `groupvar' == `l' & `sp'	
					replace `uci' =  `S_422' if `use' == 2 & `groupvar' == `l' & `sp'
				}
				
				//Weights
				local parameters "`se' `sp'"
				foreach parm of local parameters {
					sum _WT if `use' == 1 & `groupvar' == `l' & `parm'
					local groupwt = r(sum)
					replace _WT = `groupwt' if `use' == 2 & `groupvar' == `l' & `parm'
				}
			}
		}
		else {
			bys `se' : egen `grptotal' = count(`id') //# studies total
		}
		
		//Overall
		//Overall
		if "`overall'" == "" {	
			gsort  `se' `groupvar' `sortby' `id'
			bys `se' : replace `expand' = 1 + 2*(_n==_N)
			expand `expand'
			gsort  `se' `groupvar' `sortby' `id' `expand'
			bys `se' : replace `use' = 3 if _n==_N-1  //Overall
			bys `se' : replace `use' = 0 if _n==_N //blank
			bys `se' : replace `id' = `id' + 1 if _n==_N-1  //Overall
			bys `se' : replace `id' = `id' + 2 if _n==_N //blank
			//Fill in the right info
			if "`outplot'" == "abs" {
				local senrows = rowsof(`popabsoutse')
				local spnrows = rowsof(`popabsoutsp')
				local S_11 = `popabsoutse'[`senrows', 1] //se
				local S_31 = `popabsoutse'[`senrows', 4] //ll
				local S_41 = `popabsoutse'[`senrows', 5] //ul
				
				local S_12 = `popabsoutsp'[`spnrows', 1] //sp
				local S_32 = `popabsoutsp'[`spnrows', 4] //ll
				local S_42 = `popabsoutsp'[`spnrows', 5] //ul
				}
			else { 
				if "`outplot'" == "or" {
					local senrows = rowsof(`popseorout')
					local spnrows = rowsof(`popsporout')
					local S_11 = `popseorout'[`senrows', 1] //p (se)
					local S_31 = `popseorout'[`senrows', 4] //ll
					local S_41 = `popseorout'[`senrows', 5] //ul
					
					local S_12 = `popsporout'[`spnrows', 1] //p (sp)
					local S_32 = `popsporout'[`spnrows', 4] //ll
					local S_42 = `popsporout'[`spnrows', 5] //ul
				}
				else {
					local senrows = rowsof(`popserrout')
					local spnrows = rowsof(`popsprrout')
					local S_11 = `popserrout'[`senrows', 1] //p (se)
					local S_31 = `popserrout'[`senrows', 4] //ll
					local S_41 = `popserrout'[`senrows', 5] //ul
					
					local S_12 = `popsprrout'[`spnrows', 1] //p (sp)
					local S_32 = `popsprrout'[`spnrows', 4] //ll
					local S_42 = `popsprrout'[`spnrows', 5] //ul
				}
			}
			
			replace `es' = `S_11'*`se' + `S_12'*`sp' if `use' == 3	
			replace `lci' = `S_31'*`se' + `S_32'*`sp' if `use' == 3
			replace `uci' = `S_41'*`se' + `S_42'*`sp' if `use' == 3
			replace `label' = "Population Mean" if `use' == 3
		}
		
		count if `use'==1 & `se'==1
		replace `grptotal' = `=r(N)' if `use'==3
		replace `grptotal' = `=r(N)' if _n==_N
		
		replace `label' = "" if `use' == 0
		replace `es' = . if `use' == 0 | `use' == -2
		replace `lci' = . if `use' == 0 | `use' == -2
		replace `uci' = . if `use' == 0 | `use' == -2
		replace _WT = . if (`use' == 3) & ("`stratify'" != "")
		replace _WT = 100 if (`use' == 3) & ("`stratify'" == "")
		
		gsort `se' `groupvar' `sortby'  `id' 
	}
	qui {
		replace `modeles' = . if `use' != 1
		replace `modellci' = . if `use' != 1
		replace `modeluci' = . if `use' != 1
		replace _WT = . if `use' == 0 | `use' == -2 | `use' == 4
	}
	
	if "`download'" != "" {
		preserve
		qui {
			cap drop _ES _LCI _UCI _USE _LABEL _PARAMETER
			gen _ES = `es'
			gen _LCI = `lci'
			gen _UCI = `uci'
			gen _MES = `modeles'
			gen _MLCI = `modellci'
			gen _MUCI = `modeluci'
			gen _USE = `use'
			gen _LABEL = `label'
			gen _PARAMETER = `se'
			gen _ID = `id'
			
			keep `depvars' `indvars' _ES _LCI _UCI _MES _MLCI _MUCI _ESAMPLE _USE _LABEL _PARAMETER _ID _WT
		}
		di _n "Data saved"
		noi save "`download'", replace
		
		restore
	}
	qui {
		//Drop unnecessary rows
		if "`abnetwork'" == "" | ("`abnetwork'" != "" & "`outplot'" == "abs") {
			drop if (`use' == 2 | `use' == 3 ) & (`grptotal' == 1) //drop summary if 1 study
		}		
		drop if (`use' == 1 & "`summaryonly'" != "" & `grptotal' > 1) 
		
		//Remove redundant info
		replace _WT = . if `use' == 1 & (`grptotal' == 1)
		
		//remove label if summary only
		replace `label' = `label'[_n-1] if (`use' == 2 & "`summaryonly'" != "") 
		
		//Drop unnecessary rows
		drop if (`use' == 2 & "`subgroup'" != "") 
		drop if (`use' == -2 & "`summaryonly'" != "") 
		drop if (`use' == 3 & "`overall'" != "") 
		
		if "`abnetwork'" != "" & "`outplot'" != "abs" {
			drop if `use' == 1 | `use' == -2
			replace `use' = 1 if `use' == 2
		}
		
		gsort `se' `groupvar' `sortby'  `id' 
		bys `se' : replace `id' = _n 
		gsort `id' `se' 
	}
end	

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: DISPHETAB +++++++++++++++++++++++++
							Display table 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop disphetab
program define disphetab
version 14.0
#delimit ;
syntax [, isq2(name) bshet(name) dp(integer 2) bvar(name) wvar(name)  p(integer 0)] ;
	#delimit cr
	
	local rho 		= `bvar'[1, 2]/sqrt(`bvar'[1, 1]*`bvar'[2, 2])
	local tau2se 	= `bvar'[1, 1]
	local tau2sp	= `bvar'[2, 2]
	local tau2g		= (1 - (`bvar'[1, 2]/sqrt(`bvar'[1, 1]*`bvar'[2, 2]))^2)*`bvar'[1, 1]*`bvar'[2, 2]
	di as txt "Between-study heterogeneity" 
	di as txt _col(28) "covar" _cont
	di as res _n _col(28) %5.`=`dp''f `covar' 
	
	di as txt _col(28) "rho" _cont
	di as res _n _col(28) %5.`=`dp''f `rho' 
	
	di as txt  _col(28) "Tau.sq" _cont
	if `p' == 0  {
		di as txt _col(45) "I^2(%)" _cont
		local isq2b  = `isq2'[1, 1]*100
		local isq2se = `isq2'[1, 2]*100
		local isq2sp = `isq2'[1, 3]*100
	}			
	di as txt _n  "Generalized" _cont	
	di as res   _col(28) %5.`=`dp''f `tau2g' _col(45) %5.`=`dp''f `isq2b'  
	di as txt  "Sensitivity" _cont	
	di as res    _col(28) %5.`=`dp''f `tau2se' _col(45) %5.`=`dp''f `isq2se'  
	di as txt  "Specificity" _cont
	di as res    _col(28) %5.`=`dp''f `tau2sp' _col(45) %5.`=`dp''f `isq2sp'

	di as txt  _col(30) "Chi2"  _skip(8) "degrees of" _cont
	di as txt _n  _col(28) "statistic" 	_skip(6) "freedom"      _skip(8)"p-val"   _cont
	
	local chisq = `bshet'[1, 1]
	local df 	= `bshet'[1, 2]
	local pv 	= `bshet'[1, 3]	
			
	di as txt _n "LR Test: RE vs FE model" _cont
	di as res _col(25) %10.`=`dp''f `chisq' _col(45) `df' _col(52) %10.4f `pv'  
	
end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: DISPTAB +++++++++++++++++++++++++
							Display table 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop disptab
program define disptab
version 14.0
	#delimit ;
	syntax varlist, [nosubgroup nooverall level(integer 95) sumstatse(string asis) 
	sumstatsp(string asis) noitable dp(integer 2) power(integer 0) isq2(name) smooth noWT
	bghet(name) bshet(name) model(string) bvar(name) catreg(string) outplot(string) cveffect(string)
	interaction(string) se_lrtest(name) sp_lrtest(name) p(integer 0) noMC general cbnetwork comparative abnetwork]
	;
	#delimit cr
	
	tempvar id se use label es lci uci df modeles modellci modeluci
	tokenize `varlist'
	qui {
		gen `id' = `1'
		gen `se' = `2' 
		gen `use' = `3'
		gen `label' = `4'
		gen `es' = `5'
		gen `lci' = `6'
		gen `uci' = `7'
		gen `df' = `8'
		if "`smooth'" != "" {
			gen `modeles' = `9'
			gen `modellci' = `10'
			gen `modeluci' = `11'
		}
		
		if "`wt'" =="" {
			local dispwt "% Weight"
		}
		
		
		if "`outplot'" == "abs" {
			local sumstat "Absolute Estimates"
		}
		else {
			local sumstat "Relative Estimates"
		}
		
		//Find the length of the estimates
		tempvar hold holdstr slimest
		gen `hold' = `uci'*(10^`power')
		tostring `hold', gen(`holdstr') format(%10.`dp'f) force
		gen `slimest' = strlen(strltrim(`holdstr'))
		sum `slimest'
		local est_i_len = r(max)
		
		local togglesmooth = 0	
		
		if "`smooth'" !="" {
			local togglesmooth = 1
			local open " ("
			local close ")"
			local aesclose "%1s"
			local aes "%`=`est_i_len''.`=`dp''f"
		}
	}
		
	preserve
	
		tempvar tlabellen 
		//study label
		local studylb: variable label `label'
		if "`studylb'" == "" {
			local studylb "Study"
		}
		qui {
			replace `se' = `se' + 1
			widesetup `label', sid(`id') idpair(`se')
			
			gen `tlabellen' = strlen(`label'0)
			summ `tlabellen'
		}
		local nlen = r(max) + 5 
		local nlense = strlen("`sumstatse'")
		local nlensp = strlen("`sumstatsp'")
		local smoothlen = 0
		if "`smooth'" != "" {
			local smoothtext = " Observed (Model-based)"
			local smoothlen = 4
		}
		local togglewt = 1
		if "`wt'" != "" {
			local togglewt = 0
		}
		di as res  "****************************************************************************************"
		di as res "{pmore2} Study specific test accuracy: `smoothtext' `sumstat'  {p_end} "
		di as res    "****************************************************************************************" 
		
		di _n as txt _col(`nlen') "| "   _skip(`=15 - round(`nlense'/2)') "`sumstatse'" ///
				  _col(`=`nlen' + 1 + 34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') 	///
				  "| " _skip(`=15 - round(`nlensp'/2)') "`sumstatsp'" _cont
				  
		di  _n  as txt _col(2) "`studylb'" _col(`nlen') "| "   "Estimate" ///
				  _skip(`=5 - 2*`togglewt' + `est_i_len'*2*`togglesmooth'') "[`level'% Conf. Interval]"  _skip(`=3*`togglewt' + `est_i_len'*`togglesmooth'') "`dispwt'" ///
				  "| " "Estimate" ///
				  _skip(`=5 - 2*`togglewt'+ `est_i_len'*2*`togglesmooth'') "[`level'% Conf. Interval]" _skip(`=3*`togglewt'+ `est_i_len'*`togglesmooth'') "`dispwt'"
				  
		di  _dup(`=`nlen'-1') "-" "+" _dup(`=34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "-" "+" _dup(`=34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "-"
		qui count
		local N = r(N)
		
		forvalues i = 1(1)`N' {
			//Weight
			if "`wt'" =="" {
				local ww1 = _WT1[`i']
				local ww0 = _WT0[`i']
			}
			
			//Smooth estimates
			if "`smooth'" !="" {
				local mes1 "`modeles'1[`i']*(10^`power')"
				local mlci1 "`modellci'1[`i']*(10^`power')"
				local muci1 "`modeluci'1[`i']*(10^`power')"
				
				local mes0 "`modeles'0[`i']*(10^`power')"
				local mlci0 "`modellci'0[`i']*(10^`power')"
				local muci0 "`modeluci'0[`i']*(10^`power')"
			}
		
			//Group labels
			if ((`use'[`i']== -2)){ 
				di _col(2) as txt `label'0[`i'] _col(`nlen') "|  " _col(`=`nlen' + 45') "|  "
			}
			//Studies -- se
			if ((`use'[`i'] ==1)) { 
				di _col(2) as txt `label'1[`i'] _col(`nlen') "|  "  ///
				as res  %`=`est_i_len''.`=`dp''f  `es'1[`i']*(10^`power') "`open'" `aes' `mes1' "`close'"  /// 
				_skip(8) %`=`est_i_len''.`=`dp''f `lci'1[`i']*(10^`power') "`open'" `aes' `mlci1'  `aesclose' "`close'" ///
				_skip(`=5 - 2*`togglesmooth'') %`=`est_i_len''.`=`dp''f `uci'1[`i']*(10^`power') "`open'" `aes' `muci1'  `aesclose' "`close'"   _skip(`=3 + 9*`togglewt' + (`est_i_len' - 7)*`togglesmooth'') as res %5.`=`dp''f `ww1' _cont
			}
			//studies - sp
			if (`use'[`i'] ==1 )   { 
				di as txt _skip(`=4*(1-`togglesmooth')*(1-`togglewt')') "|  "  ///
				as res  %`=`est_i_len''.`=`dp''f  `es'0[`i']*(10^`power') "`open'" `aes' `mes0' "`close'"  /// 
				_col(`=`nlen' + 50 + (`est_i_len'*3 + 6)*`togglesmooth' + 9*`togglewt'') %`=`est_i_len''.`=`dp''f `lci'0[`i']*(10^`power') "`open'" `aes' `mlci0'  `aesclose' "`close'" ///
				_skip(`=5 - 2*`togglesmooth'') %`=`est_i_len''.`=`dp''f `uci'0[`i']*(10^`power') "`open'" `aes' `muci0'  `aesclose' "`close'" _skip(`=3 + 9*`togglewt'+ (`est_i_len' - 7)*`togglesmooth'') as res %5.`=`dp''f `ww0' 
			}
			//Summaries			
			if ( (`use'[`i']== 3) | (`use'[`i']== 2)){
				if ((`use'[`i']== 2) ) {
					di _col(2) as txt _col(`nlen') "|  " _col(`=`nlen' + 1 + 34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "|  "
				}
				
				di  _dup(`=`nlen'-1') "-" "+" _dup(`=34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "-" "+" _dup(`=34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "-"
				di as txt _col(`nlen') "|  " _col(`=`nlen' + 1 + 34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "|  "
				
				if ("`cveffect'" != "se" & (`use'[`i']== 2)) | (`use'[`i']== 3) {
					local cont "_cont"
				}
				//se
				if ("`cveffect'" != "sp" & (`use'[`i']== 2)) | (`use'[`i']== 3) {
					di _col(2) as txt `label'0[`i'] _col(`nlen') "|  "  ///
					_skip(`=(`est_i_len' + 2)*`togglesmooth'') as res  %`=`est_i_len''.`=`dp''f  `es'1[`i']*(10^`power') /// 
					_col(`=`nlen' + 15 + (`est_i_len'*2 + 5)*`togglesmooth'') %`=`est_i_len''.`=`dp''f `lci'1[`i']*(10^`power') ///
					_skip(`=5 +  (`est_i_len' + 1 )*`togglesmooth'') %`=`est_i_len''.`=`dp''f `uci'1[`i']*(10^`power') _skip(`=3 + 9*`togglewt' + (`est_i_len' -7 )*`togglesmooth'') as res %5.`=`dp''f `ww1' `cont'
				}
				//sp
				if ("`cveffect'" != "se" & (`use'[`i']== 2)) | (`use'[`i']== 3) {
					di as txt _col(`=`nlen' + 1 + 34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "|  " ///
					_skip(`=(`est_i_len' + 2)*`togglesmooth'') as res  %`=`est_i_len''.`=`dp''f  `es'0[`i']*(10^`power') /// 
					_col(`=`nlen' + 50 + (`est_i_len'*4 + 8)*`togglesmooth' + 9*`togglewt'') %`=`est_i_len''.`=`dp''f `lci'0[`i']*(10^`power') ///
					_skip(`=5 +  (`est_i_len' + 1 )*`togglesmooth'') %`=`est_i_len''.`=`dp''f `uci'0[`i']*(10^`power') _skip(`=3 + 7*`togglewt'+ `est_i_len'*`togglesmooth'') as res %5.`=`dp''f `ww0'
				}
			}
			//Blanks
			if (`use'[`i'] == 0 ){
				di  _dup(`=`nlen'-1') "-" "+" _dup(`=34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "-" "+" _dup(`=34 + (`est_i_len'*2 + 4)*`togglesmooth' + 9*`togglewt'') "-"
				*di as txt _col(`nlen') "|  "_col(`=`nlen' + 1 + 34 + (`est_i_len'*2 + 4)*`togglesmooth' + 5*`togglewt'') "|  "
			}
		}
	restore
end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: PRINTMAT +++++++++++++++++++++++++
							Print the outplot matrix beautifully
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop printmat
program define printmat
	version 13.1
	syntax, type(string) [cveffect(string) matrixout(name) matrixoutse(name) ///
			matrixoutsp(name) sumstat(string) dp(integer 2) p(integer 0) power(integer 0) model(string) ///
			matched cbnetwork abnetwork general comparative continuous  nsims(string) level(string)]
					
		local rownamesmaxlen = 10		
		if ("`type'" != "mc") {
			local nrows = rowsof(`matrixout')
			local ncols = colsof(`matrixout')
			local rnames : rownames `matrixout'
			
			forvalues r = 1(1)`nrows' {
				local rname : word `r' of `rnames'
				local nlen : strlen local rname
				local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
			}
		}
		if ("`type'" == "gof") {
			local rspec "--`="&"*`=`nrows'-1''-"
			di as text _n "Goodness of Fit Criterion"
			
			#delimit ;
			noi matlist `matrixout',  rowtitle(Model) 
						cspec(& %7s |   %8.`=`dp''f &  %8.`=`dp''f o2&) 
						rspec(`rspec') underscore  nodotz
			;
			#delimit cr 
		} 
			
		if (strpos("`type'", "abs") != 0){	
			if `nrows' > 3 {
				if "`cveffect'" == "sesp" {
					local rspec "--`="&"*`=`nrows'/2 - 1''-`="&"*`=`nrows'/2 - 1''-"
				}
				else if "`cveffect'" == "se" {
					local rspec "---`="&"*`=`nrows'-4''---"
				}
				else if "`cveffect'" == "sp" {
					local rspec "-----`="&"*`=`nrows'-4''-"
				}
			}
			else {
				local rspec "----"
			}
		}
		if strpos("`type'", "r") != 0 {
			local rownamesmaxlen = 20
			if "`cveffect'" == "sesp" {
				if `nrows' > 3 {
					local rspec "---`="&"*`=`nrows'/2 - 2''--`="&"*`=`nrows'/2 - 2''-"
				}
				else {
					local rspec "--&-&-"
				}
			}
			else {
				local rspec "--`="&"*`=`nrows'-1''-"
			}			
		}
				
		local nlensstat : strlen local sumstat
		local nlensstat = max(10, `nlensstat')
		
		if "`type'" == "rre" {
			local rownamesmaxlen = max(`rownamesmaxlen', 20) //Check if there is a longer name
			local rspec "--`="&"*`=`nrows'-1''-"
			di as res _n "****************************************************************************************"
			di as txt _n "Wald-type test for nonlinear hypothesis"
			di as txt _n "{phang}H0: All (log)RR equal vs. H1: Some (log)RR different {p_end}"

			#delimit ;
			noi matlist `matrixout', rowtitle(Effect) 
						cspec(& %`rownamesmaxlen's |  %8.`=`dp''f &  %8.0f &  %8.4f o2&) 
						rspec(`rspec') underscore nodotz
			;
			#delimit cr			
		}
		if ("`type'" == "logit") | ("`type'" == "abs") | ("`type'" == "rr") | ("`type'" == "or")  {
			if ("`model'" == "random") {
				local typeinf "Conditional"
			}
			else {
				local typeinf "Marginal"
			}
			
			if "`type'" == "abs" {
				local parm "Absolute proportion"
			}
			else if "`type'" == "logit" {
				local parm "Log odds"
			}
			else if "`type'" == "or" {
				local parm "Odds Ratio"
			}
			else if "`type'" == "rr" {
				local parm "Proportion Ratio"
			}
		
			di as res _n "****************************************************************************************"
			di as res "{pmore2} `typeinf' summary measures of test accuracy: `parm' {p_end}"
			di as res    "****************************************************************************************" 
			*tempname mat2print
			*mat `mat2print' = `matrixout'
			local nrows = rowsof(`matrixout')
			forvalues r = 1(1)`nrows' {
				mat `matrixout'[`r', 1] = `matrixout'[`r', 1]*10^`power'
				mat `matrixout'[`r', 5] = `matrixout'[`r', 5]*10^`power'
				mat `matrixout'[`r', 6] = `matrixout'[`r', 6]*10^`power'
						
				forvalues c = 1(1)6 {
					local cell = `matrixout'[`r', `c'] 
					if "`cell'" == "." {
						mat `matrixout'[`r', `c'] == .z
					}
				}
			}
			
			#delimit ;
			noi matlist `matrixout', rowtitle(Effect) 
						cspec(& %`rownamesmaxlen's |  %`nlensstat'.`=`dp''f &  %9.`=`dp''f &  %8.`=`dp''f &  %15.`=`dp''f &  %8.`=`dp''f &  %8.`=`dp''f o2&) 
						rspec(`rspec') underscore  nodotz
			;
			#delimit cr
		}
		if ("`type'" == "popabs") | ("`type'" == "poprr") | ("`type'" == "popor")  {
			if "`type'" == "popabs" {
				local statistic "Proportion"
			}
			else if "`type'" == "poprr" {
				local statistic "Proportion Ratio"
			}
			else if "`type'" == "popor" {
				local statistic "Odds Ratio"
			}
			di as res _n "****************************************************************************************"		
			di as res _n "Population-averaged estimates: `statistic' "
			di as res    "****************************************************************************************" 
			local nrows = rowsof(`matrixout')
			forvalues r = 1(1)`nrows' {
				mat `matrixout'[`r', 1] = `matrixout'[`r', 1]*10^`power'
				mat `matrixout'[`r', 3] = `matrixout'[`r', 3]*10^`power'
				mat `matrixout'[`r', 4] = `matrixout'[`r', 4]*10^`power'
				mat `matrixout'[`r', 5] = `matrixout'[`r', 5]*10^`power'
						
				forvalues c = 1(1)5 {
					local cell = `matrixout'[`r', `c'] 
					if "`cell'" == "." {
						mat `matrixout'[`r', `c'] == .z
					}
				}
			}
			
			#delimit ;
			noi matlist `matrixout', rowtitle(Parameter) 
						cspec(& %`rownamesmaxlen's |  %`nlensstat'.`=`dp''f &  %8.`=`dp''f &  %8.`=`dp''f &  %8.`=`dp''f &  %8.`=`dp''f o2&) 
						rspec(`rspec') underscore  nodotz
			;
			#delimit cr
		}
		if ("`type'" == "bhet") {
				di as res _n "****************************************************************************************"
				di as txt _n "Between-study heterogeneity statistics"
				
			if `nrows' > 1 {
				local rspec "-`="&"*`nrows''-"
		
				*tempname mat2print
				*mat `mat2print' = `matrixout'
				forvalues r = 1(1)`nrows' {
					forvalues c = 1(1)`ncols' {
						local cell = `matrixout'[`r', `c'] 
						if "`cell'" == "." {
							mat `matrixout'[`r', `c'] == .z
						}
					}
				}
					
				#delimit ;
				noi matlist `matrixout', 
							cspec(& %`rownamesmaxlen's |  %13.`=`dp''f `="&  %13.`=`dp''f "*`=`ncols'-1'' o2&) 
							rspec(`rspec') underscore nodotz
				;
				#delimit cr	
			}
			else {			
				local rho 		= `matrixout'[1, 8]
				local covar 	= `matrixout'[1, 7]
				local tau2se 	= `matrixout'[1, 1]
				local tau2sp	= `matrixout'[1, 3]
				local tau2g		= `matrixout'[1, 5]
				
				di as txt _col(28) "covar"  _col(45) "rho"_cont
				di as res _n _col(28) %5.`=`dp''f `covar' _col(45) %5.`=`dp''f `rho' 
				
				di as txt  _col(28) "Tau.sq" _cont
				di as txt _col(45) "I^2(%)" _cont
				if `p' == 0  {					
					local isq2g  = `matrixout'[1, 6]*100
					local isq2se = `matrixout'[1, 2]*100
					local isq2sp = `matrixout'[1, 4]*100
				}
				else {
					local isq2g  = `matrixout'[1, 6]
					local isq2se = `matrixout'[1, 2]
					local isq2sp = `matrixout'[1, 4]
				}
				di as txt _n  "Generalized" _cont	
				di as res   _col(28) %5.`=`dp''f `tau2g' _col(45) %5.`=`dp''f `isq2g'  
				di as txt  "Sensitivity" _cont	
				di as res    _col(28) %5.`=`dp''f `tau2se' _col(45) %5.`=`dp''f `isq2se'  
				di as txt  "Specificity" _cont
				di as res    _col(28) %5.`=`dp''f `tau2sp' _col(45) %5.`=`dp''f `isq2sp'
			}
		}
		if ("`type'" == "mc") {
			local semat = 1
			local spmat = 1
			if "`cveffect'" == "se" {
				*cap confirm matrix `matrixoutse'
				*local semat = _rc
				local spmat = 0
			}
			if "`cveffect'" == "sp" { 
				*cap confirm matrix `matrixoutsp'
				*local spmat = _rc
				local semat = 0
			}
	
			tempname matrixout
			if (`=`semat' + `spmat'') == 2 {
				mat `matrixout' = `matrixoutse' \ `matrixoutsp' 
				local nrowse = rowsof(`matrixoutse')
				local nrowsp = rowsof(`matrixoutsp')
				local rspec "--`="&"*`=`nrowse'-1''-`="&"*`=`nrowsp'-1''-"
			}
			else if `spmat' == 0 {
				mat `matrixout' = `matrixoutse' 
				local nrowse = rowsof(`matrixoutse')
				local rspec "--`="&"*`=`nrowse'-1''-"
			}
			else if `semat' == 0 {
				mat `matrixout' = `matrixoutsp' 
				local nrowsp = rowsof(`matrixoutsp')
				local rspec "--`="&"*`=`nrowsp'-1''-"
			}
			
			local ncols = colsof(`matrixout')
			local nrows = rowsof(`matrixout')
			local rnames : rownames `matrixout'
			mat colnames `matrixout' = chi2 df pval
			
			local rownamesmaxlen = 15
			forvalues r = 1(1)`nrows' {
				local rname : word `r' of `rnames'
				local nlen : strlen local rname
				local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
				
				forvalues c = 1(1)`ncols' {
					local cell = `matrixout'[`r', `c'] 
					if "`cell'" == "." {
						mat `matrixout'[`r', `c'] == .z
					}
				}
			}
			
			di as res _n "****************************************************************************************"
			di as txt _n "Model comparison(s): Leave-one/all-out LR Test(s)"
			#delimit ;
			noi matlist `matrixout', rowtitle(Excluded Effect) 
				cspec(& %`=`rownamesmaxlen' + 2's |  %8.`=`dp''f  `="&  %8.`=`dp''f "*`=`ncols'-1'' o2&) 
				rspec(`rspec') underscore nodotz
			;
		
			#delimit cr
			if "`interaction'" !="" {
				di as txt "*NOTE: Model with and without interaction effect(s)"
			}
			else {
				di as txt "*NOTE: Model with and without main effect(s)"
			}
		}
		if ("`type'" == "refe") {
			di as txt _n "LR Test: RE vs FE model" 
			if `nrows' > 1 {
				local rspec "--`="&"*`=`nrows'-1''-"
				*tempname mat2print
				*mat `mat2print' = `matrixout'
				forvalues r = 1(1)`nrows' {
					forvalues c = 1(1)`ncols' {
						local cell = `matrixout'[`r', `c'] 
						if "`cell'" == "." {
							mat `matrixout'[`r', `c'] == .z
						}
					}
				}
					
				#delimit ;
				noi matlist `matrixout', 
							cspec(& %`rownamesmaxlen's |  %8.0f `="&  %10.`=`dp''f "*`=`ncols'-1'' o2&) 
							rspec(`rspec') underscore nodotz
				;
				#delimit cr	
			}
			else {
				local chisq = `matrixout'[1, 1]
				local df 	= `matrixout'[1, 2]
				local pv 	= `matrixout'[1, 3]
				
				di as txt  _col(10) "Chi2"  _skip(8) "degrees of" _cont
				di as txt _n  _col(8) "statistic" 	_skip(6) "freedom"      _skip(8)"p-val"   
				di as res _col(5) %10.`=`dp''f `chisq' _col(25) `df' _col(34) %10.4f `pv' 
			}
		}
		if ("`type'" == "popabs") | ("`type'" == "poprr") | ("`type'" == "popor")  {
			di as txt "NOTE: `level'% centiles obtained from `nsims' simulations of the posterior distribution"
		}
		if ("`continuous'" != "") {
			di as txt "NOTE: For continuous variable margins are computed at their respective mean"
		} 
		if ("`type'" == "abs") {
			di as txt "NOTE: H0: P = 0.5 vs. H1: P != 0.5"
		}

		
end	

/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: LONGSETUP +++++++++++++++++++++++++
							Transform data to long format
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop longsetup
program define longsetup
version 14.0

syntax varlist, rid(name) event(name) total(name) se(name) [first(name) rowid(name) assignment(name) idpair(name) cbnetwork abnetwork general comparative ]

	qui{
		tempvar tp tn fp fn 
		tokenize `varlist'
		if "`cbnetwork'" == "" {
			local nvar = 4
		}
		else {
			local nvar = 8
		}		
		/*The four variables should contain numbers*/
		forvalue i=1(1)`nvar' {
			capture confirm numeric var ``i''
				if _rc != 0 {
					di as error "The variable ``i'' must be numeric"
					exit
				}	
		}
		if "`cbnetwork'" != "" {
			gen `tp'1 = `1'
			gen `fp'1 = `2'
			gen `fn'1 = `3'
			gen `tn'1 = `4'
			gen `tp'0 = `5'
			gen `fp'0 = `6'
			gen `fn'0 = `7'
			gen `tn'0 = `8'
			gen `rowid' = _n
			gen `assignment'1 = `9'
			gen `assignment'0 = `10'
			
			reshape long `tp' `fp' `fn'  `tn' `assignment', i(`rowid') j(`idpair')
		}
		else {
			gen `tp' = `1'
			gen `fp' = `2'
			gen `fn' = `3'
			gen `tn' = `4'
		}
		
		/*4 variables per study : TP TN FP FN*/
		gen `event'1 = `tp'  /*TP*/
		gen `event'0 = `tn'  /*TN*/
		gen `total'1 = `tp' + `fn'  /*DIS = TP + FN*/
		gen `total'0 = `tn' + `fp' /*NDIS = TN + FP*/
		
		gen `rid' = _n	
		if "`abnetwork'" != "" {
			reshape long `event' `total', i(`rid' `first') j(`se')
		}
		else {
			reshape long `event' `total', i(`rid') j(`se')
		}
	}
end

	/*++++++++++++++++	SUPPORTING FUNCTIONS: BUILDEXPRESSIONS +++++++++++++++++++++
				buildexpressions the regression and estimation expressions
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop buildregexpr
	program define buildregexpr, rclass
	version 13.1
		
		syntax varlist, [cveffect(string) interaction(string) se(name) sp(name) alphasort cbnetwork abnetwork general comparative ipair(varname) baselevel(string)]
		
		tempvar holder
		tokenize `varlist'
		if "`cbnetwork'" == "" {
			macro shift 4
			local regressors "`*'"
		}
		else {
			local index = "`9'"
			local comparator = "`10'"
			macro shift 10
			local regressors "`*'"
			
			my_ncod `holder', oldvar(`index')
			drop `index'
			rename `holder' `index'
			
			my_ncod `holder', oldvar(`comparator')
			drop `comparator'
			rename `holder' `comparator'

			my_ncod `holder', oldvar(`ipair')
			drop `ipair'
			rename `holder' `ipair'
			
		}
		local p: word count `regressors'
		
		if "`general'`comparative'" != "" {
			local seregexpression = `"`se'"'
			local spregexpression = `"`sp'"'
		}
		else if "`cbnetwork'" != "" {
			if "`interaction'" == "" {	
				local seregexpression = `"`se' i.`ipair'#c.`se' i.`index'#c.`se' "'
				local spregexpression = `"`sp'  i.`ipair'#c.`sp' i.`index'#c.`sp'"'
			}
			else {
				local seregexpression = `"`se' i.`ipair'#i.`comparator'#c.`se' i.`index'#c.`se' "'
				local spregexpression = `"`sp'  i.`ipair'#i.`comparator'#c.`sp' i.`index'#c.`sp'"'
				//nulllify
				local interaction
			}
		}
		else {
			*abnetwork
			local seregexpression
			local spregexpression
		}
	
		local catreg " "
		local contreg " "
		
		local basecode 1
		tokenize `regressors'
		forvalues i = 1(1)`p' {			
			capture confirm numeric var ``i''
			if _rc != 0 {
				if "`alphasort'" != "" {
					sort ``i''
				}
				my_ncod `holder', oldvar(``i'')
				drop ``i''
				rename `holder' ``i''
				local prefix_`i' "i"
			}
			else {
				local prefix_`i' "c"
			}
			if "`abnetwork'`comparative'" != "" & `i'==1 {
				if "`abnetwork'" != "" {
					local prefix_`i' "ibn"
				}
				if "`baselevel'" != "" {
					//Find the base level
					qui label list ``i''
					local nlevels = r(max)
					local found = 0
					local level 1
					while !`found' & `level' < `=`nlevels'+1' {
						local lab:label ``i'' `level'
						if "`lab'" == "`baselevel'" {
							local found = 1
							local basecode `level'
						}
						local ++level
					}
				}
			}
			/*Add the proper expression for regression*/
			local seregexpression = "`seregexpression' `prefix_`i''.``i''#c.`se'"
			local spregexpression = "`spregexpression' `prefix_`i''.``i''#c.`sp'"
			
			if `i' > 1 & "`interaction'" != "" {
				if "`interaction'" == "se" {
					local seregexpression = "`seregexpression' `prefix_`i''.``i''#`prefix_1'.`1'#c.`se'"
				}
				else if "`interaction'" == "sp" {
					local spregexpression = "`spregexpression' `prefix_`i''.``i''#`prefix_1'.`1'#c.`sp'"
				}
				else {
					local seregexpression = "`seregexpression' `prefix_`i''.``i''#`prefix_1'.`1'#c.`se'"
					local spregexpression = "`spregexpression' `prefix_`i''.``i''#`prefix_1'.`1'#c.`sp'"
				}
			}
			//Pick out the interactor variable
			if `i' == 1 /*& "`interaction'" != "" */{
				local varx = "``i''"
				if 	"`abnetwork'" != "" {
					local prefix_`i' = "i"
				}
				local typevarx = "`prefix_`i''"
			}
			* (`i' > 1 & "`interaction'" != "" ) |  "`interaction'" == ""  { //store the rest of  variables
			if "`prefix_`i''" == "i" {
				local catreg "`catreg' ``i''"
			}
			else {
				local contreg "`contreg' ``i''"
			}
			*}/
		}
		if "`cveffect'" == "sp" {
			local seregexpression "`se'"
		}
		else if "`cveffect'" == "se" {
			local spregexpression "`sp'"
		}
		
		return local varx = "`varx'"
		return local typevarx  = "`typevarx'"
		return local  regexpression = "`seregexpression' `spregexpression'"
		return local seregexpression =  "`seregexpression'"
		return local spregexpression  = "`spregexpression'"
		return local  catreg = "`catreg'"
		return local  contreg = "`contreg'"
		return local basecode = "`basecode'"
	end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS:  ESTP +++++++++++++++++++++++++
							estimate log odds or proportions after modelling
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/	
	cap program drop estp
	program define estp, rclass
	version 14.1
		syntax [if], estimates(string) [sumstat(string) depname(string) expit se(varname) sp(varname) DP(integer 2)) ///
			cveffect(string) interaction(string) catreg(varlist) contreg(varlist) power(integer 0) ///
			level(integer 95) tradeoff by(varname) varx(varname) typevarx(string) regexpression(string) abnetwork cbnetwork stratify general comparative ]
		
			tempname outmatrix contregmatrixout catregmatrixout bycatregmatrixout secontregmatrixout spcontregmatrixout outmatrixse ///
				outmatrixsp cutmat serow sprow outmatrixse outmatrixsp outmatrixr overallse overallsp Vmatrix byVmatrix
			
			*Nullify by			
			*if "`stratify'" != "" {
			*	local by
			*}
			marksample touse, strok
			
			tokenize `regexpression'
			
			if "`tradeoff'" != "" {
				tokenize `2', parse("#")
				tokenize `1', parse(".")
				local cutoff "`3'"
			}  
							
			if "`cbnetwork'" != "" {
				if "`interaction'" != "" {
					tokenize `2', parse("#")
					tokenize `1', parse(".")
				 }
				 else {
					tokenize `3', parse("#")
					tokenize `1', parse(".")
				 }
	
				local index "`3'"
				local catreg = "`3' `catreg'"
				local varx //nullify
			}
			
			if "`abnetwork'" != "" {
				tokenize `1', parse("#")
				tokenize `1', parse(".")
				
				local assignment "`3'"
				local catreg = "`3' `catreg'"
			}
			
			if "`interaction'" != "" & "`typevarx'" == "i" {
				local idpairconcat "#`varx'"
			}
			
			if "`typevarx'" == "i"  {
				if "`catreg'" == "" {
					local catreg = "`varx'"
				}
			}
			if "`typevarx'" == "c"  {
				if "`contreg'" == "" {
					local contreg = "`varx'"
				}
			}
			
			local marginlist
			while "`catreg'" != "" {
				tokenize `catreg'
				if ("`1'" != "`by'" & "`by'" != "") | ("`by'" =="") {
					local marginlist = `"`marginlist' `1'`idpairconcat'"'
				}
				macro shift 
				local catreg `*'
			}
			
			if "`tradeoff'" != "" {
				qui estimates restore `estimates'
				qui nlcom cutoff: (_b[`se'] - _b[`sp'])/(_b[`sp'#`cutoff'] - _b[`se'#`cutoff'] ), post level(`level') iterate(200)
				mat `cutmat' = r(b)
			}
			
			qui estimates restore `estimates'

			local byncatreg 0
			if "`by'" != "" & "`stratify'"  == "" {
				qui margin if `touse', predict(xb) over(`se' `by') level(`level')
				
				mat `bycatregmatrixout' = r(table)'
				mat `byVmatrix' = r(V)
				mat `bycatregmatrixout' = `bycatregmatrixout'[1..., 1..6]
				
				local byrnames :rownames `bycatregmatrixout'
				local byncatreg = rowsof(`bycatregmatrixout')
			}
			
			/*if "`abnetwork'`cbnetwork'" == ""  {
				local grand "grand"
				local Overall "Overall"
			}*/
			/*if "`comparative'" != "" & "`stratify'" != "" {
				local grand
				local Overall
			}
			else {*/
				local grand "grand"
				local Overall "Overall"
			*}
			
			local ncatreg 0
			qui margin `marginlist' if `touse', over(`se') predict(xb) `grand' level(`level')
						
			mat `catregmatrixout' = r(table)'
			mat `Vmatrix' = r(V)
			mat `catregmatrixout' = `catregmatrixout'[1..., 1..6]
			
			local rnames :rownames `catregmatrixout'
			local ncatreg = rowsof(`catregmatrixout')
			
			local init 1
			local ncontreg 0
			local contserownames = ""
			local contsprownames = ""
			if "`contreg'" != "" {
				foreach v of local contreg {
					summ `v' if `touse', meanonly
					local vmean = r(mean)
					qui margin if `touse', over(`se') predict(xb) at(`v'=`vmean') level(`level')
					mat `contregmatrixout' = r(table)'
					mat `contregmatrixout' = `contregmatrixout'[1..., 1..6]
					if `init' {
						local init 0
						mat `secontregmatrixout' = `contregmatrixout'[2, 1...] 
						mat `spcontregmatrixout' = `contregmatrixout'[1, 1...] 
					}
					else {
						mat `secontregmatrixout' =  `secontregmatrixout' \ `contregmatrixout'[2, 1...]
						mat `spcontregmatrixout' =  `spcontregmatrixout' \ `contregmatrixout'[1, 1...]
					}
					local contserownames = "`contserownames' `v'"
					local contsprownames = "`contsprownames' `v'"
					local ++ncontreg
				}
				mat rownames `secontregmatrixout' = `contserownames'
				mat rownames `spcontregmatrixout' = `contsprownames'
			}
			
			if "`expit'" != "" {
				forvalues r = 1(1)`byncatreg' {
					mat `bycatregmatrixout'[`r', 1] = invlogit(`bycatregmatrixout'[`r', 1])
					mat `bycatregmatrixout'[`r', 5] = invlogit(`bycatregmatrixout'[`r', 5])
					mat `bycatregmatrixout'[`r', 6] = invlogit(`bycatregmatrixout'[`r', 6])
				}
				forvalues r = 1(1)`ncatreg' {
					mat `catregmatrixout'[`r', 1] = invlogit(`catregmatrixout'[`r', 1])
					mat `catregmatrixout'[`r', 5] = invlogit(`catregmatrixout'[`r', 5])
					mat `catregmatrixout'[`r', 6] = invlogit(`catregmatrixout'[`r', 6])
				}
				forvalues r = 1(1)`ncontreg' {
					mat `secontregmatrixout'[`r', 1] = invlogit(`secontregmatrixout'[`r', 1])
					mat `secontregmatrixout'[`r', 5] = invlogit(`secontregmatrixout'[`r', 5])
					mat `secontregmatrixout'[`r', 6] = invlogit(`secontregmatrixout'[`r', 6])
					
					mat `spcontregmatrixout'[`r', 1] = invlogit(`spcontregmatrixout'[`r', 1])
					mat `spcontregmatrixout'[`r', 5] = invlogit(`spcontregmatrixout'[`r', 5])
					mat `spcontregmatrixout'[`r', 6] = invlogit(`spcontregmatrixout'[`r', 6])
				}
			}
			
			local serownames = ""
			local sprownames = ""
			
			local rownamesmaxlen = 10 /*Default*/
			
			*if "`grand'" != "" {
				local nrowss = `ncatreg' + `byncatreg' - 2 //Except the grand rows
			*}
			*else {
			local nrowscat = `ncatreg' + `byncatreg' 
			*}
			
			
			//# equations
			if "`cveffect'" == "sesp" {
				local keq 2
			}
			else {
				local keq 1
			} 
			mat `serow' = J(1, 6, .)
			mat `sprow' = J(1, 6, .)

			
			local initse 0
			local initsp 0	
			local rnames = "`byrnames' `rnames'" //attach the bynames	
			
			//Except the grand rows	
			forvalues r = 1(1)`=`byncatreg' + `ncatreg'  - `=2*`="`grand'"!=""''' {
				//Labels
				local rname`r':word `r' of `rnames'
				tokenize `rname`r'', parse("#")					
				local parm = "`1'"
				local left = "`3'"
				local right = "`5'"
				
				tokenize `left', parse(.)
				local leftv = "`3'"
				local leftlabel = "`1'"
				
				if "`right'" == "" {
					if "`leftv'" != "" {
						if strpos("`rname`r''", "1b") == 0 {
							local lab:label `leftv' `leftlabel'
						}
						else {
							local lab:label `leftv' 1
						}
						local eqlab "`leftv'"
					}
					else {
						local lab "`leftlabel'"
						local eqlab ""
					}
					local nlencovl : strlen local llab
					local nlencov = `nlencovl' + 1					
				}
				else {								
					tokenize `right', parse(.)
					local rightv = "`3'"
					local rightlabel = "`1'"
					
					if strpos("`leftlabel'", "c") == 0 {
						if strpos("`leftlabel'", "o") != 0 {
							local indexo = strlen("`leftlabel'") - 1
							local leftlabel = substr("`leftlabel'", 1, `indexo')
						}
						if strpos("`leftlabel'", "1b") == 0 {
							local llab:label `leftv' `leftlabel'
						}
						else {
							local llab:label `leftv' 1
						}
					} 
					else {
						local llab
					}
					
					if strpos("`rightlabel'", "c") == 0 {
						if strpos("`rightlabel'", "o") != 0 {
							local indexo = strlen("`rightlabel'") - 1
							local rightlabel = substr("`rightlabel'", 1, `indexo')
						}
						if strpos("`rightlabel'", "1b") == 0 {
							local rlab:label `rightv' `rightlabel'
						}
						else {
							local rlab:label `rightv' 1
						}
					} 
					else {
						local rlab
					}
					
					if (("`rlab'" != "") + ("`llab'" != "")) ==  0 {
						local lab = "`leftv'#`rightv'"
						local eqlab = ""
					}
					if (("`rlab'" != "") + ("`llab'" != "")) ==  1 {
						local lab = "`llab'`rlab'" 
						local eqlab = "`leftv'*`rightv'"
					}
					if (("`rlab'" != "") + ("`llab'" != "")) ==  2 {
						local lab = "`llab'|`rlab'" 
						local eqlab = "`leftv'*`rightv'"
					}
					local nlencovl : strlen local leftv
					local nlencovr : strlen local rightv
					local nlencov = `nlencovl' + `nlencovr' + 1
				}
				
				local lab = ustrregexra("`lab'", " ", "_")
				
				local nlenlab : strlen local lab
				if "`eqlab'" != "" {
					local nlencov = `nlencov'
				}
				else {
					local nlencov = 0
				}
				local rownamesmaxlen = max(`rownamesmaxlen', min(`=`nlenlab' + `nlencov' + 1', 32)) /*Check if there is a longer name*/
				
				if "`cbnetwork'" != "" & "`eqlab'"=="`index'" {
					local eqlab "`index'"
				}
				
				//se or sp
				local parm = substr("`parm'", 1, 1)
				if `r' > `=`byncatreg'' {
					mat `outmatrixr' = `catregmatrixout'[`=`r' - `byncatreg'', 1...] //select the r'th row
				}
				else{
					mat `outmatrixr' = `bycatregmatrixout'[`r', 1...] //select the r'th row
				}

				if `parm' == 0 {
					if `initsp' == 0 {
						mat `outmatrixsp' = `outmatrixr'
					}
					else {
						mat `outmatrixsp' = `outmatrixsp' \ `outmatrixr'
					}
					local initsp 1
					local sprownames = "`sprownames' `eqlab':`lab'"
				}
				else {
					if `initse' == 0 {
						mat `outmatrixse' = `outmatrixr'
					}
					else {
						mat `outmatrixse' = `outmatrixse' \ `outmatrixr'
					}
					local initse 1
					local serownames = "`serownames' `eqlab':`lab'"
				}
			}
			
			if `nrowscat'  > 2 {
				mat rownames `outmatrixse' = `serownames'
				mat rownames `outmatrixsp' = `sprownames'
			}			
			/*if "`interaction'" != "" {
				mat rownames `serow' = "Sensitivity--**"
				mat rownames `sprow' = "**--Specificity--**" //19 characters
			}
			else {*/
			mat rownames `serow' = "Sensitivity"
			mat rownames `sprow' = "Specificity" //19 characters
			
			local rownamesmaxlen = max(`rownamesmaxlen', 19) /*Check if there is a longer name*/
			
			*if "`grand'" != "" {
				mat `overallsp' = `catregmatrixout'[`=`ncatreg'-1', 1...]
				mat `overallse' = `catregmatrixout'[`ncatreg', 1...]
			/*}
			else {
				mat `overallsp' =
				mat `overallse' =				
			}*/
			mat rownames `overallse' = "Overall"
			mat rownames `overallsp' = "Overall"
						
			
			if `=`ncatreg' + `byncatreg' - 2' > 0 | `ncontreg' > 0 {
				if "`cveffect'" == "sesp" {
					if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' > 0) {
						mat `outmatrixse' = `serow' \ `outmatrixse' \ `secontregmatrixout' \ `overallse'
						mat `outmatrixsp' = `sprow' \ `outmatrixsp' \ `spcontregmatrixout' \ `overallsp' 
					}
					else if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' == 0) {
						mat `outmatrixse' = `serow' \ `outmatrixse' \ `overallse'
						mat `outmatrixsp' = `sprow' \ `outmatrixsp' \ `overallsp'
					}
					else if (`=`ncatreg' + `byncatreg' - 2' == 0) & (`ncontreg' > 0) {
						mat `outmatrixse' = `serow' \ `secontregmatrixout' \ `overallse'
						mat `outmatrixsp' = `sprow' \ `spcontregmatrixout' \ `overallsp'
					}
				}
				else {
					if "`cveffect'" == "se" {
						if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' > 0) {
							mat `outmatrixse' = `serow' \ `outmatrixse' \ `secontregmatrixout' \ `overallse' 
							mat `outmatrixsp' = `sprow' \ `overallsp'
						}
						else if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' == 0) {
							mat `outmatrixse' = `serow' \ `outmatrixse' \ `overallse'
							mat `outmatrixsp' = `sprow' \ `overallsp'
						}
						else if (`=`ncatreg' + `byncatreg' - 2' == 0) & (`ncontreg' > 0) {
							mat `outmatrixse' = `serow' \ `secontregmaritxout' \ `overallse' 
							mat `outmatrixsp' = `sprow' \ `overallsp'
						}
					}
					else {
						if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' > 0) { 
							mat `outmatrixse' = `serow' \ `overallse'
							mat `outmatrixsp' = `sprow' \ `outmatrixsp' \ `spcontregmatrixout' \ `overallsp'
						}
						else if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' == 0) {
							mat `outmatrixse' = `serow' \ `overallse'
							mat `outmatrixsp' =	`sprow' \ `outmatrixsp' \ `overallsp'
						}
						else if (`=`ncatreg' + `byncatreg' - 2' == 0) & (`ncontreg' > 0) {
							mat `outmatrixse' = `serow' \ `overallse'
							mat `outmatrixsp' =	`sprow' \ `spcontregmatrixout' \ `overallsp'					
						}
					}
				}
				
				
				/*if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' > 0) {
					mat `outmatrixse' = `outmatrixse' \ `secontregmatrixout' \ `overallse' 
					mat `outmatrixsp' = `outmatrixsp' \ `spcontregmatrixout' \ `overallsp' 				
				}
				else if (`=`ncatreg' + `byncatreg' - 2' > 0) & (`ncontreg' == 0) {
					mat `outmatrixse' = `outmatrixse' \ `overallse' 
					mat `outmatrixsp' = `outmatrixsp' \ `overallsp' 
				}
				else if (`=`ncatreg' + `byncatreg' - 2' == 0) & (`ncontreg' > 0) {
					mat `outmatrixse' = `secontregmatrixout' \ `overallse' 
					mat `outmatrixsp' = `spcontregmatrixout' \ `overallsp' 
				}*/
			}
			else {
				mat rownames `overallse' = "Sensitivity"
				mat rownames `overallsp' = "Specificity"
				mat `outmatrixse' =  `overallse' 
				mat `outmatrixsp' = `overallsp' 
				local rownamesmaxlen = max(`rownamesmaxlen', 12) /*Check if there is a longer name*/
			}
			
			mat `outmatrix' = `outmatrixse' \ `outmatrixsp'
			
			if "`expit'" == "" {
				mat colnames `outmatrix' = `sumstat' SE z P>|z| Lower Upper
				mat colnames `outmatrixse' = `sumstat' SE z P>|z| Lower Upper
				mat colnames `outmatrixsp' = `sumstat' SE z P>|z| Lower Upper
			}
			else {
				mat colnames `outmatrix' = `sumstat' SE(logit) z(logit) P>|z| Lower Upper
				mat colnames `outmatrixse' = `sumstat' SE(logit) z(logit) P>|z| Lower Upper
				mat colnames `outmatrixsp' = `sumstat' SE(logit) z(logit) P>|z| Lower Upper
			}
			
		if "`tradeoff'" != "" {
			return matrix cutmat = `cutmat'
		}					  
		return matrix outmatrixse = `outmatrixse'
		return matrix outmatrixsp = `outmatrixsp'		
		return matrix outmatrix = `outmatrix'
		return matrix Vmatrix = `Vmatrix'
	end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: ESTR +++++++++++++++++++++++++
							Estimate RR after modelling
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop estr
	program define estr, rclass
	version 13.1
		syntax [if], estimates(string) [catreg(varlist) varx(varname) typevarx(string) sumstat(string) comparator(varname) ///
		se(varname) level(string) dp(string) power(string) ///
		cveffect(string) general cbnetwork abnetwork stratify comparative by(varname) ///
		regexpression(string) baselevel(integer 1) refpos(string) cimethod(string)]
		
		marksample touse, strok
		
		//Approximate sampling distribution critical value
		qui estimates restore `estimates'
		local df = e(N) -  e(k)
				
		if "`cimethod'" != "wald" {
			local critvalue invttail(`df', `=(100-`level')/200')
		}
		else {
			local critvalue -invnorm((100-`level')/200)
		}
		
		tokenize `regexpression'
		if "`cbnetwork'" != "" {
			if "`interaction'" != "" {
				tokenize `2', parse("#")
				tokenize `1', parse(".")
			 }
			 else {
				tokenize `3', parse("#")
				tokenize `1', parse(".")
			 }

			local index "`3'"
			*if "`by'" == "" {
			*	local by = "`3'"
			*}
			*local varx //nullify
		}
		
		if "`abnetwork'`general'" != "" {
			//nullify
			local varx
			local typevarx
		}
		
		if "`comparative'" != "" {
			*tokenize `catreg'
			*macro shift
			*local catreg "`*'"
			local idpairconcat "#`varx'"
			local typevarx "i"
			*Nullify by			
			if "`stratify'" != "" {
				local by
			}
		}
		*if "`by'`comparator'" != ""  {
			local confounders "`by' `catreg'"
		/*}
		else {
			local confounders "`catreg'"
		}*/
		local marginlist
		while "`catreg'" != "" {
			tokenize `catreg'
			local marginlist = `"`marginlist' `1'`idpairconcat'"'
			macro shift 
			local catreg "`*'"
		}
		
		tempname serow sprow  lRRcoef lRRV RRoutmatrix RRoutmatrixse RRoutmatrixsp   ///
				RRoutmatrixr overallseRR overallspRR setestnlRR sptestnlRR  ///
				lORcoef lORV ORoutmatrix ORoutmatrixse ORoutmatrixsp   ///
				ORoutmatrixr overallseOR overallspOR setestnlOR sptestnlOR  
		
		if "`marginlist'" != "" | ("`varx'" != "" & "`by'" != "" ){
			local EstRRlnexpression  //log RR
			local EstORlnexpression  //log OR
			
			foreach c of local confounders {	
				qui label list `c'
				local nlevels = r(max)
				local sp_test_`c'
				local se_test_`c'
				
				if "`typevarx'" == "i" {
					forvalues l = 1/`nlevels' {
						if `l' == 1 {
							local sp_test_`c' = "_b[sp_`c'_`l']"
							local se_test_`c' = "_b[se_`c'_`l']"
						}
						else {
							local sp_test_`c' = "_b[sp_`c'_`l'] = `sp_test_`c''"
							local se_test_`c' = "_b[se_`c'_`l'] = `se_test_`c''"
						}
						local EstRRlnexpression = "`EstRRlnexpression' (sp_`c'_`l': ln(invlogit(_b[0.`se'#`l'.`c'#2.`varx'])) - ln(invlogit(_b[0.`se'#`l'.`c'#1.`varx'])))"	
						local EstRRlnexpression = "`EstRRlnexpression' (se_`c'_`l': ln(invlogit(_b[1.`se'#`l'.`c'#2.`varx'])) - ln(invlogit(_b[1.`se'#`l'.`c'#1.`varx'])))"	
						
						local EstORlnexpression = "`EstORlnexpression' (sp_`c'_`l': _b[0.`se'#`l'.`c'#2.`varx'] - _b[0.`se'#`l'.`c'#1.`varx'])"	
						local EstORlnexpression = "`EstORlnexpression' (se_`c'_`l': _b[1.`se'#`l'.`c'#2.`varx'] - _b[1.`se'#`l'.`c'#1.`varx'])"	
					}
				}
				else {
					local sp_test_`c' = "_b[sp_`c'_`baselevel']"
					local se_test_`c' = "_b[se_`c'_`baselevel']"
					local init 1
					
					forvalues l = 1/`nlevels' {
						if "abnetwork" !="" {
							if `l' == 1 {
								local sp_test_`c' = "_b[sp_`c'_`l']"
								local se_test_`c' = "_b[se_`c'_`l']"
							}
							else {
								local sp_test_`c' = "_b[sp_`c'_`l'] = `sp_test_`c''"
								local se_test_`c' = "_b[se_`c'_`l'] = `se_test_`c''"
							}
						}
						else {
							if `l' != `baselevel' {
								if `init' == 1 {
									local sp_test_`c' = "_b[sp_`c'_`l']"
									local se_test_`c' = "_b[se_`c'_`l']"
									local init 0
								}
								else {
									local sp_test_`c' = "_b[sp_`c'_`l'] = `sp_test_`c''"
									local se_test_`c' = "_b[se_`c'_`l'] = `se_test_`c''"
								}
							}
						}
						if "`refpos'" == "top" {
							local EstRRlnexpression = "`EstRRlnexpression' (sp_`c'_`l': ln(invlogit(_b[0.`se'#`baselevel'.`c'])) - ln(invlogit(_b[0.`se'#`l'.`c'])))"	
							local EstRRlnexpression = "`EstRRlnexpression' (se_`c'_`l': ln(invlogit(_b[1.`se'#`baselevel'.`c'])) - ln(invlogit(_b[1.`se'#`l'.`c'])))"	
							
							local EstORlnexpression = "`EstORlnexpression' (sp_`c'_`l': _b[0.`se'#`baselevel'.`c'] - _b[0.`se'#`l'.`c'])"	
							local EstORlnexpression = "`EstORlnexpression' (se_`c'_`l': _b[1.`se'#`baselevel'.`c'] - _b[1.`se'#`l'.`c'])"	
						}
						else {
							local EstRRlnexpression = "`EstRRlnexpression' (sp_`c'_`l': ln(invlogit(_b[0.`se'#`l'.`c'])) - ln(invlogit(_b[0.`se'#`baselevel'.`c'])))"	
							local EstRRlnexpression = "`EstRRlnexpression' (se_`c'_`l': ln(invlogit(_b[1.`se'#`l'.`c'])) - ln(invlogit(_b[1.`se'#`baselevel'.`c'])))"	
							
							local EstORlnexpression = "`EstORlnexpression' (sp_`c'_`l': _b[0.`se'#`l'.`c'] - _b[0.`se'#`baselevel'.`c'])"	
							local EstORlnexpression = "`EstORlnexpression' (se_`c'_`l': _b[1.`se'#`l'.`c'] - _b[1.`se'#`baselevel'.`c'])"	
						}
						
					}
				
				}
			}
			//RR
			qui estimates restore `estimates'
			if "`marginlist'" != "" {
				qui margins `marginlist' if `touse', predict(xb) over(`se') post level(`level')
			}
			if "`varx'" != "" & "`by'" != ""  {
				qui margins `varx' if `touse', predict(xb) over(`se' `by') post level(`level')
			}
			
			qui nlcom `EstRRlnexpression', post level(`level') iterate(200)
			mat `lRRcoef' = e(b)
			mat `lRRV' = e(V)
			mat `lRRV' = vecdiag(`lRRV')	
			local ncols = colsof(`lRRcoef') //length of the vector
			local rnames :colnames `lRRcoef'
			
			local rowtestnl			
			local i = 1
			
			foreach c of local confounders {
				qui label list `c'
				local nlevels = r(max)
				if (`nlevels' > 2 & "`comparative'" == "") | (`nlevels' > 1 & ("`comparative'" != "" | "`cbnetwork'" != "" )){
					qui testnl (`se_test_`c''), iterate(200)
					local se_testnl_`c'_chi2 = r(chi2)				
					local se_testnl_`c'_df = r(df)
					local se_testnl_`c'_p = r(p)
					qui testnl (`sp_test_`c'')
					local sp_testnl_`c'_chi2 = r(chi2)
					local sp_testnl_`c'_df = r(df)
					local sp_testnl_`c'_p = r(p)
					if `i'==1 {
						mat `setestnlRR' =  [`se_testnl_`c'_chi2', `se_testnl_`c'_df', `se_testnl_`c'_p']
						mat `sptestnlRR' =  [`sp_testnl_`c'_chi2', `sp_testnl_`c'_df', `sp_testnl_`c'_p']
					}
					else {
						mat `setestnlRR' = `setestnlRR' \ [`se_testnl_`c'_chi2', `se_testnl_`c'_df', `se_testnl_`c'_p']
						mat `sptestnlRR' = `sptestnlRR' \ [`sp_testnl_`c'_chi2', `sp_testnl_`c'_df', `sp_testnl_`c'_p']
					}
					 
					local ++i
					local rowtestnl = "`rowtestnl' `c' "
				}
			}
			
			//OR
			qui estimates restore `estimates'
			if "`marginlist'" != "" {
				qui margins `marginlist' if `touse', predict(xb) over(`se') post level(`level')
			}
			if "`varx'" != "" & "`by'" != ""  {
				qui margins `varx' if `touse', predict(xb) over(`se' `by') post level(`level')
			}
			
			qui nlcom `EstORlnexpression', post level(`level') iterate(200)
			mat `lORcoef' = e(b)
			mat `lORV' = e(V)
			mat `lORV' = vecdiag(`lORV')	
			local ncols = colsof(`lORcoef') //length of the vector
			local rnames :colnames `lORcoef'
			
			local i = 1
			
			foreach c of local confounders {
				qui label list `c'
				local nlevels = r(max)
				if (`nlevels' > 2 & "`comparative'" == "") | (`nlevels' > 1 & ("`comparative'" != "" | "`cbnetwork'" != "" )){
					qui testnl (`se_test_`c''), iterate(200)
					local se_testnl_`c'_chi2 = r(chi2)				
					local se_testnl_`c'_df = r(df)
					local se_testnl_`c'_p = r(p)
					qui testnl (`sp_test_`c'')
					local sp_testnl_`c'_chi2 = r(chi2)
					local sp_testnl_`c'_df = r(df)
					local sp_testnl_`c'_p = r(p)
					if `i'==1 {
						mat `setestnlOR' =  [`se_testnl_`c'_chi2', `se_testnl_`c'_df', `se_testnl_`c'_p']
						mat `sptestnlOR' =  [`sp_testnl_`c'_chi2', `sp_testnl_`c'_df', `sp_testnl_`c'_p']
					}
					else {
						mat `setestnlOR' = `setestnlOR' \ [`se_testnl_`c'_chi2', `se_testnl_`c'_df', `se_testnl_`c'_p']
						mat `sptestnlOR' = `sptestnlOR' \ [`sp_testnl_`c'_chi2', `sp_testnl_`c'_df', `sp_testnl_`c'_p']
					}
					 
					local ++i
				}
			}
						
			if `i' > 1 {
				mat rownames `setestnlRR' = `rowtestnl'
				mat rownames `sptestnlRR' = `rowtestnl'
				mat colnames `setestnlRR' = chi2 df p
				mat colnames `sptestnlRR' = chi2 df p
				
				mat roweq `setestnlRR' = Relative_Sensitivity
				mat roweq `sptestnlRR' = Relative_Specificity
				
				mat rownames `setestnlOR' = `rowtestnl'
				mat rownames `sptestnlOR' = `rowtestnl'
				mat colnames `setestnlOR' = chi2 df p
				mat colnames `sptestnlOR' = chi2 df p
				
				mat roweq `setestnlOR' = OR_Sensitivity
				mat roweq `sptestnlOR' = OR_Specificity
								
				*mat `testmat2print' =  `setestnl'  \ `sptestnl' 
				*mat colnames `testmat2print' = chi2 df p
				
				local inltest = "yes"
			}
			else {
				local inltest = "no"
			}
			
			if "`comparative'" != ""  | "`cbnetwork'" != ""  {
				mat `RRoutmatrix' = J(`=`ncols' + 2', 6, .)
				mat `ORoutmatrix' = J(`=`ncols' + 2', 6, .)
			}
			else {
				mat `RRoutmatrix' = J(`ncols', 6, .)
				mat `ORoutmatrix' = J(`ncols', 6, .)
			}
			
			forvalues r = 1(1)`ncols' {
				mat `RRoutmatrix'[`r', 1] = exp(`lRRcoef'[1,`r']) /*Estimate*/
				mat `RRoutmatrix'[`r', 2] = sqrt(`lRRV'[1, `r']) /*se in log scale, power 1*/
				mat `RRoutmatrix'[`r', 3] = `lRRcoef'[1,`r']/sqrt(`lRRV'[1, `r']) /*Z in log scale*/
				mat `RRoutmatrix'[`r', 4] =  normprob(-abs(`RRoutmatrix'[`r', 3]))*2  /*p-value*/
				mat `RRoutmatrix'[`r', 5] = exp(`lRRcoef'[1, `r'] - `critvalue' * sqrt(`lRRV'[1, `r'])) /*lower*/
				mat `RRoutmatrix'[`r', 6] = exp(`lRRcoef'[1, `r'] + `critvalue' * sqrt(`lRRV'[1, `r'])) /*upper*/
				
				mat `ORoutmatrix'[`r', 1] = exp(`lORcoef'[1,`r']) /*Estimate*/
				mat `ORoutmatrix'[`r', 2] = sqrt(`lORV'[1, `r']) /*se in log scale, power 1*/
				mat `ORoutmatrix'[`r', 3] = `lORcoef'[1,`r']/sqrt(`lORV'[1, `r']) /*Z in log scale*/
				mat `ORoutmatrix'[`r', 4] =  normprob(-abs(`ORoutmatrix'[`r', 3]))*2  /*p-value*/
				mat `ORoutmatrix'[`r', 5] = exp(`lORcoef'[1, `r'] - `critvalue' * sqrt(`lORV'[1, `r'])) /*lower*/
				mat `ORoutmatrix'[`r', 6] = exp(`lORcoef'[1, `r'] + `critvalue' * sqrt(`lORV'[1, `r'])) /*upper*/
			}
		}
		else {
			mat `RRoutmatrix' = J(2, 6, .)
			mat `ORoutmatrix' = J(2, 6, .)
			local ncols = 0
		}
		
		if "`typevarx'" == "i" {
			//RR
			qui estimates restore `estimates'
			qui margins `varx' if `touse', predict(xb) over(`se') post level(`level')
					
			//log metric
			if `baselevel' == 1 {
				local indexlevel "2"
			}
			else {
				local indexlevel "1"
			}
			if "`refpos'" == "bottom" {
				qui nlcom (sp_Overall: ln(invlogit(_b[0.`se'#`indexlevel'.`varx'])) - ln(invlogit(_b[0.`se'#`baselevel'.`varx']))) ///
					  (se_Overall: ln(invlogit(_b[1.`se'#`indexlevel'.`varx'])) - ln(invlogit(_b[1.`se'#`baselevel'.`varx']))), iterate(200)
			}
			else {
				qui nlcom (sp_Overall: ln(invlogit(_b[0.`se'#`baselevel'.`varx'])) - ln(invlogit(_b[0.`se'#`indexlevel'.`varx']))) ///
						  (se_Overall: ln(invlogit(_b[1.`se'#`baselevel'.`varx'])) - ln(invlogit(_b[1.`se'#`indexlevel'.`varx']))), iterate(200)
			}			
			mat `lRRcoef' = r(b)
			mat `lRRV' = r(V)
			mat `lRRV' = vecdiag(`lRRV')
			
			forvalues r=1(1)2 {
				mat `RRoutmatrix'[`=`ncols' + `r'', 1] = exp(`lRRcoef'[1,`r'])  //rr
				mat `RRoutmatrix'[`=`ncols' + `r'', 2] = sqrt(`lRRV'[1, `r']) //se
				mat `RRoutmatrix'[`=`ncols' + `r'', 3] = `lRRcoef'[1, `r']/sqrt(`lRRV'[1, `r']) //zvalue
								if "`cimethod'" != "wald" {
					mat `RRoutmatrix'[`=`ncols' + `r'', 4] = normprob(-abs(`RRoutmatrix'[`=`ncols' + `r'', 3]))*2 //pvalue
				}
				else {
					mat `RRoutmatrix'[`=`ncols' + `r'', 4] = ttail(`df', abs(`RRoutmatrix'[`=`ncols' + `r'', 3]))*2  //pvalue
				}
				mat `RRoutmatrix'[`=`ncols' + `r'', 5] = exp(`lRRcoef'[1, `r'] - `critvalue'*sqrt(`lRRV'[1, `r'])) //ll
				mat `RRoutmatrix'[`=`ncols' + `r'', 6] = exp(`lRRcoef'[1, `r'] + `critvalue'*sqrt(`lRRV'[1, `r'])) //ul
			}
			local rnames = "`rnames' sp_Overall se_Overall"
			
			//OR
			qui estimates restore `estimates'
			qui margins `varx' if `touse', predict(xb) over(`se') post level(`level')
					
			//log metric
			if `baselevel' == 1 {
				local indexlevel "2"
			}
			else {
				local indexlevel "1"
			}
			if "`refpos'" == "bottom" {
				qui nlcom (sp_Overall: _b[0.`se'#`indexlevel'.`varx'] - _b[0.`se'#`baselevel'.`varx']) ///
					  (se_Overall: _b[1.`se'#`indexlevel'.`varx'] - _b[1.`se'#`baselevel'.`varx']), iterate(200)
			}
			else {
				qui nlcom (sp_Overall: _b[0.`se'#`baselevel'.`varx'] - _b[0.`se'#`indexlevel'.`varx']) ///
						  (se_Overall: _b[1.`se'#`baselevel'.`varx'] - _b[1.`se'#`indexlevel'.`varx']), iterate(200)
			}			
			mat `lORcoef' = r(b)
			mat `lORV' = r(V)
			mat `lORV' = vecdiag(`lORV')
			
			forvalues r=1(1)2 {
				mat `ORoutmatrix'[`=`ncols' + `r'', 1] = exp(`lORcoef'[1,`r'])  //rr
				mat `ORoutmatrix'[`=`ncols' + `r'', 2] = sqrt(`lORV'[1, `r']) //se
				mat `ORoutmatrix'[`=`ncols' + `r'', 3] = `lORcoef'[1, `r']/sqrt(`lORV'[1, `r']) //zvalue
				if "`cimethod'" != "wald" {
					mat `ORoutmatrix'[`=`ncols' + `r'', 4] = normprob(-abs(`ORoutmatrix'[`=`ncols' + `r'', 3]))*2 //pvalue
				}
				else {
					mat `ORoutmatrix'[`=`ncols' + `r'', 4] = ttail(`df', abs(`ORoutmatrix'[`=`ncols' + `r'', 3]))*2  //pvalue
				}
				mat `ORoutmatrix'[`=`ncols' + `r'', 5] = exp(`lORcoef'[1, `r'] - `critvalue'*sqrt(`lORV'[1, `r'])) //ll
				mat `ORoutmatrix'[`=`ncols' + `r'', 6] = exp(`lORcoef'[1, `r'] + `critvalue'*sqrt(`lORV'[1, `r'])) //ul
			}

		}
		
		local sprownames = ""
		local serownames = ""
		local rspec = "-" /*draw lines or not between the rows*/
		local rownamesmaxlen = 10 /*Default*/
		
		local nrows = rowsof(`RRoutmatrix')

		local initse 0
		local initsp 0
		forvalues r = 1(1)`ncols' {
			local rname`r':word `r' of `rnames'
			tokenize `rname`r'', parse("_")					
			local parm = "`1'"
			local left = "`3'"
			local right = "`5'"
			mat `RRoutmatrixr' = `RRoutmatrix'[`r', 1...] //select the r'th row
			mat `ORoutmatrixr' = `ORoutmatrix'[`r', 1...] //select the r'th row
			if "`5'" != "" {
				local lab:label `left' `right'
				local lab = ustrregexra("`lab'", " ", "_")
				local nlen : strlen local lab
				local rownamesmaxlen = max(`rownamesmaxlen', min(`nlen', 32)) //Check if there is a longer name
				local `parm'rownames = "``parm'rownames' `left':`lab'" 
				if `init`parm'' == 0 {
					mat `RRoutmatrix`parm'' = `RRoutmatrixr'
					mat `ORoutmatrix`parm'' = `ORoutmatrixr'
				}
				else {
					mat `RRoutmatrix`parm'' = `RRoutmatrix`parm'' \ `RRoutmatrixr'
					mat `ORoutmatrix`parm'' = `ORoutmatrix`parm'' \ `ORoutmatrixr'
				}
				local init`parm' 1
			}
		}
		if `ncols' > 0 {
			mat rownames `RRoutmatrixse' = `serownames'
			mat rownames `RRoutmatrixsp' = `sprownames'
			mat rownames `ORoutmatrixse' = `serownames'
			mat rownames `ORoutmatrixsp' = `sprownames'
		}
		*mat out = `outmatrix'	
		mat `serow' = J(1, 6, .)
		mat `sprow' = J(1, 6, .)
		
		mat rownames `serow' = "Relative Sensitivity"
		mat rownames `sprow' = "Relative Specificity"  //20 chars
		local rownamesmaxlen = max(`rownamesmaxlen', 21) //Check if there is a longer name
		
		if ("`comparative'" !="" | "`cbnetwork'" != "") {
			mat `overallspRR' = `RRoutmatrix'[`=`nrows'-1', 1...]
			mat `overallseRR' = `RRoutmatrix'[`nrows', 1...]
			
			mat rownames `overallseRR' = "Overall"
			mat rownames `overallspRR' = "Overall"
			
			mat `overallspOR' = `ORoutmatrix'[`=`nrows'-1', 1...]
			mat `overallseOR' = `ORoutmatrix'[`nrows', 1...]
			
			mat rownames `overallseOR' = "Overall"
			mat rownames `overallspOR' = "Overall"
		}

		if `ncols' > 0 & ("`comparative'`cbnetwork'" != ""){
			if "`cveffect'" == "sesp" {
				*local rspec "---`="&"*`=`nrows'/2 - 1''--`="&"*`=`nrows'/2 - 1''-"
				mat `RRoutmatrix' = `serow' \ `RRoutmatrixse' \ `overallseRR'  \ `sprow' \ `RRoutmatrixsp' \ `overallspRR'
				mat `ORoutmatrix' = `serow' \ `ORoutmatrixse' \ `overallseOR'  \ `sprow' \ `ORoutmatrixsp' \ `overallspOR'
			}
			else if "`cveffect'" == "se" { 
				mat `RRoutmatrix' = `serow' \ `RRoutmatrixse' \ `overallseRR'  \ `sprow' \  `overallspRR'
				mat `ORoutmatrix' = `serow' \ `ORoutmatrixse' \ `overallseOR'  \ `sprow' \  `overallspOR'
				}
			else {
				mat `RRoutmatrix' = `serow' \ `overallseRR'  \ `sprow' \ `RRoutmatrixsp' \ `overallspRR'
				mat `ORoutmatrix' = `serow' \ `overallseOR'  \ `sprow' \ `ORoutmatrixsp' \ `overallspOR'
			}
				*local rspec "--`="&"*`=`nrows'/2''-"
			
			mat `RRoutmatrixse' = `RRoutmatrixse' \ `overallseRR' 
			mat `RRoutmatrixsp' = `RRoutmatrixsp' \ `overallspRR' 
			mat `ORoutmatrixse' = `ORoutmatrixse' \ `overallseOR' 
			mat `ORoutmatrixsp' = `ORoutmatrixsp' \ `overallspOR' 
		}
		if `ncols' > 0 &  "`cbnetwork'`comparative'" =="" {
			if "`cveffect'" == "sesp" {
				*local rspec "---`="&"*`=`nrows'/2 - 1''--`="&"*`=`nrows'/2 - 1''-"
				mat `RRoutmatrix' = `serow' \ `RRoutmatrixse'  \ `sprow' \ `RRoutmatrixsp'
				mat `RRoutmatrixse' = `serow' \ `RRoutmatrixse'
				mat `RRoutmatrixsp' = `sprow'  \ `RRoutmatrixsp'
				
				mat `ORoutmatrix' = `serow' \ `ORoutmatrixse'  \ `sprow' \ `ORoutmatrixsp'
				mat `ORoutmatrixse' = `serow' \ `ORoutmatrixse'
				mat `ORoutmatrixsp' = `sprow'  \ `ORoutmatrixsp'
			}
			else {
				if "`cveffect'" == "se" { 
					mat `RRoutmatrix' = `serow' \ `RRoutmatrixse'
					mat `RRoutmatrixse' = `serow' \ `RRoutmatrixse'
					mat `RRoutmatrixsp' = J(1, 6, .)
					
					mat `ORoutmatrix' = `serow' \ `ORoutmatrixse'
					mat `ORoutmatrixse' = `serow' \ `ORoutmatrixse'
					mat `ORoutmatrixsp' = J(1, 6, .)
				}
				else {
					mat `RRoutmatrix' = `sprow'  \ `RRoutmatrixsp'
					mat `RRoutmatrixsp' = `sprow'  \ `RRoutmatrixsp'
					mat `RRoutmatrixse' = J(1, 6, .)
					
					mat `ORoutmatrix' = `sprow'  \ `ORoutmatrixsp'
					mat `ORoutmatrixsp' = `sprow'  \ `ORoutmatrixsp'
					mat `ORoutmatrixse' = J(1, 6, .)
				}
			}		
		}
		if `ncols' == 0 {
			mat `RRoutmatrixse' =  `overallseRR' 
			mat `RRoutmatrixsp' = `overallspRR'
			mat `RRoutmatrix' = `serow' \ `RRoutmatrixse'  \ `sprow' \ `RRoutmatrixsp'
			
			mat `ORoutmatrixse' =  `overallseOR' 
			mat `ORoutmatrixsp' = `overallspOR'
			mat `ORoutmatrix' = `serow' \ `ORoutmatrixse'  \ `sprow' \ `ORoutmatrixsp'			
		}	

		mat colnames `RRoutmatrixse' = `sumstat' SE(lor) z(lor) P>|z| Lower Upper
		mat colnames `RRoutmatrixsp' = `sumstat' SE(lor) z(lor) P>|z| Lower Upper
		mat colnames `RRoutmatrix' = `sumstat' SE(lor) z(lor) P>|z| Lower Upper
		
		mat colnames `ORoutmatrixse' = `sumstat' SE(lor) z(lor) P>|z| Lower Upper
		mat colnames `ORoutmatrixsp' = `sumstat' SE(lor) z(lor) P>|z| Lower Upper
		mat colnames `ORoutmatrix' = `sumstat' SE(lor) z(lor) P>|z| Lower Upper

		if "`inltest'" == "yes" {
			return matrix setestnlRR = `setestnlRR'
			return matrix sptestnlRR = `sptestnlRR'
			
			return matrix setestnlOR = `setestnlOR'
			return matrix sptestnlOR = `sptestnlOR'
		}
		return local inltest = "`inltest'"
		return matrix RRoutmatrix = `RRoutmatrix'
		return matrix RRoutmatrixse = `RRoutmatrixse'
		return matrix RRoutmatrixsp = `RRoutmatrixsp'
		return matrix ORoutmatrix = `ORoutmatrix'
		return matrix ORoutmatrixse = `ORoutmatrixse'
		return matrix ORoutmatrixsp = `ORoutmatrixsp'
	end	
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: ESTCOVAR +++++++++++++++++++++++++
							Compose the var-cov matrix
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
cap program drop estcovar
program define estcovar, rclass
version 14.0

	syntax, matrix(name) model(string) [ bcov(string) wcov(string) abnetwork cbnetwork comparative general ]
	*matrix is colvector
	tempname matcoef BVar WVar
	mat `matcoef' = `matrix''
	local nrows = rowsof(`matcoef')
	*Initialize - Default
	mat	`BVar' = (0, 0\ ///
				0, 0)
	mat	`WVar' = (0, 0\ ///
				0, 0)
	local b = 0
	local w = 0	
	
	if "`model'" == "random" {
		*WVAR
		if "`wcov'" != "" {
			if strpos("`wcov'", "ind") != 0 {
				mat	`WVar' = (exp(`matcoef'[ `nrows' - 1, 1])^2, 0\ ///
							0, exp(`matcoef'[ `nrows', 1])^2)
				local w = 2
			}
			else if strpos("`wcov'", "id") != 0 {
				mat	`WVar' = (exp(`matcoef'[ `nrows', 1])^2, 0 \ ///
					0, exp(`matcoef'[ `nrows', 1])^2)
				local w = 1
			}
		}

		*BVAR
		if strpos("`bcov'", "uns") != 0 {
			mat	`BVar' = (exp(`matcoef'[`nrows' - 2 - `w', 1])^2, exp(`matcoef'[ `nrows' - 1 - `w', 1])*exp(`matcoef'[`nrows' - 2 - `w', 1])*tanh(`matcoef'[ `nrows' - `w', 1])\ ///
						exp(`matcoef'[ `nrows' - 1 - `w', 1])*exp(`matcoef'[`nrows' - 2 - `w', 1])*tanh(`matcoef'[ `nrows' - `w', 1]), exp(`matcoef'[ `nrows' - 1 - `w', 1])^2)
			local b = 3
		}		
		else if strpos("`bcov'", "ind") != 0 {
			mat	`BVar' = (exp(`matcoef'[ `nrows' - 1 - `w', 1])^2, 0\ ///
						0, exp(`matcoef'[ `nrows' - `w', 1])^2)
			local b = 2
		}
		else if strpos("`bcov'", "exc") != 0 {
			mat	`BVar' = (exp(`matcoef'[ `nrows' - 1 - `w', 1])^2, exp(`matcoef'[ `nrows' - 1 - `w', 1])*exp(`matcoef'[ `nrows' - 1 - `w', 1])*tanh(`matcoef'[ `nrows' - `w', 1])\ ///
						exp(`matcoef'[ `nrows' - 1 - `w', 1])*exp(`matcoef'[ `nrows' - 1 - `w', 1])*tanh(`matcoef'[ `nrows' - `w', 1]), exp(`matcoef'[ `nrows' - 1 - `w', 1])^2)
			local b = 2
		}
		else if (strpos("`bcov'", "id") != 0) {
			mat	`BVar' = (exp(`matcoef'[ `nrows' - `w', 1])^2, 0\ ///
				0, exp(`matcoef'[ `nrows' - `w', 1])^2)
				
			local b = 1
		}
		else if (strpos("`bcov'", "sp") != 0) {
			mat	`BVar' = (0, 0\ ///
				0, exp(`matcoef'[ `nrows' - `w', 1])^2)
				
			local b = 1
		}
		else if (strpos("`bcov'", "se") != 0) {
			mat	`BVar' = (exp(`matcoef'[ `nrows' - `w', 1])^2, 0\ ///
				0, 0)
				
			local b = 1
		}
	}
		
	local k = `b' + `w'
	
	return matrix WVar = `WVar' 
	return matrix BVar = `BVar' 
	return local k = `k' 
end
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: 	KOOPMANCI +++++++++++++++++++++++++
								CI for RR
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop koopmanci
	program define koopmanci
	version 14.0

		syntax varlist, R(name) lowerci(name) upperci(name) [alpha(real 0.05)]
		
		qui {	
			tokenize `varlist'
			gen `r' = . 
			gen `lowerci' = .
			gen `upperci' = .
			
			count
			forvalues i = 1/`r(N)' {
				local n1 = `1'[`i']
				local N1 = `2'[`i']
				local n2 = `3'[`i']
				local N2 = `4'[`i']

				koopmancii `n1' `N1' `n2' `N2', alpha(`alpha')
				mat ci = r(ci)
				
				if (`n1' == 0) &(`n2'==0) {
					replace `r' = 0 in `i'
				}
				else {
					replace `r' = (`n1'/`N1')/(`n2'/`N2')  in `i'	
				}
				replace `lowerci' = ci[1, 1] in `i'
				replace `upperci' = ci[1, 2] in `i'
			}
		}
	end
	
	/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS: KOOPMANCII +++++++++++++++++++++++++
								CI for RR
	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop koopmancii
	program define koopmancii, rclass
	version 14.0
		syntax anything(name=data id="data"), [alpha(real 0.05)]
		
		local len: word count `data'
		if `len' != 4 {
			di as error "Specify full data: n1 N1 n2 N2"
			exit
		}
		
		foreach num of local data {
			cap confirm integer number `num'
			if _rc != 0 {
				di as error "`num' found where integer expected"
				exit
			}
		}
		
		tokenize `data'
		cap assert ((`1' <= `2') & (`3' <= `4'))
		if _rc != 0{
			di as err "Order should be n1 N1 n2 N2"
			exit _rc
		}
		
		mata: koopman_ci((`1', `2', `3', `4'), `alpha')
		
		return matrix ci = ci
		return scalar alpha = `alpha'	

	end
/*	SUPPORTING FUNCTIONS: FPLOTCHECK ++++++++++++++++++++++++++++++++++++++++++
			Advance housekeeping for the fplot
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	capture program drop fplotcheck
	program define fplotcheck, rclass
	version 14.1	
	#delimit ;
	syntax  [,
		/*Passed from top*/
		COMParative 
		cbnetwork
		abnetwork
		general
		noWT
		/*passed via foptions*/
		AStext(integer 50) 				
		CIOpt(passthru) 
		DIAMopt(passthru) 
		DOUble 
 		LCols(varlist) 
		noOVLine 
		noSTATS
		ARRowopt(passthru) 		
		OLineopt(passthru) 
		OUTplot(string)  //or|rr
		PLOTstat(passthru) //comma seperated
		POINTopt(passthru) 
		SUBLine 
		TEXts(real 1.5) 
		XLIne(passthru)	/*silent option*/
		XLAbel(passthru) 
		XTick(passthru) 
		GRID
		GRAphsave(passthru)
		logscale
		first(varname)
		by(varname) 
		*
	  ];
	#delimit cr
	
		if `astext' > 90 | `astext' < 10 {
		di as error "Percentage of graph as text (ASTEXT) must be within 10-90%"
		di as error "Must have some space for text and graph"
		exit
	}
	if `texts' < 0 {
		di as res "Warning: Negative text size (TEXTSize) are ignored"
		local texts 1
	}	
	
	if "`outplot'" == "" {
		local outplot abs
	}
	else {
		local outplot = strlower("`outplot'")
		local rc_ = ("`outplot'" == "or") + ("`outplot'" == "rr") + ("`outplot'" == "abs")
		if `rc_' != 1 {
			di as error "Options outplot(`outplot') incorrectly specified"
			di as error "Allowed options: abs, rr, or"
			exit
		}
		/*if "`outplot'" != "abs" {
			cap assert "`abnetwork'`comparative'`cbnetwork'" != "" 
			if _rc != 0 {
				di as error "Option outplot(`outplot') only avaialable with abnetwork/comparative/cbnetwork analysis"
				di as error "Specify analysis with -abnetwork/comparative/cbnetwork- option"
				exit _rc
			}
		}
		if ("`first'" != "" & "`by'" != "") & "`comparative'" != "" {
			cap assert ("`first'" != "`by'") & ("`output'" == "abs")
			if _rc != 0 { 
					di as error "Remove the option by(`by') or specify a different by-variable"
					exit _rc
			}
		}*/
	}
	foreach var of local lcols {
		cap confirm var `var'
		if _rc!=0  {
			di in re "Variable `var' not in the dataset"
			exit _rc
		}
	}
	if "`lcols'" =="" {
		local lcols " "
	}
	if "`astext'" != "" {
		local astext "astext(`astext')"
	}
	if "`texts'" != "" {
		local texts "texts(`texts')"
	}
	local foptions `"`astext' `ciopt' `diamopt' `arrowopt' `double' `ovline' `stats' `wt' `olineopt' `plotstat' `pointopt' `subline' `texts' `xlabel' `xtick' `grid' `xline'  `logscale' `graphsave' `options'"'
	return local outplot = "`outplot'"
	return local lcols ="`lcols'"
	return local foptions = `"`foptions'"'
end

/*	SUPPORTING FUNCTIONS: FPLOT ++++++++++++++++++++++++++++++++++++++++++++++++
			The forest plot
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Some re-used code from metaprop

	capture program drop fplot
	program define fplot
	version 14.1	
	#delimit ;
	syntax varlist [if] [in] [,
	    /*Passed from top options*/
		STudyid(varname)
		POWer(integer 0)
		DP(integer 2) 
		Level(integer 95)
		smooth
		abnetwork
		cbnetwork
		comparative
		general
		cveffect(string)
		varxlabs(string)
		/*passed from within*/	
		Groupvar(varname)
		/*passed via foptions*/
		AStext(string)
		ARRowopt(string) 		
		CIOpt(string) 
		DIAMopt(string) 
		DOUble 
 		LCols(varlist) 
		noOVLine 
		noSTATS
		noWT
		noBox
		summaryonly
		OLineopt(string) 
		OUTplot(string) 
		PLOTstat(string asis) /*comma seperated*/
		POINTopt(string) 
		SUBLine 
		TEXts(real 1.5) 
		XLIne(string asis)
		XLAbel(string) 
		XTick(string)
		GRID
		GRAPHSave(string asis)
		logscale
		*
	  ];
	#delimit cr
	
	local foptions `"`options'"'
	if strpos(`"`foptions'"', "graphregion") == 0 {
			local foptions `"graphregion(color(white)) `foptions'"'
		}
	
	tempvar es lci uci modeles modellci modeluci use label tlabel id newid se  df  expand orig ///
		
	tokenize "`varlist'", parse(" ")

	qui {
		gen `es'=`1'*(10^`power')
		gen `lci'   =`2'*(10^`power')
		gen `uci'   =`3'*(10^`power')
		gen byte `use'=`4'
		gen str `label'=`5'
		gen `df' = `6'
		gen `id' = `7'
		gen `se' = `8'
		
		if "`smooth'" !="" {
			gen `modeles' 	= `9'*(10^`power')
			gen `modellci' 	= `10'*(10^`power')
			gen `modeluci' 	= `11'*(10^`power')
		}
		
		if "`plotstat'"=="" {
			local outplot = strlower("`outplot'")
			if "`outplot'" == "rr" {
				local plotstatse "Relative Sensitivity"
				local plotstatsp "Relative Specificity"
			}
			else if "`outplot'" == "or" {
				local plotstatse "Sensitivity OR"
				local plotstatsp "Specificity OR"
			}
			else {
				local plotstatse "Sensitivity"
				local plotstatsp "Specificity"
			}
		}
		else {
			tokenize "`plotstat'", parse(",")
			local plotstatse "`1'"
			local plotstatsp "`3'"
		}
							
		/*
		qui summ `id'
		gen `expand' = 1
		replace `expand' = 1 + 1*(`id'==r(min)) 
		expand `expand'
		
		replace `id' = `id' + 1 if _n>2
		replace `label' = "" if `id'==1
		replace `use' = 0 if `id'==1
		*/
		tempvar newid
		gen  `newid' = `id'
		//Add five spaces on top of the dataset and 1 space below
		*gen `expand' = 1
		gen `expand' = 1 + 5*(_n==1)  /*+ 1*(_n==_N) */
		expand `expand'
		sort `newid' `se' `use'

		replace `newid' = _n in 1/5
		replace `newid' = `newid' + 5 if _n>5
		replace `label' = "" in 1/5
		replace `use' = -2 in 1/4
		replace `use' = 0 in 5
		/*
		replace `newid' = _N  if _N==_n
		replace `use' = 0  if _N==_n
		replace `label' = "" if _N==_n
		*/
		/*
		tempvar flag
		gen `flag' = 1
		replace `flag' = 0 in 1/4
		*/
		replace `id' = `newid'
		drop `newid'
		sort `id' `se'
		
		//studylables
		if "`abnetwork'" != "" & "`outplot'" != "abs" {
			local studylb: variable label `groupvar'
			if "`studylb'" == "" {
				label var `label' "`groupvar'"
			}
			else {
				label var `label' "`studylb'"
			}
		}
		else {
			local studylb: variable label `studyid'
			if "`studylb'" == "" {
				label var `label' "`studyid'"
			}
			else {
				label var `label' "`studylb'"
			}
		}
		if "`lcols'" == "" {
			local lcols "`label'"
		}
		else {
			local lcols "`label' `lcols'"
		}
		
		/*
		egen `newid' = group(`id')
		replace `id' = `newid'
		drop `newid'
		*/

		tempvar estText modelestText estTextse estTextsp index wtText
		gen str `estText' = string(`es', "%10.`=`dp''f") + " (" + string(`lci', "%10.`=`dp''f") + ", " + string(`uci', "%10.`=`dp''f") + ")"  if (`use' == 1 | `use' == 2 | `use' == 3)
		replace `estText' = " " if `use' == 2 & ((`se' == 1 & "`cveffect'" == "sp") | (`se' == 0 & "`cveffect'" == "se"))
		
		if "`smooth'" !="" & "`summaryonly'" == ""  {
			gen str `modelestText' = string(`modeles', "%10.`=`dp''f") + " (" + string(`modellci', "%10.`=`dp''f") + ", " + string(`modeluci', "%10.`=`dp''f") + ")"  if (`use' == 1 )
			replace `modelestText' = string(`es', "%10.`=`dp''f") + " (" + string(`lci', "%10.`=`dp''f") + ", " + string(`uci', "%10.`=`dp''f") + ")"  if (`use' == 2 | `use' == 3)
			replace `modelestText' = " " if `use' == 2 & ((`se' == 1 & "`cveffect'" == "sp") | (`se' == 0 & "`cveffect'" == "se"))
			replace `estText' = " " if (`use' == 2 | `use' == 3)
		}
		if "`wt'" == "" {
			gen str `wtText' = string(_WT, "%10.`=`dp''f") if (`use' == 1 | `use' == 2 | `use' == 3) & _WT !=.
		}
		
		//RCOLS
		if "`stats'" == "" {
			tempvar estTextse estTextsp
			gen `estTextse' = `estText'  if `se' == 1
			gen `estTextsp' = `estText' if `se' == 0
			label var `estTextse' "`plotstatse' (`level'% CI)"
			label var `estTextsp' "`plotstatsp' (`level'% CI)"
			local rcolse "`estTextse'"
			local rcolsp "`estTextsp'"
		}
		if "`smooth'" !="" & "`summaryonly'" == "" {
			tempvar  modelestTextse modelestTextsp
			gen `modelestTextse' = `modelestText'  if `se'== 1
			gen `modelestTextsp' = `modelestText' if `se'== 0
			label var `modelestTextse' "Model-based Est (`level'% centile)"
			label var `modelestTextsp' "Model-based Est (`level'% centile)"
			local rcolse = "`rcolse' `modelestTextse'"
			local rcolsp = "`rcolsp' `modelestTextsp'"
		}
		if "`wt'" == "" {
			tempvar wtTextse wtTextsp
			gen `wtTextse' = `wtText' if `se'== 1
			gen `wtTextsp' = `wtText' + " " if `se'== 0
			label var `wtTextse' "% Weight"
			label var `wtTextsp' "% Weight"
			local rcolse = "`rcolse' `wtTextse'"
			local rcolsp = "`rcolsp' `wtTextsp'"
		}
		
		tempvar extra 
		gen `extra' = " "
		label var `extra' " "
		
		if "`rcolse'`rcolsp'" == "" {
			local rcolse = "`extra'" 
			local rcolsp = "`extra'"
		}
		
		// GET MIN AND MAX DISPLAY
		// SORT OUT TICKS- CODE PINCHED FROM MIKE AND FIRandomED. TURNS OUT I'VE BEEN USING SIMILAR NAMES...
		// AS SUGGESTED BY JS JUST ACCEPT ANYTHING AS TICKS AND RESPONSIBILITY IS TO USER!

		if "`logscale'" != "" {
			replace `es' 	= ln(`es')
			replace `lci' 		= ln(`lci')
			replace `uci' 		= ln(`uci')
			
			if "`smooth'" !="" { 
				replace `modeles'  	= ln(`modeles')
				replace `modellci' 		= ln(`modellci')
				replace `modeluci' 		= ln(`modeluci')
			}
		}
		qui summ `lci', detail
		local DXmin = r(min)
		qui summ `uci', detail
		local DXmax = r(max)
				
		if "`xlabel'" != "" {
			local DXmin = min(`xlabel')
			local DXmax = max(`xlabel')
			if "`logscale'" != "" {
				if "`DXmin'" == "0" {
					local DXmin = ln(`=10^(-`dp')')
				}
				else {
					local DXmin = ln(`DXmin')
				}
				if "`DXmax'" == "0" {
					local DXmax = ln(`=10^(-`dp')')
				}
				else {
					local DXmax = ln(`DXmax')
				}
			}
		}
		if "`xlabel'"=="" {
			local xlabel "0, `DXmax'"
		}

		local lblcmd ""
		tokenize "`xlabel'", parse(",")
		while "`1'" != ""{
			if "`1'" != ","{
				local lbl = string(`1',"%7.3g")
				if "`logscale'" != "" {
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
				}
				else {
					local val = `1'
				}

				local lblcmd `lblcmd' `val' "`lbl'"
			}
			mac shift
		}
		
		if "`xtick'" == ""{
			local xtick = "`xlabel'"
		}

		local xtick2 = ""
		tokenize "`xtick'", parse(",")
		while "`1'" != ""{
			if "`1'" != ","{
				if "`logscale'" != "" {
					if "`1'" == "0" {
						local val = ln(`=10^(-`dp')')
					}
					else {
						local val = ln(`1')
					}
				}
				else {
					local val = `1'
				}
				local xtick2 = "`xtick2' " + string(`val')
			}
			if "`1'" == ","{
				local xtick2 = "`xtick2'`1'"
			}
			mac shift
		}
		local xtick = "`xtick2'"

		local DXmin1= (min(`xtick',`DXmin'))
		local DXmax1= (max(`xtick',`DXmax'))
		*local DXwidth = `DXmax1'-`DXmin1'
		local plotwidth = `DXmax1' - `DXmin1'
		
	} // END QUI

	/*===============================================================================================*/
	/*==================================== COLUMNS   ================================================*/
	/*===============================================================================================*/
	
	qui {	// KEEP QUIET UNTIL AFTER DIAMONDS
		
		// DOUBLE LINE OPTION
		if "`double'" != "" /*& ("`lcols'" != "" | "`stats'" == "")*/{
			*gen `orig' = `id'
			replace `expand' = 1
			replace `expand' = 2 if `use' == 1
			expand `expand'
			sort `id' `se'
			bys `id' `se': gen `index' = _n
			sort  `se' `id' `index'
			egen `newid' = group(`id' `index')
			replace `id' = `newid'
			drop `newid'
			
			replace `use' = 1 if `index' == 2
			replace `es' = . if `index' == 2
			replace `lci' = . if `index' == 2
			replace `uci' = . if `index' == 2
			replace `estText' = "" if `index' == 2			
			/*
			replace `id' = `id' + 0.75 if `id' == `id'[_n-1] & `se' == `se'[_n-1] & (`use' == 1)
			replace `use' = 1 if mod(`id',1) != 0 
			replace `es' = .  if mod(`id',1) != 0
			replace `lci' = . if mod(`id',1) != 0
			replace `uci' = . if mod(`id',1) != 0
			replace `estText' = "" if mod(`id',1) != 0
			*/
			foreach var of varlist `lcols' `rcolse' `rcolsp'{
			   cap confirm string var `var'
			   if _rc == 0 {				
				tempvar length words tosplit splitwhere best
				gen `splitwhere' = 0
				gen `best' = .
				gen `length' = length(`var')
				summ `length', det
				gen `words' = wordcount(`var')
				gen `tosplit' = 1 if `length' > r(max)/2+1 & `words' >= 2
				summ `words', det
				local max = r(max)
				forvalues i = 1/`max'{
					replace `splitwhere' = strpos(`var',word(`var',`i')) ///
					 if abs( strpos(`var',word(`var',`i')) - length(`var')/2 ) < `best' ///
					 & `tosplit' == 1
					replace `best' = abs(strpos(`var',word(`var',`i')) - length(`var')/2) ///
					 if abs(strpos(`var',word(`var',`i')) - length(`var')/2) < `best' 
				}

				replace `var' = substr(`var',1,(`splitwhere'-1)) if (`tosplit' == 1) & (`index' == 1)
				replace `var' = substr(`var',`splitwhere',length(`var')) if (`tosplit' == 1) & (`index' == 2)
				replace `var' = "" if (`tosplit' != 1) & (`index' == 2) & (`use' == 1)
				drop `length' `words' `tosplit' `splitwhere' `best'
			   }
			   if _rc != 0{
				replace `var' = . if (`index' == 2) & (`use' == 1)
			   }
			}
		}
	
		local maxline = 1
		if "`lcols'" != "" {
			tokenize "`lcols'"
			local lcolsN = 0

			while "`1'" != "" {
				cap confirm var `1'
				if _rc!=0  {
					di in re "Variable `1' not defined"
					exit _rc
				}
				local lcolsN = `lcolsN' + 1
				tempvar left`lcolsN' leftLB`lcolsN' leftWD`lcolsN'
				cap confirm string var `1'
				if _rc == 0{
					gen str `leftLB`lcolsN'' = `1'
				}
				if _rc != 0{
					cap decode `1', gen(`leftLB`lcolsN'')
					if _rc != 0{
						local f: format `1'
						gen str `leftLB`lcolsN'' = string(`1', "`f'")
						replace `leftLB`lcolsN'' = "" if `leftLB`lcolsN'' == "."
					}
				}
				replace `leftLB`lcolsN'' = "" if (`se' != 1)
				if `lcolsN' > 1 {
					replace `leftLB`lcolsN'' = "" if (`use' != 1)
				}
				
				local colName: variable label `1'
				if "`colName'"==""{
					local colName = "`1'"
				}
				local titleln = strlen("`colName'")
				tempvar tmpln
				gen `tmpln' = strlen(`leftLB`lcolsN'')
				qui summ `tmpln' if `use' != 0
				local otherln = r(max)
				drop `tmpln'
				
				local spread = max(round(`titleln'/`otherln'), 1)
				if `spread' > 4 {
					local spread = 4
				}
				
				gettoken word1 remain : colName
				local index = 1
				foreach cword of local remain {		
					local lencword = strlen("`cword'")
					local lenword = strlen("`word`index''")
					if (`=`lenword' + `lencword' + 1' > `otherln') & (`index' < `spread') & (`lenword' > `lencword') {
						local ++index
						local word`index' = "`cword'"
					}
					else {
						local word`index' = "`word`index''" +  " " + "`cword'"
					}
				}
				
				forvalues line=1(1)`index'{
					replace `leftLB`lcolsN'' = "`word`line''" in `line'
				}
				if `spread' > `maxline'{
					local maxline = `spread'
				}
				mac shift
			}
		}

		local rcolseN = 0
		local rcolspN = 0
		local parameters "se sp"
		foreach param of local parameters  {
			if "`rcol`param''" != "" {
				tokenize "`rcol`param''"
				local rcol`param'N = 0
				while "`1'" != "" {
					
					cap confirm var `1'
					if _rc!=0  {
						di in re "Variable `1' not defined"
						exit _rc
					}
					local ++rcol`param'N 

					cap confirm string var `1'
					if _rc == 0{
						tempvar `param'rightLB`rcol`param'N'
						gen str ``param'rightLB`rcol`param'N'' = `1'
					}
					if _rc != 0 {
						local f: format `1'
						tempvar `param'rightLB`rcol`param'N'
						gen str ``param'rightLB`rcol`param'N'' = string(`1', "f")
						replace ``param'rightLB`rcol`param'N'' = "" if ``param'rightLB`rcol`param'N'' == "."
					}
					
					local colName: variable label `1'
					if "`colName'"==""{
						local colName = "`1'"
					}
					local titleln = strlen("`colName'")
					tempvar tmpln
					gen `tmpln' = strlen(``param'rightLB`rcol`param'N'')
					qui summ `tmpln' if `use' != 0
					local otherln = r(max)
					drop `tmpln'
					
					local spread = max(round(`titleln'/`otherln'), 1)
					if `spread' > 4 {
						local spread = 4
					}
										
					gettoken word1 remain : colName
					local index = 1
					foreach cword of local remain {		
						local lencword = strlen("`cword'")
						local lenword = strlen("`word`index''")
						if (`=`lenword' + `lencword' + 1' > `otherln') & (`index' < `spread') & (`lenword' > `lencword') {
							local ++index
							local word`index' = "`cword'"
						}
						else {
							local word`index' = "`word`index''" +  " " + "`cword'"
						}
					}
					
					forvalues line=1(1)`index'{
						replace ``param'rightLB`rcol`param'N'' = "`word`line''" in `line'
					}
					if `spread' > `maxline'{
						local maxline = `spread'
					}
					mac shift
				}
			}
		}
		

		// now get rid of extra title rows if they weren't used
		replace `id' = `id' - 5 + `maxline' if `id' > `maxline'
		if `maxline'==3 {
			drop in 4	
		}
		if `maxline'==2 {
			drop in 3/4 
		}
		if `maxline'==1 {
			drop in 2/4 
		}
		
		local borderline = `maxline' + 0.35
				
		local leftWDtot = 0
		local serightWDtot = 0
		local leftWDtotNoTi = 0
		
		forvalues i = 1/`lcolsN'{
			tempvar leftWD`i'	
			getlen `leftLB`i'' `leftWD`i''
			summ `leftWD`i'' if `use' != 3 	// DON'T INCLUDE OVERALL STATS AT THIS POINT
			local maxL = r(max)
			local leftWDtotNoTi = `leftWDtotNoTi' + `maxL'
			replace `leftWD`i'' = `maxL'
		}

		tempvar titleLN				// CHECK IF OVERALL LENGTH BIGGER THAN REST OF LCOLS
		getlen `leftLB1' `titleLN'	
		summ `titleLN' if `use' == 3
		local leftWDtot = max(`leftWDtotNoTi', r(max))

		forvalues i = 1/`rcolseN'{
			tempvar serightWD`i'
			getlen `serightLB`i'' `serightWD`i''
			summ `serightWD`i'' if  `use' != 3
			replace `serightWD`i'' = r(max)
			local serightWDtot = `serightWDtot' + r(max)
		}	

		local LEFT_WD = `leftWDtot'
		local RIGHT_WD = `serightWDtot'*2
		local textWD = ((2*`plotwidth')/(1 - `astext'/100) - (2*`plotwidth')) /(`leftWDtot' + `serightWDtot'*2)
		*local textWD = `astext'/100
		local leftWDtot = `leftWDtot'*`textWD'
		forvalues i = 1/`lcolsN'{
			tempvar left`i'
			gen `left`i'' = `DXmin' - `leftWDtot'
			*local leftWDtot = `leftWDtot' - `leftWD`i''
			local leftWDtot = `leftWDtot' - `leftWD`i''*`textWD'
		}

		tempvar seright1
		gen `seright1' = `DXmax1'
		forvalues i = 2/`=`rcolseN'+1'{
			local r2 = `i' - 1
			tempvar seright`i'
			gen `seright`i'' = `seright`r2'' + `serightWD`r2''*`textWD'
		}
		
		local shift =  `seright`=`rcolseN'+1'' + (`plotwidth' - `DXmax1')*(`DXmin1' < 0) - `DXmin1'*(`DXmin1' > 0)
	
		forvalues i = 1/`=`rcolseN'+1'{
			tempvar spright`i'
			*gen `spright`i'' = `seright`i'' + `shift'
			gen `spright`i'' = `seright`i'' + `shift' 
		}
		
		*tempvar spright`=`rcolspN'+1'
		*gen `spright`=`rcolspN'+1'' =  `spright`rcolspN'' + `serightWD`rcolseN''*`textWD'
	
		local AXmin = `left1'
		local AXmax = `spright`=`rcolspN'+1''
		
		*local DXmin2= `DXmin1' + `shift'
		*local DXmax2= `DXmax1' + `shift'
		
		local DXmin2= `DXmin1' + `shift' 
		local DXmax2= `DXmax1' + `shift' 
		
		// DIAMONDS 
		tempvar DIAMleftX DIAMrightX DIAMbottomX DIAMtopX DIAMleftY1 DIAMrightY1 DIAMleftY2 DIAMrightY2 DIAMbottomY DIAMtopY
		gen `DIAMleftX'   = `lci' if `use' == 2 | `use' == 3 
		gen `DIAMleftY1'  = `id' if (`use' == 2 | `use' == 3) 
		gen `DIAMleftY2'  = `id' if (`use' == 2 | `use' == 3) 
		
		gen `DIAMrightX'  = `uci' if (`use' == 2 | `use' == 3)
		gen `DIAMrightY1' = `id' if (`use' == 2 | `use' == 3)
		gen `DIAMrightY2' = `id' if (`use' == 2 | `use' == 3)
		
		gen `DIAMbottomY' = `id' - 0.4 if (`use' == 2 | `use' == 3)
		gen `DIAMtopY' 	  = `id' + 0.4 if (`use' == 2 | `use' == 3)
		gen `DIAMtopX'    = `es' if (`use' == 2 | `use' == 3)

		replace `DIAMleftX' = `DXmin1' if (`lci' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMleftX' = . if (`es' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		//If one study, no diamond
		replace `DIAMleftX' = . if (`df' < 2) & (`use' == 2 | `use' == 3) 
		
		replace `DIAMleftY1' = `id' + 0.4*(abs((`DXmin1' -`lci')/(`es'-`lci'))) if (`lci' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMleftY1' = . if (`es' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
	
		replace `DIAMleftY2' = `id' - 0.4*( abs((`DXmin1' -`lci')/(`es'-`lci')) ) if (`lci' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMleftY2' = . if (`es' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		
		replace `DIAMrightX' = `DXmax`r'' if (`uci' > `DXmax1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMrightX' = . if (`es' > `DXmax1' ) & (`use' == 2 | `use' == 3) 
		//If one study, no diamond
		replace `DIAMrightX' = . if (`df' == 1) & (`use' == 2 | `use' == 3) 
	
		replace `DIAMrightY1' = `id' + 0.4*( abs((`uci'-`DXmax1' )/(`uci'-`es')) ) if (`uci' > `DXmax1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMrightY1' = . if (`es' > `DXmax1' ) & (`use' == 2 | `use' == 3) 

		replace `DIAMrightY2' = `id' - 0.4*( abs((`uci'-`DXmax1' )/(`uci'-`es')) ) if (`uci' > `DXmax1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMrightY2' = . if (`es' > `DXmax1' ) & (`use' == 2 | `use' == 3) 

		replace `DIAMbottomY' = `id' - 0.4*( abs((`uci'-`DXmin1' )/(`uci'-`es')) ) if (`es' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMbottomY' = `id' - 0.4*( abs((`DXmax1' -`lci')/(`es'-`lci')) ) if (`es' > `DXmax1' ) & (`use' == 2 | `use' == 3) 

		replace `DIAMtopY' = `id' + 0.4*( abs((`uci'-`DXmin1' )/(`uci'-`es')) ) if (`es' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMtopY' = `id' + 0.4*( abs((`DXmax1' -`lci')/(`es'-`lci')) ) if (`es' > `DXmax1' ) & (`use' == 2 | `use' == 3) 

		replace `DIAMtopX' = `DXmin1'  if (`es' < `DXmin1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMtopX' = `DXmax1'  if (`es' > `DXmax1' ) & (`use' == 2 | `use' == 3) 
		replace `DIAMtopX' = . if ((`uci' < `DXmin1' ) | (`lci' > `DXmax1' )) & (`use' == 2 | `use' == 3) 
	
		gen `DIAMbottomX' = `DIAMtopX'

	} // END QUI

	forvalues i = 1/`lcolsN'{
		local lcolCommands`i' "(scatter `id' `left`i'', msymbol(none) mlabel(`leftLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
	}
	
	foreach param of local parameters  {
		forvalues i = 1/`rcol`param'N' {
			local `param'rcolCommands`i' "(scatter `id' ``param'right`i'', msymbol(none) mlabel(``param'rightLB`i'') mlabcolor(black) mlabpos(3) mlabsize(`texts'))"
		}
	}
	
	if `"`diamopt'"' == "" {
		local diamopt "lcolor("255 0 0")"
	}
	else {
		if strpos(`"`diamopt'"',"hor") != 0 | strpos(`"`diamopt'"',"vert") != 0 {
			di as error "Options horizontal/vertical not allowed in diamopt()"
			exit
		}
		if strpos(`"`diamopt'"',"con") != 0{
			di as error "Option connect() not allowed in diamopt()"
			exit
		}
		if strpos(`"`diamopt'"',"lp") != 0{
			di as error "Option lpattern() not allowed in diamopt()"
			exit
		}
		local diamopt `"`diamopt'"'
	}
	
	//Box options
	if "`box'" == "" {
		local iw = "[aw = _WT]"
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"msy") == 0{
			local boxopts = `"`boxopts' msymbol(square)"' 
		}
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"msi") == 0{
			local boxopts = `"`boxopts' msize(0.5)"' 
		}
		if `"`boxopts'"' != "" & strpos(`"`boxopts'"',"mco") == 0{
			local boxopts = `"`boxopts' mcolor("180 180 180")"' 
		}
		if `"`boxopts'"' == "" {
			local boxopts "msymbol(square) msize(.5) mcolor("180 180 180")"
		}
		else{
			local boxopts `"`boxopts'"'
		}
	}
	if ("`box'" != "") {
		local boxopts "msymbol(none)"
		local iw
	}
	
	//Point options
	if "`smooth'" != "" {
		local pointsymbol "msymbol(Oh)"
		local pointcolor "mcolor("128 128 128")"
		local pointsize "msize(small)"
	}
	else {
		local pointsymbol "msymbol(O)"
		local pointcolor "mcolor("0 0 0")"
		local pointsize "msize(vsmall)"
	}
	
	if `"`pointopt'"' != "" & strpos(`"`pointopt'"',"msy") == 0{
		local pointopt = `"`pointopt' `pointsymbol'"' 
	}
	if `"`pointopt'"' != "" & strpos(`"`pointopt'"',"msi") == 0{
		local pointopt = `"`pointopt' `pointsize'"' 
	}
	if `"`pointopt'"' != "" & strpos(`"`pointopt'"',"mc") == 0{
		local pointopt = `"`pointopt' `pointcolor'"' 
	}
	if `"`pointopt'"' == ""{
		local pointopt "`pointsymbol' `pointsize' `pointcolor'"
	}
	else{
		local pointopt `"`pointopt'"'
	}
	
	//Smooth Point options
	if `"`smoothpointopts'"' != "" & strpos(`"`smoothpointopts'"',"msy") == 0 {
		local smoothpointopts = `"`smoothpointopts' msymbol(D)"' 
	}
	if `"`smoothpointopts'"' != "" & strpos(`"`smoothpointopts'"',"ms") == 0 {
		local smoothpointopts = `"`smoothpointopts' msize(vsmall)"' 
	}
	if `"`smoothpointopts'"' != "" & strpos(`"`smoothpointopts'"',"mc") == 0 {
		local smoothpointopts = `"`smoothpointopts' mcolor("0 0 0")"' 
	}
	if `"`smoothpointopts'"' == ""{
		local smoothpointopts "msymbol(D) msize(vsmall) mcolor("0 0 0")"
		if "`compabs'" != "" {
			local smoothpointopts0 "msymbol(D) msize(vsmall) mcolor("0 0 0")"
			local smoothpointopts1 "msymbol(D) msize(vsmall) mcolor("255 127 0")"
		}
	}
	else{
		local smoothpointopts `"`smoothpointopts'"'
	}
	
	// CI options
	if `"`ciopt'"' == "" {
		if "`smooth'" != "" {
			local ciopt "lwidth(1.25) lcolor(gs13)"
		}
		else {
			local ciopt "lcolor("0 0 0")"
		}
	}
	else {
		if strpos(`"`ciopt'"',"hor") != 0 | strpos(`"`ciopt'"',"ver") != 0{
			di as error "Options horizontal/vertical not allowed in ciopt()"
			exit
		}
		if strpos(`"`ciopt'"',"con") != 0{
			di as error "Option connect() not allowed in ciopt()"
			exit
		}
		if strpos(`"`ciopt'"',"lp") != 0{
			di as error "Option lpattern() not allowed in ciopt()"
			exit
		}
		if `"`ciopt'"' != "" & strpos(`"`ciopt'"',"lc") == 0{
			local ciopt = `"`ciopt' lcolor("0 0 0")"' 
		}
		if "`smooth'" != "" {
			if strpos(`"`ciopt'"',"lc") == 0 {
				local ciopts = `"`ciopt' lcolor(red)"' 	
			}
			if strpos(`"`ciopt'"',"lw") == 0 {
				local ciopts = `"`ciopt' lwidth(2)"' 
			}
		}
		local ciopt `"`ciopt'"'
	}
	
	//Smooth ci
	if `"`smoothciopts'"' == "" {
		local smoothciopts "lcolor("0 0 0")"
	}
	else {
		if strpos(`"`smoothciopts'"',"hor") != 0 | strpos(`"`smoothciopts'"',"vert") != 0{
			di as error "Options horizontal/vertical not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"con") != 0{
			di as error "Option connect() not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"lp") != 0 {
			di as error "Option lpattern() not allowed in ciopts()"
			exit
		}
		if strpos(`"`smoothciopts'"',"lw") == 0 {
				local smoothciopts = `"`smoothciopts' lwidth(.5)"' 
		}

		local smoothciopts `"`smoothciopts'"'
	}
	
	// Arrow options
	if `"`arrowopt'"' == "" {
		local arrowopt "mcolor("0 0 0") lstyle(none)"
	}
	else {
		local forbidden "connect horizontal vertical lpattern lwidth lcolor lsytle"
		foreach option of local forbidden {
			if strpos(`"`arrowopt'"',"`option'")  != 0 {
				di as error "Option `option'() not allowed in arrowopt()"
				exit
			}
		}
		if `"`arrowopt'"' != "" & strpos(`"`arrowopt'"',"mc") == 0{
			local arrowopt = `"`arrowopt' mcolor("0 0 0")"' 
		}
		local arrowopt `"`arrowopt' lstyle(none)"'
	}

	// END GRAPH OPTS

	qui {
		//Generate indicator on direction of the off-scale arro
		tempvar rightarrow leftarrow biarrow noarrow rightlimit leftlimit offRhiY offRhiX offRloY offRloX offLloY offLloX offLhiY offLhiX
		gen `rightarrow' = 0
		gen `leftarrow' = 0
		gen `biarrow' = 0
		gen `noarrow' = 0
		
		replace `rightarrow' = 1 if ///
			(((round(`uci', 0.001) > round(`DXmax1' , 0.001)) & (round(`lci', 0.001) >= round(`DXmin1' , 0.001)))  |  ///	
			((round(`uci', 0.001) > round(`DXmax1' , 0.001)) & (round(`lci', 0.001) > round(`DXmax1' , 0.001))))  &  ///
			(`use' == 1) & (`uci' != .) & (`lci' != .)
			
			
		replace `leftarrow' = 1 if ///
			(((round(`lci', 0.001) < round(`DXmin1' , 0.001)) & (round(`uci', 0.001) <= round(`DXmax1' , 0.001))) | ///
			((round(`lci', 0.001) < round(`DXmin1' , 0.001)) & (round(`uci', 0.001) < round(`DXmin1' , 0.001)))) & ///
			(`use' == 1) & (`uci' != .) & (`lci' != .)
		
		replace `biarrow' = 1 if ///
			(round(`lci', 0.001) < round(`DXmin1' , 0.001)) & ///
			(round(`uci', 0.001) > round(`DXmax1' , 0.001)) & ///
			(`use' == 1) & (`uci' != .) & (`lci' != .)
			
		replace `noarrow' = 1 if ///
			(`leftarrow' != 1) & (`rightarrow' != 1) & (`biarrow' != 1) & ///
			(`use' == 1) & (`uci' != .) & (`lci' != .)	

		replace `lci' = `DXmin1'  if (round(`lci', 0.001) < round(`DXmin1' , 0.001)) & (`lci' !=.) & (`use' == 1) 
		replace `uci' = `DXmax1'  if (round(`uci', 0.001) > round(`DXmax1' , 0.001)) & (`uci' !=.) & (`use' == 1) 
		
		replace `lci' = `DXmax1' - 0.00001  if (round(`lci', 0.001) > round(`DXmax1' , 0.001)) & (`lci' !=. ) & (`use' == 1) 
		replace `uci' = `DXmin1' + 0.00001  if (round(`uci', 0.001) < round(`DXmin1' , 0.001)) & (`uci' !=. ) & (`use' == 1) 
		
		*replace `lci' = . if (round(`lci', 0.001) > round(`DXmax1' , 0.001)) & (`lci' !=. ) & (`use' == 1) 
		*replace `uci' = . if (round(`uci', 0.001) < round(`DXmin1' , 0.001)) & (`uci' !=. ) & (`use' == 1) 

		replace `es' = . if (round(`es', 0.001) < round(`DXmin1' , 0.001)) & (`use' == 1) 
		replace `es' = . if (round(`es', 0.001) > round(`DXmax1' , 0.001)) & (`use' == 1)
		
		if "`smooth'" != "" {		
			replace `modellci' = `DXmin'  if (round(`modellci', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1) 
			replace `modeluci' = `DXmax'  if (round(`modeluci', 0.001) > round(`DXmax' , 0.001)) & (`modeluci' !=.) & (`use' == 1 ) 
			
			replace `modellci' = . if (round(`modeluci', 0.001) < round(`DXmin' , 0.001)) & (`modeluci' !=. ) & (`use' == 1 ) 
			replace `modeluci' = . if (round(`modellci', 0.001) > round(`DXmax' , 0.001)) & (`modellci' !=. ) & (`use' == 1 )
			replace `modeles' = . if (round(`modeles', 0.001) < round(`DXmin' , 0.001)) & (`use' == 1 ) 
			replace `modeles' = . if (round(`modeles', 0.001) > round(`DXmax' , 0.001)) & (`use' == 1 )
		}

		summ `id'
		local xaxislineposition = r(max)

		local xaxis1 "(pci `xaxislineposition' `DXmin1' `xaxislineposition' `DXmax1', lwidth(thin) lcolor(black))"
		local xaxis2 "(pci `xaxislineposition' `DXmin2' `xaxislineposition' `DXmax2', lwidth(thin) lcolor(black))"
		
		/*Xaxis 1 title */
		local xaxistitlex1 `=(`DXmax1' + `DXmin1')*0.5'
		local xaxistitlex2 `=(`DXmax1' + `DXmin1')*0.5 + `shift''
		
		if "`comparative'" != "" & "`outplot'" != "abs"  {
			local varx :word 1 of `varxlabs'
			local indexlab :word 2 of `varxlabs'
			local baselab :word 3 of `varxlabs'
			
			local xaxistitle1  (scatteri `=`xaxislineposition' + 2.25' `xaxistitlex1' "`plotstatse' (`varx' = `indexlab' / `baselab')", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))
			local xaxistitle2  (scatteri `=`xaxislineposition' + 2.25' `xaxistitlex2' "`plotstatsp' (`varx' = `indexlab' / `baselab')", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))
		} 
		else {
			local xaxistitle1  (scatteri `=`xaxislineposition' + 2.25' `xaxistitlex1' "`plotstatse'", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))
			local xaxistitle2  (scatteri `=`xaxislineposition' + 2.25' `xaxistitlex2' "`plotstatsp'", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))
		}
		/*xticks*/
		local ticksx1
		local ticksx2
		tokenize "`xtick'", parse(",")	
		while "`1'" != "" {
			if "`1'" != "," {
				forvalues r=1(1)2 {
					local where = `1'
					if `r' == 2 {          
						local where = `1' + `DXmin2' - `DXmin1'
					}
					local ticksx`r' "`ticksx`r'' (pci `xaxislineposition'  `where' 	`=`xaxislineposition'+.25' 	`where' , lwidth(thin) lcolor(black)) "
				}
			}
			macro shift 
		}
		/*labels*/
		local xaxislabels
		tokenize `lblcmd'
		while "`1'" != ""{
			forvalues r = 1(1)2 {
				local where = `1'
				if `r' == 2 {
					*local where = `1' + `DXmax1' - `DXmin1' + 0.5*`rightWDtot'*`textWD' 
					local where = `1' + `DXmin2' - `DXmin1'
				}
				local xaxislabels`r' "`xaxislabels`r'' (scatteri `=`xaxislineposition'+1' `where' "`2'", msymbol(i) mlabcolor(black) mlabpos(0) mlabsize(`texts'))"
			}
			macro shift 2
		}
		if "`grid'" != "" {
			tempvar gridy gridxmax gridxmin
			
			gen `gridy' = `id' + 0.5
			gen `gridxmax' = `AXmax'
			gen `gridxmin' = `left1'
			local betweengrids "(pcspike `gridy' `gridxmin' `gridy' `gridxmax'  if `use' == 1 , lwidth(vvthin) lcolor(gs12))"	
		}
		
		//Shift the position of sensitivitiy plot
		replace `lci' = `lci' + `shift' if !`se' 
		replace `uci' = `uci' + `shift'  if !`se' 
		replace `es' = `es' + `shift'  if !`se' 
		
		if "`smooth'" != "" { 
			replace `modellci' = `modellci' + `shift'  if !`se' 
			replace `modeluci' = `modeluci' + `shift'  if !`se' 
			replace `modeles' = `modeles' + `shift'  if !`se' 
		}
		
		replace `DIAMbottomX' = `DIAMbottomX' + `shift'  if !`se' 
		replace `DIAMtopX' = `DIAMtopX' + `shift'  if !`se' 
		replace `DIAMleftX' = `DIAMleftX' + `shift' if !`se' 
		replace `DIAMrightX' = `DIAMrightX' + `shift'  if !`se' 			
	}	// end qui	
	
	tempvar tempOv overrallLine ovMin ovMax h0Line
	
	if `"`olineopt'"' == "" {
		local olineopt "lwidth(thin) lcolor(red) lpattern(shortdash)"
	}
	if `"`vlineopt'"' == "" {
		local vlineopt "lwidth(thin) lcolor(black) lpattern(solid)"
	}
	qui summ `id'
	local DYmin = r(min)
	local DYmax = r(max)+2
	
	forvalues r= 1(1)2 {
		qui summ `es' if `use' == 3 & `se' == `=2 - `r''
		local overall`r' = r(max)				
		if `overall1' > `DXmax1' | `overall1' < `DXmin1' | "`ovline'" != "" {	// ditch if not on graph
			local overallCommand`r' ""
		}
		else {
			local overallCommand`r' `" (pci `=`DYmax'-2' `overall`r'' `borderline' `overall`r'', `olineopt') "'
		
		}
		if "`ovline'" != "" {
			local overallCommand`r' ""
		}
		if "`subline'" != "" & "`groupvar'" != "" {
			local sublineCommand`r' ""
			
			qui label list `groupvar'
			local nlevels = r(max)
			forvalues l = 1/`nlevels' {
				qui summ `es' if `use' == 2  & `groupvar' == `l' & (`se' == `=`r' - 1')
				local tempSub`l' = r(mean)
				qui summ `id' if `use' == 1 & `groupvar' == `l'
				local subMax`l' = r(max) + 1
				local subMin`l' = r(min) - 2
				qui count if `use' == 1 & `groupvar' == `l' & (`se' == `=`r' - 1')
				if r(N) > 1 {
					local sublineCommand`r' `" `sublineCommand`r'' (pci `subMin`l'' `tempSub`l'' `subMax`l'' `tempSub`l'', `olineopt')"'
				}
			}
		}
		else {
			local sublineCommand`r' ""
		}
	}
	if `"`xline'"' != `""' {
			tokenize "`xline'", parse(",")
			if "`logscale'" != "" {
				if "`1'" == "0" {
					local xlineval = ln(`=10^(-`dp')')
				}
				else {
					local xlineval = ln(`1')
				}
			}
			else {
				local xlineval = `1'
			}
			if "`3'" != "" {
				local xlineopts = "`3'"
			}
			else {
				local xlineopts = "lcolor(black)"
			}
			local xlineval2 =  `xlineval' + `shift'
			local vlineCommand1 `" (pci `=`DYmax'-2' `xlineval' `borderline' `xlineval', `xlineopts') "'
			local vlineCommand2 `" (pci `=`DYmax'-2' `xlineval2' `borderline' `xlineval2', `xlineopts') "'
			
			if (`xlineval' > `DXmax1' | `xlineval' < `DXmin1') & "`xline'" != "" {	// ditch if not on graph
				local vlineCommand1 ""
				local vlineCommand2 ""
			}
	}
	
	
	/*===============================================================================================*/
	/*====================================  GRAPH    ================================================*/
	/*===============================================================================================*/
	if "`smooth'" != "" {
		local xboxcenter "`modeles'"
		local smoothcommands1 "(pcspike `id' `modellci' `id' `modeluci' if `use' == 1 , `smoothciopts')"
		local smoothcommands2 "(scatter `id' `modeles' if `use' == 1 , `smoothpointopts')"
	}
	else {
		local xboxcenter "`es'"
	}
	
	//Observed CI
	local cicommand "(pcspike `id' `lci' `id' `uci' if `use' == 1 , `ciopt')"
	
	//Diamond
	local diamondcommand1 "(pcspike `DIAMleftY1' `DIAMleftX' `DIAMtopY' `DIAMtopX' if (`use' == 2 | `use' == 3) , `diamopt')"
	local diamondcommand2 "(pcspike `DIAMtopY' `DIAMtopX' `DIAMrightY1' `DIAMrightX' if (`use' == 2 | `use' == 3) , `diamopt')"
	local diamondcommand3 "(pcspike `DIAMrightY2' `DIAMrightX' `DIAMbottomY' `DIAMbottomX' if (`use' == 2 | `use' == 3) , `diamopt')"
	local diamondcommand4 "(pcspike `DIAMbottomY' `DIAMbottomX' `DIAMleftY2' `DIAMleftX' if (`use' == 2 | `use' == 3) , `diamopt')"
	
	//Give name if none
	if strpos(`"`foptions'"',"name(") == 0 {
		local fname = "name(fplot, replace)"
	}
	if "$by_index_" != "" {
		local fname = "name(fplot" + "$by_index_" + ", replace)"
		noi di as res _n  "NOTE: forest plot name -> fplot$by_index_"
	}
	#delimit ;
	twoway
		`notecmd' 
		`overallCommand1' `sublineCommand1' `overallCommand2' `sublineCommand2' 
		`hetGroupCmd'  
		`xaxis1' `xaxistitle1' `ticksx1' `xaxislabels1'  
		`xaxis2' `xaxistitle2'  `ticksx2' `xaxislabels2' 
		`vlineCommand1' `vlineCommand2' 
	 /*COLUMN VARIABLES */
		`lcolCommands1' `lcolCommands2' `lcolCommands3' `lcolCommands4' `lcolCommands5' `lcolCommands6'
		`lcolCommands7' `lcolCommands8' `lcolCommands9' `lcolCommands10' `lcolCommands11' `lcolCommands12' 
		`sercolCommands1' `sercolCommands2' `sercolCommands3' `sprcolCommands1' `sprcolCommands2' `sprcolCommands3'
	 /*PLOT EMPTY POINTS AND PUT ALL THE GRAPH OPTIONS IN THERE */ 
		(scatter `id' `xboxcenter' `iw' if `use' == 1 , 
			`boxopts'			
			yscale(range(`DYmin' `DYmax') noline reverse)
			ylabel(none) ytitle("")
			xscale(range(`AXmin' `AXmax') noline)
			xlabel(none)
			yline(`borderline', lwidth(thin) lcolor(gs12))
			xtitle("") legend(off) xtick("")) 
	 /*HERE ARE GRIDS */
		`betweengrids'			
	 /*HERE ARE THE CONFIDENCE INTERVALS */
		`cicommand' `smoothcommands1'
	 /*ADD ARROWS  `ICICmd1' `ICICmd2' `ICICmd3'*/
		(pcarrow `id' `uci' `id' `lci' if `leftarrow' == 1 , `arrowopt')	
		(pcarrow `id' `lci' `id' `uci' if `rightarrow' == 1 , `arrowopt')	
		(pcbarrow `id' `lci' `id' `uci' if `biarrow' == 1 , `arrowopt')	
	 /*DIAMONDS FOR SUMMARY ESTIMATES -START FROM 9 O'CLOCK */
		
		`diamondcommand1'
		`diamondcommand2'
		`diamondcommand3'
		`diamondcommand4'
	 /*LAST OF ALL PLOT EFFECT MARKERS TO CLARIFY  */
		(scatter `id'  `es'  if `use' == 1 , `pointopt') 
		`smoothcommands2' 
		`overallCommand1' `overallCommand2'	
		,`foptions' `fname'
		;
		#delimit cr		
			

		if `"`graphsave'"' != `""' {
			di _n
			noi graph save `graphsave', replace
		}				
end

/*==================================== GETWIDTH  ================================================*/
/*===============================================================================================*/
capture program drop getlen
program define getlen
version 14.1
//From metaprop

qui{

	gen `2' = 0
	count
	local N = r(N)
	forvalues i = 1/`N'{
		local this = `1'[`i']
		local width: _length "`this'"
		replace `2' =  `width' in `i'
	}
} 

end
/*+++++++++++++++++++	SUPPORTING FUNCTIONS: SROC ++++++++++++++++++++++++++++++++++++
				   DRAW THE SROC CURVES, CROSSES, CONFIDENCE & PREDICTION REGION
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop sroc
	program define sroc
		version 14.0

		#delimit ;
		syntax varlist,	
			selogodds(name) /*Log odds for se*/
			splogodds(name) /*Log odds for sp*/
			popabsoutse(name) /*pop se*/
			popabsoutsp(name) /*pop sp*/
			v(name) /*Var-cov for log odds se & se*/
			bvar(name) /*Between study var-cov*/
			model(string) /*model*/ 
			[
			cveffect(string)
			groupvar(name) /*Grouping variable*/
			p(integer 0) /*No of parameters in the regression equation*/
			cimethod(string) /*How to compute the study-specific CI*/
			LEVel(integer 95) /*Significance level*/
			SUMMaryonly
			COLorpalette(string) /*soptions: Colors seperated by space*/
			noPOPest  /*Suppress population statistic*/
			noPREDiction  /*soptions:no Prediction region*/
			noCREG  /*no confidence region*/
			noCURve  /*soptions:no curve*/
			Bubbles /*soptions: size of study by bubbles*/
			BUBbleid /*soptions: Identify the bubbles by index*/
			SPointopt(string) /*soptions: options study points*/
			OPointopt(string) /*soptions: options the overall summary points*/
			CUrveopt(string) /*soptions: options the CI points*/
			CIopt(string) /*soptions: options the CI points*/
			PREDCIopt(string) /*soptions: options the PREDCI points*/
			BUBOpt(string) /*soptions: options the bubble points*/
			BIDopt(string) /*soptions: options the bubble ID points*/
			GRAPHSave(string asis)
			STRAtify  /*Stratify*/
			* /*soptions:Other two-way options*/
			]
			;
		#delimit cr
		tempvar se selci seuci sp splci spuci csp Ni rowid Dis tp fp NDis fn tn mu gvar
		tempname vi bvari
		
		local soptions `"`options'"'
		
		tokenize `varlist'
		gen `tp' = `1'
		gen `fp' = `2'
		gen `fn' = `3'
		gen `tn' = `4'
		gen `Dis' = (`tp' + `fn')
		gen `NDis' = (`fp' + `tn')		
		gen `Ni' = `Dis' +`NDis'
		gen `rowid' = _n
		
		if "`groupvar'" != "" {
			my_ncod `gvar', oldvar(`groupvar') //
		}
		//CI
		metadta_propci `Dis' `tp', p(`se') lowerci(`selci') upperci(`seuci') cimethod(`cimethod') level(`level')
		metadta_propci `NDis' `tn', p(`sp') lowerci(`splci') upperci(`spuci') cimethod(`cimethod') level(`level')
		
		gen `csp' = 1 - `sp'	
				
		/*If categorical variable, obtain the sroc and the drawing parameters for each level*/
		if "`groupvar'" != "" {
			qui label list `gvar'
			local nlevels = r(max)
		}
		else {
			gen `gvar' = 1
			local nlevels = 1
		}
		if "`colorpalette'" == "" {
			local colorpalette  "black forest_green cranberry blue sienna orange emerald magenta dknavy gray purple"
		}
		else {
			local kcolors : word count `colorpalette'
			if (`kcolors' < `nlevels') {
				di as error "Please specify the colours to be used for all the `nlevels' levels of `gvar'" 
				di as error "colours should be separated by space in the colorpalette() option"
				exit
			}
		}
		local index 0
		local centre
		local kross
		local sroc
		local points
		local rings
		local idbubble
		local cregion
		local pregion
		local legendlabel
		local legendorder
		
		local shiftindex 1
		
		/*Options*/
		// CI options
		if `"`ciopt'"' == "" {
			if "`model'" == "random" {
				local ciopt "lpattern(dash)"
			}
			else {
				local ciopt "lpattern(solid)"
			}
		}
		else {
			local forbidden "lcolor"
			foreach option of local forbidden {
				if strpos(`"`ciopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in ciopt()"
					exit
				}
			}
		}
		
		// PREDCI options
		if `"`predciopt'"' == "" {
			local predciopt "lpattern(-.)"
		}
		else {
			local forbidden "lcolor"
			foreach option of local forbidden {
				if strpos(`"`predciopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in predciopt()"
					exit
				}
			}
			if `"`predciopt'"' != "" & strpos(`"`predciopt'"',"lpattern") == 0{
				local predciopt = `"`predciopt' lpattern(-.)"' 
			}
		}
		
		// Overall Point options
		if `"`opointopt'"' == "" {
			local opointopt "msymbol(D)"
			local popopointopt "msymbol(Dh)"
		}
		else {
			local forbidden "mcolor"
			foreach option of local forbidden {
				if strpos(`"`opointopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in opointopt()"
					exit
				}
			}
			if `"`opointopt'"' != "" & strpos(`"`opointopt'"',"msymbol") == 0{
				local opointopt = `"`opointopt' msymbol(D)"' 
			}
		}
		
		// Study point options
		if `"`spointopt'"' == "" {
			local spointopt "msymbol(o)"
		}
		else {
			local forbidden "mcolor"
			foreach option of local forbidden {
				if strpos(`"`spointopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in spointopt()"
					exit
				}
			}
			if `"`spointopt'"' != "" & strpos(`"`spointopt'"',"msymbol") == 0{
				local spointopt = `"`spointopt' msymbol(o)"' 
			}
		}
		
		// Curve options
		if `"`curveopt'"' != "" {
			local forbidden "lcolor"
			foreach option of local forbidden {
				if strpos(`"`curveopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in curveopt()"
					exit
				}
			}
		}
		
		// Bubble options
		if `"`bubopt'"' == "" {
			local bubopt "msymbol(Oh)"
		}
		else {
			local forbidden "mcolor"
			foreach option of local forbidden {
				if strpos(`"`bubopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in bubopt()"
					exit
				}
			}
			if `"`bubopt'"' != "" & strpos(`"`bubopt'"',"msymbol") == 0{
				local bubopt = `"`bubopt' msymbol(Oh)"' 
			}
		}
		
		//Bubble ID options
		if `"`bidopt'"' == "" {
			local bidopt "mlabsize(`texts') msymbol(i) mlabel(`rowid')"
		}
		else {
			local forbidden "mcolor mlabcolor "
			foreach option of local forbidden {
				if strpos(`"`bidopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in bidopt()"
					exit
				}
			}
			if `"`bidopt'"' != "" & strpos(`"`bidopt'"',"mlabsize") == 0{
				local bidopt = `"`bidopt' mlabsize(`texts')"' 
			}
			if `"`bidopt'"' != "" & strpos(`"`bidopt'"',"msymbol") == 0{
				local bidopt = `"`bidopt' msymbol(i)"' 
			}
			if `"`bidopt'"' != "" & strpos(`"`bidopt'"',"mlabel") == 0{
				local bidopt = `"`bidopt' mlabel(`rowid')"' 
			}
		}
	
		qui {
			local already 0
			local nrowsp = rowsof(`splogodds')
			local nrowse = rowsof(`selogodds')
			local ovindex = rowsof(`v')
			local nrowsv = rowsof(`v')
			local ncolsv = colsof(`v')
			local layer = `=`nrowsv'/(`nlevels'*2)'
	
			mat `bvari' = (`bvar'[1 ,1], `bvar'[1, 7] \ `bvar'[1, 7], `bvar'[1, 3])
			
			forvalues j=1/`nlevels' {
				if `p' == 0 & `nlevels' == 1 {
					local shiftindex 0
				}
				*Take out the right matrices
				local seindex = `=`j'+`shiftindex''
				local spindex = `=`j'+`shiftindex''

				if (`p' == 1){					
					mat `vi' = (`v'[`=`nlevels' + `j'', `=`nlevels' + `j''],  `v'[`j', `=`nlevels' + `j''] \ `v'[`j', `=`nlevels' + `j''], `v'[`j', `j'])	
				}
				if (`p' != 1) {					
					mat `vi' = (`v'[`nrowsv', `ncolsv'],  `v'[`=`nrowsv'-1', `ncolsv'] \ `v'[`=`nrowsv'-1', `ncolsv'], `v'[`=`nrowsv'-1', `=`ncolsv'-1'])
				}							
				if "`stratify'" != "" {					
					mat `bvari' = (`bvar'[`j' ,1], `bvar'[`j', 7] \ `bvar'[`j', 7], `bvar'[`j', 3])
					mat `vi' = (`v'[`=`ncolsv'*`j'' , `ncolsv'] , `v'[`=`ncolsv'*`j'-1', `ncolsv'] \ `v'[`=`ncolsv'*`j'-1', `ncolsv'], `v'[`=`ncolsv'*`j'-1', `=`ncolsv'-1'])
					
					if (`p' == 0){
						local layerse = 1
						local layersp = 1
					}
					else {
						local layerse = (`=`nrowse'-1')/`nlevels'
						local layersp = (`=`nrowsp'-1')/`nlevels'
					}
					local seindex = `=`j'*`layerse' + `shiftindex''
					local spindex = `=`j'*`layersp' + `shiftindex''
				}
				
				if "`cveffect'" == "sp" | `p' > 1 {
					local seindex  = `nrowse'
				}
				if "`cveffect'" == "se" | `p' > 1{
					local spindex = `nrowsp'
				}
				
				local color:word `j' of `colorpalette'
				qui count if `gvar' == `j'
				//Centre
				if r(N) > 1 {	
					local mux`j' = 1 - invlogit(`splogodds'[`spindex', 1])
					local muy`j' = invlogit(`selogodds'[`seindex', 1])
					
					local popmux`j' = 1 - `popabsoutsp'[`spindex', 1]
					local popmuy`j' = `popabsoutse'[`seindex', 1]
					}
				else {						
					qui summ `sp' if `gvar' == `j'
					local mux`j' = 1 - r(mean)
					
					qui summ `se' if `gvar' == `j'
					local muy`j' = r(mean)
				}
				local centre `"`centre' (scatteri `muy`j'' `mux`j'', mcolor(`color') `opointopt')"'
				if `nlevels' == 1 {
					local ++index
					local legendlabel `"lab(`index' "Summary point") `legendlabel'"'
					local legendorder `"`index'  `legendorder'"'				
				}
				if "`popest'" == "" {
					local popcentre `"`popcentre' (scatteri `popmuy`j'' `popmux`j'', mcolor(`color') `popopointopt')"'
					
					if `nlevels' == 1 {
						local ++index
						local legendlabel `"lab(`index' "Population Mean") `legendlabel'"'
						local legendorder `"`index'  `legendorder'"'				
					}
				}
				//Crosses
				if "`model'" == "fixed" | (r(N) < 3 & "`model'" == "random" ){ 
					if r(N) > 1 {
						local leftX`j' = 1 - invlogit(`splogodds'[`spindex', 6])
						local leftY`j' = invlogit(`selogodds'[`seindex', 1])
						
						local rightX`j' = 1 - invlogit(`splogodds'[`spindex', 5])
						local rightY`j' = invlogit(`selogodds'[`seindex', 1])
						
						local topX`j' = 1 - invlogit(`splogodds'[`spindex', 1])
						local topY`j' = invlogit(`selogodds'[`seindex', 6])
						
						local bottomX`j' = 1 - invlogit(`splogodds'[`spindex', 1])
						local bottomY`j' = invlogit(`selogodds'[`seindex', 5])
					}
					else {						
						qui summ `splci' if `gvar' == `j'
						local leftX`j' = 1 - r(mean)
						local leftY`j' = `muy`j''
						
						qui summ `spuci' if `gvar' == `j'
						local rightX`j' = 1 - r(mean)
						local rightY`j' = `muy`j''
						
						local topX`j' = `mux`j''
						qui summ `selci' if `gvar' == `j'
						local topY`j' =  r(mean)
						
						local bottomX`j' = `mux`j''
						qui summ `seuci' if `gvar' == `j'
						local bottomY`j' =  r(mean)
					}
					local kross `"`kross' (pci `leftY`j'' `leftX`j'' `rightY`j'' `rightX`j'', lcolor(`color') `ciopt') (pci `topY`j'' `topX`j'' `bottomY`j'' `bottomX`j'', lcolor(`color') `ciopt') "'
					if `nlevels' == 1 {
						local ++index
						local legendlabel `"lab(`index' "Confidence intervals") `legendlabel'"'
						local legendorder `"`index'  `legendorder'"'
						local ++index						
					}
				}
				else {
				//Confidence & prediction ellipses
					
					if !`already' {
						qui set obs 500
						local already 1
					}
					qui summ `csp' if `gvar' == `j' 
					local max`j' = min(0.9999, r(min))
					local min`j' = max(0.0001, r(max))
					local N`j' = r(N)
					
					tempvar fpr`j' sp`j' se`j' xcregion`j' ycregion`j' xpregion`j' ypregion`j'
					
					range `fpr`j'' `min`j'' `max`j''
					gen `sp`j'' = 1 - `fpr`j''
					
					/*HsROC parameters*/
					local b = (sqrt(`bvari'[2,2])/sqrt(`bvari'[1,1]))^0.5
					local beta = ln(sqrt(`bvari'[2,2]) / sqrt(`bvari'[1,1]))
					
					local lambda = `b' * `selogodds'[`seindex', 1] + `splogodds'[`spindex', 1] / `b'
					local theta = 0.5 * (`b' * `selogodds'[`seindex', 1] -  `splogodds'[`spindex', 1]  /`b')
										
					local var_accu =  2*( sqrt(`bvari'[2,2]*`bvari'[1,1]) + `bvari'[2,1]) 
					local var_thresh = 0.5*( sqrt(`bvari'[2,2]*`bvari'[1,1]) - `bvari'[2,1]) 

					/*The curves*/
					if "`curve'" =="" {
						gen `se`j'' = invlogit(`lambda' * exp(-`beta' / 2) + exp(-`beta') * logit(`fpr`j''))
						local sroc "`sroc' (line `se`j'' `fpr`j'', lcolor(`color') `curveopt')"
						if `nlevels' == 1 {
							local ++index
							local legendlabel `"lab(`index' "SROC curve") `legendlabel'"'
							local legendorder `"`index'  `legendorder'"'					
						}
					}
						
					/*Joint confidence region*/
					if "`creg'" == "" {
						local t = sqrt(2*invF(2, `=`N`j'' - 2', `level'/100))
						local nlen = 500
						
						local rho = `vi'[1, 2]/sqrt(`vi'[1, 1]*`vi'[2, 2])
						
						tempvar a 
						range `a' 0 2*_pi `nlen'
						
						gen `xcregion`j'' = 1 - invlogit(`splogodds'[`spindex', 1] + sqrt(`vi'[2, 2]) * `t' * cos(`a' + acos(`rho')))
						gen `ycregion`j'' = invlogit(`selogodds'[`seindex', 1] +  sqrt(`vi'[1, 1]) * `t' * cos(`a'))
						
						local cregion `"`cregion' (line `ycregion`j'' `xcregion`j'', lcolor(`color') `ciopt')"'
						if `nlevels' == 1 {
							local ++index
							local legendlabel `"lab(`index' "Confidence region") `legendlabel'"'
							local legendorder `"`index'  `legendorder'"'					
						}
					}
					
					/*Joint prediction region*/
					if "`prediction'" == "" {
						local rho =  (`vi'[1, 2] + `bvari'[1,2])/ sqrt((`vi'[1, 1] + `bvari'[1, 1]) * (`vi'[2, 2] + `bvari'[2, 2]))
						
						gen `xpregion`j'' = 1 - invlogit(`splogodds'[`spindex', 1] + `t' * sqrt(`v'[2, 2] + `bvari'[2, 2]) * cos(`a' + acos(`rho')))
						gen `ypregion`j'' = invlogit(`selogodds'[`seindex', 1] +  sqrt(`vi'[1, 1] + `bvari'[1, 1]) * `t' * cos(`a'))
						
						local pregion `"`pregion' (line `ypregion`j'' `xpregion`j'', `predciopt' lcolor(`color'))"'
						if `nlevels' == 1 {
							local ++index
							local legendlabel `"lab(`index' "Prediction region") `legendlabel'"'
							local legendorder `"`index'  `legendorder'"'					
						}
					}
				}
				if "`summaryonly'" =="" {
					if "`bubbles'" != "" {
					//bubbles
						local rings `"`rings' (scatter `se' `csp' [fweight = `Ni'] if `gvar' == `j',   mcolor(`color') `bubopt')"'
						
						if "`bubbleid'" != "" {
							local idbubble `"`idbubble' (scatter `se' `csp' if `gvar' == `j',  mcolor(`color') mlabcolor(`color') `bidopt')"'
						}
					}
					else {
					//points
						local points `"`points' (scatter `se' `csp' if `gvar' == `j',  mcolor(`color') `spointopt')"'
					}
					if `nlevels' == 1 {
						local ++index
						local legendlabel `"lab(`index' "Observed data") `legendlabel'"'
						local legendorder `"`index'  `legendorder'"'					
					}
				}
				if `nlevels' > 1 {
					local lab:label `gvar' `j' /*label*/
					local legendlabel `"lab(`j' "`groupvar' = `lab'") `legendlabel'"'
					local legendorder `"`j'  `legendorder'"'
				}
			}
		}
		
		if strpos(`"`soptions'"', "legend") == 0 {
			local legendstr `"legend(order(`legendorder') `legendlabel' cols(1) ring(0) position(6))"'
		}
		if strpos(`"`soptions'"', "xscale") == 0 {
			local soptions `"xscale(range(0 1)) `soptions'"'
		}
		if strpos(`"`soptions'"', "yscale") == 0 {
			local soptions `"yscale(range(0 1)) `soptions'"'
		}
		if strpos(`"`soptions'"', "xtitle") == 0 {
			local soptions `"xtitle("1 - Specificity") `soptions'"'
		}
		if strpos(`"`soptions'"', "ytitle") == 0 {
			local soptions `"ytitle("Sensitivity") `soptions'"'
		}
		if strpos(`"`soptions'"', "xlabel") == 0 {
			local soptions `"xlabel(0(0.2)1) `soptions'"'
		}
		if strpos(`"`soptions'"', "ylabel") == 0 {
			local soptions `"ylabel(0(0.2)1, nogrid) `soptions'"'
		}
		if strpos(`"`soptions'"', "graphregion") == 0 {
			local soptions `"graphregion(color(white)) `soptions'"'
		}
		if strpos(`"`soptions'"', "plotregion") == 0 {
			local soptions `"plotregion(margin(medium)) `soptions'"'
		}
		if strpos(`"`soptions'"', "aspectratio") == 0 {
			local soptions `"aspectratio(1) `soptions'"'
		}
		if strpos(`"`soptions'"', "xsize") == 0 {
			local soptions `"xsize(5)  `soptions'"'
		}
		if strpos(`"`soptions'"', "ysize") == 0 {
			local soptions `"ysize(5)  `soptions'"'
		}
		
		//Give name if none
		if strpos(`"`soptions'"',"name(") == 0 {
			local sname = "name(sroc, replace)"
		}
		if "$by_index_" != "" {
			local sname = "name(sroc" + "$by_index_" + ", replace)"
			noi di as res _n  "NOTE: SROC plot name -> sroc$by_index_"
		}
		
		#delimit ;
		graph tw 
			`centre'
			`popcentre'
			`kross'
			`sroc'
			`cregion'
			`pregion'
			`points'
			`rings'
			`idbubble'
			,
			`legendstr' `soptions' `sname'
		;
		#delimit cr
				
		if `"`graphsave'"' != `""' {
			di _n
			noi graph save `graphsave', replace
		}
	end
	
	
/*+++++++++++++++++++	SUPPORTING FUNCTIONS: TRADEOFF ++++++++++++++++++++++++++++++++++++
				   DRAW THE TRADEOFF CURVES
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	cap program drop tradeoff
	program define tradeoff
		version 14.0

		#delimit ;
		syntax varlist,
			cutmat(name) /*cutpoint*/
			coefmat(name) /*coef*/
			sevar(string) /*se*/
			spvar(string) /*sp*/
			[
			groupvar(name) /*Grouping variable*/
			p(integer 0) /*No of parameters in the regression equation*/			
			COLorpalette(string) /*soptions: Colors seperated by space*/
			noPREDiction  /*soptions:no Prediction region*/
			noCURve  /*soptions:no curve*/
			Bubbles /*soptions: size of study by bubbles*/
			BUBbleid /*soptions: Identify the bubbles by index*/
			SPointopt(string) /*soptions: options study points*/
			OPointopt(string) /*soptions: options the overall summary points*/
			CUrveopt(string) /*soptions: options the CI points*/
			CIopt(string) /*soptions: options the CI points*/
			PREDCIopt(string) /*soptions: options the PREDCI points*/
			BUBOpt(string) /*soptions: options the bubble points*/
			BIDopt(string) /*soptions: options the bubble ID points*/
			GRAPHSave(string asis)
			STRAtify  /*Stratify*/
			* /*soptions:Other two-way options*/
			]
			;
		#delimit cr
		
		tempvar se selci seuci sp splci spuci csp Ni rowid Dis tp fp NDis fn tn mu gvar cutoff
		tempname vi bvari
		
		local soptions `"`options'"'
		
		tokenize `varlist'
		gen `tp' = `1'
		gen `fp' = `2'
		gen `fn' = `3'
		gen `tn' = `4'
		local cutoff = "`5'"
		gen `Dis' = (`tp' + `fn')
		gen `NDis' = (`fp' + `tn')		
		gen `Ni' = `Dis' +`NDis'
		gen `rowid' = _n
		
		if "`groupvar'" != "" {
			my_ncod `gvar', oldvar(`groupvar') //
		}
		
		gen `se' = 	`tp'/`Dis'
		gen `sp' =  `tn'/`NDis'				
		
		//If categorical variable, obtain the sroc and the drawing parameters for each level
		if "`groupvar'" != "" {
			qui label list `gvar'
			local nlevels = r(max)
		}
		else {
			gen `gvar' = 1
			local nlevels = 1
		}
		if "`colorpalette'" == "" {
			local colorpalette  "black forest_green cranberry blue sienna orange emerald magenta dknavy gray purple"
		}
		else {
			local kcolors : word count `colorpalette'
			if (`kcolors' < `nlevels') {
				di as error "Please specify the colours to be used for all the `nlevels' levels of `gvar'" 
				di as error "colours should be separated by space in the colorpalette() option"
				exit
			}
		}
		local index 0
		local centre
		local kross
		local sroc
		local points
		local rings
		local idbubble
		local cregion
		local pregion
		local legendlabel
		local legendorder
		//Options
		
		// Study point options
		if `"`spointopt'"' == "" {
			local spointopt "msymbol(o)"
		}
		else {
			local forbidden "mcolor"
			foreach option of local forbidden {
				if strpos(`"`spointopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in spointopt()"
					exit
				}
			}
			if `"`spointopt'"' != "" & strpos(`"`spointopt'"',"msymbol") == 0 {
				local spointopt = `"`spointopt' msymbol(o)"' 
			}
		}
		
		// Curve options
		if `"`curveopt'"' != "" {
			local forbidden "lcolor"
			foreach option of local forbidden {
				if strpos(`"`curveopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in curveopt()"
					exit
				}
			}
		}
		
		// Bubble options
		if `"`bubopt'"' == "" {
			local bubopt "msymbol(Oh)"
		}
		else {
			local forbidden "mcolor"
			foreach option of local forbidden {
				if strpos(`"`bubopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in bubopt()"
					exit
				}
			}
			if `"`bubopt'"' != "" & strpos(`"`bubopt'"',"msymbol") == 0{
				local bubopt = `"`bubopt' msymbol(Oh)"' 
			}
		}
		
		//Bubble ID options
		if `"`bidopt'"' == "" {
			local bidopt "mlabsize(`texts') msymbol(i) mlabel(`rowid')"
		}
		else {
			local forbidden "mcolor mlabcolor "
			foreach option of local forbidden {
				if strpos(`"`bidopt'"',"`option'")  != 0 {
					di as error "Option `option'() not allowed in bidopt()"
					exit
				}
			}
			if `"`bidopt'"' != "" & strpos(`"`bidopt'"',"mlabsize") == 0{
				local bidopt = `"`bidopt' mlabsize(`texts')"' 
			}
			if `"`bidopt'"' != "" & strpos(`"`bidopt'"',"msymbol") == 0{
				local bidopt = `"`bidopt' msymbol(i)"' 
			}
			if `"`bidopt'"' != "" & strpos(`"`bidopt'"',"mlabel") == 0{
				local bidopt = `"`bidopt' mlabel(`rowid')"' 
			}
		}
	
		qui {
			local already 0
			
			local coefnames :colnames `coefmat'
			local nrows = colsof(`coefmat')
			local ase 0
			local bse 0
			local asp 0
			local bsp 0
			forvalues r = 1(1)`nrows' {
				//stop if all values are found
				if `=`ase'*`bse'*`asp'*`bsp'' > 0 {
					break
				}
				
				local rname:word `r' of `coefnames' 
				
				if "`rname'" == "`sevar'" {
					local ase "`r'"
				}
				else if "`rname'" == "`spvar'" {
					local asp "`r'"
				}
				else if "`rname'" == "c.`cutoff'#c.`sevar'" {
					local bse "`r'"
				}
				else if "`rname'" == "c.`cutoff'#c.`spvar'" {
					local bsp "`r'"
				}
			}
				
			forvalues j=1/`nlevels' {
						
				local color:word `j' of `colorpalette'
				qui count if `gvar' == `j'
				if `p' > 0 {
						//Intercept
						local asp`j' = `coefmat'[1, `asp']
						local ase`j' = `coefmat'[1, `ase']
						
						//slope
						local bsp`j' = `coefmat'[1, `bsp']
						local bse`j' = `coefmat'[1, `bse']
				}

				//Curves
				if !`already' {
					qui set obs 500
					local already 1
				}
				qui summ `cutoff' if `gvar' == `j' 
				local maxcutoff = r(max)
				local mincutoff = r(min)
				local stepcutoff = round((`maxcutoff' - `mincutoff')/5, 1)
				local N`j' = r(N)
				
				tempvar x`j' sphat`j' sehat`j' xcregion`j' ycregion`j' xpregion`j' ypregion`j'
				
				range `x`j'' `mincutoff' `maxcutoff'
				
				gen `sehat`j'' = invlogit(`ase`j'' + `bse`j'' * `x`j'')
				gen `sphat`j'' = invlogit(`asp`j'' + `bsp`j''*`x`j'')
				local optimal = `cutmat'[1, 1]  
				local sroc "`sroc' (line `sehat`j'' `x`j'', lcolor(blue) yaxis(1) `curveopt')"
				local sroc "`sroc' (line `sphat`j'' `x`j'', lcolor(red) yaxis(2) `curveopt')"
							
				//points
				local points `"`points' (scatter `se' `cutoff' if `gvar' == `j',  mcolor(blue) `spointopt')"'
				local points `"`points' (scatter `sp' `cutoff' if `gvar' == `j',  mcolor(red) `spointopt')"'
				
				//label
				if `nlevels' > 1 {
					local lab:label `gvar' `j' 
					local legendlabel `"lab(`j' "`lab'") `legendlabel'"'
					local legendorder `"`j'  `legendorder'"'
				}
			}
		}
		
		
		if strpos(`"`soptions'"', "legend") == 0 {
			local legendstr `"legend(order(2 3) lab(2 "Sensitivity") lab(3 "Specificity") cols(1) ring(0) position(6))"'
		}
		if strpos(`"`soptions'"', "xscale") == 0 {
			local soptions `"xscale(range(`mincutoff' `maxcutoff')) `soptions'"'
		}
		if strpos(`"`soptions'"', "yscale") == 0 {
			local soptions `"yscale(range(0 1) axis(1)) yscale(range(0 1) axis(2)) `soptions'"'
		}
		if strpos(`"`soptions'"', "xtitle") == 0 {
			local soptions `"xtitle("Cutoff") `soptions'"'
		}
		if strpos(`"`soptions'"', "ytitle") == 0 {
			local soptions `"ytitle("Sensitivity", axis(1))  ytitle("Specificity", axis(2))`soptions'"'
		}
		if strpos(`"`soptions'"', "ylabel") == 0 {
			local soptions `"ylabel(0(0.2)1, axis(1)) ylabel(0(0.2)1, axis(2)) `soptions'"'
		}
		if strpos(`"`soptions'"', "xlabel") == 0 {
			local soptions `"xlabel(`mincutoff'(`stepcutoff')`maxcutoff', nogrid) `soptions'"'
		}
		if strpos(`"`soptions'"', "graphregion") == 0 {
			local soptions `"graphregion(color(white)) `soptions'"'
		}
		if strpos(`"`soptions'"', "plotregion") == 0 {
			local soptions `"plotregion(margin(medium)) `soptions'"'
		}
		if strpos(`"`soptions'"', "aspectratio") == 0 {
			local soptions `"aspectratio(1) `soptions'"'
		}
		if strpos(`"`soptions'"', "xsize") == 0 {
			local soptions `"xsize(5)  `soptions'"'
		}
		if strpos(`"`soptions'"', "ysize") == 0 {
			local soptions `"ysize(5)  `soptions'"'
		}
		
		#delimit ;
		graph tw 
			`centre'
			`kross'
			`sroc'
			`cregion'
			`pregion'
			`points'
			`rings'
			`idbubble'			
			, xline(`optimal') `legendstr' `soptions' name(tradeoff, replace)
		;
		#delimit cr
		if "$by_index_" != "" {
			qui graph dir
			local gnames = r(list)
			local gname: word $by_index_ of `gnames'
			tokenize `gname', parse(".")
			if "`3'" != "" {
				local ext =".`3'"
			}
			
			qui graph rename tradeoff`ext' tradeoff$by_index_`ext', replace
			
		}
		
		if `"`graphsave'"' != `""' {
			di _n
			noi graph save `graphsave', replace
		}
		
	end
/*+++++++++++++++++++++++++	SUPPORTING FUNCTIONS:  postsim +++++++++++++++++++++++++
							Simulate & summarize posterior distribution
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/	
cap program drop postsim
program define postsim, rclass
version 14.1
	#delimit ;
	syntax [if] [in], todo(string) orderid(varname) studyid(varname) estimates(name) 
	[cveffect(string) regressors(varlist) absoutse(name) absoutsp(name) serrout(name) sprrout(name) orout(name) link(string) se(varname) sp(varname)
	modeles(varname) modellci(varname) modeluci(varname) outplot(string) baselevel(integer 1) bcov(string) wcov(string)
	model(string) by(varname) level(real 95) general comparative abnetwork cbnetwork  varx(varname) nsims(string)  ]
	;
	#delimit cr 
	
	marksample touse, strok
	
	tempname betacoef rawcoef varrawcoef ///
				fullrawcoef fullvarrawcoef X beta sims ///
				popabsout popabsouti poprrout popserrout  popsprrout popserrouti popsprrouti ///
				popseorouti popsporouti popselorouti popsplorouti  ///
				poporout popseorout popsporout poporouti poplorout ///
				popselorout popsplorout poplorouti simvar ///
				popabsoutse popabsoutsei popabsoutsp popabsoutspi serow sprow
	
	tempvar feff sfeff reff sreff reffse sreffse reffsse sreffsse simsreffse simsreffsse  ///
			reffsp sreffsp reffssp sreffssp simsreffsp simsreffssp eta insample ///
			newobs idpair gid rid sid hold holdleft holdright ///
			simmu sumsehat meansehat sumsphat meansphat subset subsetid subsetid1 sumsehat1 ///
			meansehat1 gid1 modelp modelrr modelor modellor modelse ///
			sumserrhat sumsprrhat sumseorhat sumsporhat sumselorhat sumsplorhat ///
			meanserrhat meanseorhat meanselorhat meansprrhat meansporhat ///
			meansplorhat varsp varse vars fisherrho ///
			simvars simvarse simvarsp simfisherrho simrho covar simcovar ///
			simwvars simwvarse simwvarsp wvars wvarse wvarsp
			
	if "`link'" == "cloglog" {
		local invfn "invcloglog"
	}
	else if "`link'" == "loglog" {
		local invfn "1-invcloglog"
	}
	else {
		local invfn "invlogit"
	}		
	
	//Restore 
	qui {
		
		estimates restore `estimates'		
		local predcmd = e(predict)
		gen `insample' = e(sample)
		
		//Coefficients estimates and varcov
		mat `fullrawcoef' = e(b)
		mat `fullvarrawcoef' = e(V)
		
		local ncoef = colsof(`fullrawcoef')
		local rho = 0
		local w = 0	
		if "`model'" == "random" {
			if "`wcov'" != "" {
					if strpos("`wcov'", "ind") != 0 {
						local w = 2
						local wvarnames "`wvarse' `wvarsp'"
						local wsimvarnames "`simwvarse'  `simwvarsp'"
					}
					else if strpos("`wcov'", "id") != 0 {
						local w = 1
						local wvarnames "`wvars'"
						local wsimvarnames "`simwvars'"
					}
			}		
			if strpos("`bcov'", "uns") != 0 {
				local nfeff = `=`ncoef' - 3 - `w''
				
				if "`predcmd'" == "meqrlogit_p" {
					local rho = tanh(`fullrawcoef'[1, `=`ncoef'-`w''])
					local varnames "`varse' `varsp' `fisherrho' `wvarnames'"
					local simvarnames "`simvarse'  `simvarsp' `simfisherrho' `wsimvarnames'"
				}
				else {
					local rho = `fullrawcoef'[1, `=`ncoef'-`w'']/(sqrt(`fullrawcoef'[1, `=`ncoef'-1-`w'']*`fullrawcoef'[1, `=`ncoef'-2-`w'']))
					local varnames "`varse' `varsp' `covar' `wvarnames'"
					local simvarnames "`simvarse'  `simvarsp' `simcovar' `wsimvarnames'"
				}
			}
			else if "`bcov'" =="identity" {
				local nfeff = `=`ncoef' - 1 - `w''
				local varnames "`vars' `wvarnames'"
				local simvarnames "`simvars' `wsimvarnames'"
			}
			else if strpos("`bcov'", "exc") != 0 | strpos("`bcov'", "ind") != 0 {
				local nfeff = `=`ncoef' - 2 - `w''
				
				if strpos("`bcov'", "exc") != 0 {
					if "`predcmd'" == "meqrlogit_p" {
						local rho = tanh(`fullrawcoef'[1, `=`ncoef'-`w''])
						local varnames "`vars' `fisherrho' `wvarnames'"
						local simvarnames "`simvars'  `simfisherrho' `wsimvarnames'"
					}
					else {
						local rho = `fullrawcoef'[1, `=`ncoef'-`w'']/(sqrt(`fullrawcoef'[1, `=`ncoef'-1-`w'']*`fullrawcoef'[1, `=`ncoef'-2-`w'']))
						local varnames "`vars' `covar' `wvarnames'"
						local simvarnames "`simvars'  `simcovar' `wsimvarnames'"
					}
				}
				
				if strpos("`bcov'", "ind") != 0 {
					local varnames "`varse' `varsp' `wvarnames'"
					local simvarnames "`simvarse'  `simvarsp' `wsimvarnames'"
				}
			}
			else if (strpos("`bcov'", "sp") != 0) {
				local varnames "`varsp' `wvarnames'"
				local simvarnames "`simvarsp' `wsimvarnames'"
			}
			else if (strpos("`bcov'", "se") != 0) {
				local varnames "`varse' `wvarnames'"
				local simvarnames "`simvarse' `wsimvarnames'"
			}
			
		}
		else {
			local nfeff = `ncoef'
		}
		//Get the FE parameters and their covariances
		mat `betacoef' = `fullrawcoef'[1, 1..`nfeff']

		//Predict		
		//Fill data if less than 7
		count
		local nobs = r(N)
		if ((`nobs' < 7) & ("`model'" == "random")) {
			local multipler = int(ceil(7/`nobs'))
			qui expand `multipler', gen(`newobs')
		}
		
		predict `feff', xb  //FE
		predict `sfeff', stdp //se of FE
		
		if "`model'" == "random" {
			if "`wcov'" == "" {
				predict `reffse' `reffsp', reffects reses(`sreffse' `sreffsp')  //se=1  sp=2 
					
				gen `reff' = `reffse'*`se' + `reffsp'*`sp'
				gen `sreff' = sqrt(`sreffse'^2 + `sreffsp'^2)
			}
			else {
				predict `reffse' `reffsp' `reffsse' `reffssp', reffects reses(`sreffse' `sreffsp' `sreffsse' `sreffssp')  //intse=1  intsp=2 slopese=3  slopesp=4 
					
				gen `reff' = `reffse'*`se' + `reffsse'*`se'  + `reffsp'*`sp' + `reffssp'*`sp'
				gen `sreff' = sqrt(`sreffse'^2 + `sreffsp'^2 + `sreffsse'^2 + `sreffssp'^2)
			}
						
			
			gen `eta' = `feff' + `reff' //linear predictor				
			gen `modelse' = sqrt(`sreffse'^2 + `sreffsp'^2 + `sfeff'^2) if `insample'==1
			
		}
		else {
			gen `eta' = `feff' 
			gen `modelse' = `sfeff' if `insample'==1
		}
		
		//Revert to original data if filler data was generated
		if (("`model'" == "random") & (`nobs' < 7))  {
			keep if !`newobs'
		}
		
		//Smooth p estimates
		gen `modelp' = `invfn'(`eta') if `insample' == 1
		
		//identifiers
		gsort -`insample' `se' `orderid'
		*egen `rid' = seq() if `insample'==1  //rowid
		
		if "`comparative'`cbnetwork'`abnetwork'" != ""  {
			count if `insample' == 1  
			local Nobs = `=r(N)'*0.25
			
			egen `sid' = group(`studyid') //Studyid
			
			if "`regressors'" != "" {
				gsort -`insample' `se' `regressors' `orderid'
			}
			if "`cbnetwork'" != "" {
				sort `insample' `se' `varx' `orderid'
			}
			egen `gid' = seq() if `insample' == 1, f(1) t(`Nobs') b(1) 
			egen `idpair' = group(`varx') if `insample' == 1
			sort `se' `gid' `idpair'
			
			/*
			egen `gid' = group(`studyid' `by') if `insample'==1  
			gsort `gid' `se' `orderid' `varx'
			bys `gid' `se': egen `idpair' = seq()
			*/
			egen `rid' = seq() if `insample'==1  //rowid
			
			if "`abnetwork'" == "" {
				gen `modelrr' = `modelp'[_n] / `modelp'[_n-1] if (`gid'[_n]==`gid'[_n-1]) & (`se'[_n]==`se'[_n-1]) & (`idpair' == 2)
				gen `modelor' = (`modelp'[_n]/(1 - `modelp'[_n])) / (`modelp'[_n-1]/(1 - `modelp'[_n-1])) if (`gid'[_n]==`gid'[_n-1]) & (`se'[_n]==`se'[_n-1]) & (`idpair' == 2)
				gen `modellor' = ln(`modelor')
			}
		}
		else {
			egen `sid' = group(`studyid') if `insample'==1  //Studyid
			gsort `insample' `se' `orderid' 
			bys `insample' `se' : egen `gid' = seq()  //gid
			bys `insample' : egen `rid' = seq()   //rowid
		}
				
		//Generate designmatrix
		local colnames :colnames `betacoef'
		local nvars: word count `colnames'
		forvalues i=1(1)`nvars' {
			tempvar v`i' beta`i'
			
			local var`i' : word `i' of `colnames'	
			local left
			local right
			local param
			local rightleft
			local rightright
			local leftleft
			local leftright
			
			//Split the term
			if strpos("`var`i''", "#") != 0 {
				tokenize `var`i'', parse("#")
				if "`5'" != "" {
					local left = "`1'"
					local right = "`3'"
					local param = "`5'"
				}
				else {
					local right = "`1'"
					local param = "`3'"
				}
				
				tokenize `param', parse(.)
				local param  = "`3'"
				
				tokenize `right', parse(.)
				local rightleft = "`1'"
				local rightright = "`3'"
				
				if "`left'" != "" {
					tokenize `left', parse(.)
					local leftleft = "`1'"
					local leftright = "`3'"
				}
			}
			else {
				local param "`var`i''"
			}
			
			//Constants
			if "`right'" == "" {
				cap confirm var `param'
				if _rc!=0  {	
					gen `v`i'' = 0
				}
				else {				
					gen `v`i'' = `param'
				}
			}
			
			//Main effects
			if "`right'" != "" & "`left'" == ""  {
				//Continous  
				if strpos("`rightleft'", "c") != 0 {
					gen `v`i'' = `param'*`rightright'
				}
				else {
					//Categorical
					if strpos("`rightleft'", "bn") != 0 {
						local rightleft = ustrregexra("`rightleft'", "bn", "")
					}
					if strpos("`rightleft'", "b") != 0 {
						local rightleft = ustrregexra("`rightleft'", "b", "")
					}		
					gen `v`i'' = 0 +  1*`param'*(`rightright' == `rightleft')
				}
			}
			
			//Interactions
			if "`left'" != "" {	
				//continous left
				if strpos("`leftleft'", "c") == 1 {
					local factorleft 0
				}
				else {
					local factorleft 1
				}
				
				//continous right
				if strpos("`rightleft'", "c") == 1 {
					local factorright 0
				}
				else {
					local factorright 1
				}
						
				local part 1
				local prefices "`leftleft' `rightleft'"
				foreach prefix of local prefices {
					//Continous  
					if strpos("`prefix'", "c") != 0 {
						local level`part' = `part'
					}
					else {
						//Categorical
						if strpos("`prefix'", "bn") != 0 {
							local level`part' = ustrregexra("`prefix'", "bn", "")
						}
						else if strpos("`prefix'", "b") != 0 {
							local level`part' = ustrregexra("`prefix'", "b", "")
						}
						else if strpos("`prefix'", "o") != 0 {
							local level`part' = ustrregexra("`prefix'", "o", "")
						}
						else {
							local level`part' = `prefix'
						}
					}
					local ++part
					
				}
				gen `v`i'' = ((`factorleft'*(`leftright'==`level1') + !`factorleft'*`leftright') * (`factorright'*(`rightright'==`level2') + !`factorright'*`rightright'))*`param'
				
			}
			local vnamelist "`vnamelist' `v`i''"
			local bnamelist "`bnamelist' `beta`i''"
		}
		
		if "`model'" == "random" {
			//Add varnames
			local bnamelist "`bnamelist' `varnames'"
		}
		
		set matsize `nsims'

		//make matrices from the dataset
		//roweq(`idpair')
		mkmat `vnamelist' if `insample'==1, matrix(`X')  rownames(`rid')
		
		tempvar present
		gen `present' = 1	
		
		//Simulate the parameters
		if `nobs' < `nsims' {
			set obs `nsims'
		}
		
		drawnorm `bnamelist', n(`nsims') cov(`fullvarrawcoef') means(`fullrawcoef') seed(1)

		mkmat `bnamelist', matrix(`beta')
		
		//Subset the matrix
		if `ncoef' > `nfeff' {
			mat `simvar' = `beta'[1..`nsims', `=`nfeff'+1'..`ncoef']
			mat `beta' = `beta'[1..`nsims', 1..`nfeff']
		}
		
		mat `sims' = `beta'*`X''
		
		//Construct the names
		local ncols = colsof(`sims') //length of the vector
		local cnames :colnames `sims'
		
		local matcolnames
		
		forvalues c=1(1)`ncols' {
			local matrid : word `c' of `cnames'

			tempvar ferow`matrid'
			local matcolnames = "`matcolnames' `ferow`matrid''"
		}
		
		if "`model'" == "random" {
			//Append the var matrix
			mat `sims' = (`sims', `simvar')
			
			//Add varnames
			local matcolnames "`matcolnames' `simvarnames'"
		}
		
		//pass the names
		matname `sims' `matcolnames', col(.) explicit

		//Bring the matrix to the dataset
		svmat `sims', names(col)
		
		if "`model'" == "random" {
			if "`wcov'" !="" {
					if strpos("`wcov'", "ind") != 0 {
						if "`predcmd'" == "meqrlogit_p" {
							gen `simsreffsse' = exp(`simwvarse') 
							gen `simsreffssp' = exp(`simwvarsp')
						}
					}
					else if strpos("`wcov'", "id") != 0 {
						if "`predcmd'" == "meqrlogit_p" {
							gen `simsreffsse' = exp(`simwvars') 
							gen `simsreffssp' = exp(`simwvars')
						}
					}
			}
			if "`bcov'" !="" {
				if (strpos("`bcov'", "uns") == 0) & (strpos("`bcov'", "exc") == 0) {
					gen `simrho' = 0				
				}
				if (strpos("`bcov'", "exc") != 0) | (strpos("`bcov'", "id") != 0){
					gen `simvarse' = `simvars' 
					gen `simvarsp' = `simvars'
				}
				if (strpos("`bcov'", "sp") != 0) {
					gen `simvarse' = 0	
				}
				if (strpos("`bcov'", "se") != 0) {
					gen `simvarsp' = 0	
				}
				if "`predcmd'" == "melogit_p" {				
					//Truncate the values to zero
					replace `simvarse' = 0 if `simvarse' < 0
					replace `simvarsp' = 0 if `simvarsp' < 0
				}
				if (strpos("`bcov'", "uns") != 0) | strpos("`bcov'", "exc") != 0   {
					if "`predcmd'" == "meqrlogit_p" {
						gen `simrho' = tanh(`simfisherrho')
					}
					else {				
						//Truncate the values to zero
						replace `simcovar' = 0 if `simvarse' == 0 & `simvarsp' == 0	
						gen `simrho' = `simcovar'/sqrt(`simvarsp'*`simvarse')
					}
				}
				if "`predcmd'" == "meqrlogit_p" {
					gen `simsreffse' = exp(`simvarse') //Marginal se
					gen `simsreffsp' = sqrt((1 - (`simrho')^2)*(exp(`simvarsp')^2)) //Conditional se
				}
				else {			
					gen `simsreffse' = sqrt(`simvarse') //Marginal se
					gen `simsreffsp' = sqrt((1 - (`simrho')^2)*(`simvarsp')) //Conditional se
				}
			}
		}
		
		//# of obs
		count if `insample' == 1
		local nobs = r(N)
		
		//Generate the p's and r's
		forvalues j=1(1)`nobs' { 
			
			//Study index
			sum `sid' if `rid' == `j'
			local sidj = r(min)
			
			//group index
			sum `gid' if `rid' == `j'
			local index = r(min)
			
			//Parameter
			sum `se' if `rid' == `j'
			local param = r(min)
			
			if `param' == 1 {
				local prefix "se"
			}
			else {
				local prefix "sp"
			}
			
			tempvar phat`j' `prefix'hat`index' 
			
			//Figure out if new study or not
			local new`prefix'study 1
			if `j' > 1 {
				sum `sid' if (`rid' == `=`j'-1' )
				local prevsid = r(min)
				local new`prefix'study = (`sidj' != `prevsid')
			}				
				
			if "`comparative'`cbnetwork'`abnetwork'" == "" {	
				if "`model'" == "random" {
					//EB re 
					sum `reff' if `rid' == `j'
					local reff_`j' = r(mean)
					
					//re  simulate					
					if `new`prefix'study' {						
						//Generate the variables
						tempvar `prefix'study`sidj'  
						gen ``prefix'study`sidj'' = rnormal(0, `simsreff`prefix'') 
					}
					
					gen `phat`j'' = `invfn'(`reff_`j'' + ``prefix'study`sidj'' + `ferow`j'')
										
				}
				else {
					//Create variable
					gen `phat`j'' = `invfn'(`ferow`j'')
				}
			}
			
			if "`comparative'" != "" | "`cbnetwork'" != "" | "`abnetwork'" != "" {
				sum `idpair' if `rid' == `j'
				local pair = r(min)
				
				//Generate the variables
				tempvar `prefix'hat_`pair'`index' 
				
				if "`model'" == "random" {
					//EB re 
					sum `reff' if `rid' == `j' 
					local reff_`j' = r(mean)
						
					if `pair' == 1 {					
						//re  simulate
						if `new`prefix'study' {						
							//Generate the variables
							tempvar `prefix'study`sidj'  
							gen ``prefix'study`sidj'' = rnormal(0, `simsreff`prefix'') 
						}
					}
					if "`wcov'" != "" {
						tempvar s`prefix'study`j' 
						
						gen `s`prefix'study`j'' = rnormal(0, `simsreffs`prefix'')
						replace ``prefix'study`sidj'' = ``prefix'study`sidj'' + `s`prefix'study`j''
					}
					gen `phat`j'' = `invfn'(`reff_`j'' + ``prefix'study`sidj'' + `ferow`j'')				
				}
				else {
					gen `phat`j'' = `invfn'(`ferow`j'')
				}
				if "`abnetwork'" == "" {					
					//Create the pairs
					gen ``prefix'hat_`pair'`index'' = `phat`j''
					
					if `pair' == 2 {
						tempvar `prefix'rrhat`index' `prefix'orhat`index' `prefix'lorhat`index'
						
						gen ``prefix'rrhat`index'' = ``prefix'hat_2`index'' / ``prefix'hat_1`index''
						gen ``prefix'orhat`index'' = (``prefix'hat_2`index''/(1 - ``prefix'hat_2`index'')) / (``prefix'hat_1`index'' /(1 - ``prefix'hat_1`index''))
						gen ``prefix'lorhat`index'' = ln(``prefix'orhat`index'')
					}
				}
			}
		}			
				
		if "`todo'" == "p" {
			//Do twice; 
			local shortnames "se sp"
			local parameters "`se' `sp'"
			local fullnames "Sensitivity Specificity"
			
			forvalues run=1(1)2 {
			
				local prefix : word `run' of `shortnames'
				local name : word `run' of `fullnames'
				local parameter : word `run' of `parameters'
			
				//Summarize 
				local nrow`prefix' = rowsof(`absout`prefix'') //length of the vector
				local rnames :rownames `absout`prefix''
				local eqnames :roweq `absout`prefix''
				local newnrows = 0
				local mindex = 0
					
				foreach vari of local eqnames {		
					local ++mindex
					local group : word `mindex' of `rnames'
					
					//Skip if continous variable
					if ((strpos("`vari'", "_") == 1) & ("`group'" != "`name'") & ("`group'" != "Overall")) | (`nrow`prefix'' > 1 & `mindex' == 1) {
						continue
					}					
					
					cap drop `subset' `subsetid'
					
					if ("`group'" != "`name'") & ("`group'" != "Overall") {
						if strpos("`vari'", "*") == 0 {
							cap drop `hold'
							decode `vari', gen(`hold')
							cap drop `subset'
							local latentgroup = ustrregexra("`group'", "_", " ")
							gen `subset' = 1 if `hold' == "`latentgroup'" & `insample' == 1 & `parameter'
						}
						else {
							tokenize `vari', parse("*")
							local leftvar = "`1'"
							local rightvar = "`3'"
							
							tokenize `group', parse("|")
							local leftgroup = "`1'"
							local rightgroup = "`3'"
							
							cap drop `holdleft' `holdright'
							decode `leftvar', gen(`holdleft')
							decode `rightvar', gen(`holdright')
							cap drop `subset'
							local latentleftgroup = ustrregexra("`leftgroup'", "_", " ")
							local latentrightgroup = ustrregexra("`rightgroup'", "_", " ")
							gen `subset' = 1 if (`holdleft' == "`latentleftgroup'") & (`holdright' == "`latentrightgroup'") & (`insample' == 1) & `parameter'
						}					
						egen `subsetid' = seq() if `subset' == 1 
					}
					else {
						//All
						gen `subset' = 1 if `insample' == 1 & `parameter'
						egen `subsetid' = seq() if `subset' == 1
					}
					
					count if `subset' == 1 
					local nsubset = r(N)
					
					//Compute mean of simulated values
					cap drop `sum`prefix'hat' `mean`prefix'hat' 
					gen `sum`prefix'hat' = 0	
					
					forvalues j=1(1)`nsubset' {
						//Study index
						sum `rid' if `subsetid' == `j'
						local index = r(min)
											
						replace `sum`prefix'hat' = `sum`prefix'hat' + `phat`index''
					}
					
					//Obtain mean of modelled estimates
					sum `modelp' if `subset' == 1 
					local meanmodel`prefix' = r(mean)
							
					gen `mean`prefix'hat' = `sum`prefix'hat'/`nsubset'
					
					//Standard error
					sum `mean`prefix'hat'	
					local post`prefix'se = r(sd)
								
					//Obtain the quantiles
					centile `mean`prefix'hat', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
					local `prefix'median = r(c_1) //Median
					local `prefix'lowerp = r(c_2) //Lower centile
					local `prefix'upperp = r(c_3) //Upper centile
									
					mat `popabsout`prefix'i' = (`meanmodel`prefix'', `post`prefix'se', ``prefix'median', ``prefix'lowerp', ``prefix'upperp')
					mat rownames `popabsout`prefix'i' = `vari':`group'
									
					//Stack the matrices
					local ++newnrows
					if `newnrows' == 1 {
						mat `popabsout`prefix'' = `popabsout`prefix'i'	
					}
					else {
						mat `popabsout`prefix'' = `popabsout`prefix''	\  `popabsout`prefix'i'
					}
				}
			}
			
			mat `serow' = J(1, 5, .)
			mat `sprow' = J(1, 5, .)
					
			mat rownames `serow' = "Sensitivity"
			mat rownames `sprow' = "Specificity" 
			if `nrowse' > 1 {
				mat `popabsoutse' = `serow' \ `popabsoutse'
			}
			if `nrowsp' > 1 {
				mat `popabsoutsp' = `sprow' \ `popabsoutsp'
			}
			mat `popabsout' = `popabsoutse' \ `popabsoutsp'
			
			/*if `nrowse' > 1 {
					mat `popabsoutse' = `serow' \ `popabsoutse'
			}
			else {
				mat `popabsout' =  `popabsoutse'
			}
			if `nrowsp' > 1 {
				mat `popabsout' = `popabsout' \ `sprow' \ `popabsoutsp'
			}
			else {
				mat `popabsout' = `popabsout' \ `popabsoutsp'
			}*/
		}
		
			local runs = 2 //Do twice; 
			local shortnames "se sp"
			local parameters "`se' `sp'"
			local fullnames "Sensitivity Specificity"
			
		if "`todo'" == "r" {			
			forvalues run = 1(1)`runs' {
			
				local prefix : word `run' of `shortnames'
				local name : word `run' of `fullnames'
				local parameter : word `run' of `parameters'
			
				//Summarize RR
				local nrows = rowsof(``prefix'rrout') //length of the vector
				local rnames :rownames ``prefix'rrout'
				local eqnames :roweq  ``prefix'rrout'
				local newnrows 0
								
				if "`comparative'`cbnetwork'" == "" {
					local catvars : list uniq eqnames	
					foreach vari of local catvars {					
						//Skip if continous variable
						if (strpos("`vari'", "_") == 1) {
							continue
						}
						
						cap drop `hold'	
						decode `vari', gen(`hold')
						label list `vari'
						local ngroups = r(max)
						local baselab:label `vari' `baselevel'
						
						//count in basegroup
						tempvar mean`prefix'hat`baselevel' mean`prefix'rrhat`baselevel' mean`prefix'orhat`baselevel' mean`prefix'lorhat`baselevel' gid`baselevel' sum`prefix'hat`baselevel' subsetid`baselevel'
						tempname pop`prefix'rrouti`baselevel' pop`prefix'orouti`baselevel' pop`prefix'lorouti`baselevel'
						
						count if `vari' == `baselevel' & `insample' == 1 & `parameter'
						local ngroup`baselevel' = r(N)
						
						egen `subsetid`baselevel'' = group(`rid') if `vari' == `baselevel' & `insample' == 1 & `parameter'
						
						cap drop `sum`prefix'hat`baselevel'' `mean`prefix'hat`baselevel''
						//basegroup
						gen `sum`prefix'hat`baselevel'' = 0
						
						forvalues j=1(1)`ngroup`baselevel'' {
							sum `rid' if `subsetid`baselevel'' == `j'
							local index = r(min)
							
							qui replace `sum`prefix'hat`baselevel'' = `sum`prefix'hat`baselevel'' + `phat`index''
						}
						gen `mean`prefix'hat`baselevel'' = `sum`prefix'hat`baselevel''/`ngroup`baselevel''
						
						sum `modelp' if `vari' == `baselevel' & `insample' == 1 & `parameter'
						local meanmodel`prefix'`baselevel' = r(mean)
						
						mat `pop`prefix'rrouti`baselevel'' = (1, 0, 1, 1, 1)
						mat `pop`prefix'orouti`baselevel'' = (1, 0, 1, 1, 1)
						mat `pop`prefix'lorouti`baselevel'' = (1, 0, 1, 1, 1)
						
						local baselab = ustrregexra("`baselab'", " ", "_")
						mat rownames `pop`prefix'rrouti`baselevel'' = `vari':`baselab'
						mat rownames `pop`prefix'orouti`baselevel'' = `vari':`baselab'
						mat rownames `pop`prefix'lorouti`baselevel'' = `vari':`baselab'
						
						//Other groups
						forvalues g=1(1)`ngroups' {
							if `g' != `baselevel' {
								tempvar mean`prefix'hat`g' mean`prefix'rrhat`g' mean`prefix'orhat`g' mean`prefix'lorhat`g' gid`g' sum`prefix'hat`g' subsetid`g'
								tempname pop`prefix'rrouti`g' pop`prefix'orouti`g' pop`prefix'lorouti`g'
								
								local glab:label `vari' `g'
								count if `vari' == `g' & `insample' == 1 & `parameter'
								local ngroup`g' = r(N)	
								egen `subsetid`g'' = group(`rid') if `vari' == `g' & `insample' == 1 & `parameter'
								
								gen `sum`prefix'hat`g'' = 0
								
								//Group of interest
								forvalues j=1(1)`ngroup`g'' {
									sum `rid' if `subsetid`g'' == `j'
									local index = r(min)			
									qui replace `sum`prefix'hat`g'' = `sum`prefix'hat`g'' + `phat`index''
								}
								
								gen `mean`prefix'hat`g'' = `sum`prefix'hat`g''/`ngroup`g''
								
								//Generate R 
								gen `mean`prefix'rrhat`g'' = `mean`prefix'hat`g'' / `mean`prefix'hat`baselevel''
								gen `mean`prefix'orhat`g'' = (`mean`prefix'hat`g''/(1 - `mean`prefix'hat`g'')) / (`mean`prefix'hat`baselevel''/(1 - `mean`prefix'hat`baselevel''))
								gen `mean`prefix'lorhat`g'' = ln(`mean`prefix'orhat`g'')

								//Obtain mean of modelled estimates
								sum `modelp' if `vari' == `g' & `insample' == 1 & `parameter'
								local meanmodel`prefix'`g' = r(mean)
								local meanmodel`prefix'rr`g' = `meanmodel`prefix'`g'' / `meanmodel`prefix'`baselevel''
								local meanmodel`prefix'or`g' = (`meanmodel`prefix'`g''/(1 - `meanmodel`prefix'`g'')) / (`meanmodel`prefix'`baselevel''/(1 - `meanmodel`prefix'`baselevel''))
								local meanmodel`prefix'lor`g' = ln((`meanmodel`prefix'`g''/(1 - `meanmodel`prefix'`g'')) / (`meanmodel`prefix'`baselevel''/(1 - `meanmodel`prefix'`baselevel'')))
								
								//Standard error
								sum `mean`prefix'rrhat`g''
								local post`prefix'rrse = r(sd)
								
								sum `mean`prefix'orhat`g''
								local post`prefix'orse = r(sd)
								
								sum `mean`prefix'lorhat`g''
								local post`prefix'lorse = r(sd)
								
								//Obtain the quantiles
								centile `mean`prefix'rrhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
								local `prefix'medianrr = r(c_1) //Median
								local `prefix'lowerprr = r(c_2) //Lower centile
								local `prefix'upperprr = r(c_3) //Upper centile
								
								centile `mean`prefix'orhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
								local `prefix'medianor = r(c_1) //Median
								local `prefix'lowerpor = r(c_2) //Lower centile
								local `prefix'upperpor = r(c_3) //Upper centile
								
								centile `mean`prefix'lorhat`g'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
								local `prefix'medianlor = r(c_1) //Median
								local `prefix'lowerplor = r(c_2) //Lower centile
								local `prefix'upperplor = r(c_3) //Upper centile
														
								mat `pop`prefix'rrouti`g'' = (`meanmodel`prefix'rr`g'', `post`prefix'rrse', ``prefix'medianrr', ``prefix'lowerprr', ``prefix'upperprr')
								mat `pop`prefix'orouti`g'' = (`meanmodel`prefix'or`g'', `post`prefix'orse', ``prefix'medianor', ``prefix'lowerpor', ``prefix'upperpor')
								mat `pop`prefix'lorouti`g'' = (`meanmodel`prefix'lor`g'', `post`prefix'lorse', ``prefix'medianlor', ``prefix'lowerplor', ``prefix'upperplor')
								
								local glab = ustrregexra("`glab'", " ", "_")
								mat rownames `pop`prefix'rrouti`g'' = `vari':`glab'
								mat rownames `pop`prefix'orouti`g'' = `vari':`glab'
								mat rownames `pop`prefix'lorouti`g'' = `vari':`glab'
							}
							if `g' == 1 {
								mat `pop`prefix'rrouti' = `pop`prefix'rrouti`g''
								mat `pop`prefix'orouti' = `pop`prefix'orouti`g''
								mat `pop`prefix'lorouti' = `pop`prefix'lorouti`g''
							}
							else {
								//Stack the matrices
								mat `pop`prefix'rrouti' = `pop`prefix'rrouti'	\  `pop`prefix'rrouti`g''
								mat `pop`prefix'orouti' = `pop`prefix'orouti'	\  `pop`prefix'orouti`g''
								mat `pop`prefix'lorouti' = `pop`prefix'lorouti'	\  `pop`prefix'lorouti`g''
							}
						}
						//Stack the matrices
						local ++newnrows
						if `newnrows' == 1 {
							mat `pop`prefix'rrout' = `pop`prefix'rrouti'
							mat `pop`prefix'orout' = `pop`prefix'orouti'
							mat `pop`prefix'lorout' = `pop`prefix'lorouti'
						}
						else {
							mat `pop`prefix'rrout' = `pop`prefix'rrout'	\  `pop`prefix'rrouti'
							mat `pop`prefix'orout' = `pop`prefix'orout'	\  `pop`prefix'orouti'
							mat `pop`prefix'lorout' = `pop`prefix'lorout'	\  `pop`prefix'lorouti'
						}
					}
				}
				
				if "`comparative'" != "" | "`cbnetwork'" != "" {
					//Comparative R
					local mindex 0
					local newnrows 0
					foreach vari of local eqnames {		
						local ++mindex
						local group : word `mindex' of `rnames'
						
						cap drop `subset'
						if "`group'" != "Overall" {
							cap drop `hold'
							decode `vari', gen(`hold')
							
							local latentgroup = ustrregexra("`group'", "_", " ")
							gen `subset' = 1 if `hold' == "`latentgroup'"  & `insample' == 1 & `parameter'
						}
						else {
							//All
							gen `subset' = 1  if `insample' == 1 & `parameter'
						}
						cap drop `subsetid'
						egen `subsetid' = seq()	 if `subset' == 1
						
						count if `subset' == 1 & `idpair' == 2
						local nsubset = r(N)
						
						//Compute mean of simulated values
						cap drop `sum`prefix'rrhat' `sum`prefix'orhat' `sum`prefix'lorhat'
						
						gen `sum`prefix'rrhat' = 0
						gen `sum`prefix'orhat' = 0	
						gen `sum`prefix'lorhat' = 0
						forvalues j=1(1)`nobs' { 
							sum `gid' if `subsetid' == `j'
							local index = r(min)
							
							sum `idpair' if `subsetid' == `j'
							local pair = r(min)
							
							if `pair' == 2 {
								qui replace `sum`prefix'rrhat' = `sum`prefix'rrhat' + ``prefix'rrhat`index''
								qui replace `sum`prefix'orhat' = `sum`prefix'orhat' + ``prefix'orhat`index''
								qui replace `sum`prefix'lorhat' = `sum`prefix'lorhat' + ``prefix'lorhat`index''
							}
						}
						
						//Obtain mean of modelled estimates
						sum `modelrr' if `subset' == 1
						local meanmodel`prefix'rr = r(mean)
						
						sum `modelor' if `subset' == 1
						local meanmodel`prefix'or = r(mean)
						
						sum `modellor' if `subset' == 1
						local meanmodel`prefix'lor = r(mean)
						
						cap drop `mean`prefix'rrhat' `mean`prefix'orhat' `mean`prefix'lorhat'
						gen `mean`prefix'rrhat' = `sum`prefix'rrhat'/`nsubset'
						gen `mean`prefix'orhat' = `sum`prefix'orhat'/`nsubset'
						gen `mean`prefix'lorhat' = `sum`prefix'lorhat'/`nsubset'
						
						//Standard error
						sum `mean`prefix'rrhat'	
						local post`prefix'rrse = r(sd)
						
						sum `mean`prefix'orhat'	
						local post`prefix'orse = r(sd)
						
						sum `mean`prefix'lorhat'	
						local post`prefix'lorse = r(sd)
						
						//Obtain the quantiles
						centile `mean`prefix'rrhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						local median`prefix'rr = r(c_1) //Median
						local lowerp`prefix'rr = r(c_2) //Lower centile
						local upperp`prefix'rr = r(c_3) //Upper centile
						
						centile `mean`prefix'orhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						local median`prefix'or = r(c_1) //Median
						local lowerp`prefix'or = r(c_2) //Lower centile
						local upperp`prefix'or = r(c_3) //Upper centile
						
						centile `mean`prefix'lorhat' , centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						local median`prefix'lor = r(c_1) //Median
						local lowerp`prefix'lor = r(c_2) //Lower centile
						local upperp`prefix'lor = r(c_3) //Upper centile
						
						mat `pop`prefix'rrouti' = (`meanmodel`prefix'rr', `post`prefix'rrse', `median`prefix'rr', `lowerp`prefix'rr', `upperp`prefix'rr')
						mat `pop`prefix'orouti' = (`meanmodel`prefix'or', `post`prefix'orse', `median`prefix'or', `lowerp`prefix'or', `upperp`prefix'or')
						mat `pop`prefix'lorouti' = (`meanmodel`prefix'lor', `post`prefix'lorse', `median`prefix'lor', `lowerp`prefix'lor', `upperp`prefix'lor')
						mat rownames `pop`prefix'rrouti' = `vari':`group'
						mat rownames `pop`prefix'orouti' = `vari':`group'
						mat rownames `pop`prefix'lorouti' = `vari':`group'
						
						//Stack the matrices
						local ++newnrows
						if `newnrows' == 1 {
							mat `pop`prefix'rrout' = `pop`prefix'rrouti'
							mat `pop`prefix'orout' = `pop`prefix'orouti'
							mat `pop`prefix'lorout' = `pop`prefix'lorouti'
						}
						else {
							mat `pop`prefix'rrout' = `pop`prefix'rrout'	\  `pop`prefix'rrouti'
							mat `pop`prefix'orout' = `pop`prefix'orout'	\  `pop`prefix'orouti'
							mat `pop`prefix'lorout' = `pop`prefix'lorout'	\  `pop`prefix'lorouti'
						}
					}
				}
			}
			mat `serow' = J(1, 5, .)
			mat `sprow' = J(1, 5, .)
					
			mat rownames `serow' = "Relative Sensitivity"
			mat rownames `sprow' = "Relative Specificity" 
			if "`cveffect'" == "sesp" {
				mat `poprrout' = `serow' \ `popserrout' \ `sprow' \ `popsprrout'
				mat `poporout' = `serow' \ `popseorout' \ `sprow' \ `popsporout'
			}
			else if "`cveffect'" == "se" {
				mat `poprrout' = `serow' \ `popserrout' 
				mat `poporout' = `serow' \ `popseorout' 
			}
			else if "`cveffect'" == "sp" {
				mat `poprrout' =  `sprow' \ `popsprrout'
				mat `poporout' =  `sprow' \ `popsporout'
			}
		}
		if "`todo'" == "smooth" {
			replace `modeles' = `modelp' if  `insample' == 1
			//Smooth p's
			if "`outplot'" == "abs" {
				//postci
				local critvalue -invnorm((100-`level')/200)
				if "`link'" == "loglog" {
					local sign -
				}
				replace `modellci' = `invfn'(`eta' - `sign' `critvalue'*`modelse') if  `insample' == 1 //lower
				replace `modeluci' = `invfn'(`eta' +  `sign' `critvalue'*`modelse') if  `insample' == 1 //upper
			}
			
			//Smooth r's
			if "`outplot'" == "rr" | "`outplot'" == "or" {
				if "`outplot'" == "rr" { 
					replace `modeles' = `modelrr' if  `insample' == 1
				}
				else {
					replace `modeles' = `modelor' if  `insample' == 1
				}
				
				//Do twice; 
				local shortnames "se sp"
				local parameters "`se' `sp'"
				forvalues run = 1(1)2 {
					local prefix : word `run' of `shortnames'
					local parameter : word `run' of `parameters'
					
					*sum `gid' if `insample' == 1 
					count if `insample' == 1 & `idpair' == 2 & `parameter'
					local nstudies = r(N)
					
					cap drop `subsetid'
					egen `subsetid' = seq()	if `insample' == 1 & `idpair' == 2 & `parameter'
				
					forvalues j=1(1)`nstudies' {
						
						sum `gid' if `subsetid' == `j'
						local index = r(min)
						
						//Obtain the quantiles
						if "`outplot'" == "rr" {
							centile ``prefix'rrhat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}
						else {
							centile ``prefix'orhat`index'', centile(50 `=(100-`level')/2' `=100 - (100-`level')/2')
						}
						local median = r(c_1) //Median
						local lowerp = r(c_2) //Lower centile
						local upperp = r(c_3) //Upper centile
						
						replace `modellci' = `lowerp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1 & `parameter'
						replace `modeluci' = `upperp' if (`gid' == `index') & (`idpair' == 2) &  `insample' == 1 & `parameter'
					}
				}				
			}			
		}
		
		drop if `present' != 1
	}
		
	//Return matrices
	if "`todo'" =="p" {
		mat colnames `popabsoutse' = Mean SE Median Lower Upper
		mat colnames `popabsoutsp' = Mean SE Median Lower Upper
		mat colnames `popabsout' = Mean SE Median Lower Upper
		
		return matrix outmatrixse = `popabsoutse'
		return matrix outmatrixsp = `popabsoutsp'
		return matrix outmatrix = `popabsout'
	}
	if "`todo'" == "r" {
		mat colnames `poprrout' = Mean SE Median Lower Upper
		mat colnames `poporout' = Mean SE Median Lower Upper
		if "`cveffect'" != "sp" {
			mat colnames `popserrout' = Mean SE Median Lower Upper
			mat colnames `popseorout' = Mean SE Median Lower Upper
			mat colnames `popselorout' = Mean SE Median Lower Upper
		}
		
		if "`cveffect'" != "se" {
			mat colnames `popsprrout' = Mean SE Median Lower Upper
			mat colnames `popsporout' = Mean SE Median Lower Upper
			mat colnames `popsplorout' = Mean SE Median Lower Upper
		}
		
		return matrix oroutmatrix = `poporout'
		return matrix rroutmatrix = `poprrout'
		
		if "`cveffect'" != "sp" {
			return matrix rroutmatrixse = `popserrout'
			return matrix oroutmatrixse = `popseorout'
			return matrix loroutmatrixse = `popselorout'
		}
		if "`cveffect'" != "se" {
			return matrix rroutmatrixsp = `popsprrout'
			return matrix oroutmatrixsp = `popsporout'
			return matrix loroutmatrixsp = `popsplorout'
		}
	}

end

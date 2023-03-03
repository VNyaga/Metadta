cap log close
cap log using "C:\DATA\WIV\Projects\GitHub\Metadta\RSM2023\Logs\stataout.txt", text replace
set more off

*Install the commands
ssc install metadta
ssc install midas
ssc install metandi

/*========================================================================================*/
***Section 5.1 
/*========================================================================================*/
*------Table 2 (column 2), Figure 1, Figure 4 (top left)
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/telomerase.csv", delim(",")clear
metadta tp fp fn tn, studyid(study) sumtable(all)

*------Table 2 (column 3), Figure 4 (top middle)
rename study studyid
midas tp fp fn tn, plot sroc(both) res(all) nip(1)
*midas tp fp fn tn, id(studyid) bfor(dss) fordata  //forest plot

*------Table 2 (column 3)
*metandi tp fp fn tn //Non-convergence
metandi tp fp fn tn, nip(7)

*------ Figure 4 (top right)
metandiplot tp fp fn tn, graphregion(color(white)) legend(ring(0) 	col(1) bplacement(6))
/*========================================================================================*/
***Section 5.2
/*========================================================================================*/
*------Table 3 (column 2), Figure 5 (top left)
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/dementia.csv", delim(",")clear
metadta tp fp fn tn, studyid(study) sumtable(all)

*------Table 3 (column 3), Figure 5 (top middle)
midas tp fp fn tn, plot sroc(both) res(all) nip(1)  scheme(s1color) 
*midas tp fp fn tn, id(studyid) bfor(dss) fordata  //forest plot

*------Table 3 (column 3)
metandi tp fp fn tn
*------ Figure 5 (bottom right)
metandiplot tp fp fn tn, graphregion(color(white)) legend(ring(0) col(1) bplacement(6))

/*========================================================================================*/
***Section 5.3 
/*========================================================================================*/
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/ascus.csv", delim(",")clear

*------Table 4 (column 2), Figure 6 (top left) 
metadta tp fp fn tn test, studyid(studyid) nomc ///
	comparative sumtable(all) dp(2) nofplot
	
*------Figure 7 (left)	
metadta tp fp fn tn test, studyid(studyid) nomc ///
	comparative sumtable(rr) dp(2) nosroc ///
	foptions(logscale outplot(rr) texts(2)	xlabel(0.3, 1, 2))
	
*------Table 4 (column 3), Figure 6 (top right)
encode test, gen(Test)
replace Test = Test - 1
midas tp fp fn tn, reg(Test) sroc(both) res(all) 

*------Figure 7 (right)	
midas tp fp fn tn, reg(Test) plot 
	

/*========================================================================================*/
***Section 6 
/*========================================================================================*/
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/hpvcyto-wide.csv", delim(",")clear

*------Figure 8
*Install command	
ssc install netplot	
netplot test2 test1, label arrows type(circle)	

*------Stratified analysis (Comparative studies stratified by index)
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/hpvcyto-long1.csv", delim(",")clear

bys index:metadta tp fp fn tn test, nomc studyid(studyid) ///
	comparative sumtable(all) ref(HC2,top)  ///
	foptions(outplot(rr) texts(1.75) ///
	xlabel(0.80, 0.90, 1, 1.11, 2) logscale ///
	 astext(70))

*------Figure 9
graph combine fplot2 fplot1, col(1)  xcommon graphregion(col(white))

*------cb-network
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/hpvcyto-wide.csv", delim(",")clear

*------Figure 10
metadta tp2 fp2 fn2 tn2 tp1 fp1 fn1 tn1 test1 test2, studyid(studyid) ///
	cbnetwork cov(,identity) progress sumtable(all) ///
	foptions(outplot(rr) texts(1.75) ///
	xlabel(0.80, 0.90, 1, 1.11, 2) logscale ///
	lcols(age location) astext(70))

*------ab-network
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/hpvcyto-long2.csv", delim(",")clear
//Summaries & fplot
metadta tp fp fn tn test, studyid(studyid) progress ///
	abnetwork cov(,identity) ref(HC2, top) sumtable(all)  ///
	foptions(pointopt(msize(1.5)) outplot(rr) texts(2) ///
	xlabel(0.9, 1, 1.8))

//SROC
metadta tp fp fn tn test, studyid(studyid) progress ///
	abnetwork cov(,identity) ref(HC2, top) nofplot
*------Figure 11
graph combine fplot sroc, col(2) graphregion(col(white))

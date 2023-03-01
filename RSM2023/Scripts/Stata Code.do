cap log close
cap log using "C:\DATA\WIV\Projects\GitHub\Metadta\RSM2023\Logs\stataout.txt", text replace
set more off


/*****************************************************************************
metadta
******************************************************************************/
*Install the command
ssc install metadta

***Section 5.1 
*------Table 2 (column 2), Figure 1, Figure 4 (top left)
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/telomerase.csv", delim(",")clear
metadta tp fp fn tn, studyid(study) sumtable(all)

***Section 5.2
*------Table 3 (column 2), Figure 5 (top left)
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/dementia.csv", delim(",")clear
metadta tp fp fn tn, studyid(study) sumtable(all)

***Section 5.3 
*------Table 4 (column 2), Figure 6 (top left), 
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/ascus.csv", delim(",")clear
metadta tp fp fn tn test, studyid(studyid) nomc ///
	comparative sumtable(all) dp(2) nofplot
	
*------Figure 7 (right)	
metadta tp fp fn tn test, studyid(studyid) nomc ///
	comparative sumtable(rr) dp(2) nosroc ///
	foptions(logscale outplot(rr) texts(2)	xlabel(0.3, 1, 2)) 	

***Section 6 
import delimited using "https://raw.githubusercontent.com/VNyaga/Metadta/master/RSM2023/Data/hpvcyto-wide.csv", delim(",")clear

*------Figure 8
*Install command	
ssc install netplot	
netplot test2 test1, label arrows type(circle)	
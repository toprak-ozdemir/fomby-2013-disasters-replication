************************************************************
*                    "Replace_outliers.do"                 *
*    Replace outliers with . , instead of dropping them    *
************************************************************

* Remove if an observation has any outlier  
gen temp_growth = growth
gen temp_agrgrowth = agrgrowth
gen temp_noagrowth = noagrowth

replace temp_growth = . if abs(growth) > 0.5 & growth ~= .
replace temp_growth = . if droughtcont > log(0.2) 
replace temp_growth = . if erqcont > log(0.2) 
replace temp_growth = . if floodcont > log(0.2) 
replace temp_growth = . if stormcont > log(0.2) 

replace temp_agrgrowth = . if abs(agrgrowth) > 0.5 & agrgrowth ~= .
replace temp_agrgrowth = . if droughtcont > log(0.2) 
replace temp_agrgrowth = . if erqcont > log(0.2) 
replace temp_agrgrowth = . if floodcont > log(0.2) 
replace temp_agrgrowth = . if stormcont > log(0.2) 

replace temp_noagrowth = . if abs(noagrowth) > 0.5 & noagrowth ~= .
replace temp_noagrowth = . if droughtcont > log(0.2) 
replace temp_noagrowth = . if erqcont > log(0.2) 
replace temp_noagrowth = . if floodcont > log(0.2) 
replace temp_noagrowth = . if stormcont > log(0.2) 

drop growth agrgrowth noagrowth 
rename temp_growth growth
rename temp_agrgrowth agrgrowth
rename temp_noagrowth noagrowth

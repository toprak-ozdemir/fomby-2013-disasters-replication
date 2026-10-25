********************************************************
*                   "continuous.do"                    *
*      Generate new continuous disaster variables      *
********************************************************

** Make 4 positive continuous disaster variables
** for drought, earthquake, flood, and storm (dcont, econt, fcont, scont) 

/* Drought */
gen dcont = droughtcont + 20
gen dcontpos = .
replace dcontpos = dcont if dcont > 0
sum dcontpos, detail
gen th75_d = r(p75)
* gen th90_d = r(p90)
* gen th95_d = r(p95)

gen dcontsev75 = 0 
replace dcontsev75 = dcont if dcont > th75_d
* gen dcontsev90 = 0
* replace dcontsev90 = dcont if dcont > th90_d
* gen dcontsev95 = 0
* replace dcontsev95 = dcont if dcont > th95_d
drop th75_d dcontpos

/* Earthquake */
gen econt = erqcont + 20
gen econtpos = .
replace econtpos = econt if econt > 0
sum econtpos, detail
* gen th75_e = r(p75)
gen th90_e = r(p90)
* gen th95_e = r(p95)

* gen econtsev75 = 0
* replace econtsev75 = econt if econt > th75_e
gen econtsev90 = 0
replace econtsev90 = econt if econt > th90_e
* gen econtsev95 = 0
* replace econtsev95 = econt if econt > th95_e
drop th90_e econtpos

/* Flood */
gen fcont = floodcont + 20
gen fcontpos = .
replace fcontpos = fcont if fcont > 0
sum fcontpos, detail
* gen th75_f = r(p75)
gen th90_f = r(p90)
* gen th95_f = r(p95)

* gen fcontsev75 = 0
* replace fcontsev75 = fcont if fcont > th75_f
gen fcontsev90 = 0
replace fcontsev90 = fcont if fcont > th90_f
* gen fcontsev95 = 0
* replace fcontsev95 = fcont if fcont > th95_f
drop th90_f fcontpos

/* Storm */
gen scont = stormcont + 20
gen scontpos = .
replace scontpos = scont if scont > 0
sum scontpos, detail
* gen th75_s = r(p75)
gen th90_s = r(p90)
* gen th95_s = r(p95)
* gen scontsev75 = 0
* replace scontsev75 = scont if scont > th75_s
gen scontsev90 = 0
replace scontsev90 = scont if scont > th90_s
* gen scontsev95 = 0
* replace scontsev95 = scont if scont > th95_s
drop th90_s scontpos


/* Generate new severe continuous variables */
/* Severe if intensity > 0.01 (same as discrete case) */

gen dcontsevnew = 0
gen econtsevnew = 0
gen fcontsevnew = 0
gen scontsevnew = 0
replace dcontsevnew = dcont if droughtsev > 0  
replace econtsevnew = econt if erqsev > 0
replace fcontsevnew = fcont if floodsev > 0
replace scontsevnew = scont if stormsev > 0

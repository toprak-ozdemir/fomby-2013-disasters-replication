*****************************************************
*                   "Data_prep.do"                  * 
*         Data preparation for Matlab codes         *     
*                Revised in July 2011               * 
*****************************************************

** NOTE:
** (1) SELECT the sample: developing countries or developed countries
**     by removing stars in line number 48-56. 

** (2) SELECT the endogenous growth variable: GDP growth, Agr. growth, or Non-agr. growth
**     by removing stars in line number 102-116. 

** (3) CHANGE the name of file at the bottom of this do file
** 	 in accordance with (1) and (2) above.


clear all
set more off
set memory 80m

/* Set current directory */
 

use Data.dta, clear
egen id = group(ccode)

/* Select 95 countries used in the paper */
run "cnt95.do"

/* Remove 8 countries which don't have tot data */ 
drop if ccode == "FJI" | ccode == "BDI" | ccode == "JAM"| ccode== "NER" | ///
        ccode == "NPL"| ccode == "SGP"| ccode == "SLB" | ccode == "SUR"

/* Remove (replace with .) outliers */
run "Replace_outliers.do"

/* Generate continuous intensity measures */
run "continuous.do"

/* Make "intensity measures with fatalities" positive */
gen dfatal = droughtfatal + 20
gen efatal = erqfatal + 20
gen ffatal = floodfatal + 20
gen sfatal = stormfatal + 20

*********************************************
** SELECT developing or developed 

/* if developing countries */
run "developing.do" 
drop oecdarab 

/* if developed countries */
* run "ArabOECD.do" 
* drop oecdarab 
*********************************************

/* Remove top 0.5% of i_cpi */
gen i_cpi05 = . 
gsort -i_cpi
sum i_cpi, detail
replace i_cpi05 = i_cpi if  _n > round(0.005*r(N))
drop i_cpi
rename i_cpi05 i_cpi

sort ccode year
xtset id year

/* Set up lagged variables */
run "Lag.do"

order ccode year id totgr growth agrgrowth noagrowth ///
	lt lg lag lng l2t l2g l2ag l2ng l3t l3g l3ag l3ng ///
 	drought earthquake flood storm droughtsev erqsev floodsev stormsev ///
 	ld le lf ls lds les lfs lss l2d l2e l2f l2s l2ds l2es l2fs l2ss ///
	l3d l3e l3f l3s l3ds l3es l3fs l3ss ///
	dcont econt fcont scont dcontsev75 econtsev90 fcontsev90 scontsev90 ///
 	dcontsevnew econtsevnew fcontsevnew scontsevnew ///
	ldcont lecont lfcont lscont ldconts75 leconts90 lfconts90 lsconts90 ///
	ldcontsn lecontsn lfcontsn lscontsn ///
      l2dcont l2econt l2fcont l2scont l2dconts75 l2econts90 l2fconts90 l2sconts90 ///
	l2dcontsn l2econtsn l2fcontsn l2scontsn ///
	l3dcont l3econt l3fcont l3scont l3dconts75 l3econts90 l3fconts90 l3sconts90 ///
	l3dcontsn l3econtsn l3fcontsn l3scontsn ///	
	dfatal efatal ffatal sfatal capform i_cpi openness govcon privcred world_growth ///
	ldfatal lefatal lffatal lsfatal lcapform li_cpi lopenness lgovcon lprivcred lworld_growth ///
	l2dfatal l2efatal l2ffatal l2sfatal l2capform l2i_cpi l2openness l2govcon l2privcred l2world_growth ///
	l3dfatal l3efatal l3ffatal l3sfatal l3capform l3i_cpi l3openness l3govcon l3privcred l3world_growth ///
	drought2 earthquake2 flood2 storm2 ld2 le2 lf2 ls2 l2d2 l2e2 l2f2 l2s2 l3d2 l3e2 l3f2 l3s2 
	

**** Drop missings of additional variables ****

** Additional 6 vars (w/ 4 growth vars) ** 
drop if g_capform == . | i_cpi ==. | g_openness ==. | g_govcon ==.| g_privcred == . | world_growth == .  
drop if lgcapform == . | li_cpi ==. | lgopenness ==. | lggovcon ==.| lgprivcred == . | lworld_growth == .  
drop if l2gcapform == . | l2i_cpi ==. | l2gopenness ==. | l2ggovcon ==.| l2gprivcred == . | l2world_growth == . 
drop if l3gcapform == . | l3i_cpi ==. | l3gopenness ==. | l3ggovcon ==.| l3gprivcred == . | l3world_growth == . 


************************************************
** SELECT endogenous growth rate in the model 

* if GDP growth, remove stars 
drop if growth == . | lg ==. | l2g ==. | l3g ==.| totgr == . |lt ==.| l2t ==.| l3t ==. 
drop agrgrowth noagrowth lag lng l2ag l2ng l3ag l3ng 

* if agr growth, remove stars 
* drop if agrgrowth == .| lag ==. | l2ag ==. | l3ag ==.| totgr == . |lt ==.| l2t ==.| l3t ==. 
* drop growth noagrowth lg lng l2g l2ng l3g l3ng 

* if non-agr growth, remove stars 
* drop if noagrowth == . | lng ==. | l2ng ==. | l3ng ==.| totgr == . |lt ==.| l2t ==.| l3t ==. 
* drop growth agrgrowth lg lag l2g l2ag l3g l3ag 
**************************************************

drop if id ~= id[_n+1] & id ~= id[_n-1]	
egen newid = group(id)


**** Drop unnecessary variables ****
#delimit ;
drop  cnt95 ccode year id cname pop
      rgdp lrgdp tot ltot lwldgdp 
	agrgdpcondoll indgdpcondoll servgdpcondoll 
	agrpcgdp indpcgdp servpcgdp 
	droughtcont erqcont floodcont stormcont
	droughtfatal erqfatal floodfatal stormfatal lcpi  
;
#delimit cr


**** Export the series in csv format ****
** Developing 
outsheet using "dev_grow.csv", comma nonames replace

* outsheet using "dev_agr.csv", comma nonames replace
* outsheet using "dev_non.csv", comma nonames replace

** Developed
* outsheet using "rich_grow.csv", comma nonames replace
* outsheet using "rich_agr.csv", comma nonames replace
* outsheet using "rich_non.csv", comma nonames replace










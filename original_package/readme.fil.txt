Thomas B. Fomby, Yuki Ikeda, and Norman V. Loayza, "The Growth Aftermath 
of Natural Disasters", Journal of Applied Econometrics, Vol. 28, No. 3, 
2013, pp. 412-434.

I. The sources from which the data were taken

The data on natural disasters were obtained from the Emergency Disasters 
Database (EM-DAT), maintained by the Centre for Research on the 
Epidemiology of Disasters (CRED), and the data on macroeconomic 
variables were taken from the World Development Indicators (WDI) of the 
World Bank.

II. The number of observations

We use annual panel data for 84 countries over the period 1960-2007.  
The panel is unbalanced. We split the panel into two groups: developing 
countries and developed countries, classifying 60 countries into the 
first group and the other 24 counties into the second group. Table A-I 
in the Supplement of the paper (available online at the JAE homepage) 
gives the list of countries of these groups.  In our data for developing 
countries, we have 2097 observations with the growth rate of real GDP 
per capita, 1878 observations with the growth rate of real agricultural 
value added per capita, and 1863 observations with the growth rate of 
real non-agricultural value added per capita. As for developed 
countries, we have 853 observations with the growth rate of real GDP per 
capita, 740 observations with the growth rate of real agricultural value 
added per capita, and 740 observations with the growth rate of real 
non-agricultural value added per capita.  The size of the dataset in dta 
format for STATA, fil-data.dta, is 2.25 MB, and the size of the dataset 
in txt format, fil-data.txt, is 2.16 MB. The list of the variables in 
the dataset is the following.

III. The list of variables
[Varible name] [Description]

year	       	Year
ccode	     	Country code
cname	 	Country name
pop	 	Population
rgdp	 	Real GDP (constant LCU)
lrgdp	 	Log of real GDP per capita
growth	 	Growth rate of real GDP per capita
agrgdpcondoll	Agriculture, value added (constant 2000 US$)
agrpcgdp	Agriculture, per capita Value Added (constant 2000 US$)
agrgrowth	Growth rate of agricultural value added per capita
indgdpcondoll	Industry, value added (constant 2000 US$)
servgdpcondoll	Services, etc., value added (constant 2000 US$)
indpcgdp	Industry, per capita Value Added (constant 2000 US$)
servpcgdp	Services, per capita Value Added (constant 2000 US$)
noagrowth	Growth rate of non-agricultural value added per capita
tot		Terms of trade
ltot		Terms of trade (in logs)
totgr		Terms of trade growth
capform		Gross fixed capital formation/GDP (in logs)
g_capform	Growth rate of gross fixed capital formation/GDP
lcpi		Consumer price index (in logs)
i_cpi		1 + growth rate of CPI (in logs)
openness	(Exports+imports)/GDP (in logs)
g_openness	Growth rate of (Exports+imports)/GDP
govcon		General government final consumption expenditure/GDP (in logs)
g_govcon	Growth rate of general government final consumption
		expenditure/GDP
privcred	Domestic credit to private sector/GDP (in logs)
g_privcred	Growth rate of domestic credit to private sector/GDP
lwldgdp		World real GDP per capita (in logs)
world_growth	Growth rate of World real GDP per capita
drought		Drought, \sum_i 1(intensity > 0.0001){intensity_i}
earthquake	Earthquake, sum_i 1(intensity > 0.0001){intensity_i}
flood		Flood, sum_i 1(intensity > 0.0001){intensity_i}
storm		Storm, sum_i 1(intensity > 0.0001){intensity_i}
drought2	Drought2, sum_i 1(intensity2 > 0.0001){intensity2_i}
earthquake2	Earthquake2, sum_i 1(intensity2 > 0.0001){intensity2_i}
flood2		Flood2, sum_i 1(intensity2 > 0.0001){intensity2_i}
storm2		Storm2, sum_i 1(intensity2 > 0.0001){intensity2_i}
droughtsev	Severe drought, sum_i 1(intensity > 0.01){intensity_i}
erqsev		Severe earthquake, sum_i 1(intensity > 0.01){intensity_i}
floodsev	Severe flood, sum_i 1(intensity > 0.01){intensity_i}
stormsev	Severe storm, sum_i 1(intensity > 0.01){intensity_i}
droughtcont	Drought (continuous), ln(\sum_i intensity_i)
erqcont		Earthquake (continuous), ln(\sum_i intensity_i)
floodcont	Flood (continuous), ln(\sum_i intensity_i)
stormcont	Storm (continuous), ln(\sum_i intensity_i)
droughtfatal	Drought fatalities, ln(\sum_i intensity_i)
erqfatal	Earthquake fatalities, ln(\sum_i intensity_i)
floodfatal	Flood fatalities, ln(\sum_i intensity_i)
stormfatal	Storm fatalities, ln(\sum_i intensity_i)


IV. Data organization in the data files

[1] STATA files

The main Stata data file, fil-data.dta, is zipped in 
stata-data-file.zip. The same data is also stored in the ASCII file 
fil-data.txt, which is zipped in stata-text-files.zip, along with the 
following files:

- data_prep.do (the main do file)
- ArabOECD.do
- cnt95.do
- continuous.do
- developing.do
- lag.do
- replace_outliers.do

All do files except the main do file, data_prep.do, are sub-files used 
in it.

Unix/Linux users should use "unzip -a" for stata-text-files.zip but 
*not* for stata-data-file.zip.

[2] Spreadsheets

All spreadsheets are stored in the file spreadsheets-xls.zip. There are 
also ASCII text versions in spreadsheets-text.zip. Unix/Linux users 
should use "unzip -a" for spreadsheets-text.zip but not for 
spreadsheets-xls.zip.

- dev_grow.xls / dev_grow.txt
- dev_agr.xls / dev_agr.txt
- dev_non.xls / dev_non.txt
- rich_grow.xls / rich_grow.txt
- rich_agr.xls / rich_agr.txt
- rich_non.xls / rich_non.txt

[3] MATLAB files

All m files for MATLAB are stored in the file matlab-files.zip. These 
are ASCII files in DOS format, so Unix/Linux users should use "unzip 
-a". There are seven Matlab files:

- apply_BSBC.m (the main m file for the BSBC estimator)
  [Table III, IV in the paper and Table B-II, III, IV in the supplement]
- apply_severe_BSBC.m (the main m file for the BSBC estimator) [Table V in
  the paper]
- apply_LSDV.m (the main m file for the BSBC estimator) [Table B-I in the
  supplement] 
- apply_severe_LSDV.m (the main m file for the BSBC estimator)
- btst.m
- meanpsi.m
- psimat.m

The last 3 files are sub-files used in the first 4 files.


Yuki Ikeda
yi22 [AT] georgetown.edu


****************************************************
*                   "ArabOECD.do"                  *
*          Select Arab and OECD countries          *
****************************************************

gen oecdarab=0

replace oecdarab=1 if ccode=="AUS"|ccode=="AUT"|ccode=="BEL"| ///
ccode=="CAN"|ccode=="CZE"|ccode=="DNK"| ///
ccode=="FIN"|ccode=="FRA"|ccode=="DEU"| ///
ccode=="GRC"|ccode=="HUN"|ccode=="ISL"| ///
ccode=="IRL"|ccode=="ITA"|ccode=="JPN"| ///
ccode=="LUX"|ccode=="NLD"|ccode=="NZL"| ///
ccode=="NOR"|ccode=="POL"|ccode=="PRT"| ///
ccode=="SVK"|ccode=="ESP"|ccode=="SWE"| ///
ccode=="CHE"|ccode=="TUR"|ccode=="GBR"| ///
ccode=="USA"|ccode=="ISR"|ccode=="SAU"|ccode=="ARE"|ccode=="KWT"

drop if oecdarab == 0


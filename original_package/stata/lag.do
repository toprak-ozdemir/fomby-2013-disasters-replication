****************************************
*               "Lag.do"               *
*       Generate lagged variables      * 
*        Revised in October 2010       *
****************************************

/* First lags */
gen ld = L.drought
gen le = L.earthquake
gen lf = L.flood
gen ls = L.storm

gen ld2 = L.drought2
gen le2 = L.earthquake2
gen lf2 = L.flood2
gen ls2 = L.storm2

gen lds = L.droughtsev
gen les = L.erqsev
gen lfs = L.floodsev
gen lss = L.stormsev

gen lt = L.totgr
gen lg = L.growth
gen lag = L.agrgrowth
gen lng = L.noagrowth

gen ldcont = L.dcont
gen lecont = L.econt
gen lfcont = L.fcont
gen lscont = L.scont

gen ldconts75 = L.dcontsev75
gen leconts90 = L.econtsev90
gen lfconts90 = L.fcontsev90
gen lsconts90 = L.scontsev90

gen ldcontsn = L.dcontsevnew
gen lecontsn = L.econtsevnew
gen lfcontsn = L.fcontsevnew
gen lscontsn = L.scontsevnew

gen ldfatal = L.dfatal
gen lefatal = L.efatal
gen lffatal = L.ffatal
gen lsfatal = L.sfatal

gen lcapform = L.capform
gen li_cpi = L.i_cpi
* gen li_gdpdef = L.i_gdpdef
gen lopenness = L.openness
gen lgovcon = L.govcon
gen lprivcred = L.privcred
gen lworld_growth = L.world_growth

gen lgcapform = L.g_capform
gen lgopenness = L.g_openness
gen lggovcon = L.g_govcon
gen lgprivcred = L.g_privcred


/* Second lags */
gen l2d = L2.drought
gen l2e = L2.earthquake
gen l2f = L2.flood
gen l2s = L2.storm

gen l2d2 = L2.drought2
gen l2e2 = L2.earthquake2
gen l2f2 = L2.flood2
gen l2s2 = L2.storm2

gen l2ds = L2.droughtsev
gen l2es = L2.erqsev
gen l2fs = L2.floodsev
gen l2ss = L2.stormsev

gen l2t = L2.totgr
gen l2g = L2.growth
gen l2ag = L2.agrgrowth
gen l2ng = L2.noagrowth

gen l2dcont = L2.dcont
gen l2econt = L2.econt
gen l2fcont = L2.fcont
gen l2scont = L2.scont

gen l2dconts75 = L2.dcontsev75
gen l2econts90 = L2.econtsev90
gen l2fconts90 = L2.fcontsev90
gen l2sconts90 = L2.scontsev90

gen l2dcontsn = L2.dcontsevnew
gen l2econtsn = L2.econtsevnew
gen l2fcontsn = L2.fcontsevnew
gen l2scontsn = L2.scontsevnew

gen l2dfatal = L2.dfatal
gen l2efatal = L2.efatal
gen l2ffatal = L2.ffatal
gen l2sfatal = L2.sfatal

gen l2capform = L2.capform
gen l2i_cpi = L2.i_cpi
* gen l2i_gdpdef = L2.i_gdpdef
gen l2openness = L2.openness
gen l2govcon = L2.govcon
gen l2privcred = L2.privcred
gen l2world_growth = L2.world_growth

gen l2gcapform = L2.g_capform
gen l2gopenness = L2.g_openness
gen l2ggovcon = L2.g_govcon
gen l2gprivcred = L2.g_privcred



/* Third lags */
gen l3d = L3.drought
gen l3e = L3.earthquake
gen l3f = L3.flood
gen l3s = L3.storm

gen l3d2 = L3.drought2
gen l3e2 = L3.earthquake2
gen l3f2 = L3.flood2
gen l3s2 = L3.storm2

gen l3ds = L3.droughtsev
gen l3es = L3.erqsev
gen l3fs = L3.floodsev
gen l3ss = L3.stormsev

gen l3t = L3.totgr
gen l3g = L3.growth
gen l3ag = L3.agrgrowth
gen l3ng = L3.noagrowth

gen l3dcont = L3.dcont
gen l3econt = L3.econt
gen l3fcont = L3.fcont
gen l3scont = L3.scont

gen l3dconts75 = L3.dcontsev75
gen l3econts90 = L3.econtsev90
gen l3fconts90 = L3.fcontsev90
gen l3sconts90 = L3.scontsev90

gen l3dcontsn = L3.dcontsevnew
gen l3econtsn = L3.econtsevnew
gen l3fcontsn = L3.fcontsevnew
gen l3scontsn = L3.scontsevnew

gen l3dfatal = L3.dfatal
gen l3efatal = L3.efatal
gen l3ffatal = L3.ffatal
gen l3sfatal = L3.sfatal

gen l3capform = L3.capform
gen l3i_cpi = L3.i_cpi
* gen l3i_gdpdef = L3.i_gdpdef
gen l3openness = L3.openness
gen l3govcon = L3.govcon
gen l3privcred = L3.privcred
gen l3world_growth = L3.world_growth

gen l3gcapform = L3.g_capform
gen l3gopenness = L3.g_openness
gen l3ggovcon = L3.g_govcon
gen l3gprivcred = L3.g_privcred

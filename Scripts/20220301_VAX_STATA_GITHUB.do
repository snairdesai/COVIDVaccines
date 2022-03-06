u data/finalSampleWW, clear

global xVars1 npiStringency recoveredTotal capitaGDP aged_65_older
global xVars2 npiStringency recoveredTotal                         cdum*
global xVars3 npiStringency recoveredTotal
global xVars4                              capitaGDP aged_65_older             
global xVars5                                                      cdum*                                            
global xVars6 npiStringency recoveredTotal capitaGDP aged_65_older       VAXXDEPRx*

foreach k in infectedNew	{
								bys cnum: gen   L1`k' = l1.`k'
								bys cnum: replace `k' = `k' - L1`k'
							} 

gen laggedDepVar = .
foreach k in $y2vars {

                       replace `k' = `k'*10000

                       replace laggedDepVar = l.`k' 


                       eststo a_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1, vce(cluster iso_code continent)
                       eststo b_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2, vce(cluster iso_code continent)                        
                       eststo c_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3, vce(cluster iso_code continent)                       
                       eststo d_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4, vce(cluster iso_code continent)                      
                       eststo e_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5, vce(cluster iso_code continent)
                       eststo f_`k': reg `k'              l(2).Vaccinated $xVars5, vce(cluster iso_code continent)
                       eststo g_`k': reg `k' laggedDepVar l(2).vaxxDepred $xVars1, vce(cluster iso_code continent)                      
                       eststo h_`k': reg `k' laggedDepVar                 $xVars6, vce(cluster iso_code continent)
                       

   forval j = 1/3    {
                        display "`k'"
                        eststo a`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo b`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2 if vaxGrp==`j', vce(cluster iso_code continent)                       
                        eststo c`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo d`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo e`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo f`j'`k': reg `k'              l(2).Vaccinated $xVars5 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo g`j'`k': reg `k' laggedDepVar l(2).vaxxDepred $xVars1 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo h`j'`k': reg `k' laggedDepVar                 $xVars6 if vaxGrp==`j', vce(cluster iso_code continent)
                      }
                      
    foreach l in 0 1  {
                        display "oecd == `l'"
                        eststo _a`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1 if oecd==`l', vce(cluster iso_code continent)
                        eststo _b`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2 if oecd==`l', vce(cluster iso_code continent)
                        eststo _c`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3 if oecd==`l', vce(cluster iso_code continent)
                        eststo _d`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4 if oecd==`l', vce(cluster iso_code continent)
                        eststo _e`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5 if oecd==`l', vce(cluster iso_code continent)
                        eststo _f`l'`k': reg `k'              l(2).Vaccinated $xVars5 if oecd==`l', vce(cluster iso_code continent)
                      }                     
                     
    foreach m in 0 1  {
                        display "EME == `m'"
                        eststo _`m'a`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1 if EME==`m', vce(cluster iso_code continent)
                        eststo _`m'b`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2 if EME==`m', vce(cluster iso_code continent)
                        eststo _`m'c`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3 if EME==`m', vce(cluster iso_code continent)
                        eststo _`m'd`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4 if EME==`m', vce(cluster iso_code continent)
                        eststo _`m'e`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5 if EME==`m', vce(cluster iso_code continent)
                        eststo _`m'f`k': reg `k'              l(2).Vaccinated $xVars5 if EME==`m', vce(cluster iso_code continent)
                      }                     
                     }

*** Log transformation

global toLog deathsNewLIcu2W recoveredTotal capitaGDP aged_65_older npiStringency Vaccinated vaxxDepred infectedNew vaxTec11 vaxTec12 vaxTec13 vaxTec14
 
drop VACCINESxAlpha VACCINESxBeta VACCINESxGamma VACCINESxDelta VACCINESxOthers VAXXDEPRxAlpha VAXXDEPRxBeta VAXXDEPRxGamma VAXXDEPRxDelta VAXXDEPRxOthers

foreach k in $toLog {                      
                      replace `k'= ln(`k')
                    }
					
gl VAXX VACCINES VAXXDEPR
gl DOMV Alpha Beta Gamma Delta Others
foreach k in $VAXX	{
	foreach j in $DOMV	{
							g `k'x`j' = ln(`k')*`j'
						}
					}
                    
foreach k in infectedNew	{     
								replace L1`k' = l1.`k'
								replace `k' = `k' - L1`k'
							}                    
                    
foreach k in $y2vars {

                       replace `k' = `k'

                       replace laggedDepVar = l.`k' 

                       eststo lna_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1, vce(cluster iso_code continent)
                       eststo lnb_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2, vce(cluster iso_code continent)
                       eststo lnc_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3, vce(cluster iso_code continent)
                       eststo lnd_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4, vce(cluster iso_code continent)
                       eststo lne_`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5, vce(cluster iso_code continent)
                       eststo lnf_`k': reg `k'              l(2).Vaccinated $xVars5, vce(cluster iso_code continent)
                       eststo lng_`k': reg `k' laggedDepVar l(2).vaxxDepred $xVars1, vce(cluster iso_code continent)
                       eststo lnh_`k': reg `k' laggedDepVar                 $xVars6, vce(cluster iso_code continent)
                       

   forval j = 1/3    {
                        display "`k'"
                        eststo lna`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo lnb`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo lnc`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo lnd`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo lne`j'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo lnf`j'`k': reg `k'              l(2).Vaccinated $xVars5 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo lng`j'`k': reg `k' laggedDepVar l(2).vaxxDepred $xVars1 if vaxGrp==`j', vce(cluster iso_code continent)
                        eststo lnh`j'`k': reg `k' laggedDepVar                 $xVars6 if vaxGrp==`j', vce(cluster iso_code continent)
 

                      }
                      
    foreach l in 0 1  {
                        display "oecd == `l'"
                        eststo ln_a`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1 if oecd==`l', vce(cluster iso_code continent)
                        eststo ln_b`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2 if oecd==`l', vce(cluster iso_code continent)
                        eststo ln_c`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3 if oecd==`l', vce(cluster iso_code continent)
                        eststo ln_d`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4 if oecd==`l', vce(cluster iso_code continent)
                        eststo ln_e`l'`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5 if oecd==`l', vce(cluster iso_code continent)
                        eststo ln_f`l'`k': reg `k'              l(2).Vaccinated $xVars5 if oecd==`l', vce(cluster iso_code continent)                     
                        }                     
                     
    foreach m in 0 1  {
                        display "EME == `m'"
                        eststo ln_`m'a`k': reg `k' laggedDepVar l(2).Vaccinated $xVars1 if EME==`m', vce(cluster iso_code continent)
                        eststo ln_`m'b`k': reg `k' laggedDepVar l(2).Vaccinated $xVars2 if EME==`m', vce(cluster iso_code continent)
                        eststo ln_`m'c`k': reg `k' laggedDepVar l(2).Vaccinated $xVars3 if EME==`m', vce(cluster iso_code continent)
                        eststo ln_`m'd`k': reg `k' laggedDepVar l(2).Vaccinated $xVars4 if EME==`m', vce(cluster iso_code continent)
                        eststo ln_`m'e`k': reg `k' laggedDepVar l(2).Vaccinated $xVars5 if EME==`m', vce(cluster iso_code continent)
                        eststo ln_`m'f`k': reg `k'              l(2).Vaccinated $xVars5 if EME==`m', vce(cluster iso_code continent)               
                        }                     
                     }

**# Tables 1, 1a, 1b
estout a_d* a3d* a1d* lna_d*,                           drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog")                                           label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 1: Impact of vaccines on M/LI)

estout a_d* a3d* a1d* lna_d* using Table1.rtf, replace drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog") label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 1: Impact of vaccines on M/LI)

estout g_d* g3d* g1d* lng_d*,                           drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog")                                           label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 1a: Impact of vaccines, adjusted for depreciation, on M/LI)

estout g_d* g3d* g1d* lng_d* using Table1a.rtf, replace drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog") label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 1a: Impact of vaccines, adjusted for depreciation, on M/LI)

estout h_d* h3d* h1d* lnh_d*,                           drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog")                                           label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 1b: Impact of depreciated vaccines and dominant variants on M/LI)

estout h_d* h3d* h1d* lnh_d* using Table1b.rtf, replace drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog") label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(25) title(Table 1b: Impact of depreciated vaccines and dominant variants on M/LI)

**# Tables 2, 2a, 2b
estout a_i* a3i* a1i* lna_i*,                          drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog")                                           label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 2: Impact of vaccines on new infections)

estout a_i* a3i* a1i* lna_i* using Table2.rtf, replace drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog") label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 2: Impact of vaccines on new infections)

estout g_i* g3i* g1i* lng_i*,                          drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog")                                           label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 2a: Impact of vaccines, adjusted for depreciation, on new infections)

estout g_i* g3i* g1i* lng_i* using Table2a.rtf, replace drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog") label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 2a: Impact of vaccines, adjusted for depreciation, on new infections)

estout h_i* h3i* h1i* lnh_i*,                          drop()  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog")                                           label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(20) title(Table 2b: Impact of depreciated vaccines and dominant variants on new infections)

estout h_i* h3i* h1i* lnh_i* using Table2b.rtf, replace  cells(b(star fmt(%9.2f)) se(par)) stats(F r2_a N, fmt(%9.1f %9.2f %9.0g) labels(F-stats adj-R-squared Obs)) legend mlabel("Full" "HighVax" "LowVax" "FullLog") label collabels(none) varlabels(_cons Constant) modelwidth(6) varwidth(25) title(Table 2b: Impact of depreciated vaccines and dominant variants on new infections) drop()

drop NOMISS*
gl OBSCHK1 deathsNewLIcu2W recoveredTotal capitaGDP aged_65_older npiStringency Vaccinated vaxxDepred infectedNew

mark NOMISS1
markout NOMISS1 $OBSCHK1
egen NOMISS1C = group(cnum) if NOMISS1==1
su NOMISS1C
bys vaxGrp: su NOMISS1C

gl OBSCHK2 deathsNewLIcu2W recoveredTotal capitaGDP aged_65_older npiStringency Vaccinated vaxxDepred infectedNew $DOMV VACCINESx*

mark NOMISS2
markout NOMISS2 $OBSCHK2
egen NOMISS2C = group(cnum) if NOMISS2==1
su NOMISS2C
bys vaxGrp: su NOMISS2C


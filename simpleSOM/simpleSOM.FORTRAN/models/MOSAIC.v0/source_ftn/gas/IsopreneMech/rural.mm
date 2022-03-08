*Reactions* 28
{1} ISOP + OH       = ISOPO2                     : 2.54e-11*exp(410/T)        ;
{2} ISOP + O3       = 0.23 MACR + 0.2 MVK + 
                      0.92 FORM + 0.05 MVOG + 
                      0.06 MVKO2 + 0.42 CO + 
                      0.09 HCOOH + 0.16 C2O3 + 
                      0.1 PAR + 0.1 OLE + 
                      0.06 AACD + 0.27 OH + 
                      0.26 HO2 + 0.06 H2O2 + 
                      0.15 XO2 + 0.003 SECOZO    : 7.86e-15*exp(-1913/T)      ;
{3} MACR + OH       = 0.55 MACRO2 + 0.45 MACO3   : 1.86e-11*exp(175/T)        ;
{4} MACR + O3       = 0.67 MGLY + 0.34 OH + 
                      0.53 CO + 0.32 HO2 + 
                      0.28 C2O3 + 0.57 FORM +
                      0.19 H2O2 + 0.12 HCOOH + 
                      0.05 ALD2 + 0.002 SECOZO   : 4.4e-15*exp(-2500/T)       ;
{5} MVK + OH        = MVKO2                      : 4.13e-12*exp(452/T)        ;
{6} MVK + O3        = 0.67 MGLY + 0.17 OH + 
                      0.53 CO + 0.06 HO2 + 
                      0.28 C2O3 + 0.57 FORM + 
                      0.19 H2O2 + 0.12 HCOOH + 
                      0.05 ALD2 + 0.002 SECOZO   : 7.51e-16*exp(-1521/T);
{7} GLY + OH        = HO2 + 2 CO                 : 1.1e-11; 
{8} ISOPO2          = OH + MVK + FORM            : 0.41*2.38e12*exp(-10770/T) ;
{9} ISOPO2          = OH + MACR + FORM           : 0.23*1.27e12*exp(-10570/T) ;
{10} ISOPO2         = OH + 2 HO2 + 0.5 HYAT + 
                      0.5 MGLY + 0.5 GLYALD + 
                      FORM                       : 5.2e8*exp(-	7700/T)     ;
{11} MACRO2         = OH + HYAT + CO             : 0.1
{12} ISOPO2 + HO2   = ISOPOOH                    : kRO2HO2                    ;
{13} ISOPOOH + OH   = 0.4 IEPOX + 0.7 OH + 
                      0.3 ISOPO2 + 0.3 ALDX      : 1.9e-11*exp(390/T)         ;
{14} ISOPOOH + hv   = OH + 0.3542 MVK + 
                      0.308 MACR + 0.695 FORM + 
                      0.693 HO2 + 0.1078 MVKOH + 
                      0.077 C2O3 + 0.138 C5CARB + 
                      0.069 ALDX + 0.23 HO2      : photo_ISOPOOH              ;

{15} IEPOX + OH     = IEPOXO2                    : 9.1E-12                    ;
{16} IEPOXO2 + HO2  = 0.725 HYAT + 0.275 GLY + 
                      0.275 GLYALD + 
                      0.275 MGLY +
                      0.375 FORM + 0.074 HCOOH + 
                      1.125 OH + 0.825 HO2 + 
                      0.251 CO                   : kRO2HO2                    ;
{17} C5CarbO2 + HO2 = 0.41 PACD + 0.15 AACD + 
                      0.15 O3 + 0.14 HYAT + 
                      0.1452 MGLY + 
                      0.2728 GLYALD + 
                      0.25 ALDX + 
                      0.44 HO2 + 0.44 OH         : kRO2HO2                    ;
{18} MACRO2 + HO2   = 0.925 ROOH + 0.075 HYAT + 
                      0.075 CO + 0.075 HO2 + 
                      0.075 OH                   : 0.88*kRO2HO2               ;
{19} MACO3 + HO2    = MPAA                       : kAPHO2                     ;
{20} MPAA + OH      = MAE + OH                   : 1.9E-11*exp(390/T)         ;
{21} MVKO2 + HO2    = 0.89 ROOH + 0.11 MGLY + 
                      0.11 FORM + 0.11 HO2 + 
                      0.11 OH                    : 0.88*kRO2HO2               ;
{22} C5CARB + OH    = C5CarbO2                   : 1.35e-11                   ;
{23} GLYALD + OH    = 0.1 ALDX + 0.2 HCOOH + 
                      0.6 FORM + 0.5 CO + 
                      0.25 OH + 0.75 HO2         : 8.00e-12                   ;
{24} HYAT + OH      = 0.65 ALDX + 0.9 HO2 + 
                      0.1 OH + 0.1 HCOOH + 
                      0.1 AACD                   : 5.98e-12                   ;
{25} ISOPO2         = 0.2417 MVK + 0.21 MACR + 
                      0.4725 FORM + 0.4725 HO2 + 
                      0.075 MVKOH + 
                      0.0525 C2O3 + 
                      0.09315 C5Carb + 
                      0.0075 GLYALD + 
                      0.0075 MGLY + 0.0093 HYAT + 
                      0.0155 FORM + 0.0466 ALDX + 
                      0.155 HO2 + 0.125 ALDX     : RO2*1.8e-12                ;
{26} C5CarbO2       = 0.24 HYAT + 0.2475 MGLY + 
                      0.465 GLYALD + 
                      0.4425 ALDX + 
                      0.75 HO2 + 0.125 AACD + 
                      0.125 ALDX                 : RO2*1.0e-11                ;
{27} MACRO2         = 0.325 HYAT + 0.325 MGLY + 
                      0.325 FORM + 0.325 CO + 
                      0.65 HO2 + 0.175 AACD + 
                      0.175 ALDX                 : RO2*2.5e-13                ;
{28} MACO3          = 0.4225 CO + 1.3 FORM + 
                      0.175 AACD + 0.175 ALDX    : RO2*1.0e-11                ;
{29} MVKO2          = 0.504 MGLY + 0.504 FORM + 
                      0.196 GLYALD + 0.196 ALDX + 
                      0.7 HO2 + 0.15*MVKOH + 
                      0.15 ALDX                  : RO2*2.5e-13                ;



















{1} ISOP + OH       = ISOPP + .08 XO2       : ARR(2.6e-11, 409.)              ;
{2} ISOP + O3       = .6 HCHO + .65 ISOPRD +
                      .27 OH + .07 HO2 + 
                      .07 CO + .39 RCOOH + 
                      .15 ALD2 + .2 XO2 +
                      .2 C2O3               : ARR(1.2e-14, -2013.)            ;
{3} ISOP + NO3      = ISOPN                 : ARR(3.0e-12, -446.)             ;
{4} ISOPRD + hv     = .97 C2O3 + .33 HO2 + 
                      .33 CO + .7 CH3O2 +
                      .2 HCHO + .07 ALD2 +
                      .03 AONE              : photo_ISOPRD                    ;
{5} ISOPRD + OH     = .5 C2O3 + .5 ISOPO2 +
                      .2 XO2                : 3.3e-11                         ;
{6} ISOPRD + O3     = .27 OH + .1 HO2 + 
                      .11 C2O3 + .07 XO2 +
                      .05 CH3O2 + .16 CO + 
                      .15 HCHO + .02 ALD2 +
                      .09 AONE + .85 MGLY +
                      .46 RCOOH             : 7.0e-18                         ;
{7} ISOPRD + NO3    = .07 C2O3 + .07 HNO3 +
                      .64 CO + .28 HCHO + 
                      .93 ONIT + .28 ALD2 +
                      .93 HO2 + .93 XO2 +
                      1.86 PAR              : 1.0e-15                         ;
{8} ISOPP + NO      = .09 ONIT + .91 NO2 + 
                      .91 HO2 + .63 HCHO +
                      .91 ISOPRD + 0.18 PAR : 4.0e-12                         ;
{9} ISOPN + NO      = 1.2 NO2 + .8 ONIT + 
                      .8 ALD2 + .8 HO2 + 
                      .2 ISOPRD + 1.6 PAR   : 4.0e-12                         ;
{10} ISOPO2 + NO    = NO2 + HO2 + .59 CO +
                      .55 ALD2 + .25 HCHO +
                      .34 MGLY + .63 AONE   : 4.0e-12                         ;
{11} ISOPP + HO2    = ROOH                  : ARR(1.7e-13, 1300.)             ;
{12} ISOPN + HO2    = ONIT + 2 PAR          : ARR(1.7e-13, 1300.)             ;
{13} ISOPO2 + HO2   = ROOH                  : ARR(1.7e-13, 1300.)             ;
{14} ISOPP          = ISOPRD                : rk_param(jisopp)                ;
{15} ISOPN          = ALD2 + ONIT + 2 PAR   : rk_param(jisopn)                ;
{16} ISOPO2         = .5 ALD2 + .5 AONE     : rk_param(jisopo2)               ;
*end*

*Reactions* 1407
#
# photolysis reactions
#
{1} MACR = CH3C2H2O2 + CO + HO2 		: J(18) 	;
{2} MACR = MACO3 + HO2 				: J(19) 	;
{3} MVK = C3H6 + CO 				: J(23) 	;
{4} MVK = C2O3 + HCHO + CO + HO2 		: J(24) 	;
{5} NISOPOOH = NISOPO + OH 			: J(41) 	;
{6} NISOPNO3 = NISOPO + NO2 			: J(53)*2.0 	;
{7} ISOPCNO3 = ISOPCO + NO2 			: J(53) 	;
{8} NC4CHO = NOA + 2 CO + 2 HO2 		: J(18) 	;
{9} MGLY = C2O3 + CO + HO2 			: J(34) 	;
{10} HC4CCHO = C2O3 + HO2 + CO + 
                HOCH2CHO 			: J(18) 	;
{11} HC4CCHO = HC4CCO3 + HO2 			: J(19) 	;
{12} ISOPAOOH = ISOPAO + OH 			: J(41) 	;
{13} ISOPANO3 = ISOPAO + NO2 			: J(53) 	;
{14} HC4ACHO = ACETOL + 2 HO2 + 2 CO 	: J(18) 	;
{15} HC4ACHO = HC4ACO3 + HO2 			: J(19) 	;
{16} ISOPBOOH = ISOPBO + OH 			: J(41) 	;
{17} ISOPBNO3 = ISOPBO + NO2 			: J(55) 	;
{18} ISOPCOOH = ISOPCO + OH 			: J(41) 	;

{19} ISOPDOOH = ISOPDO + OH 			: J(41) 	;
{20} ISOPDNO3 = ISOPDO + NO2 			: J(54) 	;
{21} HCOC5 = C2O3 + HCHO + HOCH2CO3 	: J(24) 	;
{22} NOA = C2O3 + HCHO + NO2 			: J(57) 	;
{23} NOA = CH3COCH2O + NO2 			: J(56) 	;
{24} HOCH2CHO = HO2 + HCHO + HO2 + CO 	: J(15) 	;
{25} GLY = CO + CO + H2 			: J(31) 	;
{26} GLY = CO + CO + 2 HO2 			: J(33) 	;
{27} GLY = HCHO + CO 				: J(32) 	;
{28} MACO3H = CH3C2H2O2 + OH 			: J(41) 	;

{29} MACROOH = ACETOL + CO + HO2 + OH 	: J(17) 	;
{30} MACROOH = MACRO + OH 			: J(41) 	;
{31} MACROH = ACETOL + CO + 2 HO2 		: J(17) 	;
{32} MACROHOOH = IBUTALOH + CO + HO2  + 
                  OH 				: J(17) 	;
{33} MACROHOOH = MACROHO + OH 		: J(41) 	;
{34} MACROHNO3 = MACROHO + NO2 		: J(55) 	;
{35} MACROHNO3 = NOA + HO2 + CO + HO2 	: J(17) 	;
{36} C3MDIALOH = MGLY + HO2 + CO + HO2 	: J(17)*2 	;
{37} CH3CO3H = CH3O2 + OH 			: J(41) 	;
{38} HMVKAOOH = HMVKAO + OH 			: J(41) 	;
{39} CO2H3CHO = MGLY + CO + 2 HO2 		: J(15) 	;
{40} HO12CO3C4 = C2O3 + HOCH2CHO + HO2 	: J(22) 	;
{41} HMVKBOOH = HMVKBO + OH 			: J(41) 	;
{42} BIACETOH = C2O3 + HOCH2CO3 		: J(35) 	;
{43} MVKOOH = HCHO + OH + ACO3 		: J(24) 	;
{44} MVKOOH = MVKO + OH 			: J(41) 	;
{45} MVKOH = ALLYLOH + CO 			: J(23) 	;
{46} MVKOH = HCHO + HO2 + HOCH2CO3 + CO 	: J(24) 	;
{47} VGLYOX = HO2 + CO + ACO3 		: J(34) 	;
{48} ACETOL = C2O3 + HCHO + HO2 		: J(22) 	;
{49} NO3CH2CHO = HO2 + CO + HCHO + NO2 	: J(57) 	;
{50} NO3CH2CHO = NO2 + HCOCH2O 		: J(56) 	;
{51} MACRNO3 = ACETOL + NO2 + CO + HO2 	: J(17) 	;
{52} MVKNO3 = C2O3 + HOCH2CHO + NO2 	: J(56)+J(57) 	;
{53} INCOOH = INCO + OH 			: J(41) 	;
{54} NC4CO3H = NOA + CO + HO2 + OH 		: J(41) 	;
{55} C510OOH = C510O + OH 			: J(41) 	;
{56} IBUTALOH = AONE + 2 HO2 + CO 		: J(17) 	;
{57} PR1O2HNO3 = PRONO3AO + OH 		: J(41) 	;
{58} CHOPRNO3 = HO2 + CO + ALD2 + NO2 	: J(57) 	;
{59} CHOPRNO3 = PROPALO + NO2 		: J(56) 	;
{60} PR2O2HNO3 = PRONO3BO + OH 		: J(41) 	;
{61} HYPROPO2H = HYPROPO + OH 		: J(41) 	;
{62} IPROPOLO2H = IPROPOLO + OH 		: J(41) 	;
{63} CH3CHOHCHO = ALD2 + HO2 + CO + HO2 	: J(17) 	;
{64} CH3COCO2H = C2O3 + HO2 			: J(34) 	;
{65} CO23C3CHO = C2O3 + 2 CO + HO2 		: J(34) 	;
{66} CO23C3CHO = C2O3 + HCOCO3 		: J(35) 	;
{67} HOCH2COCHO = HOCH2CO3 + CO + HO2 	: J(34) 	;
{68} IEACHO = HMVKBO2 + CO + HO2 		: J(17) 	;
{69} INAOOH = INAO + OH 			: J(41) 	;
{70} C524OOH = C524O + OH 			: J(41) 	;
{71} C524NO3 + OH = C524O + NO2 		: J(54) 	;
{72} C524CO = HOCH2CO3 + HOCH2CO3 + HCHO 	: J(24) 	;
{73} HC4ACO3H = ACETOL + CO + HO2 + OH 	: J(41) 	;
{74} C58OOH = C58O + OH 			: J(41) 	;
{75} C58NO3 = HO2 + CO + MACRNO3 + HO2 	: J(15) 	;
{76} CONM2CHO = MGLY + NO2 + CO + HO2 	: J(56)+J(57)*2 	;
{77} INB1OOH = INB1O + OH 			: J(41) 	;
{78} INB1CO = ACETOL + NO2 + HOCH2CO3 	: J(56)+J(57) 	;
{79} INB2OOH = INB2O + OH 			: J(41) 	;
{80} IECCHO = MACRO2 + HO2 + CO 		: J(15) 	;
{81} HC4CCO3H = HOCH2CHO + C2O3 + OH 	: J(41) 	;
{82} C57OOH = C57O + OH 			: J(41) 	;
{83} CO2N3CHO = CO + HO2 + MGLY + NO2 	: J(15) 	;
{84} INDOOH = INDO + OH 			: J(41) 	;
{85} HOCH2CO3H = HCHO + HO2 + OH 		: J(41) 	;
{86} C59OOH = C59O + OH 			: J(41) 	;
{87} C59OOH = HOCH2CO3 + ACETOL + OH 	: J(22) 	;
{88} INCNCHO = MACRNB + CO + NO2 + HO2 	: J(15) 	;
{89} INCGLYOX = MACRNBCO3 + CO + HO2 	: J(34) 	;
{90} HCOCO2H = HO2 + HO2 + CO 		: J(34) 	;
{91} HCOCO3H = HO2 + CO + OH 			: J(41)+J(15) 	;
{92} CHOMOHCO3H = MGLY + OH + HO2 		: J(41)+J(17) 	;
{93} HCOCH2OOH = HCOCH2O + OH 		: J(41) 	;
{94} HCOCH2OOH = HO2 + CO + HCHO + OH 	: J(15) 	;
{95} CO2H3CO3H = C2O3 + HO2 + HCOCO3H 	: J(22) 	;
{96} CO2H3CO3H = MGLY + HO2 + OH 		: J(41) 	;
{97} ACO3H = HO2 + CO + HCHO + OH 		: J(41) 	;
{98} MVKOHAOOH = MVKOHAO + OH 		: J(41)+J(22) 	;
{99} H13CO2CHO = HOCH2CHO + CO + 2 HO2 	: J(15) 	;
{100} MVKOHAOH = HOCH2CO3 + HOCH2CHO + 
                 HO2 				: J(22) 	;
{101} MVKOHBOOH = HOCH2CHO + HOCH2CO3 + 
                  OH 				: J(22) 	;
{102} MVKOHBOOH = MVKOHBO + OH 		: J(41) 	;
{103} H14CO23C4 = HOCH2CO3 + HOCH2CO3 	: J(35) 	;
{104} ACR = 0.3 ACO3 + 0.4 ETH + 0.7 CO +	
            0.3 HCHO + 0.3 HO2		: J(18) 	;
{105} NO3CH2CO3H = HCHO + NO2 + OH 		: J(41) 	;
{106} INAHPCHO = HMVKANO3 + OH + CO + HO2 : J(17)+J(41) 	;
{107} INANCHO = HMVKANO3 + NO2 + CO + HO2 : J(17) 	;
{108} INANCO = ACETOL + NO2 + NO3CH2CO3 	: J(56)+J(57) 	;
{109} INAHCHO = HMVKANO3 + 2 HO2 + CO 	: J(17) 	;
{110} IEB1OOH = HO12CO3C4 + OH + CO + HO2 : J(17) 	;
{111} IEB1OOH = IEB1O + OH 			: J(41) 	;
{112} IEB2OOH = IEB2O + OH 			: J(41) 	;
{113} IEB2OOH = MACROH + OH + CO + HO2 	: J(15) 	;
{114} MACRNCO3H = ACETOL + NO2 + OH 	: J(41) 	;
{115} INB1HPCHO = MACRNO3 + OH + CO + HO2 : J(15) 	;
{116} INB1NACHO = MACRNO3 + NO2 + 
                  CO + HO2 			: J(15) 	;
{117} INB1NBCHO = MVKNO3 + NO2 + CO + HO2 : J(17) 	;
{118} INB1GLYOX = MACRNCO3 + CO + HO2 	: J(34) 	;
{119} C57NO3 = HO2 + CO + HO12CO3C4 + NO2 : J(17) 	;
{120} IEC1OOH = ACETOL + OH + HOCH2CO3 	: J(22) 	;
{121} IEC1OOH = IEC1O + OH 			: J(41) 	;
{122} INDHPCHO = CO + HO2 + MVKNO3 + OH 	: J(17)+J(41) 	;
{123} INDHCHO = CO + HO2 + MVKNO3 + HO2 	: J(17) 	;
{124} MACRNB = NOA + HO2 + CO + HO2 	: J(17) 	;
{125} IPRHOCO3H = AONE + HO2 + OH 		: J(41) 	;
{126} PRNO3CO3H = ALD2 + NO2 + OH 		: J(41) 	;
{127} IPROPOLPER = ALD2 + HO2 + OH 		: J(41) 	;
{128} HOCH2COCO2H = HOCH2CO3 + HO2 		: J(34) 	;
{129} H1CO23CHO = 2 CO + HOCH2CO3 + HO2 	: J(34)+J(35) 	;
{130} IEACO3H = HMVKBO2 + OH 			: J(41) 	;
{131} NC524OOH = NC524O + OH 			: J(41) 	;
{132} C525OOH = C525O + OH 			: J(41) 	;
{133} C58NO3CO3H = MACRNO3 + HO2 + OH 	: J(41) 	;
{134} CONM2CO2H = CO + HO2 + NO2 + 
                   CH3COCO2H 			: J(56)+J(57) 	;
{135} CONM2CO3H = CO + HO2 + NO2 + 
                   CH3COCO3H 			: J(56)+J(57) 	;
{136} CONM2CO3H = MGLY + NO2 + OH 		: J(41) 	;
{137} CONM2PAN = CO + HO2 + NO2 + 
                  CH3COPAN 			: J(56)+J(57) 	;

{138} C4M2ALOHNO3 = CO2H3CHO + NO2 + 
                     CO + HO2 		: J(17) 	;
{139} C4M2ALOHNO3 = CONM2CHO + HO2 + 
                     CO + HO2 		: J(15) 	;
{140} IEC2OOH = BIACETOH + OH + CO + HO2 	: J(17) 	;
{141} IEC2OOH = MGLY + OH + HOCH2CO3 	: J(22) 	;
{142} IECCO3H = MACRO2 + OH 			: J(41) 	;
{143} CO2N3CO3H = MGLY + NO2 + OH 		: J(41) 	;
{144} CO2N3PAN = MGLY + NO2 + NO3 		: J(41) 	;
{145} INCNCO3H = MACRNB + NO2 + OH 		: J(41) 	;
{146} MACRNBCO3H = NOA + HO2 + OH 		: J(41) 	;
{147} HYPERACET = C2O3 + HCHO + OH 		: J(22) 	;
{148} HYPERACET = CH3COCH2O + OH 		: J(41) 	;
{149} H13CO2CO3H = HOCH2COCHO + HO2 + OH 	: J(41) 	;
{150} HOCHOCOOH = CHOCOHCO + OH 		: J(41) 	;
{151} C32OH13CO = GLY  + 2 HO2 + CO 	: J(15)*2 	;
{152} OCCOHCOOH = OCCOHCO + OH 		: J(41) 	;
{153} C42AOH = HO2 + CO + HO2 + 
                NO3CH2CHO 			: J(15) 	;
{154} INAHPCO2H = HMVKANO3 + OH + HO2 	: J(41) 	;
{155} INAHPCO3H = HMVKANO3 + 2 OH 		: J(41)*2 	;
{156} INAHPPAN = HMVKANO3 + OH + NO3 	: J(41) 	;
{157} INANCO3H = HMVKANO3 + NO2 + OH 	: J(41) 	;
{158} INAHCO3H = HMVKANO3 + HO2 + OH 	: J(41) 	;
{159} HIEB1OOH = HIEB1O + OH 			: J(41) 	;
{160} HIEB1OOH = MVKOHAOH + CO + 
                  OH + HO2 			: J(17) 	;
{161} HIEB2OOH = HIEB2O + OH 			: J(41) 	;
{162} HIEB2OOH = HMACROH + CO + OH + HO2 	: J(15) 	;
{163} HPNC524CO = HMVKNO3 + CO + 
                   HO2 + OH 			: J(17) 	;
{164} DNC524CO = HMVKNO3 + CO + HO2 +
                  NO2 				: J(17) 	;
{165} HMVKNO3 = HOCH2CHO + NO2 + CO + 
                 CO + HO2 			: J(56)+J(57) 	;
{166} H13CO2C3 = HOCH2CO3 + HCHO + HO2 	: J(22) 	;
{167} HNC524CO = HMVKNO3 + CO + 2 HO2 	: J(17) 	;
{168} HMACO3H = HOCH2CO3 + HCHO + OH 	: J(41) 	;
{169} HMACROOH = H13CO2C3 + CO + HO2 + 
                  OH 				: J(17) 	;
{170} HMACROOH = HMACRO + OH 			: J(41) 	;
{171} HMACROH = H13CO2C3 + CO + 2 HO2 	: J(17) 	;
{172} MMALNACO2H = MGLY + NO2 + HCOCO2H 	: J(56)+J(57) 	;
{173} MMALNACO3H = CONM2CHO + HO2 + OH 	: J(41) 	;
{174} MMALNACO3H = MGLY + NO2 + HCOCO3H 	: J(56)+J(57) 	;
{175} MMALNAPAN = MGLY + NO2 + GLYPAN 	: J(56)+J(57) 	;
{176} CH3COCO3H = C2O3 + OH 			: J(41)+J(35) 	;
{177} INB1HPCO3H = MACRNO3 + 2 OH 		: J(41) 	;
{178} INB1NACO3H = MACRNO3 + NO2 + OH 	: J(41) 	;
{179} INB1NBCO3H = MVKNO3 + NO2 + OH 	: J(41) 	;
{180} C57NO3CO3H = HO12CO3C4 + NO2 + OH 	: J(41) 	;
{181} INDHPCO3H = MVKNO3 + 2 OH 		: J(41)*2 	;
{182} INDHPPAN = MVKNO3 + OH + NO3 		: J(41) 	;
{183} INDHCO3H = MVKNO3 + OH + HO2 		: J(41) 	;
{184} COHM2CO2H = HCOCO2H + CO + HO2 	: J(17) 	;
{185} COHM2CO3H = GLY  + HO2 + OH 		: J(41) 	;
{186} COHM2CO3H = HCOCO3H + CO + HO2 	: J(17) 	;
{187} COHM2PAN = GLYPAN + CO + HO2 		: J(17) 	;
{188} NMGLY = NO2 + HCHO + 2 CO + HO2 	: J(53) 	;
{189} NMGLY = NO3CH2CO3 + CO + HO2 		: J(34) 	;
{190} ETHO2HNO3 = ETHENO3O + OH 		: J(41) 	;
{191} HYETHO2H = HOCH2CH2O + OH 		: J(41) 	;
{192} INANCOCO2H = NO3CH2CO3 + 
                    CH3COCO2H + NO2 	: J(56)+J(57) 	;
{193} INANCOCO3H = NO2 + CO23C4NO3 + OH 	: J(41) 	;
{194} INANCOCO3H = NO3CH2CO3 + 
                    CH3COCO3H + NO2 	: J(56)+J(57) 	;
{195} CO23C4NO3 = BIACETO + NO2 		: J(56) 	;
{196} CO23C4NO3 = HCHO + C2O3 + NO2 + CO	: J(57) 	;
{197} CO23C4NO3 = NO3CH2CO3 + C2O3 		: J(35) 	;
{198} INANCOPAN = NO3CH2CO3 + CH3COPAN + 
                   NO2 				: J(56)+J(57) 	;
{199} HMVKNGLYOX = CO + CO + HOCH2CHO + 
                    NO2 + HO2 		: J(34) 	;
{200} MMALNBCO2H = CONM2CO2H + HO2 + 
                    CO + HO2 			: J(15) 	;
{201} MMALNBCO3H = CO2H3CHO + NO2 + OH 	: J(41) 	;
{202} MMALNBCO3H = CONM2CO3H + HO2 + 
                    CO + HO2 			: J(15) 	;
{203} MMALNBPAN + OH = CO2H3CHO + 
{204} MMALNBPAN = CONM2PAN + HO2 + 
                   CO + HO2 			: J(15) 	;
{205} C2OHOCOOH = C3DIOLO2 			: J(41) 	;
{206} HCOCOHCO3H = GLY  + HO2 + OH 		: J(41) 	;
{207} C3DIOLOOH = C3DIOLO + OH 		: J(41) 	;
# bimolecular reactions
{208} NO3 + ISOP = NISOPO2 			: 3.15D-12*EXP(-450/TEMP) 	;
{209} O3 + ISOP = 0.3 CH2OOE + 0.3 MACR +
                0.2 CH2OOE + 0.2 MVK +
                0.5 HCHO + 0.3 MACROOA +
                0.2 MVKOOA     		: 1.03D-14*EXP(-1995/TEMP) 	;
{210} OH + ISOP = 0.148 ISOPAO2 +
                0.444 ISOPBO2 +
                0.102 ISOPCO2 +
                0.306 ISOPDO2 		: 2.7D-11*EXP(390/TEMP) 	;
{211} NISOPO2 + HO2 = NISOPOOH 		: KRO2HO2*0.706 	;
{212} NISOPO2 + NO = 0.052 NISOPNO3 +
	              0.948 NISOPO + 
                    0.948 NO2			: KRO2NO 	;
{213} NISOPO2 + NO3 = NISOPO + NO2 		: KRO2NO3 	;
{214} CH2OOE = 0.22 CH2OO + 0.78 CO +
              0.27 HO2 + 0.27 OH		: KDEC 	;
{215} NO3 + MACR = MACO3 + HNO3 		: 3.4D-15 	;
{216} O3 + MACR = 0.12 HCHO + 
                 0.12 MGLYOOB + 
                 0.88 MGLY + 0.88 CH2OOG 	: 1.4D-15*EXP(-2100/TEMP) 	;
{217} OH + MACR = 0.45 MACO3 + 
                 0.47 MACRO2 + 
                 0.08 MACROHO2 		: 8.0D-12*EXP(380/TEMP) 	;
{218} O3 + MVK = 0.5 MGLOOA + 0.5 HCHO + 
                0.5 MGLY + 0.5 CH2OOB 	: 8.5D-16*EXP(-1520/TEMP) 	;
{219} OH + MVK = 0.3 HMVKAO2 + 
                 0.7 HMVKBO2			: 2.6D-12*EXP(610/TEMP) 	;
{210} MACROOA = 0.255 C3H6 + 0.255 C2O3 + 
               0.255 HCHO + 0.255 HO2 + 
               0.22 MACROO + 
               0.27 OH + 0.27 CO + 
               0.27 C2O3 + 0.27 HCHO 	: KDEC 	;
{211} MVKOOA = 0.255 C3H6 + 0.255 CH3O2 + 
              0.255 HCHO + 0.255 CO + 
              0.255 HO2 + 0.22 MVKOO + 
              0.27 MVKO2 			: KDEC 	;
{212} ISOPAO2 + HO2 = ISOPAOOH 		: KRO2HO2*0.706 	;
{213} ISOPAO2 + NO = 0.1 ISOPANO3 + 
                    0.9 ISOPAO + 0.9 NO2 	: KRO2NO 	;
{214} ISOPAO2 + NO3 = ISOPAO + NO2 		: KRO2NO3 	;
{215} ISOPBO2 + HO2 = ISOPBOOH 		: KRO2HO2*0.706 	;
{216} ISOPBO2 + NO = 0.066 ISOPBNO3 + 
                    0.934 ISOPBO + 
                    0.934 NO2 		: KRO2NO 	;
{217} ISOPBO2 + NO3 = ISOPBO + NO2 		: KRO2NO3 	;
{218} ISOPCO2 + HO2 = ISOPCOOH 		: KRO2HO2*0.706 	;
{219} ISOPCO2 + NO = 0.1 ISOPCNO3 +
                    0.9 ISOPCO + 0.9 NO2 	: KRO2NO 	;
{220} ISOPCO2 + NO3 = ISOPCO + NO2 		: KRO2NO3 	;
{221} ISOPDO2 + HO2 = ISOPDOOH 		: KRO2HO2*0.706 	;
{222} ISOPDO2 + NO = 0.134 ISOPDNO3 + 
                    0.866 ISOPDO + 
                    0.866 NO2 		: KRO2NO 	;
{223} ISOPDO2 + NO3 = ISOPDO + NO2 		: KRO2NO3 	;
{224} OH + NISOPOOH = NC4CHO + OH 		: 1.03D-10 	;
{225} OH + NISOPNO3 = NC4CHO + NO2 		: 8.55D-11 	;
{226} NISOPO = NC4CHO + HO2 			: KROPRIM*O2 	;
{227} O3 + ISOPCNO3 = 0.5 GAOOB + 0.5 NOA +
                     0.5 HOCH2CHO + 
                     0.5 NC3OOA 		: 5.30D-17 	;
{228} OH + ISOPCNO3 = INCO2 			: 6.93D-11 	;
{229} NO3 + NC4CHO = NC4CO3 + HNO3 		: KNO3AL*4.25 	;
{230} O3 + NC4CHO = 0.5 NOA + 0.5 GLYOOC
                   0.5 NOAOOA + 0.5 GLY  	: 2.40D-17 	;
{231} OH + NC4CHO = 0.52 C510O2 + 
                   0.48 NC4CO3 		: 4.16D-11	;
{232} CH2OO + CO = HCHO 			: 1.20D-15 	;
{233} CH2OO + NO = HCHO + NO2 		: 1.00D-14 	;
{234} CH2OO + NO2 = HCHO + NO3 		: 1.00D-15 	;
{235} CH2OO + SO2 = HCHO + H2SO4 		: 7.00D-14 	;
{236} CH2OO = 0.375 HCHO + 0.375 H2O2 + 
             0.625 HCOOH 			: 1.60D-17*H2O 	;
{237} CH3C2H2O2 = 0.35 C2O3 + 0.35 HCHO +	
                 0.65 HCHO + 0.65 CH3O2 + 
                 0.65 CO 			: KDEC 	;
{238} MACO3 + HO2 = 0.44 CH3C2H2O2 + 
                   0.44 OH + 
                   0.15 MACO2H + 0.15 O3 +
                   0.41 MACO3H 		: KAPHO2 	;
{239} MACO3 + NO = CH3C2H2O2 + NO2 		: 8.70D-12*EXP(290/TEMP) 	;
{240} MACO3 + NO2 = MPAN 			: KFPAN 	;
{241} MACO3 + NO3 = CH3C2H2O2 + NO2 	: KRO2NO3*1.74 	;
{242} MGLYOOB = 0.18 MGLYOO + 0.82 OH + 
                0.82 CO + 0.82 C2O3 	: KDEC 	;
{243} NO3 + MGLY = C2O3 + CO + HNO3 	: KNO3AL*2.4 	;
{244} OH + MGLY = C2O3 + CO 			: 1.9D-12*EXP(575/TEMP) 	;
{245} CH2OOG = 0.37 CH2OO + 0.47 CO + 
               0.16 HO2 + 0.16 CO + 
               0.16 OH 				: KDEC 	;
{246} MACRO2 + HO2 = MACROOH 			: KRO2HO2*0.625 	;
{247} MACRO2 + NO = MACRO + NO2 		: KRO2NO 	;
{248} MACRO2 + NO3 = MACRO + NO2 		: KRO2NO3 	;
{249} MACROHO2 + HO2 = MACROHOOH 		: KRO2HO2*0.625 	;
{250} MACROHO2 + NO = 0.017 MACROHNO3 +
                      0.983 MACROHO + 
                      0.983 NO2 		: KRO2NO 	;
{251} MACROHO2 + NO3 = MACROHO + NO2 	: KRO2NO3 	;
{252} NO3 + C3H6 = 0.35 PRONO3AO2 +	
                   0.65 PRONO3BO2 		: 4.6D-13*EXP(-1155/TEMP) 	;
{253} O3 + C3H6 = 0.5 CH2OOB + 0.5 ALD2 + 
                  0.5 CH3CHOOA + 0.5 HCHO : 5.5D-15*EXP(-1880/TEMP) 	;
{254} OH + C3H6 = 0.87 HYPROPO2 + 
                  0.13 IPROPOLO2 		: KMT16 	;
{255} MGLOOA = 0.2 ALD2 + 0.2 C2O3 + 
               0.2 HCHO + 0.2 HO2 + 
               0.24 MGLOO + 
               0.36 OH + 0.36 CO + 
               0.36 C2O3	 		: KDEC 	;
{256} CH2OOB = 0.24 CH2OO +
               0.4 CO + 0.36 HO2 + 
               0.36 CO + 0.36 OH 		: KDEC 	;
{257} HMVKAO2 + HO2 = HMVKAOOH 		: KRO2HO2*0.625 	;
{258} HMVKAO2 + NO = 0.017 HMVKANO3 + 
                     0.983 HMVKAO + 
                     0.983 NO2 		: KRO2NO 	;
{259} HMVKAO2 + NO3 = HMVKAO + NO2 		: KRO2NO3 	;
{260} HMVKBO2 + HO2 = HMVKBOOH 		: KRO2HO2*0.625 	;
{261} HMVKBO2 + NO = HMVKBO + NO2 		: KRO2NO 	;
{262} HMVKBO2 + NO3 = HMVKBO + NO2 		: KRO2NO3 	;
{263} MACROO + CO = MACR 			: 1.2D-15 	;
{264} MACROO + NO = MACR + NO2 		: 1.0D-14 	;
{265} MACROO + NO2 = MACR + NO3 		: 1.0D-15 	;
{266} MACROO + SO2 = MACR + H2SO4 		: 7.0D-14 	;
{267} MACROO = 0.625 MACO2H +	
               0.375 MACR + 0.375 H2O2 	: 1.6D-17*H2O 	;
{268} MVKOO + CO = MVK 				: 1.2D-15 	;
{269} MVKOO + NO = MVK + NO2 			: 1.0D-14 	;
{270} MVKOO + NO2 = MVK + NO3 		: 1.0D-15 	;
{271} MVKOO + SO2 = MVK + H2SO4 		: 7.0D-14 	;
{272} MVKOO = MVK + H2O2 			: 6.0D-18*H2O 	;
{273} MVKO2 + HO2 = MVKOOH 			: KRO2HO2*0.625 	;
{274} MVKO2 + NO = MVKO + NO2 		: KRO2NO 	;
{275} MVKO2 + NO3 = MVKO + NO2 		: KRO2NO3 	;
{276} OH + ISOPAOOH = 0.05 HC4ACHO +  
                      0.93 IEPOXA + 
                      0.98 OH + 
                      0.02 ISOPAO2 		: 1.54D-10 	;

{277} O3 + ISOPANO3 = 0.5 ACETOL + 
                      0.5 NC2OOA + 
                      0.5 ACLOOA + 
                      0.5 NO3CH2CHO 	: 5.30D-17 	;
{278} OH + ISOPANO3 = INAO2 			: 6.93D-11 	;
{279} ISOPAO = 0.25 C524O2 + 
               0.75 HC4CCHO + 0.75 HO2 	: KDEC 	;
{280} NO3 + HC4ACHO = HC4ACO3 + HNO3 	: KNO3AL*4.25 	;
{281} O3 + HC4ACHO = 0.5 ACETOL + 
                     0.5 GLYOOC + 
                     0.5 ACLOOA + 0.5 GLY : 2.40D-17 	;
{282} OH + HC4ACHO = 0.52 C58O2 + 
                     0.48 HC4ACO3 		: 4.52D-11 	;
{283} OH + ISOPAOH = 0.5 HC4ACHO + 
                     0.5 HO2 + 
                     0.5 HC4CCHO 		: 9.30D-11 	;
{284} OH + ISOPBOOH = 0.92 IEPOXB + 
                      0.92 OH + 
                      0.08 ISOPBO2 		: 5.00D-11 	;
{285} O3 + ISOPBNO3 = 0.5 HCHO + 
                      0.5 MACRNOOA + 
                      0.5 MACRNO3 + 
                      0.5 CH2OOB 		: 1.06D-16 	;
{286} OH + ISOPBNO3 = 0.72 INB1O2 + 
                      0.28 INB2O2 		: 1.36D-11 	;
{287} ISOPBO = MVK + HCHO + HO2 		: KDEC 	;
{288} OH + ISOPBOH = ISOPBO 			: 3.85D-11 	;
{289} OH + ISOPCOOH = 0.05 HC4CCHO + 
                      0.93 IEPOXC + 
                      0.98 OH + 
                      0.02 ISOPCO2 		: 1.54D-10 	;
{290} ISOPCO = 0.75 HC4ACHO + HO2 +	
               0.25 HC4CCHO 			: KDEC 	;
{291} NO3 + HC4CCHO = HC4CCO3 + HNO3 	: KNO3AL*4.25 	;
{292} O3 + HC4CCHO = 0.5 MGLYOOA + 
                     0.5 HOCH2CHO + 
                     0.5 MGLY + 0.5 GAOOB : 2.40D-17 	;
{293} OH + HC4CCHO = 0.52 C57O2 + 
                     0.48 HC4CCO3 		: 4.52D-11 	;
{294} OH + ISOPDOOH = 0.22 HCOC5 + 
                      0.75 IEPOXB + 
                      0.97 OH + 
                      0.03 ISOPDO2 		: 1.15D-10 	;
{295} O3 + ISOPDNO3 = 0.5 CH2OOC + 
                      0.5 MVKNO3 + 
                      0.5 HCHO + 
                      0.5 NC4OOA 		: 1.06D-16 	;
{296} OH + ISOPDNO3 = INDO2 			: 2.54D-11 	;
{297} ISOPDO = MACR + HCHO + HO2 		: KDEC 	;
{298} OH + HCOC5 = C59O2 			: 3.81D-11 	;
{299} OH + ISOPDOH = HCOC5 + HO2 		: 7.38D-11 	;
{300} GAOOB = 0.11 GAOO + 
              0.89 OH + 0.89 HO2 + 
              0.89 CO + 0.89 HCHO 		: KDEC 	;
{301} NOA + OH = MGLY + NO2 			: 1D-12 	;
{302} HOCH2CHO + NO3 = HOCH2CO3 + HNO3 	: KNO3AL 	;
{303} HOCH2CHO + OH = 0.2 GLY + 0.2 HO2 + 
                      0.8 HOCH2CO3 		: 1.00D-11 	;
{304} NC3OOA = 0.11 NC3OO + 
{305} NC3OOA = 0.89 OH + 0.89 NO2 + 
               0.89 MGLY 			: KDEC 	;
{306} INCO2 + HO2 = INCOOH 			: KRO2HO2*0.706 	;
{307} INCO2 + NO = 0.145 INCNO3 + 
                   0.855 INCO + 0.855 NO2 : KRO2NO 	;
{308} INCO2 + NO3 = INCO + NO2 		: KRO2NO3 	;
{309} NC4CO3 + HO2 = 0.15 NC4CO2H + 
                     0.15 O3 + 
                     0.41 NC4CO3H +	
                     0.44 NOA + 0.44 CO + 
                     0.44 HO2 + 0.44 OH 	: KAPHO2 	;
{310} NC4CO3 + NO = NOA + CO + HO2 + NO2 	: KAPNO 	;
{311} NC4CO3 + NO2 = C5PAN18 			: KFPAN 	;
{312} NC4CO3 + NO3 = NOA + CO + HO2 + NO2 : KRO2NO3*1.74 	;
{313} GLYOOC = 0.11 GLYOO + 
               0.89 OH + 0.89 HO2 + 
               1.78 CO 				: KDEC 	;
{314} NOAOOA = 0.11 NOAOO + 
               0.89 OH + 0.89 NO2 + 
               0.89 MGLY 			: KDEC 	;
{315} GLY + NO3 = 1.2 CO + 0.6 HO2 + 
                  0.4 HCOCO3 + HNO3 	: KNO3AL 	;
{316} GLY + OH = 1.2 CO + 0.6 HO2 +	
                 0.4 HCOCO3 			: 3.1D-12*EXP(340/TEMP) 	;
{317} C510O2 + HO2 = C510OOH 			: KRO2HO2*0.706 	;
{318} C510O2 + NO = C510O + NO2 		: KRO2NO 	;
{319} C510O2 + NO3 = C510O + NO2 		: KRO2NO3 	;
{320} OH + MACO2H = CH3C2H2O2 		: 1.51D-11 	;
{321} OH + MACO3H = MACO3 			: 1.661D-11 	;
{322} MPAN = MACO3 + NO2 			: 1.6D+16*EXP(-13500/TEMP) 	;
{323} O3 + MPAN = HCHO + C2O3 + NO3 	: 8.2D-18 	;
{324} OH + MPAN = ACETOL + CO + NO2 	: 2.9D-11 	;
{325} MGLYOO + CO = MGLY 			: 1.2D-15 	;
{326} MGLYOO + NO = MGLY + NO2 		: 1.0D-14 	;
{327} MGLYOO + NO2 = MGLY + NO3 		: 1.0D-15 	;
{328} MGLYOO + SO2 = MGLY + H2SO4 		: 7.0D-14 	;
{329} MGLYOO = MGLY + H2O2 			: 6.0D-18*H2O 	;
{330} OH + MACROOH = ACETOL + CO + OH 	: 3.77D-11 	;
{331} MACRO = ACETOL + CO + HO2 		: KDEC 	;
{332} OH + MACROH = C3MDIALOH + HO2 	: 3.42D-11 	;
{333} MACROHOOH + OH = C3MDIALOH + OH 	: 5.55D-11 	;
{334} MACROHNO3 + OH = NOA + CO + HO2 	: 2.18D-11 	;
{335} MACROHO = MGLY + HCHO + HO2 		: KDEC 	;
{336} C3MDIALOH + OH = CHOMOHCO3 		: 1.32D-10 	;
{337} PRONO3AO2 + HO2 = PR1O2HNO3 		: KRO2HO2*0.520 	;
{338} PRONO3AO2 + NO = PRONO3AO + NO2 	: KRO2NO 	;
{339} PRONO3AO2 + NO3 = PRONO3AO + NO2 	: KRO2NO3 	;
{340} PRONO3BO2 + HO2 = PR2O2HNO3 		: KRO2HO2*0.520 	;
{341} PRONO3BO2 + NO = PRONO3BO + NO2 	: KRO2NO 	;
{342} PRONO3BO2 + NO3 = PRONO3BO + NO2 	: KRO2NO3 	;
{343} CH3CHOOA = 0.24 CH3CHOO + 
                 0.36 CH3O2 + 
                 0.36 CO + 0.36 OH + 
                 0.2 CH3O2 + 0.2 HO2 + 
                 0.2 CH4 			: KDEC 	;
{344} HYPROPO2 + HO2 = HYPROPO2H 		: KRO2HO2*0.520 	;
{345} HYPROPO2 + NO3 = HYPROPO + NO2 	: KRO2NO3 	;
{346} NO + HYPROPO2 = 0.977 HYPROPO + 
                      0.977 NO2 + 
                      0.023 PROPOLNO3 	: KRO2NO 	;
{347} IPROPOLO2 + HO2 = IPROPOLO2H 		: KRO2HO2*0.520 	;
{348} IPROPOLO2 + NO = 0.991 IPROPOLO + 
                       0.991 NO2 + 
                       0.009 PROLNO3 	: KRO2NO 	;
{349} IPROPOLO2 + NO3 = IPROPOLO + NO2 	: KRO2NO3 	;
{350} CH3CO2H + OH = CH3O2 			: 8.00D-13 	;
{351} CH3CO3H + OH = C2O3 			: 3.70D-12 	;
{352} MGLOO + CO = MGLY 			: 1.2D-15 	;
{353} MGLOO + NO = MGLY + NO2 		: 1.0D-14 	;
{354} MGLOO + NO2 = MGLY + NO3 		: 1.0D-15 	;
{355} MGLOO + SO2 = MGLY + H2SO4 		: 7.0D-14 	;
{356} MGLOO = 0.625 CH3COCO2H + 
              0.375 MGLY + 0.375 H2O2 	: 1.6D-17*H2O 	;
{357} OH + HMVKAOOH = CO2H3CHO + OH 	: 5.77D-11 	;
{358} OH + HMVKANO3 = CO2H3CHO + NO2 	: 2.23D-12 	;
{359} HMVKAO = MGLY + HCHO + HO2 		: KDEC 	;
{360} NO3 + CO2H3CHO = CO2H3CO3 + HNO3 	: KNO3AL*4.0 	;
{361} OH + CO2H3CHO = CO2H3CO3 		: 2.45D-11 	;
{362} OH + HO12CO3C4 = BIACETOH + HO2 	: 1.88D-11 	;
{363} OH + HMVKBOOH = BIACETOH + OH 	: 3.95D-11 	;
{364} HMVKBO = C2O3 + HOCH2CHO 		: KDEC 	;
{365} OH + BIACETOH = CO23C3CHO + HO2 	: 2.69D-12 	;
{366} OH + MVKOOH = MVKO2 			: 1.90D-12*EXP(190/TEMP) 	;
{367} OH + MVKOOH = VGLYOX + OH 		: 2.55D-11 	;
{368} MVKO = HCHO + ACO3 			: KDEC 	;
{369} MVKOH + O3 = 0.5 HMGLOOA + 
                   0.5 HCHO +
                   0.5 HOCH2COCHO + 
                   0.5 CH2OOB 		: 7.51D-16*EXP(-1521/TEMP) 	;
{370} MVKOH + OH = 0.3 MVKOHAO2 + 
                   0.7 MVKOHBO2 		: 4.60D-12*EXP(452/TEMP) 	;
{371} NO3 + VGLYOX = CO + ACO3 + HNO3 	: KNO3AL*2.0 	;
{372} OH + VGLYOX = CO + ACO3 		: 2.95D-11 	;
{373} OH + IEPOXA = IEACHO + HO2 		: 5.23D-12 	;
{374} ACETOL + OH = MGLY + HO2 		: 1.6D-12*EXP(305/TEMP) 	;
{375} NC2OOA = 0.11 NC2OO + 
               0.89 OH + 0.89 NO2 + 
               0.89 GLY  			: KDEC 	;
{376} ACLOOA = 0.11 ACLOO + 
               0.89 OH + 0.89 HO2 + 
               0.89 MGLY 			: KDEC 	;
{377} NO3CH2CHO + NO3 = NO3CH2CO3 + HNO3 	: KNO3AL 	;
{378} NO3CH2CHO + OH = NO3CH2CO3 		: 3.40D-12 	;
{379} INAO2 + HO2 = INAOOH 			: KRO2HO2*0.706 	;
{380} INAO2 + NO = 0.072 INANO3 + 
                   0.928 INAO + 0.928 NO2 : KRO2NO 	;
{381} INAO2 + NO3 = INAO + NO2 		: KRO2NO3 	;
{382} C524O2 + HO2 = C524OOH 			: KRO2HO2*0.706 	;
{383} C524O2 + NO = 0.134 C524NO3 +	
                    0.866 C524O + 
                    0.866 NO2 		: KRO2NO 	;
{384} C524O2 + NO3 = C524O + NO2 		: KRO2NO3 	;
{385} HC4ACO3 + HO2 = 0.44 ACETOL + 
                      0.44 CO + 0.44 HO2 + 
                      0.44 OH + 
                      0.15 HC4ACO2H + 
                      0.15 O3 + 
                      0.41 HC4ACO3H 	: KAPHO2 	;
{386} HC4ACO3 + NO = ACETOL + CO + 
                     HO2 + NO2 		: KAPNO 	;
{387} HC4ACO3 + NO2 = C5PAN17 		: KFPAN 	;
{388} HC4ACO3 + NO3 = ACETOL + CO + 
                      HO2 + NO2 		: KRO2NO3*1.74 	;
{389} C58O2 + HO2 = C58OOH 			: KRO2HO2*0.706 	;
{390} C58O2 + NO = 0.019 C58NO3 + 
                   0.981 C58O + 0.981 NO2 : KRO2NO 	;
{391} C58O2 + NO3 = C58O + NO2 		: KRO2NO3 	;
{392} OH + IEPOXB = 0.5 IEB1O2 + 
                    0.5 IEB2O2 		: 9.05D-12 	;
{393} MACRNOOA = 0.36 ACETOL + 0.36 NO2 + 
                 0.36 CO + 0.36 OH + 
                 0.2 ACETOL + 0.2 NO2 + 
                 0.2 HO2 + 0.24 MACRNOO + 
                 0.2 PROPOLNO3 		: KDEC 	;
{394} MACRNO3 + OH = 0.16 CONM2CHO + 
                     0.16 HO2 + 
                     0.84 MACRNCO3 		: 4.34D-12 	;
{395} INB1O2 + HO2 = INB1OOH 			: KRO2HO2*0.706 	;
{396} INB1O2 + NO = 0.145 INB1NO3 +
                    0.855 INB1O + 
                    0.855 NO2 		: KRO2NO 	;
{397} INB1O2 + NO3 = INB1O + NO2 		: KRO2NO3 	;
{398} INB2O2 + HO2 = INB2OOH 			: KRO2HO2*0.706 	;
{399} INB2O2 + NO = 0.064 INANO3 +
                    0.936 INB2O + 
                    0.936 NO2 		: KRO2NO 	;
{400} OH + IEPOXC = 0.719 IEC1O2 +
                    0.281 IECCHO + 
                    0.281 HO2 		: 1.50D-11 	;
{401} HC4CCO3 + HO2 = 0.44 C2O3 + 
                      0.44 HOCH2CHO + 
                      0.44 OH +
                      0.15 HC4CCO2H + 
                      0.15 O3 + 
                      0.41 HC4CCO3H 	: KAPHO2 	;
{402} HC4CCO3 + NO = C2O3 + HOCH2CHO + 
                     NO2 			: KAPNO 	;
{403} HC4CCO3 + NO2 = C5PAN19 		: KFPAN 	;
{404} HC4CCO3 + NO3 = C2O3 + HOCH2CHO + 
                      NO2 			: KRO2NO3*1.74 	;
{405} MGLYOOA = 0.11 MGLYOO + 
                0.89 OH + 0.89 CO + 
                0.89 C2O3 			: KDEC 	;
{406} C57O2 + HO2 = C57OOH 			: KRO2HO2*0.706 	;
{407} C57O2 + NO = C57O + NO2 		: KRO2NO 	;
{408} C57O2 + NO3 = C57O + NO2 		: KRO2NO3 	;
{409} CH2OOC = 0.18 CH2OO + 0.82 HO2 + 
               0.82 CO + 0.82 OH 		: KDEC	;
{410} MVKNO3 + OH = 0.33 BIACETOH + 
                    0.33 NO2 + 
                    0.67 CO2N3CHO + 
                    0.67 HO2 			: 1.33D-12 	;
{411} NC4OOA = 0.18 NC4OO + 
               0.82 OH + 0.82 NO2 + 
               0.82 BIACETOH 			: KDEC 	;
{412} INDO2 + HO2 = INDOOH 			: KRO2HO2*0.706 	;
{413} INDO2 + NO = 0.072 INB1NO3 + 
                   0.928 INDO + 0.928 NO2 : KRO2NO 	;
{414} INDO2 + NO3 = INDO + NO2 		: KRO2NO3 	;
{415} HOCH2CO3 + HO2 = 0.44 HO2 + 
                       0.44 HCHO + 
                       0.44 OH +
                       0.15 HOCH2CO2H + 
                       0.15 O3 + 
                       0.41 HOCH2CO3H 	: KAPHO2 	;
{416} HOCH2CO3 + NO = NO2 + HO2 + HCHO 	: KAPNO 	;
{417} HOCH2CO3 + NO2 = PHAN 			: KFPAN 	;
{418} HOCH2CO3 + NO3 = NO2 + HO2 + HCHO 	: KRO2NO3*1.74 	;
{419} C59O2 + HO2 = C59OOH 			: KRO2HO2*0.706 	;
{420} C59O2 + NO = C59O + NO2 		: KRO2NO 	;
{421} C59O2 + NO3 = C59O + NO2 		: KRO2NO3 	;
{422} GAOO + CO = HOCH2CHO 			: 1.2D-15 	;
{423} GAOO + NO = HOCH2CHO + NO2 		: 1.0D-14 	;
{424} GAOO + NO2 = HOCH2CHO + NO3 		: 1.0D-15 	;
{425} GAOO + SO2 = HOCH2CHO + H2SO4 	: 7.0D-14 	;
{426} GAOO = 0.375 HOCH2CHO + 0.375 H2O2 +
             0.675 HOCH2CO2H 			: 1.6D-17*H2O 	;
{427} NC3OO + CO = NOA 				: 1.2D-15 	;
{428} NC3OO + NO = NOA + NO2 			: 1.0D-14 	;
{429} NC3OO + NO2 = NOA + NO3 		: 1.0D-15 	;
{430} NC3OO + SO2 = NOA + H2SO4 		: 7.0D-14 	;
{431} NC3OO = NOA + H2O2 			: 6.0D-18*H2O 	;
{432} OH + INCOOH = 0.89 INCCO + 0.89 OH +
                    0.11 INCO2 		: 3.31D-11 	;
{433} OH + INCNO3 = 0.445 INCCO +  
                    0.414 INCNCHO + 
                    0.414 HO2 +
                    0.141 NOA + 
                    0.141 HOCH2CHO + 
                    0.586 NO2 		: 1.98D-12 	;
{434} INCO = NOA + HO2 + HOCH2CHO 		: KDEC 	;
{435} OH + INCCO = INCGLYOX + HO2 		: 3.30D-12 	;
{436} OH + INCOH = INCCO + HO2 		: 1.53D-11 	;
{437} OH + NC4CO2H = NOA + CO + HO2 	: 2.16D-11 	;
{438} OH + NC4CO3H = NC4CO3 			: 2.52D-11 	;
{439} C5PAN18 = NC4CO3 + NO2 			: KBPAN 	;
{440} OH + C5PAN18 = NOA + 2 CO + NO2 	: 2.16D-11 	;
{441} GLYOO + CO = GLY  			: 1.2D-15 	;
{442} GLYOO + NO = GLY  + NO2 		: 1.0D-14 	;
{443} GLYOO + NO2 = GLY  + NO3 		: 1.0D-15 	;
{444} GLYOO + SO2 = GLY  + H2SO4 		: 7.0D-14 	;
{445} GLYOO = 0.375 GLY  + 0.375 H2O2 + 
              0.625 HCOCO2H 			: 1.6D-17*H2O 	;
{446} NOAOO + CO = NOA 				: 1.2D-15 	;
{447} NOAOO + NO = NOA + NO2 			: 1.0D-14 	;
{448} NOAOO + NO2 = NOA + NO3 		: 1.0D-15 	;
{449} NOAOO + SO2 = NOA + H2SO4 		: 7.0D-14 	;
{450} NOAOO = NOA + H2O2 			: 6.0D-18*H2O 	;
{451} HCOCO3 + HO2 = 0.15 HCOCO2H + 
                     0.15 O3 + 
                     0.41 HCOCO3H + 
                     0.44 HO2 + 0.44 CO + 
                     0.44 OH 			: KAPHO2 	;
{452} HCOCO3 + NO = HO2 + CO + NO2 		: KAPNO 	;
{453} HCOCO3 + NO2 = GLYPAN 			: KFPAN 	;
{454} HCOCO3 + NO3 = HO2 + CO + NO2 	: KRO2NO3*1.74 	;
{455} OH + C510OOH = C510O2 			: 2.81D-11 	;
{456} C510O = NOA + GLY  + HO2 		: KDEC 	;
{457} OH + C510OH = C510O 			: 2.69D-11 	;
{458} IBUTALOH + OH = IPRHOCO3 		: 1.4D-11 	;
{459} CHOMOHCO3 + HO2 = 0.56 CHOMOHCO3H + 
                        0.44 MGLY + 
                        0.44 HO2 + 
                        0.44 OH 		: KAPHO2 	;
{460} CHOMOHCO3 + NO = MGLY + NO2 + HO2 	: KAPNO 	;
{461} CHOMOHCO3 + NO2 = CHOMOHPAN 		: KFPAN 	;
{462} CHOMOHCO3 + NO3 = MGLY + NO2 + HO2 	: KRO2NO3*1.74 	;
{463} PR1O2HNO3 + OH = CHOPRNO3 + OH 	: 1.69D-12 	;
{464} PR1O2HNO3 + OH = PRONO3AO2 		: 1.90D-12*EXP(190/TEMP) 	;
{465} PRONO3AO = CHOPRNO3 + HO2 		: KROPRIM*O2 	;
{466} PRONO3AO = HCHO + ALD2 + NO2 		: 7.00D+03 	;
{467} CHOPRNO3 + NO3 = PRNO3CO3 + HNO3 	: KNO3AL*2.4 	;
{468} CHOPRNO3 + OH = PRNO3CO3 		: 3.55D-12 	;
{469} PROPOLNO3 + OH = ACETOL + NO2 	: 9.16D-13 	;
{470} PR2O2HNO3 + OH = NOA + OH 		: 3.47D-12 	;
{471} PR2O2HNO3 + OH = PRONO3BO2 		: 1.90D-12*EXP(190/TEMP) 	;
{472} PRONO3BO = ALD2 + HCHO + NO2 		: 7.00D+03 	;
{473} PRONO3BO = NOA + HO2 			: KROSEC*O2 	;
{474} PROLNO3 + OH = CH3CHOHCHO + NO2 	: 1.71D-12 	;
{475} HCOCH2O2 + HO2 = HCOCH2OOH 		: KRO2HO2*0.387 	;
{476} HCOCH2O2 + NO = NO2 + HCOCH2O 	: KRO2NO 	;
{477} HCOCH2O2 + NO3 = HCOCH2O + NO2 	: KRO2NO3 	;
{478} CH3CHOO + CO = ALD2 			: 1.20D-15 	;
{479} CH3CHOO + NO = ALD2 + NO2 		: 1.00D-14 	;
{480} CH3CHOO + NO2 = ALD2 + NO3 		: 1.00D-15 	;
{481} CH3CHOO + SO2 = ALD2 + H2SO4 		: 7.00D-14 	;
{482} CH3CHOO = 0.375 ALD2 + 0.375 H2O2 +
                0.625 CH3CO2H 		: 1.60D-17*H2O 	;
{483} HYPROPO2H + OH = ACETOL + OH 		: 2.44D-11 	;
{484} HYPROPO2H + OH = HYPROPO2 		: 1.90D-12*EXP(190/TEMP) 	;
{485} HYPROPO = ALD2 + HCHO + HO2 		: 2.00D+14*EXP(-6410/TEMP) 	;
{486} PROPGLY + OH = 0.613 ACETOL + 
                     0.387 CH3CHOHCHO +
                     HO2 		 	: 1.20D-11 	;
{487} IPROPOLO2H + OH = CH3CHOHCHO + OH 	: 1.83D-11 	;
{488} IPROPOLO2H + OH = IPROPOLO2 		: 1.90D-12*EXP(190/TEMP) 	;
{489} IPROPOLO = ALD2 + HCHO + HO2 		: 2.00D+14*EXP(-5505/TEMP) 	;
{490} CH3CHOHCHO + NO3 = CH3CHOHCO3 + 
                         HNO3 		: KNO3AL*2.4 	;
{491} CH3CHOHCHO + OH = CH3CHOHCO3 		: 1.7D-11 	;
{492} OH + CH3COCO2H = C2O3 			: 8.0D-13 	;
{493} CO2H3CO3 + HO2 = 0.56 CO2H3CO3H +
                       0.44 MGLY + 
                       0.44 HO2 + 0.44 OH : KAPHO2 	;
{494} CO2H3CO3 + NO = MGLY + HO2 + NO2 	: KAPNO 	;
{495} CO2H3CO3 + NO2 = C4PAN6 		: KFPAN 	;
{496} CO2H3CO3 + NO3 = MGLY + HO2 + NO2 	: KRO2NO3*1.74 	;
{497} NO3 + CO23C3CHO = C2O3 + 
                        2 CO + HNO3 	: KNO3AL*4.0 	;
{498} OH + CO23C3CHO = C2O3 + 2 CO  	: 1.23D-11 	;
{499} ACO3 + HO2 = 0.15 ACO2H + 0.15 O3 + 
                   0.41 ACO3H + 
                   0.44 HO2 + 0.44 CO + 
                   0.44 HCHO + 0.44 OH 	: KAPHO2 	;
{500} ACO3 + NO = HO2 + CO + HCHO + NO2 	: KAPNO 	;
{501} ACO3 + NO2 = ACRPAN 			: KFPAN 	;
{502} ACO3 + NO3 = HO2 + CO + HCHO + NO2 	: KRO2NO3*1.74 	;
{503} HMGLOOA = 0.24 HMGLOO +	
                0.2 HOCH2CHO + 
                0.2 HOCH2CO3 + 0.2 HO2 + 
                0.36 OH + 0.36 CO + 
                0.36 HOCH2CO3 		: KDEC 	;
{504} NO3 + HOCH2COCHO = HOCH2CO3 + 
                         CO + HNO3 		: KNO3AL*2.4 	;
{505} OH + HOCH2COCHO = HOCH2CO3 + CO 	: 1.44D-11 	;
{506} MVKOHAO2 + HO2 = MVKOHAOOH 		: KRO2HO2*0.625 	;
{507} MVKOHAO2 + NO = 0.017 MVKOHANO3 + 
                      0.983 MVKOHAO + 
                      0.983 NO2 		: KRO2NO 	;
{508} MVKOHAO2 + NO3 = MVKOHAO + NO2 	: KRO2NO3 	;
{509} MVKOHBO2 + HO2 = MVKOHBOOH 		: KRO2HO2*0.625 	;
{510} MVKOHBO2 + NO = MVKOHBO + NO2 	: KRO2NO 	;
{511} MVKOHBO2 + NO3 = MVKOHBO + NO2 	: KRO2NO3 	;
{512} OH + ALLYLOH = ACR + HO2 		: 2.59D-11 	;
{513} NO3 + IEACHO = IEACO3 + HNO3 		: KNO3AL*7.5 	;
{514} OH + IEACHO = IEACO3 			: 2.20D-11 	;
{515} NC2OO + CO = NO3CH2CHO 			: 1.2D-15 	;
{516} NC2OO + NO = NO3CH2CHO + NO2 		: 1.0D-14 	;
{517} NC2OO + NO2 = NO3CH2CHO + NO3 	: 1.0D-15 	;
{518} NC2OO + SO2 = NO3CH2CHO + H2SO4 	: 7.0D-14 	;
{519} NC2OO = NO3CH2CHO + H2O2 		: 6.0D-18*H2O 	;
{520} NC2OO = NO3CH2CO2H 			: 1.0D-17*H2O 	;
{521} ACLOO + CO = ACETOL 			: 1.2D-15 	;
{522} ACLOO + NO = ACETOL + NO2 		: 1.0D-14 	;
{523} ACLOO + NO2 = ACETOL + NO3 		: 1.0D-15 	;
{524} ACLOO + SO2 = ACETOL + H2SO4 		: 7.0D-14 	;
{525} ACLOO = ACETOL + H2O2 			: 6.0D-18*H2O 	;
{526} NO3CH2CO3 + HO2 = 0.44 HCHO + 
                        0.44 NO2 + 
                        0.44 OH + 
                        0.15 NO3CH2CO2H + 
                        0.15 O3 + 
                        0.41 NO3CH2CO3H 	: KAPHO2 	;
{527} NO3CH2CO3 + NO = HCHO + 2 NO2 	: KAPNO 	;
{528} NO3CH2CO3 + NO2 = NO3CH2PAN 		: KFPAN 	;
{529} NO3CH2CO3 + NO3 = HCHO + 2 NO2 	: KRO2NO3*1.74 	;
{530} HCOCH2O = HCHO + CO + HO2 		: KDEC 	;
{531} OH + INAOOH = 0.65 INAHPCHO + 
                    0.65 HO2 + 0.35 INAO2 : 1.01D-11 	;
{532} OH + INANO3 = 0.07 ACETOL + 
                    0.14 NO2 + 
                    0.07 NO3CH2CHO + 
                    0.39 C58NO3 + 
                    0.39 NO2 + 
                    0.07 HCHO + 
                    0.07 HMVKANO3 + 
                    0.33 INANCHO + 
                    0.47 HO2 + 
                    0.14 INANCO 		: 2.00D-12 	;
{533} INAO = ACETOL + NO3CH2CHO + HO2 	: KDEC 	;
{534} OH + INAOH = INAHCHO + HO2 		: 6.68D-12 	;
{535} C524OOH + OH = 0.22 C524CO + 
                     0.22 OH + 
                     0.03 C524O2 + 
                     0.75 HIEPOXB + 
                     0.75 OH 			: 1.18D-10 	;
{536} C524NO3 + OH = NC524O2 			: 2.94D-11 	;
{537} C524O = HMACR + HCHO + HO2 		: KDEC 	;
{538} C524CO + OH = C525O2 			: 4.21D-11 	;
{539} C524OH + OH = C524CO + HO2 		: 7.78D-11 	;
{540} OH + HC4ACO2H = ACETOL + CO + HO2 	: 2.52D-11 	;
{541} OH + HC4ACO3H = HC4ACO3 		: 2.88D-11 	;
{542} C5PAN17 = HC4ACO3 + NO2 		: KBPAN 	;
{543} OH + C5PAN17 = MACROH + CO + NO2 	: 2.52D-11 	;
{544} OH + C58OOH = C58O2 			: 3.16D-11 	;
{545} OH + C58NO3 = C58NO3CO3 		: 2.32D-11 	;
{546} C58O = ACETOL + GLY  + HO2 		: KDEC 	;
{547} OH + C58OH = C58O 			: 3.04D-11 	;
{548} IEB1O2 + HO2 = IEB1OOH 			: KRO2HO2*0.706 	;
{549} IEB1O2 + NO = IEB1O + NO2 		: KRO2NO 	;
{550} IEB1O2 + NO3 = IEB1O + NO2 		: KRO2NO3 	;
{551} IEB2O2 + HO2 = IEB2OOH 			: KRO2HO2*0.706 	;
{552} IEB2O2 + NO = IEB2O + NO2 		: KRO2NO 	;
{553} IEB2O2 + NO3 = IEB2O + NO2 		: KRO2NO3 	;

{554} MACRNOO + CO = MACRNO3 			: 1.2D-15 	;
{555} MACRNOO + NO = MACRNO3 + NO2 		: 1.0D-14 	;
{556} MACRNOO + NO2 = MACRNO3 + NO3 	: 1.0D-15 	;
{557} MACRNOO + SO2 = MACRNO3 + H2SO4 	: 7.0D-14 	;
{558} MACRNOO = MACRNCO2H 			: 1.0D-17*H2O 	;
{559} MACRNOO = MACRNO3 + H2O2 		: 6.0D-18*H2O 	;
{560} CONM2CHO + OH = CONM2CO3 		: 6.78D-12 	;
{561} MACRNCO3 + HO2 = 0.44 ACETOL + 
                       0.44 NO2 + 
                       0.44 OH + 
                       0.15 MACRNCO2H + 
                       0.15 O3 +
                       0.41 MACRNCO3H 	: KAPHO2 	;
{562} MACRNCO3 + NO = ACETOL + 2 NO2 	: KAPNO 	;
{563} MACRNCO3 + NO2 = MACRNPAN 		: KFPAN 	;
{564} MACRNCO3 + NO3 = ACETOL + 2 NO2 	: KRO2NO3*1.74 	;
{565} INB1OOH + OH = 0.35 INB1CO + 
                     0.35 OH + 
                     0.34 INB1HPCHO +
                     0.31 INB1O2 		: 1.27D-11 	;
{566} INB1NO3 + OH = 0.5 INB1NACHO + HO2 +
                     0.5 INB1NBCHO	 	: 1.63D-12 	;
{567} INB1O = HOCH2CHO + ACETOL + NO2 	: KDEC 	;
{568} INB1CO + OH = INB1GLYOX + HO2 	: 3.27D-12 	;
{569} INB1OH + OH = 0.71 C58NO3 + HO2 + 
                    0.29 INB1CO  		: 6.65D-12 	;
{570} INB2OOH + OH = 0.73 C58NO3 + 
                     0.73 OH + 
                     0.27 INB2O2 		: 1.59D-11 	;
{571} INB2O = C57NO3 + HO2 			: KDEC 	;
{572} IEC1O2 + HO2 = IEC1OOH 			: KRO2HO2*0.706 	;
{573} IEC1O2 + NO = IEC1O + NO2 		: KRO2NO 	;
{574} IEC1O2 + NO3 = IEC1O + NO2 		: KRO2NO3 	;
{575} NO3 + IECCHO = IECCO3 + HNO3 		: KNO3AL*7.5 	;
{576} OH + IECCHO = IECCO3 			: 2.76D-11 	;
{577} OH + HC4CCO2H = C2O3 + HOCH2CHO 	: 2.52D-11 	;
{578} OH + HC4CCO3H = HC4CCO3 		: 2.88D-11 	;
{579} C5PAN19 = HC4CCO3 + NO2 		: KBPAN 	;
{580} OH + C5PAN19 = HO12CO3C4 + CO + NO2 : 2.52D-11 	;
{581} OH + C57OOH = C57O2 			: 3.16D-11 	;
{582} C57O = MGLY + HOCH2CHO + HO2 		: KDEC 	;
{583} OH + C57OH = C57O 			: 3.04D-11 	;
{584} CO2N3CHO + OH = CO2N3CO3 		: 7.20D-12 	;
{585} NC4OO + CO = MVKNO3 			: 1.2D-15 	;
{586} NC4OO + NO = MVKNO3 + NO2 		: 1.0D-14 	;
{587} NC4OO + NO2 = MVKNO3 + NO3 		: 1.0D-15 	;
{588} NC4OO + SO2 = MVKNO3 + H2SO4 		: 7.0D-14 	;
{589} NC4OO = MVKNO3 + H2O2 			: 6.0D-18*H2O 	;
{590} INDOOH + OH = 0.61 INDHPCHO + 
                    0.61 HO2 + 
                    0.39 INDO2 		: 9.20D-12 	;
{591} INDO = ACETOL + HOCH2CHO + NO2 	: 1.80D+13*(TEMP/298)**1.7*EXP(-4733/TEMP) 	;
{592} INDO = HCHO + HO2 + MVKNO3 		: 1.80D+13*(TEMP/298)**1.7*EXP(-4079/TEMP) 	;
{593} INDOH + OH = INDHCHO + HO2 		: 5.60D-12 	;
{594} HOCH2CO2H + OH = HCHO + HO2 		: 2.73D-12 	;
{595} HOCH2CO3H + OH = HOCH2CO3 		: 6.19D-12 	;
{596} PHAN + OH = HCHO + CO + NO2 		: 1.12D-12 	;
{597} PHAN = HOCH2CO3 + NO2 			: KBPAN 	;
{598} OH + C59OOH = C59O2 			: 9.70D-12 	;
{599} C59O = ACETOL + HOCH2CO3 		: KDEC 	;
{600} INCNCHO + OH = 0.19 INCGLYOX + 
                     0.19 NO2 + 
                     0.81 INCNCO3 		: 4.52D-12 	;
{601} OH + INCGLYOX = MACRNBCO3 + CO 	: 1.34D-11 	;
{602} OH + HCOCO2H = CO + HO2 		: 1.23D-11 	;
{603} OH + HCOCO3H = HCOCO3 			: 1.58D-11 	;
{604} GLYPAN = HCOCO3 + NO2 			: KBPAN 	;
{605} OH + GLYPAN = CO + CO + NO2 		: 1.22D-11 	;
{606} IPRHOCO3 + HO2 = 0.44 AONE + 
                       0.44 HO2 + 
                       0.44 OH +
                       0.15 IPRHOCO2H + 
                       0.15 O3 + 
                       0.41 IPRHOCO3H 	: KAPHO2 	;
{607} IPRHOCO3 + NO = AONE + HO2 + NO2 	: KAPNO 	;
{608} IPRHOCO3 + NO2 = C4PAN5 		: KFPAN 	;
{609} IPRHOCO3 + NO3 = AONE + HO2 + NO2 	: KRO2NO3*1.74 	;
{610} CHOMOHCO3H + OH = CHOMOHCO3 		: 6.99D-11 	;
{611} CHOMOHPAN + OH = MGLY + CO + NO2 	: 6.64D-11 	;
{612} CHOMOHPAN = CHOMOHCO3 + NO2 		: KBPAN 	;
{613} PRNO3CO3 + HO2 = 0.44 ALD2 + 
                       0.44 NO2 + 
                       0.44 OH + 
                       0.15 PRNO3CO2H + 
                       0.15 O3 +
                       0.41 PRNO3CO3H 	: KAPHO2	;
{614} PRNO3CO3 + NO = ALD2 + 2 NO2 		: KAPNO 	;
{615} PRNO3CO3 + NO2 = PRNO3PAN 		: KFPAN 	;
{616} PRNO3CO3 + NO3 = ALD2 + 2 NO2 	: KRO2NO3*1.74 	;
{617} PROPALO = ALD2 + HO2 + CO 		: KDEC 	;
{618} HCOCH2OOH + OH = GLY + OH 		: 2.91D-11 	;
{619} HCOCH2OOH + OH = HCOCH2O2 		: 1.90D-12*EXP(190/TEMP) 	;
{620} CH3CHOHCO3 + HO2 = 0.44 ALD2 + 
                         0.44 HO2 + 
                         0.44 OH +
                         0.56 IPROPOLPER 	: KAPHO2 	;
{621} CH3CHOHCO3 + NO = ALD2 + HO2 + NO2 	: KAPNO 	;
{622} CH3CHOHCO3 + NO2 = IPROPOLPAN 	: KFPAN 	;
{623} CH3CHOHCO3 + NO3 = ALD2 + HO2 + NO2 : KRO2NO3*1.74 	;
{624} OH + CO2H3CO3H = CO2H3CO3 		: 7.34D-12 	;
{625} C4PAN6 = CO2H3CO3 + NO2 		: KBPAN 	;
{626} OH + C4PAN6 = MGLY + CO + NO2 	: 3.74D-12 	;
{627} OH + ACO2H = HO2 + CO + HCHO 		: 8.66D-12 	;
{628} OH + ACO3H = ACO3 			: 1.22D-11 	;
{629} ACRPAN = ACO3 + NO2 			: KBPAN 	;
{630} OH + ACRPAN = HOCH2CHO + CO + NO2 	: 8.63D-12 	;
{631} HMGLOO + CO = HOCH2COCHO 		: 1.20D-15 	;
{632} HMGLOO + NO = HOCH2COCHO + NO2 	: 1.00D-14 	;
{633} HMGLOO + NO2 = HOCH2COCHO + NO3 	: 1.00D-15 	;
{634} HMGLOO + SO2 = HOCH2COCHO + H2SO4 	: 7.00D-14 	;
{635} HMGLOO = HOCH2COCHO + H2O2 		: 6.00D-18*H2O 	;
{636} HMGLOO = HOCH2COCO2H 			: 1.00D-17*H2O 	;
{637} MVKOHAOOH + OH = H13CO2CHO + OH 	: 5.98D-11 	;
{638} MVKOHANO3 + OH = H13CO2CHO + NO2 	: 4.37D-12 	;
{639} MVKOHAO = HOCH2COCHO + HCHO + HO2 	: KDEC 	;
{640} NO3 + H13CO2CHO = H13CO2CO3 + HNO3 	: KNO3AL*4.0 	;
{641} OH + H13CO2CHO = H13CO2CO3 		: 2.66D-11 	;
{642} MVKOHAOH + OH = H13CO2CHO + HO2 	: 2.10D-11 	;
{643} MVKOHBOOH + OH = H14CO23C4 + OH 	: 4.39D-12 	;
{644} MVKOHBO = HOCH2CHO + HOCH2CO3 	: KDEC 	;
{645} H14CO23C4 + OH = H1CO23CHO + HO2 	: 4.44D-12 	;
{646} ACR + NO3 = ACO3 + HNO3 		: 1.72D-13*EXP(-1190/TEMP) 	;
{647} ACR + OH = 0.68 ACO3 + 
                 0.255 ACO3B + 
                 0.065 OCCOHCO2 		: 2.00E-11 	;
{648} O3 + ACR = 0.5 CH2OOB + 0.5 GLY +
                 0.5 GLYOOB + 0.5 HCHO 	: 2.9E-19 	;
{649} IEACO3 + HO2 = 0.44 HMVKBO2 + 
                     0.44 OH + 
                     0.56 IEACO3H 		: KAPHO2 	;
{650} IEACO3 + NO = HMVKBO2 + NO2 		: KAPNO 	;
{651} IEACO3 + NO2 = IEAPAN 			: KFPAN 	;
{652} IEACO3 + NO3 = HMVKBO2 + NO2 		: KRO2NO3*1.74 	;
{653} NO3CH2CO2H + OH = HCHO + NO2 		: 1.68D-13 	;
{654} NO3CH2CO3H + OH = NO3CH2CO3 		: 3.63D-12 	;
{655} NO3CH2PAN + OH = HCHO + CO + 2 NO2 	: 1.12D-14 	;
{656} NO3CH2PAN = NO3CH2CO3 + NO2 		: KBPAN 	;
{657} OH + INAHPCHO = INAHPCO3 		: 2.67D-11 	;
{658} OH + INANCHO = INANCO3 			: 4.22D-12 	;
{659} OH + INANCO = 0.56 INANCOCHO + 
                    0.56 HO2 + 
                    0.44 INB1GLYOX + 
                    0.44 NO2		 	: 1.21D-12 	;
{660} OH + INAHCHO = INAHCO3 			: 2.29D-11 	;
{661} HIEPOXB + OH = 0.667 HIEB1O2 + 
                     0.333 HIEB2O2 		: 1.31D-11 	;
{662} NC524O2 + HO2 = NC524OOH 		: KRO2HO2*0.706 	;
{663} NC524O2 + NO = 0.072 NC524NO3 +
                     0.928 NC524O + 
                     0.928 NO2 		: KRO2NO 	;
{664} NC524O2 + NO3 = NC524O + NO2 		: KRO2NO3 	;
{665} HMACR + NO3 = HMACO3 + HNO3 		: 3.40D-15 	;
{666} HMACR + O3 = CH2OOA + HOCH2COCHO 	: 7.80D-19*0.5 	;
{667} HMACR + O3 = HCHO + HMGLYOOA 		: 7.80D-19*0.5 	;
{668} HMACR + OH = 0.376 HMACO3 + 
                   0.624 HMACRO2 		: 4.83D-11 	;
{669} C525O2 + HO2 = C525OOH 			: KRO2HO2*0.706 	;
{670} C525O2 + NO = C525O + NO2 		: KRO2NO 	;
{671} C525O2 + NO3 = C525O + NO2 		: KRO2NO3 	;
{672} C58NO3CO3 + HO2 = 0.15 C58NO3CO2H + 
                        0.15 O3 + 
                        0.41 C58NO3CO3H + 
                        0.44 MACRNO3 + 
                        0.44 HO2 +
                        0.44 OH 		: KAPHO2 	;
{673} C58NO3CO3 + NO = MACRNO3 + 
                       HO2 + NO2 		: KAPNO 	;
{674} C58NO3CO3 + NO2 = C58NO3PAN 		: KFPAN 	;
{675} OH + IEB1OOH = HO12CO3C4 + OH + CO 	: 3.90D-11 	;
{676} IEB1O = MGLY + HOCH2CHO + HO2 	: KDEC 	;
{677} OH + IEB2OOH = MACROH + OH + CO 	: 5.34D-11 	;
{678} IEB2O = GLY  + ACETOL + HO2 		: KDEC 	;
{679} MACRNCO2H + OH = 0.75 ACETOL + 
                       0.75 NO2 + 
                       0.25 CONM2CO2H + 
                       0.25 HO2 		: 7.9D-13 	;
{680} CONM2CO3 + HO2 = 0.15 CONM2CO2H + 
                       0.15 O3 +
                       0.41 CONM2CO3H +
                       0.44 MGLY + 
                       0.44 NO2 + 0.44 OH : KAPHO2 	;
{681} CONM2CO3 + NO = MGLY + 2 NO2 		: KAPNO 	;
{682} CONM2CO3 + NO2 = CONM2PAN 		: KFPAN 	;
{683} CONM2CO3 + NO3 = MGLY + 2 NO2 	: KRO2NO3*1.74 	;
{684} MACRNCO3H + OH = 0.15 CONM2CO3H + 
                       0.15 HO2 + 
                       0.85 MACRNCO3	: 4.42D-12 	;
{685} MACRNPAN + OH = CONM2PAN + HO2 	: 8.21D-13 	;
{686} MACRNPAN = MACRNCO3 + NO2 		: KBPAN 	;
{687} OH + INB1HPCHO = INB1HPCO3 		: 2.41D-11 	;
{688} OH + INB1NACHO = INB1NACO3 		: 1.85D-11 	;
{689} OH + INB1NBCHO = INB1NBCO3 		: 1.85D-11 	;
{690} INB1GLYOX + OH = MACRNCO3 + CO 	: 1.35D-11 	;
{691} C57NO3 + OH = 0.54 C4M2ALOHNO3 + 
                    0.54 HO2
                    0.46 C57NO3CO3 		: 9.37D-12 	;
{692} OH + IEC1OOH = IEC1O2 			: 3.60D-12 	;
{693} OH + IEC1OOH = IEC2OOH + HO2 		: 1.57D-11 	;
{694} IEC1O = BIACETOH + HCHO + HO2 	: KDEC 	;
{695} IECCO3 + HO2 = 0.56 IECCO3H +	
                     0.44 MACRO2 + 
                     0.44 OH 			: KAPHO2 	;
{696} IECCO3 + NO = MACRO2 + NO2 		: KAPNO 	;
{697} IECCO3 + NO2 = IECPAN 			: KFPAN 	;
{698} IECCO3 + NO3 = MACRO2 + NO2 		: KRO2NO3*1.74 	;
{699} CO2N3CO3 + HO2 = 0.56 CO2N3CO3H +
                       0.44 MGLY + 
                       0.44 NO2 + 0.44 OH : KAPHO2*0.44 	;
{700} CO2N3CO3 + NO = MGLY + 2 NO2 		: KAPNO 	;
{701} CO2N3CO3 + NO2 = CO2N3PAN 		: KFPAN 	;
{702} CO2N3CO3 + NO3 = MGLY + 2 NO2 	: KRO2NO3*1.74 	;
{703} INDHPCHO + OH = INDHPCO3 		: 2.58D-11 	;
{704} INDHCHO + OH = INDHCO3 			: 2.27D-11 	;
{705} INCNCO3 + HO2 = 0.15 INCNCO2H + 
                      0.15 O3 +
                      0.41 INCNCO3H + 
                      0.44 MACRNB + 
                      0.44 NO2 + 0.44 OH 	: KAPHO2 	;
{706} INCNCO3 + NO = MACRNB + 2 NO2 	: KAPNO 	;
{707} INCNCO3 + NO2 = INCNPAN 		: KFPAN 	;
{708} INCNCO3 + NO3 = MACRNB + 2 NO2 	: KRO2NO3*1.74 	;
{709} OH + MACRNB = MACRNBCO3 		: 2.15D-11 	;
{710} MACRNBCO3 + HO2 = 0.15 MACRNBCO2H + 
                        0.15 O3 + 
                        0.41 MACRNBCO3H + 
                        0.44 NOA + 
                        0.44 HO2 + 
                        0.44 OH 		: KAPHO2 	;
{711} MACRNBCO3 + NO = NOA + HO2 + NO2 	: KAPNO 	;
{712} MACRNBCO3 + NO2 = MACRNBPAN 		: KFPAN 	;
{713} MACRNBCO3 + NO3 = NOA + HO2 + NO2 	: KRO2NO3*1.74 	;
{714} IPRHOCO2H + OH = AONE + HO2 		: 1.72D-12 	;
{715} OH + IPRHOCO3H = IPRHOCO3 		: 4.80D-12 	;
{716} C4PAN5 = IPRHOCO3 + NO2 		: KBPAN 	;
{717} OH + C4PAN5 = AONE + CO + NO2 	: 4.75D-13 	;
{718} PRNO3CO2H + OH = ALD2 + NO2 		: 3.14D-13 	;
{719} PRNO3CO3H + OH = PRNO3CO3 		: 3.77D-12 	;
{720} PRNO3PAN + OH = ALD2 + CO + 2 NO2 	: 1.43D-13 	;
{721} PRNO3PAN = PRNO3CO3 + NO2 		: KBPAN 	;
{722} IPROPOLPER + OH = CH3CHOHCO3 		: 9.34D-12 	;
{723} IPROPOLPAN + OH = ALD2 + CO + NO2 	: 2.34D-12 	;
{724} IPROPOLPAN = CH3CHOHCO3 + NO2 	: KBPAN 	;
{725} HOCH2COCO2H + OH = HOCH2CO3 		: 2.89D-12 	;
{726} H13CO2CO3 + HO2 = 0.56 H13CO2CO3H +
                        0.44 HOCH2COCHO + 
                        0.44 HO2 + 
                        0.44 OH 		: KAPHO2 	;
{727} H13CO2CO3 + NO = HOCH2COCHO + 
                       HO2 + NO2 		: KAPNO 	;
{728} H13CO2CO3 + NO2 = C4PAN10 		: KFPAN 	;
{729} H13CO2CO3 + NO3 = HOCH2COCHO + 
                        HO2 + NO2 		: KRO2NO3*1.74 	;
{730} H1CO23CHO + OH = 2 CO + HOCH2CO3 	: 1.44D-11 	;
{731} ACO3B + HO2 = HOCHOCOOH 		: KRO2HO2*0.52 	;
{732} ACO3B + NO = CHOCOHCO + NO2 		: KRO2NO 	;
{733} ACO3B + NO3 = CHOCOHCO + NO2 		: KRO2NO3 	;
{734} OCCOHCO2 + HO2 = 0.4 C32OH13CO + 
                        0.4 O3 + 
                        0.6 OCCOHCOOH 	: 0.52*KRO2HO2 	;
{735} OCCOHCO2 + NO = 0.05 C42AOH +
                       0.95 OCCOHCO + 
                       0.95 NO2 		: KRO2NO 	;
{736} OCCOHCO2 + NO3 = C42AOH + NO2 	: KRO2NO3 	;
{737} GLYOOB = 0.24 GLYOO + 0.2 HCHO + 
                0.76 HO2 + 0.92 CO +  
                0.36 OH 			: KDEC 	;
{738} OH + IEACO3H = IEACO3 			: 4.81D-12 	;
{739} IEAPAN = IEACO3 + NO2 			: KBPAN 	;
{740} OH + IEAPAN = HMVKBO2 + CO + NO2 	: 1.21D-12 	;
{741} INAHPCO3 + HO2 = 0.44 HMVKANO3 + 
                        0.88 OH + 
                        0.15 INAHPCO2H + 
                        0.15 O3 + 
                        0.41 INAHPCO3H 	: KAPHO2 	;
{742} INAHPCO3 + NO = HMVKANO3 + 
                       OH + NO2 		: KAPNO 	;
{743} INAHPCO3 + NO2 = INAHPPAN 		: KFPAN 	;
{744} INAHPCO3 + NO3 = HMVKANO3 + 
                        OH + NO2 		: KRO2NO3*1.74 	;
{745} INANCO3 + HO2 = 0.44 HMVKANO3 + 
                       0.44 NO2 + 
                       0.44 OH + 
                       0.15 INANCO2H + 
                       0.15 O3 + 
                       0.41 INANCO3H 	: KAPHO2 	;
{746} INANCO3 + NO = HMVKANO3 + 2 NO2 	: KAPNO 	;
{747} INANCO3 + NO2 = INANPAN 		: KFPAN 	;
{748} INANCO3 + NO3 = HMVKANO3 + 2 NO2 	: KRO2NO3*1.74 	;
{749} OH + INANCOCHO = INANCOCO3 		: 3.79D-12 	;
{750} INAHCO3 + HO2 = 0.44 HMVKANO3 + 
                       0.44 HO2 + 
                       0.44 OH + 
                       0.15 INAHCO2H + 
                       0.15 O3 + 
                       0.41 INAHCO3H 	: KAPHO2 	;
{751} INAHCO3 + NO = HMVKANO3 + 
                      HO2 + NO2 		: KAPNO 	;
{752} INAHCO3 + NO2 = INAHPAN 		: KFPAN 	;
{753} INAHCO3 + NO3 = HMVKANO3 + HO2 + 
                       NO2 			: KRO2NO3*1.74 	;
{754} HIEB1O2 + HO2 = HIEB1OOH 		: KRO2HO2*0.706 	;
{755} HIEB1O2 + NO = HIEB1O + NO2 		: KRO2NO 	;
{756} HIEB1O2 + NO3 = HIEB1O + NO2 		: KRO2NO3 	;
{757} HIEB2O2 + HO2 = HIEB2OOH 		: KRO2HO2*0.706 	;
{758} HIEB2O2 + NO = HIEB2O + NO2 		: KRO2NO 	;
{759} HIEB2O2 + NO3 = HIEB2O + NO2		: KRO2NO3 	;
{760} NC524OOH + OH = 0.728 HPNC524CO + 
                       0.728 HO2 + 	
                       0.272 NC524O2 	: 1.32D-11 	;
{761} NC524NO3 + OH = DNC524CO + HO2 	: 2.43D-12 	;
{762} NC524O = HMVKNO3 + HCHO + HO2 	: 1.80D+13*(TEMP/298)**1.7*EXP(-4079/TEMP) 	;
{763} NC524O = HOCH2CHO + NO2 + H13CO2C3	: 1.80D+13*(TEMP/298)**1.7*EXP(-4733/TEMP) 	;
{764} NC524OH + OH = HNC524CO + HO2 	: 9.60D-12 	;
{765} HMACO3 + HO2 = 0.15 HMACO2H + 
                      0.15 O3 + 
                      0.41 HMACO3H + 
                      0.44 HOCH2CO3 + 
                      0.44 HCHO + 0.44 OH : KAPHO2 	;
{766} HMACO3 + NO = HOCH2CO3 + 
                     HCHO + NO2 		: KAPNO 	;
{767} HMACO3 + NO2 = HMPAN 			: KFPAN 	;
{768} HMACO3 + NO3 = HOCH2CO3 + 
                      HCHO + NO2 		: KRO2NO3*1.74 	;
{769} CH2OOA = 0.37 CH2OO + 0.5 CO +
                0.13 HO2 + 0.13 CO + 
                0.13 OH 			: KDEC 	;
{770} HMGLYOOA = 0.18 HMGLYOO + 
                  0.82 HOCH2CO3 + 
                  0.82 CO + 0.82 HO2 	: KDEC 	;
{771} HMACRO2 + HO2 = HMACROOH 		: KRO2HO2*0.625 	;
{772} HMACRO2 + NO = HMACRO + NO2 		: KRO2NO 	;
{773} HMACRO2 + NO3 = HMACRO + NO2 		: KRO2NO3 	;
{774} C525OOH + OH = C525O2 			: 1.37D-11 	;
{775} C525O = HOCH2CO3 + H13CO2C3 		: KDEC 	;
{776} OH + C58NO3CO2H = MMALNACO2H + HO2	: 2.49D-12 	;
{777} OH + C58NO3CO3H = 0.68 C58NO3CO3 +
                         0.32 MMALNACO3H + 
                         0.32 HO2 		: 5.57D-12 	;
{778} C58NO3PAN = C58NO3CO3 + NO2 		: KBPAN 	;
{779} OH + C58NO3PAN = MMALNAPAN + HO2 	: 1.97D-12 	;
{780} CONM2CO2H + OH = MGLY + NO2 		: 3.70D-12 	;
{781} CONM2CO3H + OH = CONM2CO3 		: 6.78D-12 	;
{782} CONM2PAN + OH = MGLY + NO2 + NO3 	: 3.18D-12 	;
{783} CONM2PAN = CONM2CO3 + NO2 		: KBPAN 	;
{784} INB1HPCO3 + HO2 = 0.15 INB1HPCO2H + 
                         0.15 O3 + 
                         0.41 INB1HPCO3H + 
                         0.44 MACRNO3 + 
                         0.88 OH  		: KAPHO2 	;
{785} INB1HPCO3 + NO = MACRNO3 + 
                        OH + NO2 		: KAPNO 	;
{786} INB1HPCO3 + NO2 = INB1HPPAN 		: KFPAN 	;
{787} INB1HPCO3 + NO3 = MACRNO3 + 
                         OH + NO2 		: KRO2NO3*1.74 	;
{788} INB1NACO3 + HO2 = 0.15 INB1NACO2H + 
                         0.15 O3 +
                         0.41 INB1NACO3H +
                         0.44 MACRNO3 + 
                         0.44 NO2 + 
                         0.44 OH 		: KAPHO2 	;
{789} INB1NACO3 + NO = MACRNO3 + 2 NO2 	: KAPNO 	;
{790} INB1NACO3 + NO2 = INB1NAPAN 		: KFPAN 	;
{791} INB1NACO3 + NO3 = MACRNO3 + 2 NO2 	: KRO2NO3*1.74 	;
{792} INB1NBCO3 + HO2 = 0.15 INB1NBCO2H + 
                         0.15 O3 + 
                         0.41 INB1NBCO3H +
                         0.44 MVKNO3 + 
                         0.44 NO2 + 
                         0.44 OH 		: KAPHO2 	;
{793} INB1NBCO3 + NO = MVKNO3 + 2 NO2 	: KAPNO 	;
{794} INB1NBCO3 + NO2 = INB1NBPAN 		: KFPAN 	;
{795} INB1NBCO3 + NO3 = MVKNO3 + 2 NO2 	: KRO2NO3*1.74 	;
{796} C4M2ALOHNO3 + OH = 0.86 MMALNACO3 + 
                          0.14 MMALNBCO3 	: 2.53D-11 	;
{797} C57NO3CO3 + HO2 = 0.15 C57NO3CO2H +
                         0.15 O3 +
                         0.41 C57NO3CO3H +
                         0.44 HO12CO3C4 + 
                         0.44 NO2 + 
                         0.44 OH 		: KAPHO2 	;
{798} C57NO3CO3 + NO = HO12CO3C4 + 
                        2 NO2 		: KAPNO 	;
{799} C57NO3CO3 + NO2 = C57NO3PAN 		: KFPAN 	;
{800} C57NO3CO3 + NO3 = HO12CO3C4 + 
                         2 NO2 		: KRO2NO3*1.74 	;
{801} C57NO3CO3 + NO3 = MACRNO3 + HO2 + 
                         NO2 			: KRO2NO3*1.74 	;
{802} OH + IEC2OOH = BIACETOH + OH + CO 	: 2.73D-11 	;
{803} OH + IECCO3H = IECCO3 			: 1.04D-11 	;
{804} IECPAN = IECCO3 + NO2 			: KBPAN 	;
{805} OH + IECPAN = MACRO2 + CO + NO2 	: 6.80D-12 	;
{806} CO2N3CO3H + OH = CO2N3CO3 		: 4.11D-12 	;
{807} CO2N3PAN + OH = CO2N3CO3 		: 5.10D-13 	;
{808} CO2N3PAN = CO2N3CO3 + NO2 		: KBPAN 	;
{809} INDHPCO3 + HO2 = 0.56 INDHPCO3H + 
                        0.44 MVKNO3 + 
                        0.88 OH 		: KAPHO2	;
{810} INDHPCO3 + NO = MVKNO3 + OH + NO2 	: KAPNO 	;
{811} INDHPCO3 + NO2 = INDHPPAN 		: KFPAN 	;
{812} INDHPCO3 + NO3 = MVKNO3 + OH + NO2	: KRO2NO3*1.74 	;
{813} INDHCO3 + HO2 = 0.56 INDHCO3H + 
                       0.44 MVKNO3 + 
                       0.44 HO2 + 0.44 OH	: KAPHO2 	;
{814} INDHCO3 + NO = MVKNO3 + HO2 + NO2 	: KAPNO 	;
{815} INDHCO3 + NO2 = INDHPAN 		: KFPAN 	;
{816} INDHCO3 + NO3 = MVKNO3 + HO2 + NO2	: KRO2NO3*1.74 	;
{817} OH + INCNCO2H = MACRNB + NO2 		: 1.66D-12 	;
{818} OH + INCNCO3H = INCNCO3 		: 4.74D-12 	;
{819} INCNPAN = INCNCO3 + NO2 		: KBPAN 	;
{820} OH + INCNPAN = MACRNB + NO2 + NO3 	: 1.14D-12 	;
{821} MACRNBCO2H + OH = COHM2CO2H + NO2 	: 1.23D-12*0.32 	;
{822} MACRNBCO2H + OH = NOA + HO2 		: 1.23D-12*0.68 	;
{823} MACRNBCO3H + OH = COHM2CO3H + NO2 	: 4.31D-12*0.09 	;
{824} MACRNBCO3H + OH = MACRNBCO3 		: 4.31D-12*0.91 	;
{825} MACRNBPAN + OH = COHM2PAN + NO2 	: 7.10D-13 	;
{826} MACRNBPAN = MACRNBCO3 + NO2 		: KBPAN 	;
{827} HYPERACET + OH = ANO2 			: 1.90D-12*EXP(190/TEMP) 	;
{828} HYPERACET + OH = MGLY + OH 		: 8.39D-12 	;
{829} OH + H13CO2CO3H = H13CO2CO3 		: 9.43D-12 	;
{830} C4PAN10 = H13CO2CO3 + NO2 		: KBPAN 	;
{831} OH + C4PAN10 = HOCH2COCHO + 
                     CO + NO2 		: 5.83D-12 	;
{832} HOCHOCOOH + OH = HOCH2COCHO + OH 	: 4.77D-11 	;
{833} CHOCOHCO = HOCH2CHO + HO2 + CO 	: KDEC 	;
{834} OCCOHCOH + OH = A2PANOO 		: 6.22D-11 	;
{835} C32OH13CO + OH = HCOCOHCO3 		: 1.36D-10 	;
{836} OCCOHCOOH + OH = OCCOHCO2 		: 9.258E-11 	;
{837} OH + C42AOH = NMGLY + HO2 		: 2.92D-11 	;
{838} OCCOHCO = HCHO + GLY  + HO2 		: KDEC 	;
{839} ETHENO3O2 + HO2 = ETHO2HNO3 		: KRO2HO2*0.387 	;
{840} ETHENO3O2 + NO = ETHENO3O + NO2 	: KRO2NO 	;
{841} ETHENO3O2 + NO3 = ETHENO3O + NO2 	: KRO2NO3 	;
{842} HOCH2CH2O2 + HO2 = HYETHO2H 		: 1.53D-13*EXP(1300/TEMP) 	;
{843} HOCH2CH2O2 + NO = ETHOHNO3 		: KRO2NO*0.005 	;
{848} HOCH2CH2O2 + NO = HOCH2CH2O + NO2 	: KRO2NO*0.995 	;
{849} HOCH2CH2O2 + NO3 = HOCH2CH2O + NO2	: KRO2NO3 	;
{850} OH + INAHPCO2H = HMVKANO3 + OH 	: 6.50D-12 	;
{851} OH + INAHPCO3H = INAHPCO3 		: 9.58D-12 	;
{852} INAHPPAN = INAHPCO3 + NO2 		: KBPAN 	;
{853} OH + INAHPPAN = HMVKANO3 + OH + 
                       CO + NO2 		: 5.98D-12 	;
{854} OH + INANCO2H = HMVKANO3 + NO2 	: 1.36D-12 	;
{855} OH + INANCO3H = INANCO3 		: 4.08D-12 	;
{856} INANPAN = INANCO3 + NO2 		: KBPAN 	;
{857} OH + INANPAN = HMVKANO3 + NO2 + 
                      CO + NO2 		: 4.85D-13 	;
{858} INANCOCO3 + HO2 = 0.15 INANCOCO2H + 
                         0.15 O3 +
                         0.41 INANCOCO3H +
                         0.44 CO23C4NO3 +
                         0.44 NO2 + 
                         0.44 OH 		: KAPHO2 	;
{859} INANCOCO3 + NO = CO23C4NO3 + 2 NO2	: KAPNO 	;
{860} INANCOCO3 + NO2 = INANCOPAN 		: KFPAN 	;
{861} INANCOCO3 + NO3 = CO23C4NO3 +
                         2 NO2		: KRO2NO3*1.74 	;
{862} OH + INAHCO2H = HMVKANO3 + HO2 	: 3.04D-12 	;
{863} OH + INAHCO3H = INAHCO3 		: 6.12D-12 	;
{864} INAHPAN = INAHPCO3 + NO2 		: KBPAN 	;
{865} OH + INAHPAN = HMVKANO3 + HO2 + 
                      CO + NO2 		: 2.52D-12 	;
{866} HIEB1OOH + OH = MVKOHAOH + CO + OH	: 4.30D-11 	;
{867} HIEB1O = HOCH2COCHO + HOCH2CHO + 
                HO2 				: KDEC 	;
{868} HIEB2OOH + OH = HMACROH + CO + OH 	: 5.74D-11 	;
{869} HIEB2O = H13CO2C3 + GLY  + HO2 	: KDEC 	;
{870} HPNC524CO + OH = HMVKNO3 + CO + OH	: 2.98D-11 	;
{871} DNC524CO + OH = HMVKNO3 + CO + NO2	: 1.93D-11 	;
{872} HMVKNO3 + OH = HMVKNGLYOX + HO2 	: 3.85D-12 	;
{873} H13CO2C3 + OH = HOCH2COCHO + HO2 	: 5.25D-12 	;
{874} HNC524CO + OH = HMVKNO3 + CO + HO2	: 2.67D-11 	;
{875} HMACO2H + OH = HOCH2CO3 + HCHO 	: 1.84D-11 	;
{876} HMACO3H + OH = HMACO3 			: 1.99D-11 	;
{877} HMPAN + OH = H13CO2C3 + HCHO + 
                    CO + NO2 			: 2.90D-11 	;
{878} HMPAN = HMACO3 + NO2 			: KBPAN 	;
{879} HMGLYOO + CO = HOCH2COCHO 		: 1.20D-15 	;
{880} HMGLYOO + NO = HOCH2COCHO + NO2 	: 1.00D-14 	;
{881} HMGLYOO + NO2 = HOCH2COCHO + NO3 	: 1.00D-15 	;
{882} HMGLYOO + SO2 = HOCH2COCHO + H2SO4	: 7.00D-14 	;
{883} HMGLYOO = HOCH2COCHO + H2O2 		: 6.00D-18*H2O 	;
{884} HMACROOH + OH = HMACRO2 		: 4.17D-11 	;
{885} HMACRO = H13CO2C3 + CO + HO2 		: KDEC 	;
{886} HMACROH + OH = HMACRO 			: 3.82D-11 	;
{887} MMALNACO2H + OH = CONM2CHO + HO2 	: 4.93D-12 	;
{888} MMALNACO3H + OH = MMALNACO3 		: 8.01D-12 	;
{889} MMALNAPAN + OH = CONM2CHO + HO2 + 
                        NO3 			: 4.41D-12 	;
{890} MMALNAPAN = MMALNACO3 + NO2 		: KBPAN 	;
{891} CH3COCO3H + OH = CH3COCO3 		: 3.69D-12 	;
{892} CH3COPAN + OH = HCHO + CO + 
                       CO + NO2 		: 1.02D-13 	;
{893} CH3COPAN = CH3COCO3 + NO2 		: KBPAN 	;
{894} INB1HPCO2H + OH = MACRNO3 + OH 	: 7.40D-12 	;
{895} INB1HPCO3H + OH = INB1HPCO3 		: 1.09D-11 	;
{896} INB1HPPAN + OH = MACRNO3 + CO + 
                        OH + NO2 		: 7.26D-12 	;
{897} INB1HPPAN = INB1HPCO3 + NO2 		: KBPAN 	;
{898} INB1NACO2H + OH = MACRNO3 + NO2 	: 1.72D-12 	;
{899} INB1NACO3H + OH = INB1NACO3 		: 5.18D-12 	;
{900} INB1NAPAN + OH = MACRNO3 + CO + 
                        2 NO2 		: 1.58D-12 	;
{901} INB1NAPAN = INB1NACO3 + NO2 		: KBPAN 	;
{902} INB1NBCO2H + OH = MVKNO3 + NO2 	: 1.73D-12 	;
{903} INB1NBCO3H + OH = INB1NBCO3 		: 5.19D-12 	;
{904} INB1NBPAN + OH = MVKNO3 + CO + 
                        2 NO2 		: 1.59D-12 	;
{905} INB1NBPAN = INB1NBCO3 + NO2 		: KBPAN 	;
{906} MMALNACO3 + HO2 = 0.44 CONM2CHO + 
                         0.44 HO2 + 
                         0.44 OH +
                         0.15 MMALNACO2H + 
                         0.15 O3 +
                         0.41 MMALNACO3H 	: KAPHO2 	;
{907} MMALNACO3 + NO = CONM2CHO + HO2 + 
                        NO2 			: KAPNO 	;
{908} MMALNACO3 + NO2 = MMALNAPAN 		: KFPAN 	;
{909} MMALNACO3 + NO3 = CONM2CHO + HO2 +
                         NO2 			: KRO2NO3*1.74 	;
{910} MMALNBCO3 + HO2 = 0.44 CO2H3CHO + 
                         0.44 NO2 + 
                         0.44 OH +
                         0.15 MMALNBCO2H + 
                         0.15 O3 + 
                         0.41 MMALNBCO3H 	: KAPHO2 	;
{911} MMALNBCO3 + NO = CO2H3CHO + 2 NO2	: KAPNO 	;
{912} MMALNBCO3 + NO2 = MMALNBPAN 		: KFPAN 	;
{913} MMALNBCO3 + NO3 = CO2H3CHO + 2 NO2	: KRO2NO3*1.74 	;
{914} C57NO3CO2H + OH = MMALNBCO2H + HO2	: 6.52D-12 	;
{915} C57NO3CO3H + OH = 0.39 C57NO3CO3 +
                         0.61 MMALNBCO3H + 
                         0.61 HO2         : 6.52D-12 	;
{916} C57NO3PAN + OH = MMALNBPAN + HO2 	: 6.00D-12 	;
{917} C57NO3PAN = C57NO3CO3 + NO2 		: KBPAN 	;
{918} INDHPCO3H + OH = INDHPCO3 		: 8.64D-12 	;
{919} INDHPPAN + OH = MVKNO3 + NO3 		: 5.04D-12 	;
{920} INDHPPAN = INDHPCO3 + NO2 		: KBPAN 	;
{921} INDHCO3H + OH = INDHCO3 		: 5.66D-12 	;
{922} INDHPAN + OH = MVKNO3 + NO3 		: 1.96D-12 	;
{923} INDHPAN = INDHCO3 + NO2 		: KBPAN 	;
{924} COHM2CO2H + OH = GLY  + HO2 		: 2.16D-11 	;
{925} COHM2CO3H + OH = COHM2CO3 		: 2.47D-11 	;
{926} COHM2PAN + OH = GLY  + NO3 		: 2.11D-11 	;
{927} COHM2PAN = COHM2CO3 + NO2 		: KBPAN 	;
{928} A2PANOO + HO2 = 0.44 A2PANO + 
                       0.44 OH + 
                       0.15 C2OHOCO2H + 
                       0.15 O3 +
                       0.41 C2OHOCOOH 	: KAPHO2 	;
{929} A2PANOO + NO = A2PANO + NO2 		: KAPNO 	;
{930} A2PANOO + NO2 = A2PAN 			: KFPAN 	;
{931} A2PANOO + NO3 =  A2PANO + NO2 	: KRO2NO3*1.74 	;
{932} HCOCOHCO3 + HO2 = 0.44 GLY + 
                         0.44 HO2 + 
                         0.44 OH +
                         0.56 HCOCOHCO3H	: KAPHO2 	;
{933} HCOCOHCO3 + NO = GLY  + HO2 + NO2 	: KAPNO 	;
{934} HCOCOHCO3 + NO2 = HCOCOHPAN 		: KFPAN 	;
{935} HCOCOHCO3 + NO3 = GLY  + HO2 + NO2	: KRO2NO3*1.74 	;
{936} NO3 + NMGLY = NO3CH2CO3 + CO + 
                     HNO3 			: KNO3AL*2.4 	;
{937} OH + NMGLY = NO3CH2CO3 + CO 		: 1.24D-11 	;
{938} ETHO2HNO3 + OH = ETHENO3O2 		: 1.90D-12*EXP(190/TEMP) 	;
{939} ETHO2HNO3 + OH = NO3CH2CHO + OH 	: 1.62D-12 	;
{940} ETHENO3O = NO2 + HCHO + HCHO 		: 7.00D+03 	;
{941} ETHENO3O = NO3CH2CHO + HO2 		: KROPRIM*O2 	;
{942} ETHOHNO3 + OH = HOCH2CHO + NO2 	: 8.40D-13 	;
{943} HYETHO2H + OH = HOCH2CH2O2 		: 1.90D-12*EXP(190/TEMP) 	;
{944} HYETHO2H + OH = HOCH2CHO + OH 	: 1.38D-11 	;
{945} HOCH2CH2O = HO2 + HCHO + HCHO 	: 9.50D+13*EXP(-5988/TEMP) 	;
{946} HOCH2CH2O = HO2 + HOCH2CHO 		: KROPRIM*O2 	;
{947} ETHGLY + OH = HOCH2CHO + HO2 		: 1.45D-11 	;
{948} OH + INANCOCO2H = NO2 + CO23C4NO3 	: 9.35D-13 	;
{949} OH + INANCOCO3H = INANCOCO3 		: 4.02D-12 	;
{950} CO23C4NO3 + OH = CO23C3CHO + NO2 	: 1.30D-13 	;
{951} INANCOPAN = INANCOCO3 + NO2 		: KBPAN 	;
{952} OH + INANCOPAN = NO3 + CO23C4NO3 	: 4.15D-13 	;
{953} HMVKNGLYOX + OH = CO + CO + 
                         HOCH2CHO + NO2 	: 1.36D-11 	;
{954} CH3COCO3 + HO2 = 0.44 C2O3 + 
                        0.44 OH +
                        0.56 CH3COCO3H	: KAPHO2 	;
{955} CH3COCO3 + NO = C2O3 + NO2 		: KAPNO 	;
{956} CH3COCO3 + NO2 = CH3COPAN 		: KFPAN 	;
{957} CH3COCO3 + NO3 = C2O3 + NO2 		: KRO2NO3*1.74 	;
{958} MMALNBCO2H + OH = CO2H3CHO + NO2 	: 2.23D-11 	;
{959} MMALNBCO3H + OH = MMALNBCO3 		: 2.59D-11 	;
                   NO2 + NO3 			: 2.18D-11 	;
{960} MMALNBPAN = MMALNBCO3 + NO2 		: KBPAN 	;
{961} COHM2CO3 + HO2 = 0.15 COHM2CO2H + 
                        0.15 O3 +	 	
                        0.41 COHM2CO3H +	
                        0.44 GLY + 
                        0.44 HO2 + 
                        0.44 OH 		: KAPHO2 	;
{962} COHM2CO3 + NO = GLY  + HO2 + NO2 	: KAPNO 	;
{963} COHM2CO3 + NO2 = COHM2PAN 		: KFPAN 	;
{964} COHM2CO3 + NO3 = GLY  + HO2 + NO2 	: KRO2NO3*1.74 	;
{965} A2PANO = HOCH2CHO + HO2 		: KDEC 	;
{966} C2OHOCO2H + OH = C3DIOLO2 		: 1.867E-11 	;
{967} C2OHOCOOH + OH = A2PANOO 		: 1.513E-11 	;
{968} A2PAN + OH = HOCH2CHO + CO + NO2 	: 1.865E-11 	;
{969} A2PAN = A2PANOO + NO2 			: KBPAN 	;
{970} HCOCOHCO3H + OH = HCOCOHCO3 		: 7.33D-11 	;
{971} HCOCOHPAN + OH = GLY  + CO + NO2 	: 6.97D-11 	;
{972} HCOCOHPAN = HCOCOHCO3 + NO2 		: KBPAN 	;
{973} BIACETO = C2O3 + HCHO + CO 		: KDEC 	;
{974} C3DIOLO2 + HO2 = C3DIOLOOH 		: KRO2HO2*0.52 	;
{975} C3DIOLO2 + NO = C3DIOLO + NO2 	: KRO2NO 	;
{976} C3DIOLO2 + NO3 = C3DIOLO + NO2 	: KRO2NO3 	;
{977} C3DIOLOOH + OH = C3DIOLO2 		: 2.78D-11 	;
{978} C3DIOLO = HOCH2CHO + HCHO + HO2 	: KDEC 	;
#
# RO2 + RO2 reactions
#
{14} NISOPO2 = 0.2 ISOPCNO3 +
               0.2 NC4CHO +
               0.6 NISOPO 			: 1.30D-12*RO2 	;
{50} ISOPAO2 = 0.1 HC4ACHO + 
               0.8 ISOPAO + 
               0.1 ISOPAOH 			: 2.40D-12*RO2 	;
{57} ISOPBO2 = 0.8 ISOPBO + 0.2 ISOPBOH 	: 8.00D-13*RO2 	;
{63} ISOPCO2 = 0.1 HC4CCHO + 
               0.1 ISOPAOH + 
               0.8 ISOPCO 			: 2.00D-12*RO2 	;
{70} ISOPDO2 = 0.1 HCOC5 + 0.8 ISOPDO +
               0.1 ISOPDOH 			: 2.90D-12*RO2 	;
{102} MACO3 = 0.7 CH3C2H2O2 +	
              0.3 MACO2H 			: 1.00D-11*RO2 	;
{115} MACRO2 = 0.7 MACRO + 0.3 MACROH 	: 9.20D-14*RO2 	;
{121} MACROHO2 = 0.2 C3MDIALOH + 
                 0.2 MACROH +
                 0.6 MACROHO 			: 1.4D-12*RO2 	;
{149} HMVKAO2 = 0.2 CO2H3CHO + 
                0.6 HMVKAO + 
                0.2 HO12CO3C4 		: 2.00D-12*RO2 	;
{155} HMVKBO2 = 0.2 BIACETOH + 
                0.6 HMVKBO + 
                0.2 HO12CO3C4 		: 8.80D-13*RO2 	;
{181} MVKO2 = 0.6 MVKO + 0.2 MVKOH + 
              0.2 VGLYOX 			: 2.00D-12*RO2 	;
{253} INCO2 = 0.1 INCCO + 
              0.8 INCO + 
              0.1 INCOH 			: 2.90D-12*RO2 	;
{262} NC4CO3 = 0.3 NC4CO2H + 
               0.7 NOA + 0.7 HO2 + 0.7 CO : 1.00D-11*RO2 	;
{278} C510O2 = 0.7 C510O + 
               0.3 C510OH 			: 9.20D-14*RO2 	;
{310} PRONO3AO2 = 0.2 CHOPRNO3 + 
                  0.6 PRONO3AO + 
                  0.2 PROPOLNO3 		: 6.00D-13*RO2 	;
{316} PRONO3BO2 = 0.2 NOA + 0.2 PROLNO3 + 
                  0.6 PRONO3BO 		: 4.00D-14*RO2 	;
{329} HYPROPO2 = 0.2 ACETOL +
                 0.6 HYPROPO + 
                 0.2 PROPGLY 			: 8.80D-13*RO2 	;
{338} IPROPOLO2 = 0.2 CH3CHOHCHO + 
                  0.6 IPROPOLO + 
                  0.2 PROPGLY 		: 2.00D-12*RO2 	;
{403} INAO2 = 0.8 INAO + 0.2 INAOH 		: 8.00D-13*RO2 	;
{409} C524O2 = 0.1 C524CO + 
               0.8 C524O + 
               0.1 C524OH 			: 2.90D-12*RO2 	;
{418} HC4ACO3 = 0.7 ACETOL + 0.7 HO2 + 
                0.7 CO + 0.3 HC4ACO2H 	: 1.00D-11*RO2 	;
{424} C58O2 = 0.7 C58O + 0.3 C58OH 		: 9.20D-14*RO2 	;
{439} INB1O2 = 0.1 INB1CO + 
               0.8 INB1O + 
               0.1 INB1OH 			: 2.90D-12*RO2 	;
{445} INB2O2 = 0.2 C58NO3 +
               0.2 INB1OH +
               0.6 INB2O 			: 8.80D-13*RO2 	;
{456} HC4CCO3 = 0.7 C2O3 + 0.7 HOCH2CHO +	
                0.3 HC4CCO2H 			: 1.00D-11*RO2 	;
{463} C57O2 = 0.7 C57O + 0.3 C57OH 		: 9.20D-14*RO2 	;
{476} INDO2 = 0.8 INDO + 0.2 INDOH		: 8.00D-13*RO2 	;
{484} HOCH2CO3 = 0.7 HCHO + 0.7 HO2 + 
                 0.3 HOCH2CO2H 		: 1.00D-11*RO2 	;
{489} C59O2 = C59O 				: 9.20D-14*RO2 	;
{533} HCOCO3 = 0.7 CO + 0.7 HO2 + 
               0.3 HCOCO2H 			: 1.00D-11*RO2 	;
{546} CHOMOHCO3 = MGLY + HO2 			: 1.00D-11*RO2 	;
{566} HCOCH2O2 = 0.2 GLY + 0.6 HCOCH2O + 
                 0.2 HOCH2CHO 		: 2.00D-12*RO2 	;
{597} CO2H3CO3 = MGLY + HO2 			: 1.00D-11*RO2 	;
{608} ACO3 = 0.3 ACO2H + 
             0.7 HO2 + 0.7 CO + 0.7 HCHO 	: 1.00D-11*RO2 	;
{621} MVKOHAO2 = 0.2 H13CO2CHO + 
                 0.6 MVKOHAO + 
                 0.2 MVKOHAOH 		: 2.00D-12*RO2 	;
{627} MVKOHBO2 = 0.2 H14CO23C4 + 
                 0.2 MVKOHAOH + 
                 0.6 MVKOHBO 			: 8.80D-13*RO2 	;
{651} NO3CH2CO3 = 0.7 HCHO + 0.7 NO2 + 
                  0.3 NO3CH2CO2H 		: 1.00D-11*RO2 	;
{688} IEB1O2 = IEB1O 				: 9.20D-14*RO2 	;
{692} IEB2O2 = IEB2O 				: 8.80D-13*RO2 	;
{707} MACRNCO3 = 0.7 ACETOL + 0.7 NO2 + 
                 0.3 MACRNCO2H 		: 1.00D-11*RO2 	;
{727} IEC1O2 = IEC1O 				: 9.20D-14*RO2 	;
{779} IPRHOCO3 = 0.7 AONE + 0.7 HO2 + 
                 0.3 IPRHOCO2H 		: 1.00D-11*RO2 	;
{793} PRNO3CO3 = 0.7 ALD2 + 0.7 NO2 + 
                 0.3 PRNO3CO2H 		: 1.00D-11*RO2 	;
{805} CH3CHOHCO3 = ALD2 + HO2 		: 1.00D-11*RO2 	;
{851} IEACO3 = HMVKBO2 				: 1.00D-11*RO2 	;
{872} NC524O2 = 0.8 NC524O + 
                0.2 NC524OH 			: 8.00D-13*RO2 	;
{882} C525O2 = C525O 				: 9.20D-14*RO2 	;
{888} C58NO3CO3 = 0.3 C58NO3CO2H +
                  0.7 MACRNO3 + 0.7 HO2 	: 1.00D-11*RO2 	;
{906} CONM2CO3 = 0.3 CONM2CO2H + 
                 0.7 MGLY + 0.7 NO2 	: 1.00D-11*RO2 	;
{934} IECCO3 = MACRO2 				: 1.00D-11*RO2 	;
{940} CO2N3CO3 = MGLY + NO2 			: 1.00D-11*RO2 	;
{951} INCNCO3 = 0.3 INCNCO2H +
                0.7 MACRNB + 0.7 NO2 	: 1.00D-11*RO2 	;
{961} MACRNBCO3 = 0.3 MACRNBCO2H +	
                  0.7 NOA + 0.7 HO2 	: 1.00D-11*RO2 	;
{991} H13CO2CO3 = HOCH2COCHO + HO2 		: 1.00D-11*RO2 	;
{997} ACO3B = 0.6 CHOCOHCO + 	
              0.2 HOCH2COCHO + 
              0.2 OCCOHCOH 			: 8.8D-13*RO2 	;
{1005} OCCOHCO2 = 0.2 C32OH13CO + 
                  0.6 OCCOHCO + 
                  0.2 OCCOHCOH 		: 2.0D-12*RO2 	;
{1025} INAHPCO3 = 0.7 HMVKANO3 + 0.7 OH +
                  0.3 INAHPCO2H 		: 1.00D-11*RO2 	;
{1033} INANCO3 = 0.7 HMVKANO3 + 0.7 NO2 + 
                 0.3 INANCO2H 		: 1.00D-11*RO2 	;
{1042} INAHCO3 = 0.7 HMVKANO3 + 0.7 HO2 +	
                 0.3 INAHCO2H 		: 1.00D-11*RO2 	;
{1047} HIEB1O2 = HIEB1O 			: 9.20D-14*RO2 	;
{1051} HIEB2O2 = HIEB2O 			: 8.80D-13*RO2 	;
{1065} HMACO3 = 0.3 HMACO2H +
                0.7 HOCH2CO3 + 0.7 HCHO 	: 1.00D-11*RO2 	;
{1075} HMACRO2 = 0.7 HMACRO +	
                 0.3 HMACROH 			: 9.20D-14*RO2 	;
{1100} INB1HPCO3 = 0.3 INB1HPCO2H +	
                   0.7 MACRNO3 + 0.7 OH 	: 1.00D-11*RO2 	;
{1108} INB1NACO3 = 0.3 INB1NACO2H + 
                   0.7 MACRNO3 + 0.7 NO2 	: 1.00D-11*RO2 	;
{1116} INB1NBCO3 = 0.3 INB1NBCO2H + 
                   0.7 MVKNO3 + 0.7 NO2 	: 1.00D-11*RO2 	;
{1129} C57NO3CO3 = 0.3 C57NO3CO2H +
                   0.7 HO12CO3C4 + 
                   0.7 NO2 			: 1.00D-11*RO2 	;
{1148} INDHPCO3 = MVKNO3 + OH 		: 1.00D-11*RO2 	;
{1154} INDHCO3 = MVKNO3 + HO2 		: 1.00D-11*RO2 	;
{1189} ETHENO3O2 = 0.6 ETHENO3O + 
                   0.2 ETHOHNO3 + 
                   0.2 NO3CH2CHO 		: 6.00D-13*RO2 	;
{1196} HOCH2CH2O2 = 0.2 ETHGLY + 
                    0.6 HOCH2CH2O +
                    0.2 HOCH2CHO 		: 2*(KCH3O2*7.8D-14*EXP(1000/TEMP))**0.5*RO2 	;
{1217} INANCOCO3 = 0.3 INANCOCO2H + 
                   0.7 NO2 + 
                   0.7 CO23C4NO3 		: 1.00D-11*RO2 	;
{1291} MMALNACO3 = 0.7 CONM2CHO + 
                   0.7 HO2 +
                   0.3 MMALNACO2H 		: 1.00D-11*RO2 	;
{1299} MMALNBCO3 = 0.7 CO2H3CHO + 
                   0.7 NO2 + 
                   0.3 MMALNBCO2H 		: 1.00D-11*RO2 	;
{1330} A2PANOO = 0.7 A2PANO + 
                 0.3 C2OHOCO2H 		: 1.00E-11*RO2 	;
{1337} HCOCOHCO3 = GLY  + HO2 		: 1.00D-11*RO2 	;
{1373} CH3COCO3 = C2O3 				: 1.00D-11*RO2 	;
{1388} COHM2CO3 = 0.3 COHM2CO2H +
                  0.7 GLY  + 0.7 HO2 	: 1.00D-11*RO2 	;
{1404} C3DIOLO2 = C3DIOLO 			: 2.00D-12*RO2 	;
*end*

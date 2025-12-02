/*=======================================

* Author: Zaria Roller
* Date: 12/5/25

* Extension of "Childcare, Labor Supply, and 
Business Development: Experimental Evidence 
from Uganda" 

* This file runs the regressions for my extension

========================================*/

local Attrition_bounds = 1

*set globals
global T T_2 T_3 T_4
global strat1 "district_? strat_younger_sibling strat_school occup_? strat_rel_tchild"

global clv	", cl(village_id)"
global rob ", robust"

**************************
* HETEROGENEITY ANALYSIS *
**************************

use "$constructed/baseline_outcomes.dta", clear

* create indicator for target child having an older sister
gen has_older_sis = sib_nb_female!=0 & sib_nb_female!=.
label var has_older_sis "Older sister"

* set global with heterogneous variables
global het_new "has_older_sis muslim"

foreach z of global het_new {
	foreach t in 2 3 4 {
		gen T`t'X`z' = (T_`t'==1 & `z'==1)
	}
}

save "$constructed/baseline_outcomes_cleaned.dta", replace

foreach z of global het_new  { 
	
	use "$constructed/baseline_outcomes_cleaned.dta", clear
	
	*intensive margin
		foreach y of var Y1_A_coo Y1_C Y1_B_coo Y2_A_pro_w Y2_C_w Y2_B_pro_w Y1_A_coo_pt Y1_C_pt Y1_B_coo_pt Y2_A_pro_pt_w Y2_C_pt_w Y2_B_pro_pt_w  {

		*regressions
		reg `y' $T `z' T?X`z' $strat1 `y'_bl_r `y'_bl_miss if wave==4 | wave==5 $rob
		estimates store `z'_`y'						
		
		*p-values
		su `y' if T_1==1 & e(sample)==1 
		estadd scalar be=r(mean)
		su `y' if T_1==1 & e(sample)==1 & `z'==1
		estadd scalar be1=r(mean)
		lincom T_2 + T2X`z'
		estadd scalar c2=r(estimate)
		estadd scalar se2=r(se)
		estadd scalar pe2=r(p)
		lincom T_3 + T3X`z'
		estadd scalar c3=r(estimate)
		estadd scalar se3=r(se)
		estadd scalar pe3=r(p)
		lincom T_4+T4X`z'
		estadd scalar c4=r(estimate)
		estadd scalar se4=r(se)
		estadd scalar pe4=r(p)
		test T_2=T_3
		estadd scalar p23=r(p)
		test T_2=T_4
		estadd scalar p24=r(p)
		test T_3=T_4
		estadd scalar p34=r(p)
		test T_2+T_3=T_4
		estadd scalar p234=r(p)
		test T_2 + T2X`z'=T_3 + T3X`z'
		estadd scalar p23b=r(p)
		test T_2 + T2X`z'=T_4 + T4X`z'
		estadd scalar p24b=r(p)
		test T_3 + T3X`z'=T_4 + T4X`z'
		estadd scalar p34b=r(p)
		test T_2 + T2X`z' + T_3 + T3X`z'=T_4 + T4X`z'
		estadd scalar p234b=r(p)	
		su `z' if e(sample)==1 
		estadd scalar be2=r(mean)
		
		parmest, label saving(`z'_`y', replace) stars(0.1 0.05 0.01) idstr("`y'_x_`z'") escal(be be1  c2 se2 pe2 c3 se3 pe3 c4 se4 pe4 p23b p24b p34b p234b N be2) idnum(1) ylab
		}
	
	
	*extensive margin

		foreach y of var Y1_A_coo_D Y1_C_D Y1_B_coo_D Y1_A_coo_pt_D Y1_C_pt_D Y1_B_coo_pt_D {		
		
		*regressions
		reg `y' $T `z' T?X`z' $strat1 `y'_bl_r `y'_bl_miss if wave==4 | wave==5 $rob
		
		*p-values
		su `y' if T_1==1 & e(sample)==1 
		estadd scalar be=r(mean)
		su `y' if T_1==1 & e(sample)==1 & `z'==1
		estadd scalar be1=r(mean)
		lincom T_2 + T2X`z'
		estadd scalar c2=r(estimate)
		estadd scalar se2=r(se)
		estadd scalar pe2=r(p)
		lincom T_3 + T3X`z'
		estadd scalar c3=r(estimate)
		estadd scalar se3=r(se)
		estadd scalar pe3=r(p)
		lincom T_4 + T4X`z'
		estadd scalar c4=r(estimate)
		estadd scalar se4=r(se)
		estadd scalar pe4=r(p)
		test T_2=T_3
		estadd scalar p23=r(p)
		test T_2=T_4
		estadd scalar p24=r(p)
		test T_3=T_4
		estadd scalar p34=r(p)
		test T_2+T_3=T_4
		estadd scalar p234=r(p)
		test T_2 + T2X`z'=T_3 + T3X`z'
		estadd scalar p23b=r(p)
		test T_2 + T2X`z'=T_4 + T4X`z'
		estadd scalar p24b=r(p)
		test T_3 + T3X`z'=T_4 + T4X`z'
		estadd scalar p34b=r(p)
		test T_2 + T2X`z' + T_3 + T3X`z'=T_4 + T4X`z'
		estadd scalar p234b=r(p)
		su `z' if e(sample)==1 
		estadd scalar be2=r(mean)
		
		parmest, label saving(`z'_`y', replace) stars(0.1 0.05 0.01) idstr("`y'_x_`z'") escal(be be1 c2 se2 pe2 c3 se3 pe3 c4 se4 pe4 p23b p24b p34b p234b N be2) idnum(1) ylab

		}
	
	}	

// append result files
	clear 
	foreach z of global het_new { 
		foreach i in Y1_A_coo Y1_C Y1_B_coo Y1_A_coo_pt Y1_C_pt Y1_B_coo_pt	///
		Y2_A_pro_w Y2_C_w Y2_B_pro_w Y2_A_pro_pt_w Y2_C_pt_w  Y2_B_pro_pt_w  ///
		Y1_A_coo_D Y1_C_D Y1_B_coo_D Y1_A_coo_pt_D Y1_C_pt_D Y1_B_coo_pt_D ///
			{
		append using `z'_`i'
		cap erase `z'_`i'.dta 
		}
}
	
* Keep the output of interest
	keep if parm =="T_2"|parm =="T_3"|parm =="T_4" | parm =="has_older_sis" | parm =="muslim"  ///
		| strpos(parm,"X") 
	rename es_1 mean_control_0
	rename es_2 mean_control_1
	rename es_3 T2_1_b	
	rename es_4 T2_1_se
	rename es_5 T2_1_p
	rename es_6 T3_1_b
	rename es_7 T3_1_se 	
	rename es_8 T3_1_p
	rename es_9 T4_1_b
	rename es_10 T4_1_se	
	rename es_11 T4_1_p
	rename es_12 p23
	rename es_13 p24
	rename es_14 p34
	rename es_15 p234
	rename es_16 obs
	rename es_17 mean_hetvar

save "$constructed/estimates_het_new", replace

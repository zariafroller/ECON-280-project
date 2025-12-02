/*=======================================

* Author: Zaria Roller
* Date: 12/5/25

* Extension of "Childcare, Labor Supply, and 
Business Development: Experimental Evidence 
from Uganda" 

* This file creates the tables for my extension. 

========================================*/

*****************
* CREATE TABLES *
*****************

use "$constructed/estimates_het_new", clear 

replace label = "Childcare $\times$ older sister" if parm=="T2Xhas_older_sis"
replace label = "Cash $\times$ older sister" if parm=="T3Xhas_older_sis"
replace label = "Childcare \& cash $\times$ older sister" if parm=="T4Xhas_older_sis"

replace label = "Childcare $\times$ Muslim" if parm=="T2Xmuslim"
replace label = "Cash $\times$ Muslim" if parm=="T3Xmuslim"
replace label = "Childcare \& cash $\times$ Muslim" if parm=="T4Xmuslim"

* a) define families for mht
* --------------------------

	split idstr, p("_x_")
	drop idstr 
	rename idstr1 idstr

	* main
	gen 	family = 1  if idstr == "Y1_A_coo" | idstr == "Y1_A_coo_D"| idstr =="Y1_C" | idstr == "Y1_C_D"  |idstr == "Y1_B_coo"  | idstr == "Y1_B_coo_D"

	replace family = 2 if idstr == "Y1_A_coo_pt" | idstr == "Y1_A_coo_pt_D"| idstr == "Y1_C_pt" | idstr == "Y1_C_pt_D"  |idstr == "Y1_B_coo_pt"  | idstr == "Y1_B_coo_pt_D"
	
	replace family = 3 if idstr == "Y2_A_pro_w" | idstr == "Y2_B_pro_w" |idstr == "Y2_C_w" 
	
	replace family = 4 if idstr == "Y2_A_pro_pt_w" | idstr == "Y2_B_pro_pt_w" | idstr == "Y2_C_pt_w" 
	
	* only keep current families
	drop if family ==.
	
* b) mht 
* ------

	gen 	fam_x = 1 if parm =="T_2" | strpos(parm,"T2X")
	replace fam_x = 2 if parm =="T_3" | strpos(parm,"T3X")
	replace fam_x = 3 if parm =="T_4" | strpos(parm,"T4X")
	replace fam_x = 4 if strpos(parm,"T2X")
	replace fam_x = 5 if strpos(parm,"T3X")
	replace fam_x = 6 if strpos(parm,"T4X")

	egen familyx = group(family fam_x idstr2) 
	
	* loop over families	
	sum familyx,d 
	forval f=1/`r(max)' {
	preserve
	keep if familyx ==`f'
	gen pval = p
			
	** using the code "fdr_sharpened_qvalues.do" from Michael Anderson http://are.berkeley.edu/~mlanderson/downloads/fdr_sharpened_qvalues.do.zip
			
	* Collect the total number of p-values tested

	quietly sum pval
	local totalpvals = r(N)

	* Sort the p-values in ascending order and generate a variable that codes each p-value's rank

	quietly gen int original_sorting_order = _n
	quietly sort pval
	quietly gen int rank = _n if pval != .

	* Set the initial counter to 1 

	local qval = 1

	* Generate the variable that will contain the BKY (2006) sharpened q-values

	gen bky06_qval = 1 if pval != .

	* Set up a loop that begins by checking which hypotheses are rejected at q = 1.000, then checks which hypotheses are rejected at q = 0.999, then checks which hypotheses are rejected at q = 0.998, etc.  The loop ends by checking which hypotheses are rejected at q = 0.001.


	while `qval' > 0 {
		* First Stage
		* Generate the adjusted first stage q level we are testing: q' = q/1+q
		local qval_adj = `qval' / ( 1 + `qval' )
		* Generate value q'*r/M
		gen fdr_temp1 = `qval_adj' * rank / `totalpvals'
		* Generate binary variable checking condition p(r) <= q'*r/M
		gen reject_temp1 = ( fdr_temp1 >= pval ) if pval != .
		* Generate variable containing p-value ranks for all p-values that meet above condition
		gen reject_rank1 = reject_temp1 * rank
		* Record the rank of the largest p-value that meets above condition
		egen total_rejected1 = max( reject_rank1 )

		* Second Stage
		* Generate the second stage q level that accounts for hypotheses rejected in first stage: q_2st = q'*(M/m0)
		local qval_2st = `qval_adj' * ( `totalpvals' / ( `totalpvals' - total_rejected1[1] ) )
		* Generate value q_2st*r/M
		gen fdr_temp2 = `qval_2st' * rank / `totalpvals'
		* Generate binary variable checking condition p(r) <= q_2st*r/M
		gen reject_temp2 = ( fdr_temp2 >= pval ) if pval != .
		* Generate variable containing p-value ranks for all p-values that meet above condition
		gen reject_rank2 = reject_temp2 * rank
		* Record the rank of the largest p-value that meets above condition
		egen total_rejected2 = max( reject_rank2 )

		* A p-value has been rejected at level q if its rank is less than or equal to the rank of the max p-value that meets the above condition
		replace bky06_qval = `qval' if rank <= total_rejected2 & rank != .
		* Reduce q by 0.001 and repeat loop
		drop fdr_temp* reject_temp* reject_rank* total_rejected*
		local qval = `qval' - .001
	}


	quietly sort original_sorting_order

	display "Code has completed."
	display "Benjamini Krieger Yekutieli (2006) sharpened q-vals are in variable 'bky06_qval'"
	display	"Sorting order is the same as the original vector of p-values"
	keep familyx idstr idstr2 parm bky06_qval
	save "$results/mht06_family_`f'" ,replace
	restore

	** T2, T3 and T4 effect in subgroup (interacted variable=1)
	forval t = 2/4 {
		preserve
		keep if familyx ==`f'
		gen pval = T`t'_1_p
				
		** using the code "fdr_sharpened_qvalues.do" from Michael Anderson http://are.berkeley.edu/~mlanderson/downloads/fdr_sharpened_qvalues.do.zip
				
		* Collect the total number of p-values tested

		quietly sum pval
		local totalpvals = r(N)

		* Sort the p-values in ascending order and generate a variable that codes each p-value's rank

		quietly gen int original_sorting_order = _n
		quietly sort pval
		quietly gen int rank = _n if pval != .

		* Set the initial counter to 1 

		local qval = 1

		* Generate the variable that will contain the BKY (2006) sharpened q-values

		gen bky06_qval = 1 if pval != .

		* Set up a loop that begins by checking which hypotheses are rejected at q = 1.000, then checks which hypotheses are rejected at q = 0.999, then checks which hypotheses are rejected at q = 0.998, etc.  The loop ends by checking which hypotheses are rejected at q = 0.001.


		while `qval' > 0 {
			* First Stage
			* Generate the adjusted first stage q level we are testing: q' = q/1+q
			local qval_adj = `qval' / ( 1 + `qval' )
			* Generate value q'*r/M
			gen fdr_temp1 = `qval_adj' * rank / `totalpvals'
			* Generate binary variable checking condition p(r) <= q'*r/M
			gen reject_temp1 = ( fdr_temp1 >= pval ) if pval != .
			* Generate variable containing p-value ranks for all p-values that meet above condition
			gen reject_rank1 = reject_temp1 * rank
			* Record the rank of the largest p-value that meets above condition
			egen total_rejected1 = max( reject_rank1 )

			* Second Stage
			* Generate the second stage q level that accounts for hypotheses rejected in first stage: q_2st = q'*(M/m0)
			local qval_2st = `qval_adj' * ( `totalpvals' / ( `totalpvals' - total_rejected1[1] ) )
			* Generate value q_2st*r/M
			gen fdr_temp2 = `qval_2st' * rank / `totalpvals'
			* Generate binary variable checking condition p(r) <= q_2st*r/M
			gen reject_temp2 = ( fdr_temp2 >= pval ) if pval != .
			* Generate variable containing p-value ranks for all p-values that meet above condition
			gen reject_rank2 = reject_temp2 * rank
			* Record the rank of the largest p-value that meets above condition
			egen total_rejected2 = max( reject_rank2 )

			* A p-value has been rejected at level q if its rank is less than or equal to the rank of the max p-value that meets the above condition
			replace bky06_qval = `qval' if rank <= total_rejected2 & rank != .
			* Reduce q by 0.001 and repeat loop
			drop fdr_temp* reject_temp* reject_rank* total_rejected*
			local qval = `qval' - .001
		}


		quietly sort original_sorting_order

		display "Code has completed."
		display "Benjamini Krieger Yekutieli (2006) sharpened q-vals are in variable 'bky06_qval'"
		display	"Sorting order is the same as the original vector of p-values"
		keep familyx idstr idstr2 parm bky06_qval
		rename bky06_qval T`t'_1_bky06_qval
		save "$results/mht06_family_`f'_T`t'" ,replace
		restore
		}
	}

	* append the families and their qvalues 
	sum familyx,d 
	gen T2_1_bky06_qval = .
	gen T3_1_bky06_qval = .
	gen T4_1_bky06_qval = .
	forval f=1/`r(max)' {
		merge 1:1 familyx idstr idstr2 parm using  "$results/mht06_family_`f'", update 
		erase "$results/mht06_family_`f'.dta"
		drop _merge
		forval t = 2/4 {
			merge 1:1 familyx idstr idstr2 parm using  "$results/mht06_family_`f'_T`t'", update  keepusing(T`t'_1_bky06_qval)
			erase "$results/mht06_family_`f'_T`t'.dta"
			drop _merge
		}
	}

	* need to fill in missing values because of how this is reshaped when exporting the table
	sort   idstr idstr2 parm	
	forval t = 2/4 {
			replace T`t'_1_bky06_qval= T`t'_1_bky06_qval[_n-1] if  T`t'_1_bky06_qval==.
		}

	* Create stars 
		gen q = "\star"	 if bky06_qval <= 0.1
			replace q = "\star\star" if bky06_qval <= 0.05
			replace q = "\star\star\star" if bky06_qval <= 0.01
		gen temp =  "\ast" if stars=="*"
			replace temp = "\ast\ast" if stars=="**"
			replace temp = "\ast\ast\ast" if stars=="***"
		gen star_q = "$\overset{"+temp+"}{\scriptstyle{"+q+"}}$"
		forval t = 2/4 {
			gen q_T`t' = "\star"	 if T`t'_1_bky06_qval <= 0.1
				replace q_T`t' = "\star\star" if T`t'_1_bky06_qval <= 0.05
				replace q_T`t' = "\star\star\star" if T`t'_1_bky06_qval <= 0.01
			gen temp_T`t' =  "\ast" if T`t'_1_p<= 0.1
				replace temp_T`t' = "\ast\ast" if  T`t'_1_p<= 0.05
				replace temp_T`t' = "\ast\ast\ast" if  T`t'_1_p<= 0.01
			gen star_q_T`t' = "$\overset{"+temp_T`t'+"}{\scriptstyle{"+q_T`t'+"}}$"
		}
		sort familyx idstr parm	

	save "$constructed/estimates_mht_new", replace	
	
	
	use "$constructed/estimates_mht_new", clear
	
	
* c) Make the tables
* ------------------

	gen table = "mother" if idstr == "Y1_A_coo" | idstr == "Y1_A_coo_D" | idstr == "Y1_C" | idstr == "Y1_C_D"  |idstr == "Y1_B_coo"  | idstr == "Y1_B_coo_D" | idstr == "Y2_A_pro_w" | idstr == "Y2_C_w"|  idstr == "Y2_B_pro_w"
	replace table = "father" if idstr == "Y1_A_coo_pt" | idstr == "Y1_A_coo_pt_D" | idstr == "Y1_C_pt" | idstr == "Y1_C_pt_D"  |idstr == "Y1_B_coo_pt"  | idstr == "Y1_B_coo_pt_D"  |idstr == "Y2_A_pro_pt_w"| idstr == "Y2_C_pt_w"| idstr == "Y2_B_pro_pt_w"
	replace table = "consumption" if idstr == "cons_tot_w" |  idstr == "cons_food_w" |  idstr == "cons_nonfood_w"  | idstr == "Y5_A_pro_w"
	
rename stderr se
foreach x of var se estimate {
		replace `x' = round(`x', .01) if table=="mother" |table=="father" 
		format `x' %9.2f
		sdecode `x', gen(`x'_st)	
		}
		replace se_st = "("+se_st+")"

	
// order the rows (that will become table columns)
		
	* mother 
	replace idstr = "A2" if idstr == "Y2_B_pro_w"
	replace idstr = "A4" if idstr == "Y2_C_w"
	replace idstr = "A6" if idstr == "Y2_A_pro_w"
	replace idstr = "A8" if idstr == "Y1_B_coo_D"
	replace idstr = "A10" if idstr == "Y1_B_coo"
	replace idstr = "A12" if idstr == "Y1_C_D"
	replace idstr = "A14" if idstr == "Y1_C"
	replace idstr = "A16" if idstr == "Y1_A_coo_D"
	replace idstr = "A18" if idstr == "Y1_A_coo"

	* father
	replace idstr = "A2" if idstr == "Y2_B_pro_pt_w"
	replace idstr = "A4" if idstr == "Y2_C_pt_w"
	replace idstr = "A6" if idstr == "Y2_A_pro_pt_w"
	replace idstr = "A8" if idstr == "Y1_B_coo_pt_D"
	replace idstr = "A10" if idstr == "Y1_B_coo_pt"
	replace idstr = "A12" if idstr == "Y1_C_pt_D"
	replace idstr = "A14" if idstr == "Y1_C_pt"
	replace idstr = "A16" if idstr == "Y1_A_coo_pt_D"
	replace idstr = "A18" if idstr == "Y1_A_coo_pt" 

	drop mean_control_0

* d) Export the Tables 
* --------------------
	
drop if table ==""
replace table = table+"_"+idstr2 
levelsof table, local(tab)

foreach T in `tab' {
	preserve
		keep if table == "`T'"	
		levelsof idstr, local(deps)
		keep estimate_st se_st obs mean_hetvar mean_control_1 label idstr idstr2 star_q T2_1_b T2_1_se star_q_T2 T3_1_b T3_1_se star_q_T3 T4_1_b T4_1_se star_q_T4 p23 p24 p34 p234  // 
		reshape wide estimate_st se_st T2_1_b T2_1_se star_q_T2 T3_1_b T3_1_se star_q_T3 T4_1_b T4_1_se star_q_T4 p23 p24 p34 p234  obs  mean_control_1 mean_hetvar star_q, i(label) j(idstr) string
				dis in red "`T'"
		expand 2, generate(expanded)
		foreach x in `deps' {
			gen col1_`x' = estimate_st`x' if expanded ==0
				replace col1_`x'  = se_st`x' if expanded ==1
			gen col2_`x' = star_q`x'	if expanded ==0
			}
		* add lines:
		local new = _N + 17
			set obs `new'
		replace label = "Mean het. variable" if _n == _N-16
		replace label = "p-value (equal treatment effects)" if _n ==_N-15
		replace label = "\midrule1" if _n ==_N-14
		replace label = "Obs." if _n ==_N-13
		replace label = "lincom" if _n == _N-12
		replace label = "Mean Control" if _n ==_N-11
		replace label = "\midrule2" if _n ==_N-10
		replace label = "Childcare effect" if _n ==_N-9
		replace label = "se2" if _n ==_N-8
		replace label = "Cash effect" if _n ==_N-7
		replace label = "se3" if _n ==_N-6
		replace label = "Childcare \& cash effect" if _n ==_N-5
		replace label = "se4" if _n ==_N-4
		replace label = "Childcare = cash" if _n ==_N-3
		replace label = "Childcare = childcare \& cash" if _n ==_N-2
		replace label = "Cash = childcare \& cash" if _n ==_N-1
		replace label = "Childcare \& cash = childcare + cash" if _n ==_N

		carryforward obs* mean* T2_* T3_* T4_* p* idstr2 star_q_T4* star_q_T3* star_q_T2*, replace
			dis `deps'
		foreach x in `deps'  { // A10 A12 A14 A16 A18 
			replace mean_control_1`x' = round(mean_control_1`x') if mean_control_1`x'>=1 | mean_control_1`x'<=-1  //, .001
			replace mean_control_1`x' = round(mean_control_1`x', .01) if mean_control_1`x'<1 & mean_control_1`x'>-1  //
			replace mean_hetvar`x' = round(mean_hetvar`x') if mean_hetvar`x'>=1 | mean_hetvar`x'<=-1  //, .001
			replace mean_hetvar`x' = round(mean_hetvar`x', .01) if mean_hetvar`x'<1 & mean_hetvar`x'>-1  //
			replace p23`x' = round(p23`x', .001)
			replace p24`x' = round(p24`x', .001)
			replace p34`x' = round(p34`x', .001)
			replace p234`x' = round(p234`x', .001)
			
			replace T2_1_b`x' = round(T2_1_b`x', .01)
			replace T2_1_se`x' = round(T2_1_se`x', .01)
			replace T3_1_b`x' = round(T3_1_b`x', .01)
			replace T3_1_se`x' = round(T3_1_se`x', .01)
			replace T4_1_b`x' = round(T4_1_b`x', .01)
			replace T4_1_se`x' = round(T4_1_se`x', .01)
			
			tostring obs`x', gen(N`x'_st)	force	
			tostring mean_control_1`x', gen(M1`x'_st)	force	
			tostring mean_hetvar`x', gen(M2`x'_st)	force	

			tostring p23`x', gen(p23`x'_st)	force	
			tostring p24`x', gen(p24`x'_st)	force	
			tostring p34`x', gen(p34`x'_st)	force	
			tostring p234`x', gen(p234`x'_st)	force
			
			tostring T2_1_b`x', gen(T2_1_b`x'_st)	force
			tostring T2_1_se`x', gen(T2_1_se`x'_st)	force
				replace T2_1_se`x'_st = "("+T2_1_se`x'_st+")"
			tostring T3_1_b`x', gen(T3_1_b`x'_st)	force
			tostring T3_1_se`x', gen(T3_1_se`x'_st)	force
				replace T3_1_se`x'_st = "("+T3_1_se`x'_st+")"
			tostring T4_1_b`x', gen(T4_1_b`x'_st)	force
			tostring T4_1_se`x', gen(T4_1_se`x'_st)	force
				replace T4_1_se`x'_st = "("+T4_1_se`x'_st+")"

			replace col1_`x' = N`x'_st if label ==  "Obs." 
			replace col1_`x' = M1`x'_st if label ==  "Mean Control" 
			replace col1_`x' = M2`x'_st if label ==  "Mean het. variable" 

			replace col1_`x' = T2_1_b`x'_st if label ==  "Childcare effect" 
			replace col2_`x' = star_q_T2`x' if label == "Childcare effect"
			
			replace col1_`x' = T2_1_se`x'_st if label ==  "se2" 
			replace col1_`x' = T3_1_b`x'_st if label ==  "Cash effect" 
				replace col2_`x' = star_q_T3`x' if label ==  "Cash effect" 
				
			replace col1_`x' = T3_1_se`x'_st if label ==  "se3" 
			replace col1_`x' = T4_1_b`x'_st if label ==  "Childcare \& cash effect" 
				replace col2_`x' = star_q_T4`x' if label == "Childcare \& cash effect"
				
			replace col1_`x' = T4_1_se`x'_st if label ==  "se4" 
			replace col1_`x' = p23`x'_st if label ==  "Childcare = cash" 
			replace col1_`x' = p24`x'_st if label ==  "Childcare = childcare \& cash" 
			replace col1_`x' = p34`x'_st if label ==  "Cash = childcare \& cash" 
			replace col1_`x' = p234`x'_st if label ==  "Childcare \& cash = childcare + cash" 

			 gen `x'_col12 = col1_`x' + col2_`x'
			}
			
		* Order the lines	
		gen order = 1 if label == "Childcare"	
			replace order = 2 if label == "Cash"
			replace order = 3 if label == "Childcare \& cash"
			replace order = 4 if label == "Older sister" |  label == "Mother's religion is Islam" 
			replace order = 5 if strpos(label,"Childcare $\times$")
			replace order = 6 if strpos(label,"Cash $\times$")
			replace order = 7 if strpos(label,"Childcare \& cash $\times$")
			replace order = 8 if label == "\midrule1"	
			replace order = 9 if strpos(label,"Impact for households")
			replace order = 10 if label == "lincom"
			replace order = 11 if label == "Childcare effect"	
			replace order = 12 if label == "se2"	
			replace order = 13 if label == "Cash effect"	
			replace order = 14 if label == "se3"	
			replace order = 15 if label == "Childcare \& cash effect"	
			replace order = 16 if label == "se4"	
			replace order = 17 if label == "p-value (equal treatment effects)"
			replace order = 18 if label == "Childcare = cash"	
			replace order = 19 if label == "Childcare = childcare \& cash"	
			replace order = 20 if label == "Cash = childcare \& cash"	
			replace order = 21 if label == "Childcare \& cash = childcare + cash"	
			replace order = 22 if label == "\midrule2"	
			replace order = 23 if label == "Mean Control"
			replace order = 24 if label == "Mean het. variable"
			replace order = 25 if label == "Obs." 
			replace order = 26 if label == "\bottomrule"
			
			replace label = "Impact with older sister" if label == "lincom" & idstr2 == "has_older_sis"

			replace label = "Impact for Muslim mothers" if label == "lincom" & idstr2 == "muslim"

			
			replace label = "Mean Control (with older sister)" if label == "Mean Control" & idstr2 == "has_older_sis"

			replace label = "Mean Control (Muslim mothers)" if label == "Mean Control" & idstr2 == "muslim"
			
			sort order expanded	
		keep label *col12 //col1* col2*	
		replace label = "\midrule" if strpos(label, "\midrule")
		replace label = "" if label[_n] == label[_n-1]
		replace label = "" if inlist(label,"se2","se3","se4") // ADDED
		replace label = "Childcare" if label == "Childcare effect"
		replace label = "Cash" if label == "Cash effect"
		replace label = "Childcare \& cash" if label == "Childcare \& cash effect"
		gen end = "\\" if label!= "\midrule" & label!= "\bottomrule"
		qui desc // r(k) then contains the number of variables
		local n = `r(k)'-2 // need (number of vars - 2) columns of "&" for latex
		forvalues i = 1/`n'{
			local j=2*`i'-1
			gen A`j' = "&" if label!= "\midrule" & label!= "\bottomrule"
			}
		order _all, seq
		order label, first
		order end, last 
		dis in red "`T'"
		export delimited using "$results/het_new/tab_`T'.tex", delimiter(tab) novarnames nolabel replace
	 restore
	}	
	
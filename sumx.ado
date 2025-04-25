///=============================================================================
/*
 * Project: Programs for corporate finance empirical studies
 * Author: Hao Zhao
 * Version
	- sumx: 1.3.4 (19apr2025)
 */
///=============================================================================
/* sumx -> Generate summary statistics tables */
///=============================================================================

capture program drop sumx
program sumx, sortpreserve

	syntax anything [if] [in], deci(numlist) edir(string) [category(string) tgroup(string) tvar(string) ///
	ttitle(string) addn(string) RLABEL TTEST ONLYT PRESENT]
	
	marksample touse

	/* Table title: "`table_title'" */
	if ("`ttitle'" == "") {
		local table_title = "Summary statistics `addn'"
	}
	else {
		local table_title = "`ttitle'"
	}
	
	local var_len : word count `anything'
	
	/* tgroup */
	if ("`tgroup'" != "") {
		if ("`category'" != "") {
			display in red "[Error] category cannot be used with tgroup"
			exit
		}
		if ("`tvar'" != "") {
			display in red "[Error] tvar cannot be used with tgroup"
			exit
		}
		eststo clear
		ereturn clear
		estimates clear
		
		levelsof `tgroup', local(tdummies)
		local n_dummy : word count `tdummies'
		if (`n_dummy' != 2) {
			display in red "[Error] input of tgroup must be a binary variable"
		}
		else {
			local grpdumidx = 1
			foreach group in `tdummies' {
				matrix sum_stat = J(`var_len', 8, 0)
				matrix colnames sum_stat = "Obs" "Mean" "Std.Dev." "Min" "p25" "Median" "p75" "Max"
				local rowname_li = `""'
				foreach sum_v in `anything' {
					if ("`rlabel'" != "") {
						local label_v : var label `sum_v'
						if ("`label_v'" == "") {
							local label_v = "`sum_v'"
						}
						if (strlen("`label_v'") > 32) {
							local label_v = substr("`label_v'", 1, 32)
						}
						local rowname_li = `"`rowname_li' "`label_v'" "'
					}
					else {
						local rowname_li = `"`rowname_li' "`sum_v'" "'
					}
				}
				matrix rownames sum_stat = `rowname_li'
				local i = 1
				foreach var in `anything' {
					
					local numeric_value = real("`group'")
					if !missing(`numeric_value') {
						quietly: summ `var' if `touse' & `tgroup'==`group', detail
					}
					else {
						quietly: summ `var' if `touse' & `tgroup'=="`group'", detail
					}
					matrix sum_stat[`i',1] = r(N)
					matrix sum_stat[`i',2] = round(r(mean), `deci')
					matrix sum_stat[`i',3] = round(r(sd), `deci')
					matrix sum_stat[`i',4] = round(r(min), `deci')
					matrix sum_stat[`i',5] = round(r(p25), `deci')
					matrix sum_stat[`i',6] = round(r(p50), `deci')
					matrix sum_stat[`i',7] = round(r(p75), `deci')
					matrix sum_stat[`i',8] = round(r(max), `deci')
					local i = `i'+ 1
				}
				matrix summatrix`grpdumidx' = sum_stat
				if "`onlyt'" == "" {
					mat2txt, matrix(sum_stat) saving("`edir'") title("`table_title' (`tgroup'=`group')") append
				}
				local grpdumidx = `grpdumidx' + 1
			}
			estimates clear
			ereturn clear
			
			estpost ttest `anything' if `touse', by(`tgroup')

			if "`present'" == "" {
				matrix ttest_stat = J(`var_len', 4, 0)
				matrix colnames ttest_stat = "Diff" "Std.Err." "P-value"  "N"
				matrix rownames ttest_stat = `rowname_li'
				local i = 1
				foreach var in `anything' {
					local e_b = e(b)[1, "`var'"]
					matrix ttest_stat[`i',1] = round(`e_b', `deci')
					local e_se = e(se)[1, "`var'"]
					matrix ttest_stat[`i',2] = round(`e_se', `deci')
					local e_p = e(p)[1, "`var'"]
					matrix ttest_stat[`i',3] = round(`e_p', `deci')
					local e_count = e(count)[1, "`var'"]
					matrix ttest_stat[`i',4] = round(`e_count', `deci')
					local i = `i'+ 1
				}
				mat2txt, matrix(ttest_stat) saving("`edir'") title("T-test for (`tgroup') `addn'") append
			}
			else {
				matrix ttest_stat = J(`var_len', 10, 0)
				matrix colnames ttest_stat = "Mean1" "Median1" "N1" "Mean2" "Median2" "N2" "Diff-mean" "P-value" "Diff-median" "N"
				matrix rownames ttest_stat = `rowname_li'
				local i = 1
				foreach var in `anything' {
					matrix ttest_stat[`i',1] = summatrix1[`i',2]
					matrix ttest_stat[`i',2] = summatrix1[`i',6]
					matrix ttest_stat[`i',3] = summatrix1[`i',1]
					matrix ttest_stat[`i',4] = summatrix2[`i',2]
					matrix ttest_stat[`i',5] = summatrix2[`i',6]
					matrix ttest_stat[`i',6] = summatrix2[`i',1]
					local e_b = e(b)[1, "`var'"]
					matrix ttest_stat[`i',7] = round(`e_b', `deci')
					local e_p = e(p)[1, "`var'"]
					matrix ttest_stat[`i',8] = round(`e_p', `deci')
					local mediandiff = summatrix1[`i',6] - summatrix2[`i',6]
					matrix ttest_stat[`i',9] = round(`mediandiff', `deci')	
					local e_count = e(count)[1, "`var'"]
					matrix ttest_stat[`i',10] = round(`e_count', `deci')					
					local i = `i'+ 1
				}
				mat2txt, matrix(ttest_stat) saving("`edir'") title("T-test for (`tgroup') `addn'") append
			}

			estimates clear
			ereturn clear
		}
	} 
	
	/* tvar */
	if ("`tvar'" != "") {
		if ("`category'" != "") {
			display in red "[Error] category cannot be used with tvar"
			exit
		}
		local n_tvar : word count `tvar'
		if (`n_tvar' != `var_len') {
			display in red "[Error] The number of paired variables should equal to the number of variables"
			exit
		}
		else {
			local vpair_idx = 1
			foreach var_pair in anything tvar {
				matrix sum_stat = J(`var_len', 8, 0)
				matrix colnames sum_stat = "Obs" "Mean" "Std.Dev." "Min" "p25" "Median" "p75" "Max"
				local rowname_li = `""'
				foreach sum_v in ``var_pair'' {
					if ("`rlabel'" != "") {
						local label_v : var label `sum_v'
						if ("`label_v'" == "") {
							local label_v = "`sum_v'"
						}
						if (strlen("`label_v'") > 32) {
							local label_v = substr("`label_v'", 1, 32)
						}
						local rowname_li = `"`rowname_li' "`label_v'" "'
					}
					else {
						local rowname_li = `"`rowname_li' "`sum_v'" "'
					}
				}
				matrix rownames sum_stat = `rowname_li'
				local i = 1
				foreach x in ``var_pair'' {
					quiet: summ `x' if `touse', detail
					matrix sum_stat[`i',1] = r(N)
					matrix sum_stat[`i',2] = round(r(mean), `deci')
					matrix sum_stat[`i',3] = round(r(sd), `deci')
					matrix sum_stat[`i',4] = round(r(min), `deci')
					matrix sum_stat[`i',5] = round(r(p25), `deci')
					matrix sum_stat[`i',6] = round(r(p50), `deci')
					matrix sum_stat[`i',7] = round(r(p75), `deci')
					matrix sum_stat[`i',8] = round(r(max), `deci')
					local i = `i'+ 1
				}
				mat2txt, matrix(sum_stat) saving("`edir'") title("`table_title' (variable list `vpair_idx')") append
				local vpair_idx = `vpair_idx' + 1
			}

			matrix ttest_stat = J(`var_len', 4, 0)
			matrix colnames ttest_stat = "Diff" "Std.Err." "P-value"  "N"
			matrix rownames ttest_stat = `rowname_li'
			forval var_idx = 1/`var_len' {	
				local pair_v1 : word `var_idx' of `anything'
				local pair_v2 : word `var_idx' of `tvar'
				capture {
					ttest `pair_v1' == `pair_v2' if `touse'
					local e_b = r(mu_1) - r(mu_2)
					matrix ttest_stat[`var_idx', 1] = round(`e_b', `deci')
					local e_se = r(se)
					matrix ttest_stat[`var_idx', 2] = round(`e_se', `deci')
					local e_p = r(p)
					matrix ttest_stat[`var_idx', 3] = round(`e_p', `deci')
					local e_count = r(N_1)
					matrix ttest_stat[`var_idx', 4] = round(`e_count', `deci')
				}
			}
			mat2txt, matrix(ttest_stat) saving("`edir'") title("T-test for (`: word 1 of `anything''==`: word 1 of `tvar'')") append
		}
	}
	
	if ("`category'" == "" & "`tgroup'" == "" & "`tvar'" == "") {
		matrix sum_stat = J(`var_len', 8, 0)
		matrix colnames sum_stat = "Obs" "Mean" "Std.Dev." "Min" "p25" "Median" "p75" "Max"
		local rowname_li = `""'
		foreach sum_v in `anything' {
			if ("`rlabel'" != "") {
				local label_v : var label `sum_v'
				if ("`label_v'" == "") {
					local label_v = "`sum_v'"
				}
				if (strlen("`label_v'") > 32) {
					local label_v = substr("`label_v'", 1, 32)
				}
				local rowname_li = `"`rowname_li' "`label_v'" "'
			}
			else {
				local rowname_li = `"`rowname_li' "`sum_v'" "'
			}
		}
		matrix rownames sum_stat = `rowname_li'
		local i = 1
		foreach x in `anything' {
			quiet: summ `x' if `touse', detail
			matrix sum_stat[`i',1] = r(N)
			matrix sum_stat[`i',2] = round(r(mean), `deci')
			matrix sum_stat[`i',3] = round(r(sd), `deci')
			matrix sum_stat[`i',4] = round(r(min), `deci')
			matrix sum_stat[`i',5] = round(r(p25), `deci')
			matrix sum_stat[`i',6] = round(r(p50), `deci')
			matrix sum_stat[`i',7] = round(r(p75), `deci')
			matrix sum_stat[`i',8] = round(r(max), `deci')
			local i = `i'+ 1
		}
		mat2txt, matrix(sum_stat) saving("`edir'") title("`table_title'") append
		
		if ("`ttest'" != "") {
			matrix ttest_stat = J(`var_len', 4, 0)
			matrix colnames ttest_stat = "Diff" "Std.Err." "P-value"  "N"
			matrix rownames ttest_stat = `rowname_li'
			forval var_idx = 1/`var_len' {
				local pair_v1 : word `var_idx' of `anything'
				capture {
					ttest `pair_v1' == 0 if `touse'
					local e_b = r(mu_1) - 0
					matrix ttest_stat[`var_idx', 1] = round(`e_b', `deci')
					local e_se = r(se)
					matrix ttest_stat[`var_idx', 2] = round(`e_se', `deci')
					local e_p = r(p)
					matrix ttest_stat[`var_idx', 3] = round(`e_p', `deci')
					local e_count = r(N_1)
					matrix ttest_stat[`var_idx', 4] = round(`e_count', `deci')
				}
			}
			mat2txt, matrix(ttest_stat) saving("`edir'") title("T-test for (`: word 1 of `anything''== 0)") append
		}
	}
	else if ("`category'" != "") {
		if `var_len' != 1 {
			display in red "[ERROR] Only one input variable allowed"
			exit
		}
		else {
			levelsof `category' if `touse', local(groups)
			local n_group : word count `groups'
			matrix sum_stat = J(`n_group', 8, 0)
			matrix colnames sum_stat = "Obs" "Mean" "Std.Dev." "Min" "p25" "Median" "p75" "Max"
			local rowname_li = `""'
			foreach group in `groups' {
				if (strlen("`group'") > 32) {
					local r_grp = substr("`group'", 1, 32)
				}
				else {
					local r_grp = "`group'"
				}
				local rowname_li = `"`rowname_li' "`r_grp'" "'
			}
			matrix rownames sum_stat = `rowname_li'
			local i = 1
			foreach group in `groups' {
				local numeric_value = real("`group'")
				if !missing(`numeric_value') {
					quietly: summ `anything' if `touse' & `category'==`group', detail
				}
				else {
					quietly: summ `anything' if `touse' & `category'=="`group'", detail
				}
				matrix sum_stat[`i',1] = r(N)
				matrix sum_stat[`i',2] = round(r(mean), `deci')
				matrix sum_stat[`i',3] = round(r(sd), `deci')
				matrix sum_stat[`i',4] = round(r(min), `deci')
				matrix sum_stat[`i',5] = round(r(p25), `deci')
				matrix sum_stat[`i',6] = round(r(p50), `deci')
				matrix sum_stat[`i',7] = round(r(p75), `deci')
				matrix sum_stat[`i',8] = round(r(max), `deci')
				local i = `i'+ 1
			}
			mat2txt, matrix(sum_stat) saving("`edir'") title("`table_title'") append
		}		
	}
	
end


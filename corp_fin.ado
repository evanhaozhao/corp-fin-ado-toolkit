///=============================================================================
/*
 * Project: Programs for corporate finance empirical studies
 * Author: Hao Zhao
 * Created: August 19, 2023
 * Modified: January 8, 2024
 */
///=============================================================================
/* regx -> regression to output tables 

 * Identifying controls and FEs, supported by `reghdfe`, `reg`, and `xtreg`
 * Sending variables to models:
	- Baseline: y ~ x
	- Interaction: y ~ x1 * x2
		- x1#x2
	- Subsample: y ~ x `if'
 * Exporting tables, supported by `esttab`:
	- Ordering variables
	- Dropping variables
	- Assigning columns
	- Adding notes
 * In order to make the code in minimal style, unless specified, the default setting is:
	- Export file directory: ${dir_table_flow}
 * Dynamic regression: with time indicator [dyn(x, y=, b=, p=, c=, s=)]
 */
///=============================================================================
capture program drop regx
program regx

	syntax anything [if] [in] , indep(namelist) [ctrl(string) absr(string) xtp(string) ///
	inte(string) clust(namelist) dyn(string) tobit(string) ///
	tnote(string) ttitle(string) addn(string) edir(string) keepvar(string) stosuf(string) REPORT DISPLAY CHISTORE]
	
	marksample touse
	
	/* Default export file: ${dir_table_flow} */
	/* Export file priority: `edir' > ${dir_table_flow} */
	if (`"`edir'"'!="") {
		local exportfile = "`edir'"
	}
	else {
		if ("${dir_table_flow}"=="") {
			di "[ERROR] Need to set a directory in [edir] or global ${dir_table_flow}"
			exit
		}
		else {
			local exportfile = "${dir_table_flow}"
		}
	}
	
	/* `suest` does not support xtreg random effect; the estimation of xtreg fe is done by reg..absorb */
	if ("`chistore'"!="") {
		if ("`xtp'"!="") {
			di "[WARNING] `suest` does not support xtreg"
		}
		if ("`display'"=="") {
			di "[ERROR] [chistore] must be used with [display]"
			exit
		}
		if ("`absr'"!="") {
			local cate_absr = ""
			foreach absrvar in `absr' {
				local cate_absr = "`cate_absr' i.`absrvar'"
			}
		}
	}
	
	/* Tobit model: cannot use [absr] */
	if ("`tobit'"!="") {
		if ("`absr'"!="" | "`xtp'"!="") {
			di "[ERROR] [tobit] cannot be used with [xtp] or [absr]"
			exit
		}
		else {
			local tobit_list : subinstr local tobit " " "", all
			local tobit_list : subinstr local tobit_list "," " ", all
			
			local tobit_bound = ""
			foreach tobit_i in `tobit_list' {
				if (strpos("`tobit_i'", "ll=")>0) {
					local tobit_add : subinstr local tobit_i "ll=" "", all
					local tobit_bound = "`tobit_bound' ll(`tobit_add')"
				}
				if (strpos("`tobit_i'", "ul=")>0) {
					local tobit_add : subinstr local tobit_i "ul=" "", all
					local tobit_bound = "`tobit_bound' ul(`tobit_add')"
				}
			}
		}
	}
	
	/* Standard error: vce(`cse') */
	/* Priority: `clust' > `xtp' > $clustervar > robust */
	/* (1) If cluster is specified */
	if ("`clust'"!="") {
		local cse = "cl `clust'"
	}
	else {
		/* (2) If panel is claimed & no cluster specified, cluster SE at id level */
		if ("`xtp'"!="") {
			local cse = "cl `: word 1 of `xtp''"
		}
		else {
			/* (3) If both panel and cluster are not specified */
			if ("${clustervar}"!="") {
				local cse = "cl ${clustervar}"
			}
			else {
				local cse = "robust"
				local extra_addn = "[No cluster SE specified]"
			}
		}
	}
	
	/* Control variables */
	if ("`ctrl'"!="") {
		/* drop if no "REPORT" specified */
		local drop_ctrl: subinstr local ctrl "i." "*", all
	}
	
	/* default ordering in tables */
	local orderinte = ""
		
	/* Interaction term: `inte_reglist`indepidx'' */
	if ("`inte'"!="") {
		local orderinte = ""
		local indepidx = 1
		foreach indepvar in `indep' {
			local expandvars = ""		
			foreach term in `inte' {
				local new_ints : subinstr local term "#" " ", all
				local new_sets = "`indepvar' `new_ints'"
				
				/* check for duplicates: for table row ordering */
				foreach nint in `new_ints' {
					if (strpos("`orderinte'", "`nint'")==0) {
						local orderinte = "`orderinte' `nint'"
					}
				}
				
				/* check for duplicates */
				foreach nint in `new_sets' {
					if (strpos("`expandvars'", "`nint'")==0) {
						local expandvars = "`expandvars' `nint'"
					}
				}
			}
			/* split: for ordering */
			local inteidx = 1
			foreach term in `inte' {
				local new_ints : subinstr local term "#" " ", all
				local new_sets = "`indepvar' `new_ints'"
				
				/* if interaction is below triple, use self-define combinations */
				if (`: word count `new_sets''==2) {
					local ns_idx = 1
					foreach ns_item in `new_sets' {
						local tuple`ns_idx' = "`ns_item'"
						local ns_idx = `ns_idx' + 1
					}
					local tuple`ns_idx' = "`new_sets'"
				}
				else if (`: word count `new_sets''==3) {
					local ns_idx = 1
					foreach ns_item in `new_sets' {
						local tuple`ns_idx' = "`ns_item'"
						local ns_idx = `ns_idx' + 1
					}
					/* double interaction */
					forval ns_i = 1/2 {
						forval ns_j = 2/3 {
							if `ns_i'!=`ns_j' {
								local tuple`ns_idx' = "`: word `ns_i' of `new_sets'' `: word `ns_j' of `new_sets''" 
								local ns_idx = `ns_idx' + 1
							}
						}
					}
					/* triple interaction */
					local tuple`ns_idx' = "`new_sets'"
				}
				else {
					cap which tuples
					if _rc==0 {
						tuples `new_sets'
					}
					else {
						di "[ERROR] Need to ssc install `tuples`"
						exit
					}
				}
				
				local matinit = `: word count `new_sets''+1
				local matsize = 2 ^ `: word count `new_sets'' - 1
				
				local intlist = ""
				forval tu_i = `matinit'/`matsize' {
					local new_comb : subinstr local tuple`tu_i' " " "#c.", all
					local intlist = "`intlist' c.`new_comb'"
				}
				
				/* check for duplicates */
				local orderinte`inteidx' = "`intlist'"
				
				foreach nint in `intlist' {
					if (strpos("`expandvars'", "`nint'")==0) {
						local expandvars = "`expandvars' `nint'"
					}
				}
				local inteidx = `inteidx' + 1
			}
			
			/* table variable ordering */
			local len_inte = ""
			local inteidx = `inteidx' - 1
			forval i_inte = 1/`inteidx' {
				if (`i_inte'==1) {
					local len_inte = "`: word count `orderinte`i_inte'''"
				}
				else {
					local len_inte = "`len_inte', `: word count `orderinte`i_inte'''"
				}
			}
			if (`inteidx'<2) {
				local max_len_inte = 1
			}
			else {
				local max_len_inte = max(`len_inte')
			}
			forval idx = 1/`max_len_inte' {
				forval i_inte = 1/`inteidx' {
					if (strpos("`orderinte'", "`: word `idx' of `orderinte`i_inte'''")==0) {
						local orderinte = "`orderinte' `: word `idx' of `orderinte`i_inte'''"
					}
				}
			}

			/* store interaction list for regression */
			local inte_reglist`indepidx' = "`expandvars'"
			local indepidx = `indepidx' + 1
		}
	}
	
	/* dynamic regression */
	if ("`dyn'"!="") {
		if ("`inte'"!="") {
			di "[ERROR] Cannot use option [dyn] and [inte] at the same time"
			exit
		}
		else {
			local orderinte = ""
			local dynoption : subinstr local dyn " " "", all
			local dynoption : subinstr local dynoption "," " ", all
			
			/* default setting: year: [-2, 2]; benchmark: [-2]; plot: no; CI: 95% */
			local dyn_year = "-2/2"
			local dyn_bench = "0"
			local dyn_benchname = "t0"
			local dyn_plot = "Time period relative to the event"
			local dyn_ci = "95"
			local if_dyn_plot = 0
			local dyn_pdir = regexr("`exportfile'", "\/([^\/]+)$", "")
			
			local dyn_kw = "dyn"
			local dyn_idx = 1
			foreach dynopt in `dynoption' {
				if (`dyn_idx'==1) {
					local dyn_prefix = "`dynopt'"
				}
				if (strpos("`dynopt'", "y=")) {
					local dyn_year : subinstr local dynopt "y=" "", all
				}
				if (strpos("`dynopt'", "b=")) {
					local dyn_bench : subinstr local dynopt "b=" "", all
					local dyn_benchname : subinstr local dyn_bench "-" "m", all
					local dyn_benchname = "t`dyn_benchname'"
				}
				if (strpos("`dynopt'", "p=")) {
					local dyn_plot : subinstr local dynopt "p=" "", all
					local dyn_plot : subinstr local dyn_plot "_" " ", all
					local if_dyn_plot = 1
				}
				if (strpos("`dynopt'", "c=")) {
					local dyn_ci : subinstr local dynopt "c=" "", all
				}
				if (strpos("`dynopt'", "s=")) {
					local dyn_folder : subinstr local dynopt "s=" "", all
					if (strpos("`dyn_folder'", "~")) {
						local dyn_kw : subinstr local dyn_folder "~" " ", all
						local dyn_pdir = "`dyn_pdir'/`: word 1 of `dyn_kw''"
						local dyn_kw : word 2 of `dyn_kw'
					}
					else {
						local dyn_pdir = "`dyn_pdir'/`dyn_folder'"
					}
				}
				local dyn_idx = `dyn_idx' + 1
			}
			local dyn_start = min(`: subinstr local dyn_year "/" ", ", all')
			local dyn_end = max(`: subinstr local dyn_year "/" ", ", all')
			if (`dyn_bench'<`dyn_start' | `dyn_bench'>`dyn_end') {
				di "[ERROR] Benchmark (`dyn_bench') is out of range [`dyn_year']"
				exit
			}
			/* set interaction list for regression: `dyn_vars' (full) `dyn_regvars' (drop bench) */
			local dyn_vars = ""
			local dyn_bvars = ""
			local dyn_regvars = ""
			local dyn_rename = `""'
			local dyn_reorder = `""'
			forval i_dyn = `dyn_year' {
				local dyn_var = "`dyn_prefix'_t`i_dyn'"
				local dyn_var : subinstr local dyn_var "-" "m", all
				local dyn_vars = "`dyn_vars' `dyn_var'"
				if (strpos("`dyn_var'", "`dyn_benchname'")==0) {
					local dyn_regvars = "`dyn_regvars' `dyn_var'"
				}
				else {
					local dyn_bvars = "`dyn_bvars' `dyn_var'"
				}
				local dyn_suf : subinstr local i_dyn "-" "m", all
				local dyn_suf = "t`dyn_suf'"
				local dyn_rename = `"`dyn_rename' .*`dyn_suf'="`i_dyn'""'
				local dyn_reorder = `"`dyn_reorder' "`i_dyn'""'
			}
			/* check the existence of dynamic indicators */
			foreach dyn_v in `dyn_regvars' {
				capture confirm variable `dyn_v', exact
				if (_rc!=0) {
					di "[ERROR] Need to generate dynamic indicator `dyn_v'"
					exit
				}
			}
			/* construct interaction list: dyn_plotlist`indepidx' and dyn_reglist`indepidx' */
			local indepidx = 1
			foreach indepvar in `indep' {
				local dynint_bvars = ""
				local dynint_regvars = ""
				foreach dyn_name in `dyn_vars' {
					if (strpos("`dyn_name'", "`dyn_benchname'")==0) {
						local dynint_regvars = "`dynint_regvars' c.`indepvar'#c.`dyn_name'"
					}
					else {
						local dynint_bvars = "`dynint_bvars' c.`indepvar'#c.`dyn_name'"
					}
				}
				local dyn_plotlist`indepidx' = "`indepvar' `dyn_regvars' `dynint_regvars' `dyn_bvars' `dynint_bvars'"
				local dyn_reglist`indepidx' = "`indepvar' `dyn_regvars' `dynint_regvars'"
				local orderinte = "`dyn_regvars' `dynint_regvars'"
				local indepidx = `indepidx' + 1
			}
		}
	}
	
	eststo clear
	
	/* Model type: y ~ x/interactions + [controls] + [FE], [id-FE] cluster SE */
	/* Priority: `xtp' (fe>re) > `absr' > reg */
	/* (1) xtreg: xtset `xtp' */
	if ("`xtp'"!="") {
		if ("`absr'"!="") {
			/* no absorb allowed in xtreg setting */
			di "[ERROR] No [absorb] option allowed for xtreg"
			exit
		}
		else {
			/* _r_effect: random effect */
			if (strpos("`xtp'", " _r_effect")>0) {
				if ("`chistore'"!="") {
					di "[ERROR] `suest` does not support xtreg random effect"
				}
				local xteffect = "re"
				local xtpanel: subinstr local xtp " _r_effect" "", all
			}
			else {
				/* default: fixed effect */
				local xteffect = "fe"
				local xtpanel = "`xtp'"
			}
			
			/* xtreg */
			xtset `xtpanel'
			local depidx = 1
			foreach depvar in `anything' {
				local indepidx = 1
				foreach indepvar in `indep' {
					if ("`inte'"!="") {
						if ("`chistore'"!="") {
							reg `depvar' `inte_reglist`indepidx'' `ctrl' if `touse', absorb(`: word 1 of `xtp'')
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							eststo: xtreg `depvar' `inte_reglist`indepidx'' `ctrl' if `touse', `xteffect' vce(`cse')
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' if `touse', absorb(`: word 1 of `xtp'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								eststo: xtreg `depvar' `dyn_reglist`indepidx'' `ctrl' if `touse', `xteffect' vce(`cse')
								if (`if_dyn_plot'==1) {
									quietly {
										xtreg `depvar' `dyn_plotlist`indepidx'' `ctrl' if `touse', `xteffect' vce(`cse')
										coefplot, keep(c.*#c.*) vertical omitted base levels(`dyn_ci') ///
										rename(`"`dyn_rename'"', regex) order(`"`dyn_reorder'"') ///
										mcolor("0 191 255") ciopts(lcolor("0 97 154") recast(rcap)) ///
										yline(0, lcolor("179 134 0") lpattern(dash)) ///
										graphregion(fcolor(white) ifcolor(white) ilcolor(white)) ///
										xtitle(`dyn_plot') ytitle("Coefficient estimates") ///
										ysize(5) xsize(8)
										graph export "`dyn_pdir'/`indepidx'_`depvar'_`dyn_kw'.pdf", replace
										graph close
									}
								}
							}
						}
						else {
							if ("`chistore'"!="") {
								reg `depvar' `indepvar' `ctrl' if `touse', absorb(`: word 1 of `xtp'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								eststo: xtreg `depvar' `indepvar' `ctrl' if `touse', `xteffect' vce(`cse')
							}
						}
					}
					local indepidx = `indepidx' + 1
				}
				local depidx = `depidx' + 1
			}
		}
	}
	else {
		if ("`absr'"!="") {
			local depidx = 1
			foreach depvar in `anything' {
				local indepidx = 1
				foreach indepvar in `indep' {
					if ("`inte'"!="") {
						if ("`chistore'"!="") {
							reg `depvar' `inte_reglist`indepidx'' `ctrl' `cate_absr' if `touse'
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							/* Interaction model */
							eststo: reghdfe `depvar' `inte_reglist`indepidx'' `ctrl' if `touse', absorb(`absr') vce(`cse')
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' `cate_absr' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								/* dynamic effect model & plot */
								eststo: reghdfe `depvar' `dyn_reglist`indepidx'' `ctrl' if `touse', absorb(`absr') vce(`cse')
								if (`if_dyn_plot'==1) {
									quietly {
										reghdfe `depvar' `dyn_plotlist`indepidx'' `ctrl' if `touse', absorb(`absr') vce(`cse')
										coefplot, keep(c.*#c.*) vertical omitted base levels(`dyn_ci') ///
										rename(`dyn_rename', regex) order(`dyn_reorder') ///
										mcolor("0 191 255") ciopts(lcolor("0 97 154") recast(rcap)) ///
										yline(0, lcolor("179 134 0") lpattern(dash)) ///
										graphregion(fcolor(white) ifcolor(white) ilcolor(white)) ///
										xtitle(`dyn_plot') ytitle("Coefficient estimates") ///
										ysize(5) xsize(8)
										graph export "`dyn_pdir'/`indepidx'_`depvar'_`dyn_kw'.pdf", replace
										graph close
									}
								}
							}
						}
						else {
							if ("`chistore'"!="") {
								reg `depvar' `indepvar' `ctrl' `cate_absr' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								/* standard model with FEs absorbed */
								eststo: reghdfe `depvar' `indepvar' `ctrl' if `touse', absorb(`absr') vce(`cse')
							}
						}
					}
					local indepidx = `indepidx' + 1
				}
				local depidx = `depidx' + 1
			}
		}
		else {
			local depidx = 1
			foreach depvar in `anything' {
				local indepidx = 1
				foreach indepvar in `indep' {
					if ("`inte'"!="") {
						if ("`chistore'"!="") {
							reg `depvar' `inte_reglist`indepidx'' `ctrl' if `touse'
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							if ("`tobit'"!="") {
								eststo: tobit `depvar' `inte_reglist`indepidx'' `ctrl' if `touse', `tobit_bound' vce(`cse')
							}
							else {
								eststo: reg `depvar' `inte_reglist`indepidx'' `ctrl' if `touse', vce(`cse')
							}
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								if ("`tobit'"!="") {
									eststo: tobit `depvar' `dyn_reglist`indepidx'' `ctrl' if `touse', `tobit_bound' vce(`cse')
								}
								else {
									eststo: reg `depvar' `dyn_reglist`indepidx'' `ctrl' if `touse', vce(`cse')
								}
								if (`if_dyn_plot'==1) {
									quietly {
										if ("`tobit'"!="") {
											tobit `depvar' `dyn_plotlist`indepidx'' `ctrl' if `touse', `tobit_bound' vce(`cse')
										}
										else {
											reg `depvar' `dyn_plotlist`indepidx'' `ctrl' if `touse', vce(`cse')
										}
										coefplot, keep(c.*#c.*) vertical omitted base levels(`dyn_ci') ///
										rename(`"`dyn_rename'"', regex) order(`"`dyn_reorder'"') ///
										mcolor("0 191 255") ciopts(lcolor("0 97 154") recast(rcap)) ///
										yline(0, lcolor("179 134 0") lpattern(dash)) ///
										graphregion(fcolor(white) ifcolor(white) ilcolor(white)) ///
										xtitle(`dyn_plot') ytitle("Coefficient estimates") ///
										ysize(5) xsize(8)
										graph export "`dyn_pdir'/`indepidx'_`depvar'_`dyn_kw'.pdf", replace
										graph close
									}
								}
							}
						}
						else {
							if ("`chistore'"!="") {
								reg `depvar' `indepvar' `ctrl' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								if ("`tobit'"!="") {
									eststo: tobit `depvar' `indepvar' `ctrl' if `touse', `tobit_bound' vce(`cse')
								}
								else {
									eststo: reg `depvar' `indepvar' `ctrl' if `touse', vce(`cse')
								}
							}
						}
					}
					local indepidx = `indepidx' + 1
				}
				local depidx = `depidx' + 1
			}
		}
	}
	
	/* Store results */
	if ("`display'"=="") {
	
		/* Question: how many indep vars? How many columns? */
		/* (1) Order */
		if ("`inte'"!="") {
			local var_order = "`indep' `orderinte' c.*"
		}
		else {
			if ("`dyn'"!="") {
				local var_order = "`indep' `orderinte'"
			}
			local var_order = "`indep'"
		}
		
		/* (2) Column names: `" `colname' "' */
		local jnum : word count `anything'
		local inum : word count `indep'
		local colnum = `jnum' * `inum'
		
		local colname = ""
		local dep_idx = 1
		forval i = 1/`colnum' {
			
			if (mod(`i'+`inum'-1, `inum')==0) {
				local labname = "`: var label `: word `dep_idx' of `anything'''"
				if ("`labname'"=="") {
					local labname = "`: word `dep_idx' of `anything''"
				}
				local colname = `" `colname' "`labname'" "'
				local dep_idx = `dep_idx' + 1
			}
			else {
				local colname = `" `colname' "" "'
			}
		}
		
		/* (3) Table title: "`table_title'" */
		if ("`ttitle'"=="") {
			local table_title = "Dep [`: word 1 of `anything''] ~ indep [`: word 1 of `indep''] `addn' `extra_addn'"
		}
		else {
			local table_title = "`ttitle' `extra_addn'"
		}
		
		/* (4) Notes: `table_note' */
		if ("`tnote'"=="") {
			if ("`ctrl'"!="" & "`absr'"!="") {
				local table_note = "Controls: [`ctrl']; FE: [`absr']"
			}
			if ("`ctrl'"=="" & "`absr'"!="") {
				local table_note = "No extra controls; FE: [`absr']"
			}
			if ("`ctrl'"!="" & "`absr'"=="") {
				local table_note = "Controls: [`ctrl']; No FEs"
			}
			if ("`ctrl'"=="" & "`absr'"=="") {
				local table_note = "No extra controls; No FEs"
			}
		}
		else {
			local table_note = "`tnote'"
		}
		
		/* (5) Additional labels: `add_stat' `add_label' */
		if ("`tobit'"!="") {
			local add_stat = "N_lc N_unc N_rc"
			local add_label = `" "L Censored" "Uncensored" "R Censored""'
		}
		
		/* If export full variable reports */
		if ("`report'"!="") {
			/* no drop */
			esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
			stats(N r2_a `add_stat', labels("Observations" "Adjusted R-squared" `add_label')) rename(_cons "Constant") ///
			order(`var_order' `ctrl' "Constant") mtitles(`colname') title("`table_title'") note("`table_note'")
		}
		else {
			if ("`keepvar'"!="") {
				foreach kpv in `keepvar' {
					local drop_ctrl : subinstr local drop_ctrl "`kpv'" "", all
				}
				/* keep selected variables */
				esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
				drop(`drop_ctrl') stats(N r2_a `add_stat', labels("Observations" "Adjusted R-squared" `add_label')) rename(_cons "Constant") ///
				order(`var_order' `keepvar' "Constant") mtitles(`colname') title("`table_title'") note("`table_note'")			
			}
			else {
				/* drop control */
				esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
				drop(`drop_ctrl') stats(N r2_a `add_stat', labels("Observations" "Adjusted R-squared" `add_label')) rename(_cons "Constant") ///
				order(`var_order' "Constant") mtitles(`colname') title("`table_title'") note("`table_note'")
			}
		}
	}
	
	eststo clear
	
end

///=============================================================================
/* eqx -> Chi-square coefficient equality test
 * Within model test: eqx y, ... eqt(x1==x2)
 * Between category test (models separated by a categorical variable): (1) eqx y, ... eqt(x | dummy), (2) eqx y, ... eqt(~ | dummy)
 * Between model test (models separated by two sets of dependent variables): (1) eqx y1, ... eqt(x @) dep2(y2), (2) eqx y1, ... eqt(~ @) dep2(y2)
 */
///=============================================================================

capture program drop eqx
program eqx

	syntax anything [if] [in] , indep(namelist) eqt(string) [ctrl(string) absr(string) xtp(string) ///
	inte(string) clust(namelist) dyn(string) dep2(string) ///
	tnote(string) ttitle(string) addn(string) edir(string) REPORT EXPORT]
	
	marksample touse
	
	/* restore arguments for regx */
	local fullopt_args = ""
	foreach opt in ctrl absr xtp inte clust dyn tnote ttitle edir report {
		if ("``opt''"!="") {
			local `opt'_arg = "`opt'(``opt'')"
			local fullopt_args = "`fullopt_args' ``opt'_arg'"
		}
	}
	
	local deplen : word count `anything'
	local indeplen : word count `indep'
	local colidxs = `deplen' * `indeplen'
	
	/* (1) If cluster is specified */
	if ("`clust'"!="") {
		local cse = "cl `clust'"
	}
	else {
		/* (2) If panel is claimed & no cluster specified, cluster SE at id level */
		if ("`xtp'"!="") {
			local cse = "cl `: word 1 of `xtp''"
		}
		else {
			/* (3) If both panel and cluster are not specified */
			if ("${clustervar}"!="") {
				local cse = "cl ${clustervar}"
			}
			else {
				local cse = "robust"
				local extra_addn = "[No cluster SE specified]"
			}
		}
	}
	/* Default export file: ${dir_table_flow} */
	/* Export file priority: `edir' > ${dir_table_flow} */
	if (`"`edir'"'!="") {
		local exportfile = "`edir'"
	}
	else {
		if ("${dir_table_flow}"=="") {
			di "[ERROR] Need to set a directory in [edir] or global ${dir_table_flow}"
			exit
		}
		else {
			local exportfile = "${dir_table_flow}"
		}
	}
	
	/* Situation 1: comparing the coefficients of two variables */
	if (strpos("`eqt'", "==") & !strpos("`eqt'", "|")) {
		local eqvlist : subinstr local eqt " " "", all
		local eqvlist : subinstr local eqvlist "==" " ", all
		local lhsvar : word 1 of `eqvlist'
		local rhsvar : word 2 of `eqvlist'
		
		/* export the results at the same time */
		if ("`export'"!="") {
			regx `anything' if `touse', indep(`indep') `fullopt_args' keepvar(`lhsvar' `rhsvar')
		}
	
		/* estimates store results */
		quietly: regx `anything' if `touse', indep(`indep') `fullopt_args' display chistore
		
		local cellnames = ""
		local sigidx = 0
		local colidx = 1
		forval depidx = 1/`deplen' {
			forval indepidx = 1/`indeplen' {
				matrix c`colidx' = J(1, 3, 0)
				matrix colnames c`colidx' = "Chi-sq" "P-value" "Sig(0/*/**/***)"
				
				quietly: suest y`depidx'_x`indepidx', vce(`cse')
				test [mean]`lhsvar' = [mean]`rhsvar'
				
				matrix c`colidx'[1, 1] = r(chi2)
				matrix c`colidx'[1, 2] = r(p)
				if (r(p)<=0.01) {
					local sigstar = 3.3333
					local sigidx = `sigidx' + 1
				}
				else if (r(p)>0.01 & r(p)<=0.05) {
					local sigstar = 2.2222
					local sigidx = `sigidx' + 1
				}
				else if (r(p)>0.05 & r(p)<=0.1) {
					local sigstar = 1.1111
					local sigidx = `sigidx' + 1
				}
				else {
					local sigstar = 0.0000
				}
				matrix c`colidx'[1, 3] = `sigstar'
				
				local cellnames = "`cellnames' c`colidx'(fmt(%9.4f))"
				local colidx = `colidx' + 1
			}
		}
		
		forval colidx = 1/`colidxs' {
			estadd matrix c`colidx'
		}
		esttab using "`exportfile'", cells("`cellnames'") append nolines not se compress nogaps noobs plain ///
		star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`lhsvar']==[`rhsvar'] `extra_addn'") ///
		mtitle("[`sigidx' out of `colidxs' columns have significant difference]")
	}
	/* Situation 2: comparing the coefficients of independent variables between subsamples */
	else if (strpos("`eqt'", "|") & !strpos("`eqt'", "==")) {
		local eqvlist : subinstr local eqt " " "", all
		local eqvlist : subinstr local eqvlist "|" " ", all
		local eqtvar : word 1 of `eqvlist'
		local subsvar : word 2 of `eqvlist'
		levelsof `subsvar', local(subsvar_list)
		
		local subvidx = 1
		foreach subv in `subsvar_list' {
			local numeric_value = real("`subv'")
			/* export tables */
			if ("`export'"!="") {
				if !missing(`numeric_value') {
					regx `anything' if `subsvar'==`subv' & `touse', indep(`indep') `fullopt_args' addn("[`subsvar'==`subv']")
				}
				else {
					regx `anything' if `subsvar'=="`subv'" & `touse', indep(`indep') `fullopt_args' addn("[`subsvar'==`subv']")
				}
			}
			
			/* estimates store results */
			if !missing(`numeric_value') {
				quietly: regx `anything' if `subsvar'==`subv' & `touse', indep(`indep') `fullopt_args' stosuf("_s`subvidx'") display chistore
			}
			else {
				quietly: regx `anything' if `subsvar'=="`subv'" & `touse', indep(`indep') `fullopt_args' stosuf("_s`subvidx'") display chistore
			}
			mata: st_local("_kw_`subv'", "`subvidx'")
			local subvidx = `subvidx' + 1
		}
		
		/* subsample combination */
		local pair_subs = ""
		local len_subsvar : word count `subsvar_list'
		local len_subsvarm = `len_subsvar' - 1
		local subvi = 1
		
		forval subvi = 1/`len_subsvarm' {
			local subvi_next = `subvi' + 1
			forval subvj = `subvi_next'/`len_subsvar' {
				local pair_subs = "`pair_subs' `: word `subvi' of `subsvar_list''&`: word `subvj' of `subsvar_list''"
			}
		}
		/* for each indepvar */
		if ("`eqtvar'"=="~") {
			/* (1) each indepvar * first word of interaction */
			if ("`inte'"!="") {
				local first_inte : word 1 of `inte'
				local first_inte : subinstr local first_inte "#" "#c.", all
				forval indepidx = 1/`indeplen' {
					local eqs_var_x`indepidx' = "c.`: word `indepidx' of `indep''#c.`first_inte'"
				}
			}
			/* (2) each indepvar alone */
			else {
				forval indepidx = 1/`indeplen' {
					local eqs_var_x`indepidx' = "`: word `indepidx' of `indep''"
				}				
			}
			/* comparing each pair of categories */
			foreach pairsub in `pair_subs' {
				local subsvar_grp : subinstr local pairsub "&" " ", all
				local prsub1 : word 1 of `subsvar_grp'
				local prsub2 : word 2 of `subsvar_grp'
				
				local cellnames = ""
				local sigidx = 0
				local colidx = 1
				forval depidx = 1/`deplen' {
					forval indepidx = 1/`indeplen' {
						local eqs_var = "c.`: word `indepidx' of `indep''#c.`first_inte'"
						
						matrix c`colidx' = J(1, 3, 0)
						matrix colnames c`colidx' = "Chi-sq" "P-value" "Sig(0/*/**/***)"
						
						quietly: suest y`depidx'_x`indepidx'_s`_kw_`prsub1'' y`depidx'_x`indepidx'_s`_kw_`prsub2'', vce(`cse')
						test [y`depidx'_x`indepidx'_s`_kw_`prsub1''_mean]`eqs_var_x`indepidx'' = [y`depidx'_x`indepidx'_s`_kw_`prsub2''_mean]`eqs_var_x`indepidx''
						
						matrix c`colidx'[1, 1] = r(chi2)
						matrix c`colidx'[1, 2] = r(p)
						if (r(p)<=0.01) {
							local sigstar = 3.3333
							local sigidx = `sigidx' + 1
						}
						else if (r(p)>0.01 & r(p)<=0.05) {
							local sigstar = 2.2222
							local sigidx = `sigidx' + 1
						}
						else if (r(p)>0.05 & r(p)<=0.1) {
							local sigstar = 1.1111
							local sigidx = `sigidx' + 1
						}
						else {
							local sigstar = 0.0000
						}
						matrix c`colidx'[1, 3] = `sigstar'
						
						local cellnames = "`cellnames' c`colidx'(fmt(%9.4f))"
						local colidx = `colidx' + 1		
					}
				}
				forval colidx = 1/`colidxs' {
					estadd matrix c`colidx'
				}
				esttab using "`exportfile'", cells("`cellnames'") append nolines not se compress nogaps noobs plain ///
				star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqs_var_x`indeplen''], `subsvar' [=`prsub1'] vs [=`prsub2'] `extra_addn'") ///
				mtitle("[`sigidx' out of `colidxs' columns have significant difference]")
			}
		}
		/* for a fixed indepvar */
		else if ("`eqtvar'"!="~" & "`eqtvar'"!="" & `: word count `eqtvar''==1) {
			/* comparing each pair of categories */
			foreach pairsub in `pair_subs' {
				local subsvar_grp : subinstr local pairsub "&" " ", all
				local prsub1 : word 1 of `subsvar_grp'
				local prsub2 : word 2 of `subsvar_grp'
				
				local cellnames = ""
				local sigidx = 0
				local colidx = 1
				forval depidx = 1/`deplen' {
					forval indepidx = 1/`indeplen' {
						local eqs_var = "c.`: word `indepidx' of `indep''#c.`first_inte'"
						
						matrix c`colidx' = J(1, 3, 0)
						matrix colnames c`colidx' = "Chi-sq" "P-value" "Sig(0/*/**/***)"
						
						quietly: suest y`depidx'_x`indepidx'_s`_kw_`prsub1'' y`depidx'_x`indepidx'_s`_kw_`prsub2'', vce(`cse')
						test [y`depidx'_x`indepidx'_s`_kw_`prsub1''_mean]`eqtvar' = [y`depidx'_x`indepidx'_s`_kw_`prsub2''_mean]`eqtvar'
						
						matrix c`colidx'[1, 1] = r(chi2)
						matrix c`colidx'[1, 2] = r(p)
						if (r(p)<=0.01) {
							local sigstar = 3.3333
							local sigidx = `sigidx' + 1
						}
						else if (r(p)>0.01 & r(p)<=0.05) {
							local sigstar = 2.2222
							local sigidx = `sigidx' + 1
						}
						else if (r(p)>0.05 & r(p)<=0.1) {
							local sigstar = 1.1111
							local sigidx = `sigidx' + 1
						}
						else {
							local sigstar = 0.0000
						}
						matrix c`colidx'[1, 3] = `sigstar'
						
						local cellnames = "`cellnames' c`colidx'(fmt(%9.4f))"
						local colidx = `colidx' + 1		
					}
				}
				forval colidx = 1/`colidxs' {
					estadd matrix c`colidx'
				}
				esttab using "`exportfile'", cells("`cellnames'") append nolines not se compress nogaps noobs plain ///
				star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqtvar'], `subsvar' [=`prsub1'] vs [=`prsub2'] `extra_addn'") ///
				mtitle("[`sigidx' out of `colidxs' columns have significant difference]")
			}
		}
	}
	/* Situation 3: comparing the coefficients of independent variables between regressions with different dependent variables */
    else if (strpos("`eqt'", "@") & "`dep2'"!="" & `deplen'==`: word count `dep2'') {
        local eqvlist : subinstr local eqt " " "", all
		local eqvlist : subinstr local eqvlist "@" " ", all
        local eqtvar : word 1 of `eqvlist'

        /* export the results at the same time */
        if ("`export'"!="") {
			regx `anything' if `touse', indep(`indep') `fullopt_args' keepvar(`eqtvar')
            regx `dep2' if `touse', indep(`indep') `fullopt_args' keepvar(`eqtvar')
		}

        /* estimates store results */
        quietly: regx `anything' if `touse', indep(`indep') `fullopt_args' stosuf("_s1") display chistore
        quietly: regx `dep2' if `touse', indep(`indep') `fullopt_args' stosuf("_s2") display chistore

        /* for each indepvar */
        if ("`eqtvar'"=="~") {
			/* (1) each indepvar * first word of interaction */
			if ("`inte'"!="") {
				local first_inte : word 1 of `inte'
				local first_inte : subinstr local first_inte "#" "#c.", all
				forval indepidx = 1/`indeplen' {
					local eqs_var_x`indepidx' = "c.`: word `indepidx' of `indep''#c.`first_inte'"
				}
			}
			/* (2) each indepvar alone */
			else {
				forval indepidx = 1/`indeplen' {
					local eqs_var_x`indepidx' = "`: word `indepidx' of `indep''"
				}				
			}
            local cellnames = ""
            local sigidx = 0
            local colidx = 1
            forval depidx = 1/`deplen' {
                forval indepidx = 1/`indeplen' {
                    local eqs_var = "c.`: word `indepidx' of `indep''#c.`first_inte'"

                    matrix c`colidx' = J(1, 3, 0)
                    matrix colnames c`colidx' = "Chi-sq" "P-value" "Sig(0/*/**/***)"

                    quietly: suest y`depidx'_x`indepidx'_s1 y`depidx'_x`indepidx'_s2, vce(`cse')
                    test [y`depidx'_x`indepidx'_s1_mean]`eqs_var_x`indepidx'' = [y`depidx'_x`indepidx'_s2_mean]`eqs_var_x`indepidx''

                    matrix c`colidx'[1, 1] = r(chi2)
                    matrix c`colidx'[1, 2] = r(p)

                    if (r(p)<=0.01) {
                        local sigstar = 3.3333
                        local sigidx = `sigidx' + 1
                    }
                    else if (r(p)>0.01 & r(p)<=0.05) {
                        local sigstar = 2.2222
                        local sigidx = `sigidx' + 1
                    }
                    else if (r(p)>0.05 & r(p)<=0.1) {
                        local sigstar = 1.1111
                        local sigidx = `sigidx' + 1
                    }
                    else {
                        local sigstar = 0.0000
                    }
                    matrix c`colidx'[1, 3] = `sigstar'
                    
                    local cellnames = "`cellnames' c`colidx'(fmt(%9.4f))"
                    local colidx = `colidx' + 1		
                }
            }
            forval colidx = 1/`colidxs' {
                estadd matrix c`colidx'
            }
            esttab using "`exportfile'", cells("`cellnames'") append nolines not se compress nogaps noobs plain ///
            star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqs_var_x`indeplen''], Dep Var [=`: word 1 of `anything''] vs [=`: word 1 of `dep2''] `extra_addn'") ///
            mtitle("[`sigidx' out of `colidxs' columns have significant difference]")
        }
        else if ("`eqtvar'"!="~" & "`eqtvar'"!="" & `: word count `eqtvar''==1) {
            local cellnames = ""
            local sigidx = 0
            local colidx = 1
            forval depidx = 1/`deplen' {
                forval indepidx = 1/`indeplen' {
                    local eqs_var = "c.`: word `indepidx' of `indep''#c.`first_inte'"

                    matrix c`colidx' = J(1, 3, 0)
                    matrix colnames c`colidx' = "Chi-sq" "P-value" "Sig(0/*/**/***)"

                    quietly: suest y`depidx'_x`indepidx'_s1 y`depidx'_x`indepidx'_s2, vce(`cse')
                    test [y`depidx'_x`indepidx'_s1_mean]`eqtvar' = [y`depidx'_x`indepidx'_s2_mean]`eqtvar'

                    matrix c`colidx'[1, 1] = r(chi2)
                    matrix c`colidx'[1, 2] = r(p)

                    if (r(p)<=0.01) {
                        local sigstar = 3.3333
                        local sigidx = `sigidx' + 1
                    }
                    else if (r(p)>0.01 & r(p)<=0.05) {
                        local sigstar = 2.2222
                        local sigidx = `sigidx' + 1
                    }
                    else if (r(p)>0.05 & r(p)<=0.1) {
                        local sigstar = 1.1111
                        local sigidx = `sigidx' + 1
                    }
                    else {
                        local sigstar = 0.0000
                    }
                    matrix c`colidx'[1, 3] = `sigstar'
                    
                    local cellnames = "`cellnames' c`colidx'(fmt(%9.4f))"
                    local colidx = `colidx' + 1		
                }
            }
            forval colidx = 1/`colidxs' {
                estadd matrix c`colidx'
            }
            esttab using "`exportfile'", cells("`cellnames'") append nolines not se compress nogaps noobs plain ///
            star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqtvar'], Dep Var [=`: word 1 of `anything''] vs [=`: word 1 of `dep2''] `extra_addn'") ///
            mtitle("[`sigidx' out of `colidxs' columns have significant difference]")
        }
    }
	else {
		di "[ERROR] Wrong input"
		exit
	}

	ereturn clear
	estimates clear
end


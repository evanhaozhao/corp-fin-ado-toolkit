///=============================================================================
/*
 * Project: Programs for corporate finance empirical studies
 * Author: Hao Zhao
 * Created: August 19, 2023
 * Modified: February 19, 2024
 * Version
 	- regx: 1.4.1 (4feb2024)
	- eqx: 2.1.1 (5feb2024)
	- sumx: 1.2.0 (19feb2024)
 */
///=============================================================================
/* regx -> regressions to output tables */
///=============================================================================
capture program drop regx
program regx, rclass sortpreserve

	syntax anything [if] [in] , indep(namelist) [ctrl(string) absr(string) xtp(string) ///
	inte(string) clust(namelist) dyn(string) rotctrl(string) tobit(string) ///
	tnote(string) ttitle(string) addn(string) edir(string) keepvar(string) ///
	stosuf(string) sigout(string) sigkw(string) SIGMAT REPORT DISPLAY CHISTORE]
	
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
	/* Rotating control variables */
	if ("`rotctrl'"!="") {
		if (`: word count `rotctrl''!=`: word count `indep'') {
			di "[ERROR] [rotctrl] number of rotating control vars must equal number of indep vars"
			exit
		}
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
			
			/* default setting: year: [-2, 2]; benchmark: [0]; plot: no; CI: 95% */
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
							reg `depvar' `inte_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', absorb(`: word 1 of `xtp'')
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							eststo: xtreg `depvar' `inte_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `xteffect' vce(`cse')
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', absorb(`: word 1 of `xtp'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								eststo: xtreg `depvar' `dyn_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `xteffect' vce(`cse')
								if (`if_dyn_plot'==1) {
									quietly {
										xtreg `depvar' `dyn_plotlist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `xteffect' vce(`cse')
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
								reg `depvar' `indepvar' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', absorb(`: word 1 of `xtp'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								eststo: xtreg `depvar' `indepvar' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `xteffect' vce(`cse')
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
							reg `depvar' `inte_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' `cate_absr' if `touse'
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							/* Interaction model */
							eststo: reghdfe `depvar' `inte_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', absorb(`absr') vce(`cse')
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' `cate_absr' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								/* dynamic effect model & plot */
								eststo: reghdfe `depvar' `dyn_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', absorb(`absr') vce(`cse')
								if (`if_dyn_plot'==1) {
									quietly {
										reghdfe `depvar' `dyn_plotlist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', absorb(`absr') vce(`cse')
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
								reg `depvar' `indepvar' `ctrl' `: word `indepidx' of `rotctrl'' `cate_absr' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								/* standard model with FEs absorbed */
								eststo: reghdfe `depvar' `indepvar' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', absorb(`absr') vce(`cse')
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
							reg `depvar' `inte_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse'
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							if ("`tobit'"!="") {
								eststo: tobit `depvar' `inte_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `tobit_bound' vce(`cse')
							}
							else {
								eststo: reg `depvar' `inte_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', vce(`cse')
							}
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								if ("`tobit'"!="") {
									eststo: tobit `depvar' `dyn_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `tobit_bound' vce(`cse')
								}
								else {
									eststo: reg `depvar' `dyn_reglist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', vce(`cse')
								}
								if (`if_dyn_plot'==1) {
									quietly {
										if ("`tobit'"!="") {
											tobit `depvar' `dyn_plotlist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `tobit_bound' vce(`cse')
										}
										else {
											reg `depvar' `dyn_plotlist`indepidx'' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', vce(`cse')
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
								reg `depvar' `indepvar' `ctrl' `: word `indepidx' of `rotctrl'' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								if ("`tobit'"!="") {
									eststo: tobit `depvar' `indepvar' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', `tobit_bound' vce(`cse')
								}
								else {
									eststo: reg `depvar' `indepvar' `ctrl' `: word `indepidx' of `rotctrl'' if `touse', vce(`cse')
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
			order(`var_order' `rotctrl' `ctrl' "Constant") mtitles(`colname') title("`table_title'") note("`table_note'")
		}
		else {
			if ("`keepvar'"!="") {
				foreach kpv in `keepvar' {
					local drop_ctrl : subinstr local drop_ctrl "`kpv'" "", all
				}
				/* keep selected variables */
				esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
				drop(`drop_ctrl') stats(N r2_a `add_stat', labels("Observations" "Adjusted R-squared" `add_label')) rename(_cons "Constant") ///
				order(`var_order' `rotctrl' `keepvar' "Constant") mtitles(`colname') title("`table_title'") note("`table_note'")			
			}
			else {
				/* drop control */
				esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
				drop(`drop_ctrl') stats(N r2_a `add_stat', labels("Observations" "Adjusted R-squared" `add_label')) rename(_cons "Constant") ///
				order(`var_order' `rotctrl' "Constant") mtitles(`colname') title("`table_title'") note("`table_note'")
			}
		}
	}
		
	/* Summarize coefficient matrix */
	if ("`sigmat'"!="" | "`sigout'"!="") {
		quietly: esttab
		matrix A = r(coefs)
		local deplen : word count `anything'
		local indeplen : word count `indep'
		
		/* For defining the matrix width and colnames (using the first regression) */
		local col_items_f = ""
		forval depidx = 1/1 {
			forval indepidx = 1/1 {
				if ("`inte'"!="") {
					local col_items_f = "`inte_reglist`indepidx'' `: word `indepidx' of `rotctrl'' `keepvar'"
				}
				else {
					if ("`dyn'"!="") {
						local col_items_f = "`dyn_reglist`indepidx'' `: word `indepidx' of `rotctrl'' `keepvar'"
					}
					else {
						local col_items_f = "`: word `indepidx' of `indep'' `: word `indepidx' of `rotctrl'' `keepvar'"
					}
				}
				local sigmat_col_num : word count `col_items_f'
				local sigmat_col_len = `sigmat_col_num' * 2
			}
		}
		matrix B = J(1, `sigmat_col_len', 0)
		local b_col_name = `""'
		foreach col in `col_items_f' {
			local b_col_name = `"`b_col_name' "`col'" "(sig)""'
		}
		matrix colnames B = `b_col_name'
		
		forval b_col_idx = 1/`sigmat_col_num' {
			local overall_sign = 0
			local overall_sig = 0
			local c_sig_cnt = 0
			local agg_sign = 0
			local col_items = ""
			forval depidx = 1/`deplen' {
				local r_sig_cnt = 0
				forval indepidx = 1/`indeplen' {
					local yb_col = 3 * (`indepidx'-`indeplen'+`indeplen'*`depidx')-2
					local yp_col = 3 * (`indepidx'-`indeplen'+`indeplen'*`depidx')
					if ("`inte'"!="") {
						local col_items = "`inte_reglist`indepidx'' `: word `indepidx' of `rotctrl'' `keepvar'"
					}
					else {
						if ("`dyn'"!="") {
							local col_items = "`dyn_reglist`indepidx'' `: word `indepidx' of `rotctrl'' `keepvar'"
						}
						else {
							local col_items = "`: word `indepidx' of `indep'' `: word `indepidx' of `rotctrl'' `keepvar'"
						}
					}
					local cell_word : word `b_col_idx' of `col_items'
					if ("`tobit'"!="") {
						local cell_coef = A["main:`cell_word'", `yb_col']
					}
					else {
						local cell_coef = A["`cell_word'", `yb_col']
					}
					if (`cell_coef' < 0) {
						local agg_sign = `agg_sign' - 1
					}
					else if (`cell_coef' > 0) {
						local agg_sign = `agg_sign' + 1
					}
					if ("`tobit'"!="") {
						local cell_sig_cnt = A["main:`cell_word'", `yp_col']
					}
					else {
						local cell_sig_cnt = A["`cell_word'", `yp_col']
					}
					if (`cell_sig_cnt' <= 0.1) {
						local r_sig_cnt = `r_sig_cnt' + 1
					}
				}
				local c_sig_cnt = `c_sig_cnt' + `r_sig_cnt'/`indeplen'
			}
			if (`agg_sign' == `deplen'*`indeplen') {
				local overall_sign = 1
			}
			else if (`agg_sign' == -1*`deplen'*`indeplen') {
				local overall_sign = -1
			}
			else {
				local overall_sign = 0
			}
			local c_sig_score = `c_sig_cnt'/`deplen'
			if (`c_sig_score'<1/3) {
				local overall_sig = 0
			}
			else if (`c_sig_score'>2/3) {
				local overall_sig = 3
			}
			else {
				local overall_sig = 2
			}
			matrix B[1, 2*`b_col_idx'-1] = `overall_sign'
			matrix B[1, 2*`b_col_idx'] = `overall_sig'
		}
		
		/* If export the significance summary into a file */
 		if ("`sigout'"!="") {
 			eststo clear
 			ereturn clear
 			estimates clear
			
 			local sigsummaryfile : subinstr local exportfile ".csv" "_ss.csv", all
 			local sigsummaryfile : subinstr local sigsummaryfile ".rtf" "_ss.rtf", all

			local cellname_li = ""
			forval col_idx = 1/`sigmat_col_len' {
				matrix B_`col_idx' = B[1, `col_idx']
				matrix colnames B_`col_idx' = "model `sigout'"
				estadd matrix B_`col_idx'
				local cellname_li = "`cellname_li' B_`col_idx'(fmt(%9.0f))"
			}
			if ("`sigout'"=="1" | "`sigout'"=="baseline" | "`sigout'"=="`sigkw'") {
				if ("`ttitle'"=="") {
					local table_title = "Dep [`: word 1 of `anything''] ~ indep [`: word 1 of `indep''] `addn' `extra_addn'"
				}
				else {
					local table_title = "`ttitle' `extra_addn'"
				}
				esttab using "`sigsummaryfile'", cells("`cellname_li'") ///
				collabels(`b_col_name') mlabels(,none) eqlab(,none) title(`table_title') ///
				append nolines not se compress nogaps noobs plain
			}
			else {
				esttab using "`sigsummaryfile'", cells("`cellname_li'") ///
				collabels(,none) mlabels(,none) eqlab(,none) ///
				append nolines not se compress nogaps noobs plain
			}
 		}
		return matrix sigmat = B
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
	inte(string) clust(namelist) dyn(string) rotctrl(string) dep2(string) ///
	tnote(string) ttitle(string) addn(string) edir(string) sigout(string) sigkw(string) REPORT EXPORT]
	
	marksample touse
	
	/* restore arguments for regx */
	local fullopt_args = ""
	foreach opt in ctrl absr xtp inte clust dyn rotctrl tnote ttitle edir report {
		if ("``opt''"!="") {
			local `opt'_arg = "`opt'(``opt'')"
			local fullopt_args = "`fullopt_args' ``opt'_arg'"
		}
	}
 	/* significant output options */
	local sigopt_args = ""
	foreach opt in sigout {
		if ("``opt''"!="") {
			local `opt'_arg = "`opt'(``opt'')"
			local sigopt_args = "`sigopt_args' ``opt'_arg'"
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
			if ("`sigkw'"!="") {
				local sigopt_args = "`sigopt_args' sigkw(`sigkw')"
			}
			regx `anything' if `touse', indep(`indep') `fullopt_args' `sigopt_args' keepvar(`lhsvar' `rhsvar')
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

		/* Chi-sq significance table */
		if ("`sigopt_args'"!="") {
			eststo clear
			ereturn clear
			estimates clear

			local sigsummaryfile : subinstr local exportfile ".csv" "_ss.csv", all
			local sigsummaryfile : subinstr local sigsummaryfile ".rtf" "_ss.rtf", all

			local cellname_li = ""
			foreach sigcnt in sigidx colidxs {
				matrix mat`sigcnt' = J(1, 1, 0)
				matrix colnames mat`sigcnt' = "Diff"
				matrix mat`sigcnt'[1, 1] = ``sigcnt''
				estadd matrix mat`sigcnt'
				local cellname_li = "`cellname_li' mat`sigcnt'(fmt(%9.0f))"
			}
			matrix matperc = J(1, 1, 0)
			matrix colnames matperc = "Diff"
			matrix matperc[1, 1] = round(`sigidx'/`colidxs', 0.001)
			estadd matrix matperc
			local cellname_li = "`cellname_li' matperc(fmt(%9.3f))"

			esttab using "`sigsummaryfile'", cells("`cellname_li'") ///
			collabels(,none) mlabels(,none) eqlab(,none) ///
			append nolines not se compress nogaps noobs plain
		}
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
				if ("`sigkw'"!="") {
					if (`subvidx'==1) {
						local sigopt_args = "sigout(`sigout') sigkw(`sigkw')"
					}
					else {
						local sigopt_args = "sigout(`sigout'(`subv')) sigkw(`sigkw')"
					}
				}
				else {
					if (`subvidx'==1) {
						local sigopt_args = "sigout(`sigout')"
					}
					else {
						local sigopt_args = "sigout(`sigout'(`subv'))"
					}
				}
				if !missing(`numeric_value') {
					regx `anything' if `subsvar'==`subv' & `touse', indep(`indep') `fullopt_args' `sigopt_args' addn("[`subsvar'==`subv']")
				}
				else {
					regx `anything' if `subsvar'=="`subv'" & `touse', indep(`indep') `fullopt_args' `sigopt_args' addn("[`subsvar'==`subv']")
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

				/* Chi-sq significance table */
				if ("`sigopt_args'"!="") {
					eststo clear
					ereturn clear
					estimates clear

					local sigsummaryfile : subinstr local exportfile ".csv" "_ss.csv", all
					local sigsummaryfile : subinstr local sigsummaryfile ".rtf" "_ss.rtf", all

					local cellname_li = ""
					foreach sigcnt in sigidx colidxs {
						matrix mat`sigcnt' = J(1, 1, 0)
						matrix colnames mat`sigcnt' = "Diff"
						matrix mat`sigcnt'[1, 1] = ``sigcnt''
						estadd matrix mat`sigcnt'
						local cellname_li = "`cellname_li' mat`sigcnt'(fmt(%9.0f))"
					}
					matrix matperc = J(1, 1, 0)
					matrix colnames matperc = "Diff"
					matrix matperc[1, 1] = round(`sigidx'/`colidxs', 0.001)
					estadd matrix matperc
					local cellname_li = "`cellname_li' matperc(fmt(%9.3f))"

					esttab using "`sigsummaryfile'", cells("`cellname_li'") ///
					collabels(,none) mlabels(,none) eqlab(,none) ///
					append nolines not se compress nogaps noobs plain
				}
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

				/* Chi-sq significance table */
				if ("`sigopt_args'"!="") {
					eststo clear
					ereturn clear
					estimates clear

					local sigsummaryfile : subinstr local exportfile ".csv" "_ss.csv", all
					local sigsummaryfile : subinstr local sigsummaryfile ".rtf" "_ss.rtf", all

					local cellname_li = ""
					foreach sigcnt in sigidx colidxs {
						matrix mat`sigcnt' = J(1, 1, 0)
						matrix colnames mat`sigcnt' = "Diff"
						matrix mat`sigcnt'[1, 1] = ``sigcnt''
						estadd matrix mat`sigcnt'
						local cellname_li = "`cellname_li' mat`sigcnt'(fmt(%9.0f))"
					}
					matrix matperc = J(1, 1, 0)
					matrix colnames matperc = "Diff"
					matrix matperc[1, 1] = round(`sigidx'/`colidxs', 0.001)
					estadd matrix matperc
					local cellname_li = "`cellname_li' matperc(fmt(%9.3f))"

					esttab using "`sigsummaryfile'", cells("`cellname_li'") ///
					collabels(,none) mlabels(,none) eqlab(,none) ///
					append nolines not se compress nogaps noobs plain
				}
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
			if ("`sigkw'"!="") {
				local sigopt_args_a = "sigout(`sigout') sigkw(`sigkw')"
				local sigopt_args_b = "sigout(`sigout'(pair)) sigkw(`sigkw')"
			}
			else {
				local sigopt_args_a = "sigout(`sigout')"
				local sigopt_args_b = "sigout(`sigout'(pair))"
			}
			if ("`eqtvar'"=="~") {
				regx `anything' if `touse', indep(`indep') `fullopt_args' `sigopt_args_a'
				regx `dep2' if `touse', indep(`indep') `fullopt_args' `sigopt_args_b'
			}
			else {
				regx `anything' if `touse', indep(`indep') `fullopt_args' `sigopt_args_a' keepvar(`eqtvar')
				regx `dep2' if `touse', indep(`indep') `fullopt_args' `sigopt_args_b' keepvar(`eqtvar')
			}
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

			/* Chi-sq significance table */
			if ("`sigopt_args'"!="") {
				eststo clear
				ereturn clear
				estimates clear

				local sigsummaryfile : subinstr local exportfile ".csv" "_ss.csv", all
				local sigsummaryfile : subinstr local sigsummaryfile ".rtf" "_ss.rtf", all

				local cellname_li = ""
				foreach sigcnt in sigidx colidxs {
					matrix mat`sigcnt' = J(1, 1, 0)
					matrix colnames mat`sigcnt' = "Diff"
					matrix mat`sigcnt'[1, 1] = ``sigcnt''
					estadd matrix mat`sigcnt'
					local cellname_li = "`cellname_li' mat`sigcnt'(fmt(%9.0f))"
				}
				matrix matperc = J(1, 1, 0)
				matrix colnames matperc = "Diff"
				matrix matperc[1, 1] = round(`sigidx'/`colidxs', 0.001)
				estadd matrix matperc
				local cellname_li = "`cellname_li' matperc(fmt(%9.3f))"

				esttab using "`sigsummaryfile'", cells("`cellname_li'") ///
				collabels(,none) mlabels(,none) eqlab(,none) ///
				append nolines not se compress nogaps noobs plain
			}
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

			/* Chi-sq significance table */
			if ("`sigopt_args'"!="") {
				eststo clear
				ereturn clear
				estimates clear

				local sigsummaryfile : subinstr local exportfile ".csv" "_ss.csv", all
				local sigsummaryfile : subinstr local sigsummaryfile ".rtf" "_ss.rtf", all

				local cellname_li = ""
				foreach sigcnt in sigidx colidxs {
					matrix mat`sigcnt' = J(1, 1, 0)
					matrix colnames mat`sigcnt' = "Diff"
					matrix mat`sigcnt'[1, 1] = ``sigcnt''
					estadd matrix mat`sigcnt'
					local cellname_li = "`cellname_li' mat`sigcnt'(fmt(%9.0f))"
				}
				matrix matperc = J(1, 1, 0)
				matrix colnames matperc = "Diff"
				matrix matperc[1, 1] = round(`sigidx'/`colidxs', 0.001)
				estadd matrix matperc
				local cellname_li = "`cellname_li' matperc(fmt(%9.3f))"

				esttab using "`sigsummaryfile'", cells("`cellname_li'") ///
				collabels(,none) mlabels(,none) eqlab(,none) ///
				append nolines not se compress nogaps noobs plain
			}
        }
    }
	else {
		di "[ERROR] Wrong input"
		exit
	}

	ereturn clear
	estimates clear
end

///=============================================================================
/* sumx -> Generate summary statistics tables */
///=============================================================================

capture program drop sumx
program sumx, sortpreserve

	syntax anything [if] [in], deci(numlist) edir(string) [category(string) ttitle(string) addn(string) RLABEL]
	
	marksample touse

	/* Table title: "`table_title'" */
	if ("`ttitle'" == "") {
		local table_title = "Summary statistics `addn'"
	}
	else {
		local table_title = "`ttitle'"
	}
	
	local var_len : word count `anything'
	
	if ("`category'" == "") {
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
	}
	else if ("`category'" != "") {
		if `var_len' != 1 {
			di "[ERROR] Only one input variable allowed"
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
				quietly: summ `anything' if `touse' & `category'=="`group'", detail
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


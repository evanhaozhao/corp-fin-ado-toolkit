///=============================================================================
/*
 * Project: Programs for corporate finance empirical studies
 * Author: Hao Zhao
 * Created: August 19, 2023
 * Modified: April 19, 2025
 * Version
 	- regx: 1.8.1 (27jun2025)
	- eqx: 2.2.6 (1may2025)
	- sumx: 1.3.4 (19apr2025)
 */
///=============================================================================
/* regx -> regressions to output tables */
///=============================================================================
capture program drop regx
program regx, rclass sortpreserve

	syntax anything [if] [in] , indep(namelist) [ctrl(string) absr(string) xtp(string) ///
	inte(string) clust(namelist) dyn(string) rotctrl(string) tobit(string) ///
	tnote(string) ttitle(string) addn(string) edir(string) keepvar(string) ///
	stosuf(string) sigout(string) sigkw(string) rcoefidx(string) rename(string) parmsetout(string) ///
	SIGMAT RCOEF ROTINTE REPORT DISPLAY CHISTORE NOSINGLETON POISSON NOROUNDDECI ROUNDTOFOUR RWITHIN]
	
	marksample touse

	/* Default export file: ${dir_table_flow} */
	/* Export file priority: `edir' > ${dir_table_flow} */
	if (`"`edir'"'!="") {
		local exportfile = "`edir'"
	}
	else {
		if ("${dir_table_flow}"=="") {
			display in red "[ERROR] Need to set a directory in [edir] or global ${dir_table_flow}"
			exit
		}
		else {
			local exportfile = "${dir_table_flow}"
		}
	}
	
	/* No rounding decimals to 3-digit */
	if "`norounddeci'" == "" {
		if "`roundtofour'" == "" {
			local bdeciopt = `" b(%8.3f) "'
		}
		else {
			local bdeciopt = `" b(%8.4f) "'
		}
	}

	/* singleton setting for reghdfe */
	if ("`nosingleton'"!="") {
		local singletonset = `""'
	}
	else {
		local singletonset = `" keepsingleton "'
	}

	/* `suest` does not support xtreg random effect; the estimation of xtreg fe is done by reg..absorb */
	if ("`chistore'"!="") {
		if ("`xtp'"!="") {
			display in red "[WARNING] `suest` does not support xtreg"
		}
		if ("`display'"=="") {
			display in red "[ERROR] [chistore] must be used with [display]"
			exit
		}
		if ("`absr'"!="") {
			local cate_absr = ""
			local cate_absr_idx = 1
			foreach absrvar in `absr' {
				if (`cate_absr_idx'!=1) {
					local cate_absr = "`cate_absr' ibn.`absrvar'"
				}
				local cate_absr_idx = `cate_absr_idx' + 1
			}
		}
	}
	
	/* Tobit model: cannot use [absr] */
	if ("`tobit'"!="") {
		if ("`absr'"!="" | "`xtp'"!="") {
			display in red "[ERROR] [tobit] cannot be used with [xtp] or [absr]"
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
		local cse = "cluster `clust'"
	}
	else {
		/* (2) If panel is claimed & no cluster specified, cluster SE at id level */
		if ("`xtp'"!="") {
			local cse = "cluster `: word 1 of `xtp''"
		}
		else {
			/* (3) If both panel and cluster are not specified */
			if ("${clustervar}"!="") {
				local cse = "cluster ${clustervar}"
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
			display in red "[ERROR] [rotctrl] number of rotating control vars must equal number of indep vars"
			exit
		}
		else {
			local rctrlidx = 1
			foreach rctrlvar in `rotctrl' {
				local rotctrlvar`rctrlidx' = "`rctrlvar'"
				local rctrlidx = `rctrlidx' + 1
			}
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
						display in red "[ERROR] Need to ssc install `tuples`"
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

		/* interaction term for rotating control variables */
		if ("`rotctrl'"!="" & "`rotinte'"!="") {
			local indepidx = 1
			foreach indepvar in `indep' {
				local rotctrlvar`indepidx' = ""
				foreach inte_item in `inte_reglist`indepidx'' {
					if (strpos("`inte_item'", "`indepvar'")) {
						local ctrl_inte_item : subinstr local inte_item "`indepvar'" "`: word `indepidx' of `rotctrl''", all
						local rotctrlvar`indepidx' = "`rotctrlvar`indepidx'' `ctrl_inte_item'"
					}
				}
				local indepidx = `indepidx' + 1
			}
		}
	}
	
	/* dynamic regression */
	if ("`dyn'"!="") {
		if ("`inte'"!="") {
			display in red "[ERROR] Cannot use option [dyn] and [inte] at the same time"
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
				display in red "[ERROR] Benchmark (`dyn_bench') is out of range [`dyn_year']"
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
					display in red "[ERROR] Need to generate dynamic indicator `dyn_v'"
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
			display in red "[ERROR] No [absorb] option allowed for xtreg"
			exit
		}
		else {
			/* _r_effect: random effect */
			if (strpos("`xtp'", " _r_effect")>0) {
				if ("`chistore'"!="") {
					display in red "[ERROR] `suest` does not support xtreg random effect"
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
							reg `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`: word 1 of `xtp'')
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							eststo: xtreg `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', `xteffect' vce(`cse')
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`: word 1 of `xtp'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								eststo: xtreg `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', `xteffect' vce(`cse')
								if (`if_dyn_plot'==1) {
									quietly {
										xtreg `depvar' `dyn_plotlist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', `xteffect' vce(`cse')
										coefplot, keep(c.*#c.*) vertical omitted base levels(`dyn_ci') ///
										rename(`dyn_rename', regex) order(`dyn_reorder') ///
										mcolor("90 106 115") ciopts(lcolor("90 106 115") recast(rcap)) ///
										yline(0, lcolor("130 0 0") lpattern(dash)) ///
										graphregion(fcolor(white) ifcolor(white) ilcolor(white)) ///
										xtitle(`dyn_plot') ytitle("Coefficient estimates") ///
										ysize(4) xsize(8)
										graph export "`dyn_pdir'/`indepidx'_`depvar'_`dyn_kw'.pdf", replace
										graph close
									}
								}
							}
						}
						else {
							if ("`chistore'"!="") {
								reg `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`: word 1 of `xtp'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								eststo: xtreg `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse', `xteffect' vce(`cse')
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
							/* reg `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' `cate_absr' if `touse' */
							reg `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' `cate_absr' if `touse', absorb(`: word 1 of `absr'')
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							/* Interaction model */
							if ("`poisson'"!="") {
								eststo: ppmlhdfe `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
							}
							else {
								eststo: reghdfe `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
							}
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								/* reg `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' `cate_absr' if `touse' */
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' `cate_absr' if `touse', absorb(`: word 1 of `absr'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								/* dynamic effect model & plot */
								if ("`poisson'"!="") {
									eststo: ppmlhdfe `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
								}
								else {
									eststo: reghdfe `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
								}
								if (`if_dyn_plot'==1) {
									quietly {
										if ("`poisson'"!="") {
											ppmlhdfe `depvar' `dyn_plotlist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
										}
										else {
											reghdfe `depvar' `dyn_plotlist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
										}
										coefplot, keep(c.*#c.*) vertical omitted base levels(`dyn_ci') ///
										rename(`dyn_rename', regex) order(`dyn_reorder') ///
										mcolor("90 106 115") ciopts(lcolor("90 106 115") recast(rcap)) ///
										yline(0, lcolor("130 0 0") lpattern(dash)) ///
										graphregion(fcolor(white) ifcolor(white) ilcolor(white)) ///
										xtitle(`dyn_plot') ytitle("Coefficient estimates") ///
										ysize(4) xsize(8)
										graph export "`dyn_pdir'/`indepidx'_`depvar'_`dyn_kw'.pdf", replace
										graph close
									}
								}
							}
						}
						else {
							if ("`chistore'"!="") {
								/* reg `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' `cate_absr' if `touse' */
								reg `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' `cate_absr' if `touse', absorb(`: word 1 of `absr'')
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								/* standard model with FEs absorbed */
								if ("`poisson'"!="") {
									eststo: ppmlhdfe `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
								}
								else {
									eststo: reghdfe `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse', absorb(`absr') vce(`cse') `singletonset'
								}
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
							reg `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse'
							estimates store y`depidx'_x`indepidx'`stosuf'
						}
						else {
							if ("`tobit'"!="") {
								eststo: tobit `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', `tobit_bound' vce(`cse')
							}
							else if ("`poisson'"!="") {
								eststo: poisson `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
							}
							else {
								eststo: reg `depvar' `inte_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
							}
						}
					}
					else {
						if ("`dyn'"!="") {
							if ("`chistore'"!="") {
								reg `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								if ("`tobit'"!="") {
									eststo: tobit `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', `tobit_bound' vce(`cse')
								}
								else if ("`poisson'"!="") {
									eststo: poisson `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
								}
								else {
									eststo: reg `depvar' `dyn_reglist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
								}
								if (`if_dyn_plot'==1) {
									quietly {
										if ("`tobit'"!="") {
											tobit `depvar' `dyn_plotlist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', `tobit_bound' vce(`cse')
										}
										else if ("`poisson'"!="") {
											poisson `depvar' `dyn_plotlist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
										}
										else {
											reg `depvar' `dyn_plotlist`indepidx'' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
										}
										coefplot, keep(c.*#c.*) vertical omitted base levels(`dyn_ci') ///
										rename(`dyn_rename', regex) order(`dyn_reorder') ///
										mcolor("90 106 115") ciopts(lcolor("90 106 115") recast(rcap)) ///
										yline(0, lcolor("130 0 0") lpattern(dash)) ///
										graphregion(fcolor(white) ifcolor(white) ilcolor(white)) ///
										xtitle(`dyn_plot') ytitle("Coefficient estimates") ///
										ysize(4) xsize(8)
										graph export "`dyn_pdir'/`indepidx'_`depvar'_`dyn_kw'.pdf", replace
										graph close
									}
								}
							}
						}
						else {
							if ("`chistore'"!="") {
								reg `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse'
								estimates store y`depidx'_x`indepidx'`stosuf'
							}
							else {
								if ("`tobit'"!="") {
									eststo: tobit `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse', `tobit_bound' vce(`cse')
								}
								else if ("`poisson'"!="") {
									eststo: poisson `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
								}
								else {
									eststo: reg `depvar' `indepvar' `ctrl' `rotctrlvar`indepidx'' if `touse', vce(`cse')
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
			local var_order = "`indep' `orderinte' `rotctrl' c.*"
		}
		else {
			if ("`dyn'"!="") {
				local var_order = "`indep' `orderinte' `rotctrl'"
			}
			else {
				local var_order = "`indep' `rotctrl'"
			}
		}

		/* Rename explanatory variables to the same name, for horizontal alignment reason */
		if ("`rename'" != "") {
			foreach ordered_var in `var_order' {
				foreach indepvar in `indep' {
					if (strpos("`ordered_var'", "`indepvar'")) {
						local new_name : subinstr local ordered_var "`indepvar'" "`rename'", all
						local rename_str `rename_str' `ordered_var' `new_name'
					}
				}
			}
			/* Adjust ordering */
			if ("`orderinte'" != "") {
				local reorderinte = "`orderinte'"
				foreach reorderinte_var in `reorderinte' {
					foreach indepvar in `indep' {
						if (strpos("`reorderinte_var'", "c.`indepvar'") & !strpos("`reorderinte'", "c.`rename'")) {
							local reorderinte : subinstr local reorderinte "c.`indepvar'" "c.`rename'", all
						}
						else if (strpos("`reorderinte_var'", "c.`indepvar'") & strpos("`reorderinte'", "c.`rename'")) {
							local reorderinte : subinstr local reorderinte "`reorderinte_var'" "", all
						}
					}
				}
			}
			if ("`inte'"!="") {
				local var_order = "`rename' `reorderinte' `rotctrl' c.*"
			}
			else {
				if ("`dyn'"!="") {
					local var_order = "`rename' `reorderinte' `rotctrl'"
				}
				else {
					local var_order = "`rename' `rotctrl'"
				}
			}			
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
			local add_fmt = "%9.3f %9.3f %9.3f"
		}
		
		/* If export full variable reports */
		if ("`rwithin'"=="") {
			local r2setting = "r2_a"
			local r2desc = "Adjusted R-squared"
		}
		else {
			local r2setting = "r2_a_within"
			local r2desc = "Adjusted R-sq (within)"
		}

		if ("`report'"!="") {
			/* no drop */
			esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
			stats(N `r2setting' `add_stat', labels("Observations" "`r2desc'" `add_label') fmt(0 %9.3f `add_fmt')) rename(_cons "Constant" `rename_str') ///
			order(`var_order' `ctrl') mtitles(`colname') title("`table_title'") note("`table_note'") `bdeciopt'
		}
		else {
			if ("`keepvar'"!="") {
				foreach kpv in `keepvar' {
					local drop_ctrl : subinstr local drop_ctrl "`kpv'" "", all
				}
				/* keep selected variables */
				esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
				drop(`drop_ctrl') stats(N `r2setting' `add_stat', labels("Observations" "`r2desc'" `add_label') fmt(0 %9.3f `add_fmt')) rename(_cons "Constant" `rename_str') ///
				order(`var_order' `keepvar') mtitles(`colname') title("`table_title'") note("`table_note'") `bdeciopt'	
			}
			else {
				/* drop control */
				esttab using "`exportfile'", append nolines not se star(* 0.1 ** 0.05 *** 0.01) compress nogaps ///
				drop(`drop_ctrl') stats(N `r2setting' `add_stat', labels("Observations" "`r2desc'" `add_label') fmt(0 %9.3f `add_fmt')) rename(_cons "Constant" `rename_str') ///
				order(`var_order') mtitles(`colname') title("`table_title'") note("`table_note'") `bdeciopt'
			}
		}
	}
	/* if display is enabled */
	else {
		/* Store the last estimates in dta via parmest */
		if "`parmsetout'" != "" {
			di in red "Saving coefficient estimates into path `parmsetout'"
			parmest, saving("`parmsetout'", replace) format(parm %20s estimate %8.3f std %8.3f t %8.3f min95 %8.3f max95 %8.3f p %8.2f)
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
					local col_items_f = "`inte_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
				}
				else {
					if ("`dyn'"!="") {
						local col_items_f = "`dyn_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
					}
					else {
						local col_items_f = "`: word `indepidx' of `indep'' `rotctrlvar`indepidx'' `keepvar'"
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
						local col_items = "`inte_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
					}
					else {
						if ("`dyn'"!="") {
							local col_items = "`dyn_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
						}
						else {
							local col_items = "`: word `indepidx' of `indep'' `rotctrlvar`indepidx'' `keepvar'"
						}
					}
					local cell_word : word `b_col_idx' of `col_items'
					if ("`tobit'"!="" | ("`poisson'"!="" & "`absr'"=="")) {
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
					if ("`tobit'"!="" | ("`poisson'"!="" & "`absr'"=="")) {
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
			if ("`sigout'"=="1" | "`sigout'"=="baseline" | "`sigout'"=="`sigkw'" | "`sigkw'"=="_printtitle_") {
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

	/* Export regression b, t, and p */
	if ("`rcoef'" != "") {
		if ("`sigmat'"!="" | "`sigout'"!="") {
			display in red "Cannot use option sigmat/sigout with option rcoef"
			exit
		}
		else {
			quietly: esttab
			matrix A = r(coefs)
			local deplen : word count `anything'
			local indeplen : word count `indep'
			local rcoef_row_len = `deplen' * `indeplen'

			/* For defining the matrix width and colnames (using the first regression) */
			local col_items_f = ""
			forval depidx = 1/1 {
				forval indepidx = 1/1 {
					if ("`inte'"!="") {
						local col_items_f = "`inte_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
					}
					else {
						if ("`dyn'"!="") {
							local col_items_f = "`dyn_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
						}
						else {
							local col_items_f = "`: word `indepidx' of `indep'' `rotctrlvar`indepidx'' `keepvar'"
						}
					}
					local rcoef_col_num : word count `col_items_f'
					local rcoef_col_len = `rcoef_col_num' * 3
				}
			}

			/* Storing b, t, and p; matrix size: nrow = # of regressions, ncol = # of vars */
			matrix C = J(`rcoef_row_len', `rcoef_col_len', 0)
			local c_col_name = `""'
			foreach col in `col_items_f' {
				local c_col_name = `"`c_col_name' "`col'" "(t)" "(p)""'
			}
			matrix colnames C = `c_col_name'
	
			forval b_col_idx = 1/`rcoef_col_num' {
				local col_items = ""
				local reg_idx = 1
				forval depidx = 1/`deplen' {
					forval indepidx = 1/`indeplen' {
						local yb_col = 3 * (`indepidx'-`indeplen'+`indeplen'*`depidx')-2
						local yt_col = 3 * (`indepidx'-`indeplen'+`indeplen'*`depidx')-1
						local yp_col = 3 * (`indepidx'-`indeplen'+`indeplen'*`depidx')
						if ("`inte'"!="") {
							local col_items = "`inte_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
						}
						else {
							if ("`dyn'"!="") {
								local col_items = "`dyn_reglist`indepidx'' `rotctrlvar`indepidx'' `keepvar'"
							}
							else {
								local col_items = "`: word `indepidx' of `indep'' `rotctrlvar`indepidx'' `keepvar'"
							}
						}
						local cell_word : word `b_col_idx' of `col_items'
						if ("`tobit'"!="" | ("`poisson'"!="" & "`absr'"=="")) {
							local cell_b = A["main:`cell_word'", `yb_col']
							local cell_t = A["main:`cell_word'", `yt_col']
							local cell_p = A["main:`cell_word'", `yp_col']
						}
						else {
							local cell_b = A["`cell_word'", `yb_col']
							local cell_t = A["`cell_word'", `yt_col']
							local cell_p = A["`cell_word'", `yp_col']
						}
						matrix C[`reg_idx', 3*`b_col_idx'-2] = `cell_b'
						matrix C[`reg_idx', 3*`b_col_idx'-1] = `cell_t'
						matrix C[`reg_idx', 3*`b_col_idx'] = `cell_p'

						local reg_idx = `reg_idx' + 1
					}
				}
			}

			/* Output results */
			eststo clear
 			ereturn clear
 			estimates clear

			local btpsummaryfile : subinstr local exportfile ".csv" "_btp.csv", all
 			local btpsummaryfile : subinstr local btpsummaryfile ".rtf" "_btp.rtf", all

			/* For each regression (mat C row) */
			forval reg_idx = 1/`rcoef_row_len' {
				local cellname_li = ""
				forval col_idx = 1/`rcoef_col_len' {
					matrix C_`reg_idx'_`col_idx' = C[`reg_idx', `col_idx']
					matrix colnames C_`reg_idx'_`col_idx' = "reg `reg_idx'"
					quietly: estadd matrix C_`reg_idx'_`col_idx'
					local cellname_li = "`cellname_li' C_`reg_idx'_`col_idx'(fmt(%12.0g))"
				}
				if ("`rcoefidx'" != "") {
					matrix C_init_`reg_idx' = J(1, 1, `rcoefidx')
					matrix colnames C_init_`reg_idx' = "reg `reg_idx'"
					quietly: estadd matrix C_init_`reg_idx'
					local cellname_li = "`cellname_li' C_init_`reg_idx'(fmt(%12.0g))"
				}
				/* If `reg_idx' is 1 & addn contains rcoef: export title */
				if (`reg_idx' == 1 & strpos("`addn'", "(rcoef)")) {
					if ("`rcoefidx'" != "") {
						local btp_c_col_name = `"`c_col_name' "rcoef_idx""'
					}
					else {
						local btp_c_col_name = `"`c_col_name'"'
					}
					if ("`ttitle'"=="") {
						local table_title = "Dep [`: word 1 of `anything''] ~ indep [`: word 1 of `indep''] `addn' `extra_addn'"
					}
					else {
						local table_title = "`ttitle' `extra_addn'"
					}
					esttab using "`btpsummaryfile'", cells("`cellname_li'") ///
					collabels(`btp_c_col_name') mlabels(,none) eqlab(,none) title(`table_title') ///
					append nolines not se compress nogaps noobs plain					
				}
				else {
					esttab using "`btpsummaryfile'", cells("`cellname_li'") ///
					collabels(,none) mlabels(,none) eqlab(,none) ///
					append nolines not se compress nogaps noobs plain
				}
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
	inte(string) clust(namelist) dyn(string) rotctrl(string) dep2(string) ///
	tnote(string) ttitle(string) addn(string) edir(string) sigout(string) sigkw(string) rename(string) ROTINTE REPORT EXPORT NOROUNDDECI]
	
	marksample touse
	
	/* restore arguments for regx */
	local fullopt_args = ""
	foreach opt in ctrl absr xtp inte clust dyn rotctrl tnote ttitle edir report rotinte rename norounddeci {
		if ("``opt''"!="") {
			if ("``opt''"!="report" & "``opt''"!="rotinte" & "``opt''"!="norounddeci") {
				local `opt'_arg = "`opt'(``opt'')"
				local fullopt_args = "`fullopt_args' ``opt'_arg'"
			}
			else {
				local fullopt_args = "`fullopt_args' ``opt''"
			}
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
		local cse = "cluster `clust'"
	}
	else {
		/* (2) If panel is claimed & no cluster specified, cluster SE at id level */
		if ("`xtp'"!="") {
			local cse = "cluster `: word 1 of `xtp''"
		}
		else {
			/* (3) If both panel and cluster are not specified */
			if ("${clustervar}"!="") {
				local cse = "cluster ${clustervar}"
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
			display in red "[ERROR] Need to set a directory in [edir] or global ${dir_table_flow}"
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
			if ("`sigkw'"!="" & "`sigout'"!="") {
				local sigopt_args = "`sigopt_args' sigkw(`sigkw')"
			}
			regx `anything' if `touse', indep(`indep') `fullopt_args' `sigopt_args' keepvar(`lhsvar' `rhsvar') addn("`addn'")
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
				
				/* Temporary solution for suest not being able to process multi-level SE clustering */
				if (`: word count `cse'' >= 3 & strpos("`cse'", "cluster ")) {
					quietly: suest y`depidx'_x`indepidx', vce(robust)
				}
				else {
					quietly: suest y`depidx'_x`indepidx', vce(`cse')
				}
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
		star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`lhsvar']==[`rhsvar'] `addn' `extra_addn'") ///
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
				if ("`sigkw'"!="" & "`sigout'"!="") {
					if (`subvidx'==1) {
						local sigopt_args = "sigout(`sigout') sigkw(`sigkw')"
					}
					else {
						local sigopt_args = "sigout(`sigout'(`subv')) sigkw(`sigkw')"
					}
				}
				else if ("`sigkw'"=="" & "`sigout'"!="") {
					if (`subvidx'==1) {
						local sigopt_args = "sigout(`sigout')"
					}
					else {
						local sigopt_args = "sigout(`sigout'(`subv'))"
					}
				}
				if !missing(`numeric_value') {
					regx `anything' if `subsvar'==`subv' & `touse', indep(`indep') `fullopt_args' `sigopt_args' addn("[`subsvar'==`subv'] `addn'")
				}
				else {
					regx `anything' if `subsvar'=="`subv'" & `touse', indep(`indep') `fullopt_args' `sigopt_args' addn("[`subsvar'==`subv'] `addn'")
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
						
						if (`: word count `cse'' >= 3 & strpos("`cse'", "cluster ")) {
							quietly: suest y`depidx'_x`indepidx'_s`_kw_`prsub1'' y`depidx'_x`indepidx'_s`_kw_`prsub2'', vce(robust)
						}
						else {
							quietly: suest y`depidx'_x`indepidx'_s`_kw_`prsub1'' y`depidx'_x`indepidx'_s`_kw_`prsub2'', vce(`cse')
						}
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
				star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqs_var_x`indeplen''], `subsvar' [=`prsub1'] vs [=`prsub2'] `addn' `extra_addn'") ///
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
						
						if (`: word count `cse'' >= 3 & strpos("`cse'", "cluster ")) {
							quietly: suest y`depidx'_x`indepidx'_s`_kw_`prsub1'' y`depidx'_x`indepidx'_s`_kw_`prsub2'', vce(robust)
						}
						else {
							quietly: suest y`depidx'_x`indepidx'_s`_kw_`prsub1'' y`depidx'_x`indepidx'_s`_kw_`prsub2'', vce(`cse')
						}
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
				star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqtvar'], `subsvar' [=`prsub1'] vs [=`prsub2'] `addn' `extra_addn'") ///
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
			if ("`sigkw'"!="" & "`sigout'"!="") {
				local sigopt_args_a = "sigout(`sigout') sigkw(`sigkw')"
				local sigopt_args_b = "sigout(`sigout'(pair)) sigkw(`sigkw')"
			}
			else if ("`sigkw'"=="" & "`sigout'"!="") {
				local sigopt_args_a = "sigout(`sigout')"
				local sigopt_args_b = "sigout(`sigout'(pair))"
			}
			if ("`eqtvar'"=="~") {
				regx `anything' if `touse', indep(`indep') `fullopt_args' `sigopt_args_a' addn("`addn'")
				regx `dep2' if `touse', indep(`indep') `fullopt_args' `sigopt_args_b' addn("`addn'")
			}
			else {
				regx `anything' if `touse', indep(`indep') `fullopt_args' `sigopt_args_a' keepvar(`eqtvar') addn("`addn'")
				regx `dep2' if `touse', indep(`indep') `fullopt_args' `sigopt_args_b' keepvar(`eqtvar') addn("`addn'")
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

					if (`: word count `cse'' >= 3 & strpos("`cse'", "cluster ")) {
						quietly: suest y`depidx'_x`indepidx'_s1 y`depidx'_x`indepidx'_s2, vce(robust)
					}
					else {
                    	quietly: suest y`depidx'_x`indepidx'_s1 y`depidx'_x`indepidx'_s2, vce(`cse')
					}
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
            star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqs_var_x`indeplen''], Dep Var [=`: word 1 of `anything''] vs [=`: word 1 of `dep2''] `addn' `extra_addn'") ///
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

					if (`: word count `cse'' >= 3 & strpos("`cse'", "cluster ")) {
						quietly: suest y`depidx'_x`indepidx'_s1 y`depidx'_x`indepidx'_s2, vce(robust)
					}
					else {
                    	quietly: suest y`depidx'_x`indepidx'_s1 y`depidx'_x`indepidx'_s2, vce(`cse')
					}
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
            star(* 0.1 ** 0.05 *** 0.01) title("Chi-square: [`eqtvar'], Dep Var [=`: word 1 of `anything''] vs [=`: word 1 of `dep2''] `addn' `extra_addn'") ///
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
		display in red "[ERROR] Wrong input"
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

	/* 
	 * FULLPERC: report more percentile values other than min, p25, p50, p75, and max
	 */

	syntax anything [if] [in], deci(numlist) edir(string) [category(string) tgroup(string) tvar(string) ///
	ttitle(string) addn(string) FULLPERC RLABEL TTEST ONLYT PRESENT]
	
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
		if ("`fullperc'" == "") {
			matrix sum_stat = J(`var_len', 8, 0)
			matrix colnames sum_stat = "Obs" "Mean" "Std.Dev." "Min" "p25" "Median" "p75" "Max"
		}
		else {
			matrix sum_stat = J(`var_len', 14, 0)
			matrix colnames sum_stat = "Obs" "Mean" "Std.Dev." "Min" "p1" "p5" "p10" "p25" "Median" "p75" "p90" "p95" "p99" "Max"			
		}
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
		if ("`fullperc'" == "") {
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
		}
		else {
			local i = 1
			foreach x in `anything' {
				quiet: summ `x' if `touse', detail
				matrix sum_stat[`i',1] = r(N)
				matrix sum_stat[`i',2] = round(r(mean), `deci')
				matrix sum_stat[`i',3] = round(r(sd), `deci')
				matrix sum_stat[`i',4] = round(r(min), `deci')
				matrix sum_stat[`i',5] = round(r(p1), `deci')
				matrix sum_stat[`i',6] = round(r(p5), `deci')
				matrix sum_stat[`i',7] = round(r(p10), `deci')
				matrix sum_stat[`i',8] = round(r(p25), `deci')
				matrix sum_stat[`i',9] = round(r(p50), `deci')
				matrix sum_stat[`i',10] = round(r(p75), `deci')
				matrix sum_stat[`i',11] = round(r(p90), `deci')
				matrix sum_stat[`i',12] = round(r(p95), `deci')
				matrix sum_stat[`i',13] = round(r(p99), `deci')
				matrix sum_stat[`i',14] = round(r(max), `deci')
				local i = `i'+ 1
			}			
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


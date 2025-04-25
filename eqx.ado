///=============================================================================
/*
 * Project: Programs for corporate finance empirical studies
 * Author: Hao Zhao
 * Version
	- eqx: 2.2.5 (16jan2025)
 */
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
	tnote(string) ttitle(string) addn(string) edir(string) sigout(string) sigkw(string) rename(string) ROTINTE REPORT EXPORT]
	
	marksample touse
	
	/* restore arguments for regx */
	local fullopt_args = ""
	foreach opt in ctrl absr xtp inte clust dyn rotctrl tnote ttitle edir report rotinte rename {
		if ("``opt''"!="") {
			if ("``opt''"!="report" & "``opt''"!="rotinte") {
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

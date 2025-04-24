{smcl}
{* *! version 1.7.6 24apr2025}{...}
{* Date: 24apr2025}
{title:Help for {cmd:regx}}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:regx} {hline 2}}Regressions to output tables, supported by {help regress:regress}, {help xtreg:xtreg}, {help reghdfe:reghdfe}, {help tobit:tobit}, and {help esttab:esttab}{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:regx}
{depvars} {ifin}
{cmd:,} {opth indep(indepvars)}
[{opth ctrl(controlvars)} {opth absr(absorbvars)} {opth rotctrl(controlvar_list)} {opth xtp(i t [_r_effect])}
{opth inte(interactlist)} {opth clust(clustervar)} {opth dyn(event [, y=, b=, p=, c=, s=])} {opth rotctrl(controlvar_list)} {opth tobit(ll=, ul=)}
{opth tnote(string)} {opth ttitle(string)} {opth addn(string)} {opth edir(directory)} {opth keepvar(vars)}
{opth stosuf(string)} {opth sigout(string)} {opth sigkw(string)} {opth rcoefidx(string)} 
{opth rename(string)} {cmd:RCOEF} {cmd:SIGMAT} {cmd:REPORT} {cmd:DISPLAY} 
{cmd:ROTINTE} {cmd:CHISTORE} {cmd:NOSINGLETON} {cmd:POISSON} {cmd:NOROUNDDECI}
{help regx##options_table:options}]
{p_end}

{marker options_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{syntab:Variables}
{synopt: {opth indep(indepvars)}}Explanatory variables{p_end}
{synopt : }- Size of output matrix: # of {it:depvars} * # of {it:indepvars}{p_end}
{synopt: {opth ctrl(controlvars)}}Control variables{p_end}
{synopt : }- Adding prefix {it:i.} if control variable is categorical{p_end}
{synopt: {opth rotctrl(controlvars)}}The number of rotating control variables should be equivalent to the number of main independent variables (1 to 1 control){p_end}

{syntab:Model settings}
{synopt: {opth absr(absorbvars)}}Categorical variables representing the fixed effects to be absorbed{p_end}
{synopt : }- High-dimensional FEs, supported by {help reghdfe:reghdfe}{p_end}
{synopt: {opth xtp(i t [_r_effect])}}Regression on the id-time panel{p_end}
{synopt : }- Default setting: id FE{p_end}
{synopt : }- Optional setting: {it:_r_effect} for id RE{p_end}
{synopt : }- Note: cannot be used with {cmd:absr}{p_end}
{synopt: {opth tobit(ll=, ul=)}}Tobit regression{p_end}
{synopt : }- ll=lower bound number{p_end}
{synopt : }- ul=upper bound number{p_end}
{synopt: {opth clust(clustervar)}}Variable at which level standard error clustered{p_end}
{synopt : }- Default setting: global ${clustervar} for a project{p_end}
{synopt : }- Priority: {cmd:clust} > {cmd:xtp} > ${clustervar} > {it:robust}{p_end}
{synopt: {cmd:NOSINGLETON}}Drop singleton for{cmd:reghdfe}{p_end}
{synopt: {cmd:POISSON}}Using poisson regressions for the estimation{p_end}

{syntab:Interaction analysis}
{synopt: {opth inte(interactlist)}}Interaction term{p_end}
{synopt : }- Pairwise interaction with each variable of {cmd:indep}{p_end}
{synopt : }- One interaction term: {opth inte(D1)}{p_end}
{synopt : }- More than one interaction term: {opth inte(D1 D2 ...)}{p_end}
{synopt : }- Multiple interaction: {opth inte(D1#Z1#C1#...)}{p_end}
{synopt : }- More than one multiple interaction: {opth inte(D1#Z1#C1#... D2#Z2#C2#... ...)}{p_end}
{synopt: {opth dyn(event [, y=, b=, p=, c=, s=])}}Dynamic interaction{p_end}
{synopt : }- Pairwise interaction with each variable of {cmd:indep}{p_end}
{synopt : }- Default dynamic time indicator: {it:event_tm1}, {it:event_t0}, {it:event_t1}, etc.{p_end}
{synopt : }- y=: time range {it: start/end}. Default range: {it:-2/2}{p_end}
{synopt : }- b=: benchmark time. Default benchmark: {it:0}{p_end}
{synopt : }- p=: if not specified, no dynamic plot; if specified, the input will be the title of the plot, with underscore replaced by spacing{p_end}
{synopt : }- c=: confidence interval. Default: 95{p_end}
{synopt : }- s=: (1) no input: no directory specified, the default output directory is {it:parent dir of ${dir_table_flow} / indepvaridx + depvar name}{p_end}
{synopt : }- s={it:foldername}: (2) no keyword specified: default output directory is {it:parent dir of ${dir_table_flow} / foldername / indepvaridx + depvar name}{p_end}
{synopt : }- s={it:foldername ~ keyword} (3) specified directory and keyword: default output directory is {it:parent dir of ${dir_table_flow} / foldername / indepvaridx + depvar name + keyword}{p_end}
{synopt: {cmd:ROTINTE}}Rotating interaction controls; using together with {opth rotctrl(controlvar_list)} and {opth inte(controlvar_list)}{p_end}

{syntab:Export settings}
{synopt: {opth tnote(string)}}Notes under the table, e.g., "Controls/FEs: Firm[Y] Industry[Y] Province[Y] Year[Y]"{p_end}
{synopt: {opth ttitle(string)}}Table title{p_end}
{synopt : }- If not specified, the default title is {it:Dep [`: word 1 of `anything''] ~ indep [`: word 1 of `indep''] `addn'}{p_end}
{synopt: {opth addn(string)}}Additional notes following the default title{p_end}
{synopt: {opth edir(directory)}}Export directory{p_end}
{synopt : }- If not specified, setting global ${dir_table_flow} as the default export directory{p_end}
{synopt: {opth keepvar(vars)}}Extra variables (e.g., control variables) to export{p_end}
{synopt: {opth stosuf(string)}}The suffix of estimate stored results, for {cmd:eqx}{p_end}
{synopt: {opth sigout(string)}}Export sign and significance matrix, with row name as {it: Model + `sigout']}{p_end}
{synopt : }- Default table title: "Dep [`: word 1 of `anything''] ~ indep [`: word 1 of `indep''] `addn'"{p_end}
{synopt : }- If {cmd:ttitle} is specified, the title will be {cmd:ttitle}{p_end}
{synopt : }- The default name of export file is {it:export directory / export file name + "_ss" + file suffix}{p_end}
{synopt : }- Can show extra variables if specifying {cmd:keepvar} and {cmd:rotctrl} {p_end}
{synopt: {opth sigkw(string)}}Keyword for table formatting{p_end}
{synopt : }- If "`sigout'"=="1" | "`sigout'"=="baseline" | "`sigout'"=="`sigkw'", starting to append results on a new table{p_end}
{synopt: {opth rcoefidx(string)}}Input number index, and export index as a column for {cmd:RCOEF} output{p_end}
{synopt: {opth rename(string)}}Consistently rename all the independent variables to the same output name{p_end}
{synopt: {cmd:RCOEF}}Export b, t, and p for each set of regressions. If title and column names are needed for export tables, {opth addn(string)} must include "(rcoef)"{p_end}
{synopt : }- {cmd:RCOEF} cannot be used with {opth sigout(string)} or {cmd:SIGMAT}{p_end}
{synopt: {cmd:SIGMAT}}Display and return sign and significance matrix{p_end}
{synopt: {cmd:REPORT}}Reporting the coefficients for all variables{p_end}
{synopt: {cmd:DISPLAY}}Display regression results only, without exporting tables{p_end}
{synopt: {cmd:CHISTORE}}Storing estimated results for {cmd:eqx}{p_end}
{synopt: {cmd:NOROUNDDECI}}If specified, output tables will keep 3-digit non-zero numbers, instead of rounding to 3-digit after the decimal point{p_end}

{marker contact}{...}
{title:Author}

{pstd}Hao Zhao{break}
Durham University{break}
Email: {browse "mailto:hao.zhao@durham.ac.uk":hao.zhao@durham.ac.uk}
{p_end}



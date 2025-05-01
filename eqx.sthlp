{smcl}
{* *! version 2.2.5 16jan2025}{...}
{* Date: 16jan2025}
{title:Help for {cmd:eqx}}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:eqx} {hline 2}}Regressions to output tables, supported by {help suest:suest}, {help regress:regress}, and {help esttab:esttab}{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:eqx}
{depvars} {ifin}
{cmd:,} {opth indep(indepvars)} {opth eqt([var1==var2]/[~/var|subvar]/[~/var@])}
[{opth ctrl(controlvars)} {opth absr(absorbvars)} {opth rotctrl(controlvar_list)} {opth xtp(i t [_r_effect])}
{opth inte(interactlist)} {opth clust(clustervar)} {opth dyn(event [, y=, b=, p=, c=, s=])} {opth dep2(depvars2)}
{opth tnote(string)} {opth ttitle(string)} {opth addn(string)} {opth edir(directory)} {opth sigout(string)} {opth sigkw(string)} 
{opth rename(string)} {cmd:ROTINTE} {cmd:REPORT} {cmd:EXPORT} {cmd:NOROUNDDECI}
{help eqx##options_table:options}]
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

{syntab:Chi-square test}
{synopt: {opth eqt([var1==var2]/[~/var|subvar]/[~/var@])}} {p_end}
{synopt : }- [var1==var2]: within model chi-square test between coefficients, e.g., eqx y, ... eqt(x1==x2){p_end}
{synopt : }- [~/var|subvar]: between category coefficient chi-square test (models separated by a categorical variable), e.g., (1) eqx y, ... eqt(x | dummy), (2) eqx y, ... eqt(~ | dummy){p_end}
{synopt : }- [~/var@]: between model coefficient chi-square test (models separated by two sets of dependent variables; {cmd:dep2} must be specified): (1) eqx y1, ... eqt(x @) dep2(y2), (2) eqx y1, ... eqt(~ @) dep2(y2){p_end}
{synopt: {opth dep2(depvars2)}}The second dependent variable set{p_end}
{synopt : }- For between model coefficient chi-square test{p_end}
{synopt : }- Condition: # of depvars2 == # of depvars{p_end}

{syntab:Model settings}
{synopt: {opth absr(absorbvars)}}Categorical variables representing the fixed effects to be absorbed{p_end}
{synopt : }- High-dimensional FEs, supported by {help reghdfe:reghdfe}{p_end}
{synopt: {opth xtp(i t [_r_effect])}}Regression on the id-time panel{p_end}
{synopt : }- Default setting: id FE{p_end}
{synopt : }- Optional setting: {it:_r_effect} for id RE{p_end}
{synopt : }- Note: cannot be used with {cmd:absr}{p_end}
{synopt: {opth clust(clustervar)}}Variable at which level standard error clustered{p_end}
{synopt : }- Default setting: global ${clustervar} for a project{p_end}
{synopt : }- Priority: {cmd:clust} > {cmd:xtp} > ${clustervar} > {it:robust}{p_end}

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
{synopt: {opth rename(string)}}Consistently rename all the independent variables to the same output name{p_end}
{synopt : }- If not specified, setting global ${dir_table_flow} as the default export directory{p_end}
{synopt: {cmd:REPORT}}Reporting the coefficients for all variables{p_end}
{synopt: {cmd:EXPORT}}Exporting the regression tables at the same time using {cmd:regx}{p_end}
{synopt: {cmd:NOROUNDDECI}}If specified, output tables will keep 3-digit non-zero numbers, instead of rounding to 3-digit after the decimal point{p_end}

{marker contact}{...}
{title:Author}

{pstd}Hao Zhao{break}
Durham University{break}
Email: {browse "mailto:hao.zhao@durham.ac.uk":hao.zhao@durham.ac.uk}
{p_end}
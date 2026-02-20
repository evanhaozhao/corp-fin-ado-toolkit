{smcl}
{* *! version 1.3.7 20feb2026}{...}
{* Date: 19apr2025}
{title:Help for {cmd:sumx}}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:sumx} {hline 2}}Summary and t-tests to output tables{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:sumx}
{vars} {ifin}
{cmd:,} {opth deci(numlist)} {opth edir(directory)}
[{opth category(var)} {opth tgroup(var)} {opth tvar(var)}
{opth ttitle(string)} {opth addn(string)} {cmd:RLABEL} {cmd:TTEST}
{cmd:ONLYT} {cmd:PRESENT}
{help sumx##options_table:options}]
{p_end}

{marker options_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{synopt: {opth deci(numlist)}}The number of digits after the decimal point indicates the level of precision. For example, 0.001 has three digits after the decimal point{p_end}

{synopt:{opth category(string)}}A categorical variable used to group the data. Summary statistics will be reported for each group defined by this variable.{p_end}

{synopt:{opth tgroup(string)}}A binary (two-group) variable used for two-sample t-tests.{p_end}

{synopt:{opth tvar(string)}}A second variable used in the two-variable comparison for t-tests.{p_end}

{synopt:{opth edir(string)}}Export directory to save the output files.{p_end}

{synopt: {opth ttitle(string)}}Table title{p_end}

{synopt: {opth addn(string)}}Additional notes following the default title{p_end}

{synopt:{opt RLABEL}}Export results using variable labels instead of variable names.{p_end}

{synopt:{opt TTEST}}Conduct a one-sample t-test testing whether the variable's mean is different from zero.{p_end}

{synopt:{opt ONLYT}}Restrict output to the two-sample t-test results only; summary statistics for the two groups will not be reported.{p_end}

{synopt:{opt PRESENT}}Format the two-sample t-test output for presentation. Group 1 and Group 2 will each display: mean, median, and number of observations. The differences in means and medians, along with p-values and total observations, will also be shown.{p_end}

{marker contact}{...}
{title:Author}

{pstd}Hao Zhao{break}
Durham University{break}
Email: {browse "mailto:hao.zhao@durham.ac.uk":hao.zhao@durham.ac.uk}
{p_end}
/*
	 Subgroup Analysis
	1) Transport subgroup analysis from iPrEx to synthetic pop of Latino Young
	MSM Study sites in Chicago and SF
	2) Estimate NNTs for TGW vs MSM in each target pop
*/

clear
set more off
use "iprex_fortp.dta" //iprex data trimmed and formatted to match target variables
merge 1:1 ptid using "iprex futime for stset.dta", keep(master match) nogen
drop if years==0 //only include those people who contributed followup time

/*
	estimate risk difference and NNT at one year
	- include offset for follow-up time, as we have some LTFU
	- data are limited to 1 year: only infections that had occurred by 1 yr are counted
	& max follow-up time is set to 1 year
*/

/**********************************************************
	iPrEx: Observed in iPrEx
**********************************************************/

/*
	Overall
*/
	quietly glm iprexsc i.tx, fam(bin) link(log) exposure(years) robust
	margins i.tx, post
	nlcom (_b[1.tx]-_b[0.tx])
	nlcom (1/(_b[1.tx]-_b[0.tx]))

/*
	Gender
*/
	quietly glm iprexsc i.tx##i.trans, fam(bin) link(log) exposure(years) robust
	margins i.tx#i.trans, post
	nlcom _b[1.tx#1.trans]-_b[0.tx#1.trans]
	nlcom _b[1.tx#0.trans]-_b[0.tx#0.trans]

	nlcom (1/(_b[1.tx#0.trans]-_b[0.tx#0.trans]))
	nlcom (1/(_b[1.tx#1.trans]-_b[0.tx#1.trans]))

/*
	NCRAI
*/
	quietly glm iprexsc i.tx##i.ncrai, fam(bin) link(log) exposure(years) robust
	margins i.tx##i.ncrai, post
	nlcom _b[1.tx#1.ncrai]-_b[0.tx#1.ncrai]
	nlcom _b[1.tx#0.ncrai]-_b[0.tx#0.ncrai]
	nlcom (1/(_b[1.tx#0.ncrai]-_b[0.tx#0.ncrai]))
	nlcom (1/(_b[1.tx#1.ncrai]-_b[0.tx#1.ncrai]))


/*
	Sex role
*/

quietly glm iprexsc i.tx##i.role, fam(bin) link(log) exposure(years) robust
margins i.tx##i.role, post
nlcom (_b[1.tx#1.role]-_b[0.tx#1.role]) //top
nlcom (_b[1.tx#2.role]-_b[0.tx#2.role]) //bottom
nlcom (_b[1.tx#3.role]-_b[0.tx#3.role]) //versatile

nlcom (1/(_b[1.tx#1.role]-_b[0.tx#1.role])) //top
nlcom (1/(_b[1.tx#2.role]-_b[0.tx#2.role])) //bottom
nlcom (1/(_b[1.tx#3.role]-_b[0.tx#3.role])) //versatile

/*
	Cocaine Use
*/

quietly glm iprexsc i.tx##i.coke, fam(bin) link(log) exposure(years) robust
margins i.tx##i.coke, post
nlcom (_b[1.tx#1.coke]-_b[0.tx#1.coke])
nlcom (1/(_b[1.tx#0.coke]-_b[0.tx#0.coke]))
nlcom (1/(_b[1.tx#1.coke]-_b[0.tx#1.coke]))


/***************************************
	iPrEx to SF
****************************************/
clear
use "iprex_fortp.dta" 
append using "LMCI_combined.dta"
lab def sf 1"SF" 0"Chicago"
lab val sf sf
keep if sf==1 | iprex==1
merge m:1 ptid using "iprex futime for stset.dta", keep(master match) nogen
drop if years==0 & iprex==1

/*
	Overall
*/
*Generate weights*
capture drop ps num wt
quietly logit iprex c.partners##i.ncrai i.alcohol i.age i.edu i.coke 
predict ps, p
quietly logit iprex
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0

quietly glm iprexsc i.tx [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx, post 
nlcom (_b[1.tx]-_b[0.tx])
nlcom (1/(_b[1.tx]-_b[0.tx]))



/*
	TRANS
*/
*Generate weights*
capture drop ps num wt
quietly logit iprex c.partners##i.role i.ncrai i.alcohol i.trans i.coke i.age i.edu
predict ps, p
quietly logit iprex trans
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0

quietly glm iprexsc i.tx##i.trans [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.trans, post
nlcom _b[1.tx#1.trans]-_b[0.tx#1.trans] //tgw
nlcom _b[1.tx#0.trans]-_b[0.tx#0.trans] //msm

nlcom (1/(_b[1.tx#0.trans]-_b[0.tx#0.trans])) //msm
nlcom (1/(_b[1.tx#1.trans]-_b[0.tx#1.trans])) //tgw




/*
	NCRAI
*/
drop ps num wt
*Generate weights*
quietly logit iprex c.partners i.alcohol i.age i.edu i.coke##i.ncrai
predict ps, p
quietly logit iprex ncrai
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0


quietly glm iprexsc i.tx##i.ncrai [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.ncrai, post
nlcom _b[1.tx#1.ncrai]-_b[0.tx#1.ncrai] //ncrai

nlcom (1/(_b[1.tx#0.ncrai]-_b[0.tx#0.ncrai]))
nlcom (1/(_b[1.tx#1.ncrai]-_b[0.tx#1.ncrai])) //ncrai


/*
	Sex role
*/
capture drop ps num wt
*Generate weights*
quietly logit iprex c.partners i.role##i.coke i.ncrai i.alcohol c.age i.edu
predict ps, p
quietly logit iprex role
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0

quietly glm iprexsc i.tx##i.role [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.role, post
nlcom (_b[1.tx#1.role]-_b[0.tx#1.role]) //top
nlcom (_b[1.tx#2.role]-_b[0.tx#2.role]) //bottom
nlcom (_b[1.tx#3.role]-_b[0.tx#3.role]) //versatile

nlcom (1/(_b[1.tx#1.role]-_b[0.tx#1.role])) //top
nlcom (1/(_b[1.tx#2.role]-_b[0.tx#2.role])) //bottom
nlcom (1/(_b[1.tx#3.role]-_b[0.tx#3.role])) //versatile

/*
	COKE
*/
capture drop ps num wt
*Generate weights*
quietly logit iprex c.partners##i.ncrai i.alcohol i.age i.edu i.coke
predict ps, p
quietly logit iprex coke
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0


quietly glm iprexsc i.tx##i.coke [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.coke, post
nlcom _b[1.tx#1.coke]-_b[0.tx#1.coke] //coke
nlcom (1/(_b[1.tx#0.coke]-_b[0.tx#0.coke]))
nlcom (1/(_b[1.tx#1.coke]-_b[0.tx#1.coke])) //coke


*BOOTSTRAP*

quietly run bsprog.do
foreach var in trans ncrai role coke {
bootstrap rd=r(`var'), reps(2000) seed(1) nodrop saving(bstpres_sf, replace): bstp
bootstrap nnt=r(nnt`var'), reps(2000) seed(1) nodrop saving(bsnnt_sf, replace): bsnnt
}

/****************************
	iPrEx to Chicago
*****************************/

clear
use "iprex_fortp.dta"
append using "LMCI_combined.dta"
lab def sf 1"SF" 0"Chicago"
lab val sf sf
keep if sf==0 | iprex==1

merge m:1 ptid using "iprex futime for stset.dta", keep(master match) nogen
drop if years==0 & iprex==1

/*
	Overall
*/
*Generate weights*
capture drop ps num wt
quietly logit iprex c.partners##i.ncrai i.alcohol i.age i.edu i.coke 
predict ps, p
quietly logit iprex
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0

quietly glm iprexsc i.tx [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx, post
nlcom (_b[1.tx]-_b[0.tx])
nlcom (1/(_b[1.tx]-_b[0.tx]))



/*
	TRANS
*/
*Generate weights*
capture drop ps num wt
quietly logit iprex c.partners##i.role i.ncrai i.alcohol i.trans i.coke i.age i.edu
predict ps, p
quietly logit iprex trans
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0

quietly glm iprexsc i.tx##i.trans [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.trans, post
nlcom _b[1.tx#1.trans]-_b[0.tx#1.trans] //tgw
nlcom _b[1.tx#0.trans]-_b[0.tx#0.trans] //msm

nlcom (1/(_b[1.tx#0.trans]-_b[0.tx#0.trans])) //msm
nlcom (1/(_b[1.tx#1.trans]-_b[0.tx#1.trans])) //tgw




/*
	NCRAI
*/
drop ps num wt
*Generate weights*
quietly logit iprex c.partners i.alcohol i.age i.edu i.coke##i.ncrai
predict ps, p
quietly logit iprex ncrai
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0


quietly glm iprexsc i.tx##i.ncrai [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.ncrai, post
nlcom _b[1.tx#1.ncrai]-_b[0.tx#1.ncrai] //ncrai

nlcom (1/(_b[1.tx#0.ncrai]-_b[0.tx#0.ncrai]))
nlcom (1/(_b[1.tx#1.ncrai]-_b[0.tx#1.ncrai])) //ncrai


/*
	Sex role
*/
capture drop ps num wt
*Generate weights*
quietly logit iprex c.partners i.role##i.coke i.ncrai i.alcohol c.age i.edu
predict ps, p
quietly logit iprex role
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0

quietly glm iprexsc i.tx##i.role [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.role, post
nlcom (_b[1.tx#1.role]-_b[0.tx#1.role]) //top
nlcom (_b[1.tx#2.role]-_b[0.tx#2.role]) //bottom
nlcom (_b[1.tx#3.role]-_b[0.tx#3.role]) //versatile

nlcom (1/(_b[1.tx#1.role]-_b[0.tx#1.role])) //top
nlcom (1/(_b[1.tx#2.role]-_b[0.tx#2.role])) //bottom
nlcom (1/(_b[1.tx#3.role]-_b[0.tx#3.role])) //versatile

/*
	COKE
*/
capture drop ps num wt
*Generate weights*
quietly logit iprex c.partners##i.ncrai i.alcohol i.age i.edu i.coke
predict ps, p
quietly logit iprex coke
predict num, p

gen wt=((1-ps)/(ps))*(num/(1-num))
replace wt=0 if iprex==0


quietly glm iprexsc i.tx##i.coke [pweight=wt], fam(bin) link(log) robust exposure(years)
margins i.tx##i.coke, post
nlcom _b[1.tx#1.coke]-_b[0.tx#1.coke] //coke
nlcom (1/(_b[1.tx#0.coke]-_b[0.tx#0.coke]))
nlcom (1/(_b[1.tx#1.coke]-_b[0.tx#1.coke])) //coke


*BOOTSTRAP*

quietly run bsprog.do
foreach var in trans ncrai role coke {
bootstrap rd=r(`var'), reps(2000) seed(1) nodrop saving(bstpres_chic, replace): bstp
bootstrap nnt=r(nnt`var'), reps(2000) seed(1) nodrop saving(bsnnt_chic, replace): bsnnt
}

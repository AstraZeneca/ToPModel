**Simulated data**

Uncompress simdata.zip to obtain simdata.csv. This provides simulated participant level data to fit a ToP model similar to the manuscript.

Variables in simdata.csv:

- id = participant ID

- z = 1 if received treatment, else 0
- start, stop = daily at-risk intervals relative to time since dose

- ev = 1 if participant had an event between [start, stop] days after dose, else 0

- conc = daily serum mAb concentration (μg/mL) predicted from pop-PK

- wtg = prevalence-adjusted IC50 (ng/mL)

- wtg_nab = prevalence-adjusted nAb titre, calculated as (conc*1000)/wtg

- wtg_lnab = log10(wtg_nab+1)

**ToPModel.sas**

SAS code to fit the ToP model to simdata.csv and generate efficacy point estimates and bootstrapped two-sided 90% confidence intervals. Developed using SAS 9.4. Set macro variables on lines 1-2 as described. In the below example NBOOT is set to 10.

	%let DPATH = <SET TO LOCAL FILE PATH WHERE SIMDATA.CSV IS SAVED.>;
 	%let NBOOT = <SET TO NUMBER OF BOOTSTRAP REPLICATES. CAN BE SET TO SMALLER NUMBER (E.G. 10) FOR QUICK RUN TIME.>;
	goptions device=png;
	ods graphics on / reset=all noborder width=10in height=6in;
	proc datasets kill nolist mt=data lib=work;
	run;
	quit;

	*Read in simulated TTE data;
	proc import datafile="&DPATH.simdata.csv" out=simdata dbms=csv replace;
	proc sort;
		by z id start;
	run;

	*Fit ToP Cox Model on full data to get efficacy point estimates;
	ods exclude hazardratios;
	proc phreg data = simdata covout outest=mat;
		model (start, stop)*ev(0) = z z*wtg_lnab / rl ties=efron;
		hazardratio z / at(wtg_lnab = 0.001 to %sysevalf(%sysfunc(ceil(%sysfunc(log10(1e5+1))*1000))/1000) by 0.001);
		ods output hazardratios = hrfull;
	run;

![image](https://github.com/user-attachments/assets/e7bc9c75-86dc-4634-a973-14f409eaa29b)

	data efffull;
		set hrfull;
		*Derive efficacy as 1 - HR;
		eff = 100*(1 - hazardratio);
		*Get covariate value HR assessed at;
		wtg_lnab = input(scan(scan(description,3,'='),1,' '),best.);
		*Backmap to nAbs;
		nab = (10**wtg_lnab)-1;
		keep eff: nab;
	proc sort;
		by nab;
	run;

	*Bootstrapped confidence interval;
	proc sql noprint;
		select count(distinct id) into :n trimmed from simdata;
	quit;
	ods select none;
	*Resample with replacementwith each size equivalence to the total sample size from the original data;
	proc surveyselect data=simdata out=bootsamp method=urs reps=&NBOOT seed=7442 sampsize=&n outhits;
		cluster z id;
	run;
	*Fit model and derive efficacy (similar to above) for each replicate;
	proc phreg data = bootsamp;
		by replicate;
		model (start, stop)*ev(0) = z z*wtg_lnab / rl ties=efron;
		hazardratio z / at(wtg_lnab = 0.001 to %sysevalf(%sysfunc(ceil(%sysfunc(log10(1e5+1))*1000))/1000) by 0.001) alpha=0.1;
		ods output hazardratios = hrboot;
	run;
	data effboot;
		set hrboot;
		eff = 100*(1 - hazardratio);
		wtg_lnab = input(scan(scan(description,3,'='),1,' '),best.);
		nab = (10**wtg_lnab)-1;
		keep replicate eff: nab;
	proc sort;
		by nab;
	run;
	*Calculate the two-side 90 percent CI using the quantile approach;
	proc univariate data=effboot;
		by nab;
		var eff;
		output out=booteff(rename=(eff5 = eff_l eff95=eff_u)) pctlpre=eff pctlpts=5, 95;
	run;
	ods select all;

	*Combine point estimates from full data with bootstrapped CIs;
	data ToPCox;
		merge efffull(in=a) booteff;
		by nab;
		if a;
	run;

	*Plot threshold curve;
	proc sgplot data=ToPCox noautolegend noborder nowall;
		series y=eff x=nab / lineattrs=(color=cx830051);
		band lower=eff_l upper=eff_u x=nab / transparency=.8 fillattrs=(color=cx830051);
		xaxis type=log logbase=10 values=(1e1 1e2 1e3 1e4 1e5) label='Prevalence-adjusted nAb titre';
		yaxis values=(0 to 100 by 10) label='Efficacy (%)';
	run;

 ![image](https://github.com/user-attachments/assets/d6ec5dab-2e6b-4435-9e11-95a5cab65f96)



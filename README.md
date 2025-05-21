simdata.zip

	Uncompress zip to obtain simdata.csv. Provides simulated participant level data to fit a ToP model similar to the manuscript.

Variables in simdata.csv

		id = participant ID
		z = 1 if received treatment, else 0
		start, stop = daily at-risk intervals
		ev = 1 if participant had an event between [start, stop] days after dose, else 0
		conc = daily serum mAb concentration (μg/mL) predicted from pop-PK
		wtg = prevalence-adjusted IC50 (ng/mL)
		wtg_nab = prevalence-adjusted nAb titre, calculated as (conc*1000)/wtg
		wtg_lnab = log10(wtg_nab+1)

ToPModel.sas

	SAS code to fit the ToP model to simdata.csv and generate efficacy point estimates and
	bootstrapped two-sided 90% confidence intervals. Developed using SAS 9.4. Set macro variables on lines
	1-2 as described.

	

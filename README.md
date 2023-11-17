# u-pro-crm-simulations
This repository contains the code used to simulate results for the paper **U-PRO-CRM: Designing patient-centred dose-finding trials using Patient-reported Outcomes**. 

## Background 
This project evaluates the operating characteristics of the novel U-PRO-CRM design under simulation scenarios and provides the code to evaluate the design against the publushed PRO-CRM design and benchmark. 

## Description of R files
* **functions.R** - code for Section 4.2; functions required to generate trials and identify the MTD for 5,000 simulations in parallel. 
  
* **u_pro_crm_run.R** - code for Section 4.2: code required to run U-PRO-CRM and PRO-CRM simulation studies. Functions required to run this code are defined in `functions.R`.

* **benchmark_run.R** - code for Section 4.3: functions and code to evaluate a benchmark[1] for the U-PRO-CRM design

[1] Cheung YK. Simple benchmark for complex dose finding studies. Biometrics. 2014 Jun;70(2):389-97. doi: 10.1111/biom.12158. Epub 2014 Feb 25. PMID: 24571185; PMCID: PMC4061271.

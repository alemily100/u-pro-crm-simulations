# u-pro-crm-simulations
This repository contains the code used to simulate results for the paper **U-PRO-CRM: Designing patient-centred dose-finding trials using Patient-reported Outcomes**. 

## Background 
This project evaluates the operating characteristics of the novel U-PRO-CRM design under simulation scenarios where it is compared to the published PRO-CRM[1] design and a benchmark[2]. 

## Description of R files
* **functions.R** - code for Section 4.2: functions required to generate trials and identify the MTD for 5,000 simulations in parallel. 
  
* **u_pro_crm_run.R** - code for Section 4.2: code required to run U-PRO-CRM and PRO-CRM simulation studies. Functions required to run this code are defined in `functions.R`.

* **benchmark_run.R** - code for Section 4.3: functions and code to evaluate a benchmark[2] for the U-PRO-CRM design

[1] Lee SM, Lu X, Cheng B. Incorporating patient-reported outcomes in dose-finding clinical trials. Stat Med. 2020 Feb 10;39(3):310-325. doi: 10.1002/sim.8402. Epub 2019 Dec 3. PMID: 31797421; PMCID: PMC8411935.

[2] Cheung YK. Simple benchmark for complex dose finding studies. Biometrics. 2014 Jun;70(2):389-97. doi: 10.1111/biom.12158. Epub 2014 Feb 25. PMID: 24571185; PMCID: PMC4061271.

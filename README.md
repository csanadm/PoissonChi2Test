# Poisson $\chi^2$ test
Testing Poisson uncertainties and $\chi^2$ definitions. In the below formulas,
- $\chi^2_i$ is the chisquare for one bin,
- $O_i$ is the observed count in the bin,
- its variance is also $O_i$ in a Poisson case,
- $E_i$ is the expected (theoretical) value, coming from the fit.

## Default

$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}}$

## Yates' Correction

$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{O_{i}}$ 
(where the numerator is zero if $|O_{i}-E_{i}|<0.5$)

## Pearson's Chisquare

$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{E_{i}}$

## Yates & Pearson

$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{E_{i}}$
(where the numerator is zero if $|O_{i}-E_{i}|<0.5$)

## Nagy, Csanád, et al.

$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}+0.5}$

## Plots

### Summary table 
Values for a 400k bin histogram, filled with uniformly distributed random values, ordered by decreasing $\chi^2$

|	Case	|	Total Hits	|	Expected Value	|	Fitted Value	|	$\chi^2$	|	Prob (CL)	|
|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|
|	Default	|	40M	|	100	|	98.99±0.02	|	403752.41	|	0.001%	|
|		|	80M	|	200	|	198.99±0.02	|	403065.34	|	0.03%	|
|		|	160M	|	400	|	399.00±0.03	|	401828.24	|	2.0%	|
|	Nagy, Csanád, et al.	|	40M	|	100	|	99.00±0.02	|	401681.68	|	3.0%	|
|		|	80M	|	200	|	198.99±0.02	|	402044.43	|	1.1%	|
|		|	160M	|	400	|	399.00±0.03	|	401322.71	|	7.0%	|
|	Pearson	|	40M	|	100	|	100.50±0.02	|	398785.86	|	91.3%	|
|		|	80M	|	200	|	200.50±0.02	|	400412.49	|	32.2%	|
|		|	160M	|	400	|	400.50±0.03	|	400526.10	|	27.8%	|
|	Yates & Pearson |	40M	|	100	|	99.00±0.02	|	372521.66	|	100%	|
|		|	80M	|	200	|	199.00±0.02	|	380866.03	|	100%	|
|		|	160M	|	400	|	399.00±0.03	|	386072.57	|	100%	|
|	Yates	|	40M	|	100	|	100.49±0.02	|	367947.61	|	100%	|
|		|	80M	|	200	|	200.49±0.02	|	378357.42	|	100%	|
|		|	160M	|	400	|	400.49±0.03	|	384823.02	|	100%	|

### Default

<img height="250" alt="poisson_vs_gauss_test_40M" src="https://github.com/user-attachments/assets/40289048-bb72-42ce-aaed-72fedce94eb9" />
<img height="250" alt="poisson_vs_gauss_test_80M" src="https://github.com/user-attachments/assets/5393ecf4-4f7d-4388-b0c1-b66b160c0eaf" />
<img height="250" alt="poisson_vs_gauss_test_160M" src="https://github.com/user-attachments/assets/d8aa7ecb-5a9a-4d7a-bfc9-da2cd539efc4" />

### Yates' Correction

<img height="250" alt="poisson_vs_gauss_test_40M_Yates" src="https://github.com/user-attachments/assets/1c95f79f-1135-4d23-93be-5d3785dd6d5e" />
<img height="250" alt="poisson_vs_gauss_test_80M_Yates" src="https://github.com/user-attachments/assets/f039fd59-e48e-4b17-b953-346764027210" />
<img height="250" alt="poisson_vs_gauss_test_160M_Yates" src="https://github.com/user-attachments/assets/ce4e11a2-d848-4097-88a5-5aeee89b352f" />

### Pearson's Chisquare

<img height="250" alt="poisson_vs_gauss_test_40M_Pearson" src="https://github.com/user-attachments/assets/1a1b6b22-96d8-4fbb-b13f-9b1f82b89632" />
<img height="250" alt="poisson_vs_gauss_test_80M_Pearson" src="https://github.com/user-attachments/assets/e0ca88e5-a820-4b22-84b0-bafe263cef9b" />
<img height="250" alt="poisson_vs_gauss_test_160M_Pearson" src="https://github.com/user-attachments/assets/1772a22f-38a7-48d4-a5e2-2d993b8e6fe7" />

### Yates & Pearson

<img height="250" alt="poisson_vs_gauss_test_40M_YatesMod" src="https://github.com/user-attachments/assets/9848681f-b7af-4555-9f25-8ff317b844f1" />
<img height="250" alt="poisson_vs_gauss_test_80M_YatesMod" src="https://github.com/user-attachments/assets/72452c0a-13fc-434e-b118-cba6a1197751" />
<img height="250" alt="poisson_vs_gauss_test_160M_YatesMod" src="https://github.com/user-attachments/assets/db418faa-bab3-4c01-b1d4-6891fd7fee11" />

### Nagy, Csanád, et al.

<img height="250" alt="poisson_vs_gauss_test_40M_corr" src="https://github.com/user-attachments/assets/d1757ca3-4d63-427f-9cf5-d8d3fdf4921d" />
<img height="250" alt="poisson_vs_gauss_test_80M_corr" src="https://github.com/user-attachments/assets/6542f7c5-65fd-4dc0-b83f-0419822dd2ba" />
<img height="250" alt="poisson_vs_gauss_test_160M_corr" src="https://github.com/user-attachments/assets/6bff0395-5509-45a1-bdbd-7a4560b09729" />




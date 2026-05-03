# Poisson $\chi^2$ test
Testing Poisson uncertainties and $\chi^2$ definitions. In the below formulas,
- $\chi^2_i$ is the chisquare for the $i$-th bin,
- $O_i$ is the observed count in the $i$-th bin,
- $\sigma^2_i = (\Delta O_i)^2 = O_i$ is the variance in the $i$-th (assuming a Poisson distribution)
- $E_i$ is the expected (theoretical) value, coming from the fit,
- and an alternative variance estimator (à la Pearson) is then $E_i$

## Default
The classical $\chi^2$ definition is:<br>
$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}}$<br>
This generally results in an underestimation of the bin contents, due to the so-called Neyman bias.

## Yates' Correction
One can correct for the Neyman bias by subtracting 1/2 from the numerator:<br>
$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{O_{i}}$ <br>
(where the numerator is zero if $|O_{i}-E_{i}|<0.5$)

## Pearson's Chisquare
Alternatively, one can use $E_i$ for the variance, which also removes the Neyman bias:<br>
$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{E_{i}}$

## Yates & Pearson
One can combine the above two:<br>
$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{E_{i}}$<br>
(where the numerator is zero if $|O_{i}-E_{i}|<0.5$)

## Nagy, Csanád, et al.
In addition, one can investigate which Gaussian distribution is closest to a given Poisson distribution with parameter $\lambda$. It turns out, that this is the Gaussian with a variance shifted by 1/2, and one can include this in the chisquare definition:<br>
$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}+0.5}$

## Log likelihood
Finally, to remedy the problems associated with assuming a Gaussian uncertainty distribution, log likelihood optimization can be utilized. In this case, one defines chisquare as:<br>
$\chi^2_i = 2\left[ E_{i} - O_{i} + O_{i} \log(O_{i}/E_{i}) \right]$


## Testing the $\chi^2$ options

In the attached code [**poisson_vs_gauss_test.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test.C) we create a histogram with $N_{\rm bins}=400 000$ bins, and then fill it with $N_{\rm hits}$ uniformly distributed random values. In principle, if $\lambda = $N_{\rm hits}/N_{\rm bins}$, then we expect $\lambda$ hits in every bin. Furthermore, if we fit the resulting histogram with a constant, based on one of the $\chi^2$ definitions, then the fitted parameter should also be $\lambda$, within uncertainties. Finally, the $\chi^2$ should be roughly equal to $N_{\rm bins}$, and the p-value (probability, confidence level, C.L.) should have an expectation value of 50%, but certainly above 0.1% in 99.9% of the cases. We test these hypotheses below.

### Summary table 
Values for a 400k bin histogram, filled with uniformly distributed random values. The C.L. is calculated based on the number of degrees of freedom (NDF) being 399999 (400k bins minus one fitted parameter)

|	Case	|	Total Hits	|	Expected Value	|	Fitted Value	|	$\chi^2$	|	$\chi^2$/NDF	|	Prob (C.L.)	|
|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|
|	Default	|	40M	|	100	|	98.99±0.02	|	403752.41	|	1.0094	|	0.001%	|
|		|	80M	|	200	|	198.99±0.02	|	403065.34	|	1.0077	|	0.03%	|
|		|	160M	|	400	|	399.00±0.03	|	401828.24	|	1.0046	|	2.0%	|
|	Nagy, Csanád, et al.	|	40M	|	100	|	99.00±0.02	|	401681.68	|	1.0042	|	3.0%	|
|		|	80M	|	200	|	198.99±0.02	|	402044.43	|	1.0051	|	1.1%	|
|		|	160M	|	400	|	399.00±0.03	|	401322.71	|	1.0033	|	7.0%	|
|	Pearson	|	40M	|	100	|	100.50±0.02	|	398785.86	|	0.9970	|	91.3%	|
|		|	80M	|	200	|	200.50±0.02	|	400412.49	|	1.0010	|	32.2%	|
|		|	160M	|	400	|	400.50±0.03	|	400526.10	|	1.0013	|	27.8%	|
|	Yates & Pearson	|	40M	|	100	|	100.49±0.02	|	367947.61	|	0.9199	|	100%	|
|		|	80M	|	200	|	200.49±0.02	|	378357.42	|	0.9459	|	100%	|
|		|	160M	|	400	|	400.49±0.03	|	384823.02	|	0.9621	|	100%	|
|	Yates	|	40M	|	100	|	99.00±0.02	|	372521.66	|	0.9313	|	100%	|
|		|	80M	|	200	|	199.00±0.02	|	380866.03	|	0.9522	|	100%	|
|		|	160M	|	400	|	399.00±0.03	|	386072.57	|	0.9652	|	100%	|
|	Likelihood	|	40M	|	100	|	100.00±0.02	|	400405.99	|	1.0010	|	32.4%	|
|		|	80M	|	200	|	200.00±0.02	|	401284.44	|	1.0032	|	7.5%	|
|		|	160M	|	400	|	400.00±0.03	|	400956.03	|	1.0024	|	14.2%	|


## Plots

### Default

<img height="250" alt="poisson_vs_gauss_test_40M_default" src="https://github.com/user-attachments/assets/7d434414-60d9-475c-b2ac-b98a41210051" />
<img height="250" alt="poisson_vs_gauss_test_80M_default" src="https://github.com/user-attachments/assets/16f7453e-1153-4d00-86aa-2babaf487adb" />
<img height="250" alt="poisson_vs_gauss_test_160M_default" src="https://github.com/user-attachments/assets/39bdf653-c126-48ce-a840-e85b4f0a8c8a" />

### Yates' Correction

<img height="250" alt="poisson_vs_gauss_test_40M_Yates" src="https://github.com/user-attachments/assets/b15a6211-0360-420f-b304-64cccd395e88" />
<img height="250" alt="poisson_vs_gauss_test_80M_Yates" src="https://github.com/user-attachments/assets/54bde82e-904b-4661-8090-d36324d9507a" />
<img height="250" alt="poisson_vs_gauss_test_160M_Yates" src="https://github.com/user-attachments/assets/fe1cf2df-7355-4455-8087-e40c186973d5" />

### Pearson's Chisquare

<img height="250" alt="poisson_vs_gauss_test_40M_Pearson" src="https://github.com/user-attachments/assets/bd6f8fc8-0f87-4721-92c1-197e929d4682" />
<img height="250" alt="poisson_vs_gauss_test_80M_Pearson" src="https://github.com/user-attachments/assets/969bc044-5b1e-43ec-867c-3bc25b4b1ee1" />
<img height="250" alt="poisson_vs_gauss_test_160M_Pearson" src="https://github.com/user-attachments/assets/00a20372-c097-470b-8314-799281fb62b4" />

### Yates & Pearson

<img height="250" alt="poisson_vs_gauss_test_40M_YatesMod" src="https://github.com/user-attachments/assets/d1cad37b-f00d-42b8-8e96-34356d1fe8dc" />
<img height="250" alt="poisson_vs_gauss_test_80M_YatesMod" src="https://github.com/user-attachments/assets/4b54df50-c21c-4a91-91cc-ec4063213979" />
<img height="250" alt="poisson_vs_gauss_test_160M_YatesMod" src="https://github.com/user-attachments/assets/9886a03e-ffca-4392-9505-92189ace05c3" />

### Nagy, Csanád, et al.

<img height="250" alt="poisson_vs_gauss_test_40M_corr" src="https://github.com/user-attachments/assets/b9015830-1386-4759-a8a8-8ff7dadcbdf2" />
<img height="250" alt="poisson_vs_gauss_test_80M_corr" src="https://github.com/user-attachments/assets/2b68e0cd-9093-4eba-9797-a45f2ac7bf9a" />
<img height="250" alt="poisson_vs_gauss_test_160M_corr" src="https://github.com/user-attachments/assets/172c9c5e-28bc-4263-990a-8f63c0a34de3" />

### Likelihood

<img height="250" alt="poisson_vs_gauss_test_40M_LL" src="https://github.com/user-attachments/assets/eef8b14e-240a-43a5-941f-f420f3520f66" />
<img height="250" alt="poisson_vs_gauss_test_80M_LL" src="https://github.com/user-attachments/assets/039a546e-b61d-4c77-9b90-0bc910519063" />
<img height="250" alt="poisson_vs_gauss_test_160M_LL" src="https://github.com/user-attachments/assets/f98a3569-76ff-4cb2-a3d3-6bdf6661bb95" />





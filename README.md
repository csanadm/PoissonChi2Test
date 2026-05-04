# Poisson $\chi^2$ test

## File contents
- [**poisson_vs_gauss_test.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test.C): Fitting a histogram with various $\chi^2$ definitions
- [**poisson_vs_gauss_test_ratio.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test_ratio.C): Fitting a ratio of histograms with various $\chi^2$ definitions

## Introduction
This page shows a test of Poisson uncertainties and $\chi^2$ definitions.

When fitting a histogram with bin contents representing (Poisson distributed) counts:
- $\chi^2_i$ is the chisquare for the $i$-th bin,
- $O_i$ is the observed count in the $i$-th bin,
- $\sigma^2_i = (\Delta O_i)^2 = O_i$ is the variance in the $i$-th bin (assuming a Poisson distribution)
- $E_i$ is the expected (theoretical) value, coming from the fit,
- and an alternative variance estimator (à la Pearson) is then $E_i$

When fitting a ratio of (Poisson distributed) histograms, then
- $\chi^2_i$ is still the chisquare for the $i$-th bin,
- $A_i$ and $B_i$ are the counts in the two histograms, while $C_i = A_i/B_i$,
- $\sigma^2_i = (\Delta C_i)^2 = C_i\cdot(1/A_i + 1/B_i) = A_i\cdot(1+A_i/B_i)$ is the variance in the $i$-th bin
- $F_i$ is the fitted value of $C_i$, hence $F_i\cdot B_i$ is the expected count, to be compared to $A_i$

## Formulas

### Default
The classical $\chi^2$ definition is:<br>
$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}}$<br>
For fitting a ratio, an equivalent expression is:<br>
$\chi^2_i = \frac{(A_{i}-F_{i}\cdot B_{i})^{2}}{A_{i}\cdot(1-A_i/B_i)}$<br>
This generally results in an underestimation of the bin contents, due to the so-called Neyman bias.

### Yates' Correction
One can correct for the Neyman bias by subtracting 1/2 from the numerator:<br>
$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{O_{i}}$ <br>
(where the numerator is zero if $|O_{i}-E_{i}|<0.5$)<br>
For fitting a ratio, an equivalent expression is similarly:<br>
$\chi^2_i = \frac{(|A_{i}-F_{i}\cdot B_{i}|-0.5)^{2}}{A_{i}\cdot(1-A_i/B_i)}$

### Pearson's Chisquare
Alternatively, one can use $E_i$ for the variance, which also removes the Neyman bias:<br>
$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{E_{i}}$<br>
For the ratio fit:<br>
$\chi^2_i = \frac{(A_{i}-F_{i}\cdot B_{i})^{2}}{F_i\cdot B_i\cdot(1-A_i/B_i)}$

### Yates & Pearson
One can combine the above two:<br>
$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{E_{i}}$<br>
(where the numerator is zero if $|O_{i}-E_{i}|<0.5$)<br>
For fitting a ratio, the equivalent expression is:<br>
$\chi^2_i = \frac{(|A_{i}-F_{i}\cdot B_{i}|-0.5)^{2}}{F_i\cdot B_i\cdot(1-A_i/B_i)}$

### Nagy, Csanád, et al.
In addition, one can investigate which Gaussian distribution is closest to a given Poisson distribution with parameter $\lambda$. It turns out, that this is the Gaussian with a variance shifted by 1/2, and one can include this in the chisquare definition:<br>
$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}+0.5}$<br>
For fitting a ratio, this becomes:<br>
$\chi^2_i = \frac{(A_{i}-F_{i}\cdot B_{i})^{2}}{A_{i}\cdot(1-A_i/B_i)+0.5}$

### Log likelihood
Finally, to remedy the problems associated with assuming a Gaussian uncertainty distribution, log likelihood optimization can be utilized. In this case, one defines chisquare as:<br>
$\chi^2_i = 2\left[ E_{i} - O_{i} + O_{i} \log(O_{i}/E_{i}) \right]$<br>
For fitting a ratio:<br>
$\chi^2_i = \left[ F_{i}\cdot B_i - A_{i} + A_{i} \log(A_{i}/(F_i\cdot B_{i})) \right]$


## Testing the $\chi^2$ options

In the attached code [**poisson_vs_gauss_test.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test.C) we create a histogram with $N_{\rm bins}=400 000$ bins, and then fill it with $N_{\rm hits}$ uniformly distributed random values. In principle, if $\lambda = $N_{\rm hits}/N_{\rm bins}$, then we expect $\lambda$ hits in every bin. Furthermore, if we fit the resulting histogram with a constant, based on one of the $\chi^2$ definitions, then the fitted parameter should also be $\lambda$, within uncertainties. Finally, the $\chi^2$ should be roughly equal to $N_{\rm bins}$, and the p-value (probability, confidence level, C.L.) should have an expectation value of 50%, but certainly above 0.1% in 99.9% of the cases. We test these hypotheses below.

For the ratio version, a separate [**poisson_vs_gauss_test_ratio.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test_ratio.C) code is written, and the results are analyzed separately. Note that here the expected value is 1.0, as the ratio is of two histograms with the same expected contents.

### Summary table 
Values for a 400k bin histogram, filled with uniformly distributed random values. The C.L. is calculated based on the number of degrees of freedom (NDF) being 399999 (400k bins minus one fitted parameter)

|	Case	|	Total Hits	|	Expected Value	|	Fitted Value	|	$\chi^2$	|	$\chi^2$/NDF	|	Prob (C.L.)	|
|	-----	|	-----------	|	--------------	|	-------------	|	---------	|	-------------	|	-----------	|
| Default | 40M | 100 | 98.991±0.016 | 403752.41 | 1.0094 | 0.001% |
|   | 80M | 200 | 198.992±0.022 | 403065.34 | 1.0077 | 0.03% |
|   | 160M | 400 | 398.995±0.032 | 401828.24 | 1.0046 | 2.05% |
| Nagy, Csanád, et al. | 40M | 100 | 98.996±0.016 | 401681.68 | 1.0042 | 3.01% |
|   | 80M | 200 | 198.995±0.022 | 402044.43 | 1.0051 | 1.12% |
|   | 160M | 400 | 398.997±0.032 | 401322.71 | 1.0033 | 6.96% |
| Pearson | 40M | 100 | 100.498±0.016 | 398785.86 | 0.9970 | 91.26% |
|   | 80M | 200 | 200.501±0.022 | 400412.49 | 1.0010 | 32.17% |
|   | 160M | 400 | 400.501±0.032 | 400526.10 | 1.0013 | 27.77% |
| Yates & Pearson | 40M | 100 | 100.486±0.016 | 367947.61 | 0.9199 | 100.00% |
|   | 80M | 200 | 200.492±0.023 | 378357.42 | 0.9459 | 100.00% |
|   | 160M | 400 | 400.494±0.032 | 384823.02 | 0.9621 | 100.00% |
| Yates | 40M | 100 | 98.998±0.016 | 372521.66 | 0.9313 | 100.00% |
|   | 80M | 200 | 198.997±0.023 | 380866.03 | 0.9522 | 100.00% |
|   | 160M | 400 | 398.998±0.032 | 386072.57 | 0.9652 | 100.00% |
| Likelihood | 40M | 100 | 100.000±0.016 | 400405.99 | 1.0010 | 32.43% |
|   | 80M | 200 | 200.000±0.022 | 401284.44 | 1.0032 | 7.54% |
|   | 160M | 400 | 400.000±0.032 | 400956.03 | 1.0024 | 14.23% |

The same table for fitting a ratio of histograms:

|	Case	|	Total Hits	|	Fitted Value	|	$\chi^2$	|	$\chi^2$/NDF	|	Prob (C.L.)	|
|	-----	|	-----------	|	-------------	|	---------	|	-------------	|	-----------	|
| Default | 40M | 0.970±0.000 | 395655.84 | 0.9891 | 100.00% |
|   | 80M | 0.985±0.000 | 396850.29 | 0.9921 | 99.98% |
|   | 160M | 0.993±0.000 | 398584.53 | 0.9965 | 94.32% |
| Nagy, Csanád, et al. | 40M | 0.970±0.000 | 394630.73 | 0.9866 | 100.00% |
|   | 80M | 0.985±0.000 | 396345.75 | 0.9909 | 100.00% |
|   | 160M | 0.993±0.000 | 398333.16 | 0.9958 | 96.89% |
| Pearson | 40M | 1.000±0.000 | 401110.93 | 1.0028 | 10.70% |
|   | 80M | 1.000±0.000 | 399787.85 | 0.9995 | 59.30% |
|   | 160M | 1.000±0.000 | 399818.00 | 0.9995 | 57.99% |
| Yates & Pearson | 40M | 1.000±0.000 | 378947.96 | 0.9474 | 100.00% |
|   | 80M | 1.000±0.000 | 384061.88 | 0.9602 | 100.00% |
|   | 160M | 1.000±0.000 | 388657.51 | 0.9716 | 100.00% |
| Yates | 40M | 0.970±0.000 | 373538.20 | 0.9338 | 100.00% |
|   | 80M | 0.985±0.000 | 381142.64 | 0.9529 | 100.00% |
|   | 160M | 0.993±0.000 | 387426.69 | 0.9686 | 100.00% |
| Likelihood | 40M | 1.000±0.000 | 403109.30 | 1.0078 | 0.03% |
|   | 80M | 1.000±0.000 | 400796.91 | 1.0020 | 18.61% |
|   | 160M | 1.000±0.000 | 400279.85 | 1.0007 | 37.65% |


## Plots for histogram fitting

### Default

<img height="250" alt="poisson_vs_gauss_test_40M_Classical" src="https://github.com/user-attachments/assets/15ccc13f-57c3-4da8-9e92-e739afe3e76b" />
<img height="250" alt="poisson_vs_gauss_test_80M_Classical" src="https://github.com/user-attachments/assets/e8abe8c2-cfda-4881-9e1b-62646d25c8bb" />
<img height="250" alt="poisson_vs_gauss_test_160M_Classical" src="https://github.com/user-attachments/assets/6d564f97-28a7-4536-9a20-ea13cc6946c8" />

### Yates' Correction

<img height="250" alt="poisson_vs_gauss_test_40M_Yates" src="https://github.com/user-attachments/assets/36325922-ef9e-45e7-9861-281f387848a2" />
<img height="250" alt="poisson_vs_gauss_test_80M_Yates" src="https://github.com/user-attachments/assets/efee6942-dd44-4754-a3d2-e9ec716942f5" />
<img height="250" alt="poisson_vs_gauss_test_160M_Yates" src="https://github.com/user-attachments/assets/d6d7e87d-5ca4-45c2-a64f-a188d0ceda0e" />

### Pearson's Chisquare

<img height="250" alt="poisson_vs_gauss_test_40M_Pearson" src="https://github.com/user-attachments/assets/f63c1a97-28ac-4735-8741-1caf51931fdb" />
<img height="250" alt="poisson_vs_gauss_test_80M_Pearson" src="https://github.com/user-attachments/assets/1e5acd58-e651-4a8b-8a36-89d4271d2ef0" />
<img height="250" alt="poisson_vs_gauss_test_160M_Pearson" src="https://github.com/user-attachments/assets/5234dbb7-e193-4507-880e-a6f4668c8d7c" />

### Yates & Pearson

<img height="250" alt="poisson_vs_gauss_test_40M_Yates_Pearson" src="https://github.com/user-attachments/assets/0f52919d-b91f-4bd3-9435-2cf7b2e35049" />
<img height="250" alt="poisson_vs_gauss_test_80M_Yates_Pearson" src="https://github.com/user-attachments/assets/d182d22e-ec61-4e7d-a869-b06fce4c0ba7" />
<img height="250" alt="poisson_vs_gauss_test_160M_Yates_Pearson" src="https://github.com/user-attachments/assets/3b57a22c-0576-4cb5-9b70-973cccba0e22" />

### Nagy, Csanád, et al.

<img height="250" alt="poisson_vs_gauss_test_40M_Corrected" src="https://github.com/user-attachments/assets/c260486a-9d8a-4e51-9288-f96ca9991272" />
<img height="250" alt="poisson_vs_gauss_test_80M_Corrected" src="https://github.com/user-attachments/assets/455bb7ca-343f-4b13-9590-47ade468d216" />
<img height="250" alt="poisson_vs_gauss_test_160M_Corrected" src="https://github.com/user-attachments/assets/95264171-0927-4612-8c1f-9a249d24ae3a" />

### Likelihood

<img height="250" alt="poisson_vs_gauss_test_40M_LogLikelihood" src="https://github.com/user-attachments/assets/62ba6b68-e146-451b-b519-dd75bc12f373" />
<img height="250" alt="poisson_vs_gauss_test_80M_LogLikelihood" src="https://github.com/user-attachments/assets/aa37b14d-715c-49bd-86be-091ed40c737b" />
<img height="250" alt="poisson_vs_gauss_test_160M_LogLikelihood" src="https://github.com/user-attachments/assets/25429dcc-8948-495f-aa6e-3d7ab31ad972" />

## Plots for histogram fitting

### Default

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Classical" src="https://github.com/user-attachments/assets/814f0bab-96bd-4345-89da-cbea4ee5b8b5" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Classical" src="https://github.com/user-attachments/assets/87fb13a9-d8a7-4be6-8661-29395828345f" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Classical" src="https://github.com/user-attachments/assets/cb1788d6-ccef-4366-9c97-d5a20cb62385" />

### Yates' Correction

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Yates" src="https://github.com/user-attachments/assets/d7ae8d01-6a9c-49fe-b1a4-b5c5972664b7" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Yates" src="https://github.com/user-attachments/assets/995fc145-05af-44db-8186-86cce5470490" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Yates" src="https://github.com/user-attachments/assets/8d4f7d52-fb12-4bea-a031-e9581893e184" />

### Pearson's Chisquare

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Pearson" src="https://github.com/user-attachments/assets/74736547-54c7-4d8d-81f7-aae4ecb30659" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Pearson" src="https://github.com/user-attachments/assets/8b63147d-bb06-441a-a85b-4bb614a83b96" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Pearson" src="https://github.com/user-attachments/assets/08380dcc-ac48-4e77-a8c8-498c9e4044cb" />

### Yates & Pearson

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Yates_Pearson" src="https://github.com/user-attachments/assets/dc98f877-6c63-42f4-8291-4cd2b3a6660e" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Yates_Pearson" src="https://github.com/user-attachments/assets/6c53b461-d4c5-4971-8a8b-6776885e28df" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Yates_Pearson" src="https://github.com/user-attachments/assets/ca2d1110-d631-479f-b7ad-16a3132693f6" />

### Nagy, Csanád, et al.

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Corrected" src="https://github.com/user-attachments/assets/dc3b9934-e43a-473a-9464-36865d0f5d7d" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Corrected" src="https://github.com/user-attachments/assets/3f204206-7202-44bf-840e-a558e8490ae5" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Corrected" src="https://github.com/user-attachments/assets/f00d9ca7-fa9f-4aef-bcf0-7885a5c02ddf" />

### Likelihood

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_LogLikelihood" src="https://github.com/user-attachments/assets/854cfbca-4abf-4c6c-9136-a15e76f3af63" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_LogLikelihood" src="https://github.com/user-attachments/assets/0f5768a4-d1ce-4c1c-be08-1c09aa3e2c47" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_LogLikelihood" src="https://github.com/user-attachments/assets/3c5ea285-4c24-489c-b32e-bfc56b7886ee" />








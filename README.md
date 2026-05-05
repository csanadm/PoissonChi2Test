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
To remedy the problems associated with assuming a Gaussian uncertainty distribution, log likelihood optimization can be utilized. In this case, one defines chisquare as:<br>
$\chi^2_i = 2\left[ E_{i} - O_{i} + O_{i} \log(O_{i}/E_{i}) \right]$<br>
For fitting a ratio:<br>
$\chi^2_i = \left[ F_{i}\cdot B_i - A_{i} + A_{i} \log(A_{i}/(F_i\cdot B_{i})) \right]$

### E802 likelihood
Finally, for fitting a ratio, an additional definition is tested, based on nucl-ex/0204001 from the E802 Collaboration:<br>
$\chi^2_i = -2 \left[ A_i\cdot \log\left(\frac{C_i (A_i + B_i)}{A_i (C_i + 1)}\right) + B_i\cdot \log\left(\frac{A_i + B_i}{B_i (C_i + 1)}\right) \right]$


## Testing the $\chi^2$ options

In the attached code [**poisson_vs_gauss_test.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test.C) we create a histogram with $N_{\rm bins}=400 000$ bins, and then fill it with $N_{\rm hits}$ uniformly distributed random values. In principle, if $\lambda = $N_{\rm hits}/N_{\rm bins}$, then we expect $\lambda$ hits in every bin. Furthermore, if we fit the resulting histogram with a constant, based on one of the $\chi^2$ definitions, then the fitted parameter should also be $\lambda$, within uncertainties. Finally, the $\chi^2$ should be roughly equal to $N_{\rm bins}$, and the p-value (probability, confidence level, C.L.) should have an expectation value of 50%, but certainly above 0.1% in 99.9% of the cases. We test these hypotheses below.

For the ratio version, a separate [**poisson_vs_gauss_test_ratio.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test_ratio.C) code is written, and the results are analyzed separately. Note that here the expected value is 1.0, as the ratio is of two histograms with the same expected contents.

### Summary table 
Values for a 400k bin histogram, filled with uniformly distributed random values. The C.L. is calculated based on the number of degrees of freedom (NDF) being 399999 (400k bins minus one fitted parameter)

|	Case                 |	Total Hits	|	Expected Value	|	Fitted Value	|	$\chi^2$	|	$\chi^2$/NDF	|	Prob (C.L.)	|
|	-------------------- |	-----------	|	--------------	|	-------------	|	---------	|	-------------	|	-----------	|
| Default              | 40M  | 100 |  98.991±0.016 | 403752.41 | 1.0094 | 0.001% |
|                      | 80M  | 200 | 198.992±0.022 | 403065.34 | 1.0077 | 0.03% |
|                      | 160M | 400 | 398.995±0.032 | 401828.24 | 1.0046 | 2.05% |
| Nagy, Csanád, et al. | 40M  | 100 |  98.996±0.016 | 401681.68 | 1.0042 | 3.01% |
|                      | 80M  | 200 | 198.995±0.022 | 402044.43 | 1.0051 | 1.12% |
|                      | 160M | 400 | 398.997±0.032 | 401322.71 | 1.0033 | 6.96% |
| Pearson              | 40M  | 100 | 100.498±0.016 | 398785.86 | 0.9970 | 91.26% |
|                      | 80M  | 200 | 200.501±0.022 | 400412.49 | 1.0010 | 32.17% |
|                      | 160M | 400 | 400.501±0.032 | 400526.10 | 1.0013 | 27.77% |
| Yates & Pearson      | 40M  | 100 | 100.486±0.016 | 367947.61 | 0.9199 | 100.00% |
|                      | 80M  | 200 | 200.492±0.023 | 378357.42 | 0.9459 | 100.00% |
|                      | 160M | 400 | 400.494±0.032 | 384823.02 | 0.9621 | 100.00% |
| Yates                | 40M  | 100 |  98.998±0.016 | 372521.66 | 0.9313 | 100.00% |
|                      | 80M  | 200 | 198.997±0.023 | 380866.03 | 0.9522 | 100.00% |
|                      | 160M | 400 | 398.998±0.032 | 386072.57 | 0.9652 | 100.00% |
| Likelihood           | 40M  | 100 | 100.000±0.016 | 400405.99 | 1.0010 | 32.43% |
|                      | 80M  | 200 | 200.000±0.022 | 401284.44 | 1.0032 | 7.54% |
|                      | 160M | 400 | 400.000±0.032 | 400956.03 | 1.0024 | 14.23% |

The same table for fitting a ratio of histograms:

|	Case            |	Total Hits	|	Fitted Value	|	$\chi^2$	|	$\chi^2$/NDF	|	Prob (C.L.)	|
|	----------------|	-----------	|	-------------	|	---------	|	-------------	|	-----------	|
| Default              | 40M  | 0.9703±0.0002 | 395655.84 | 0.9891 | 100.00% |
|                      | 80M  | 0.9851±0.0002 | 396850.29 | 0.9921 | 99.98% |
|                      | 160M | 0.9925±0.0001 | 398584.53 | 0.9965 | 94.32% |
| Nagy, Csanád, et al. | 40M  | 0.9704±0.0002 | 394630.73 | 0.9866 | 100.00% |
|                      | 80M  | 0.9851±0.0002 | 396345.75 | 0.9909 | 100.00% |
|                      | 160M | 0.9925±0.0001 | 398333.16 | 0.9958 | 96.89% |
| Pearson              | 40M  | 1.0000±0.0002 | 401110.93 | 1.0028 | 10.70% |
|                      | 80M  | 1.0000±0.0002 | 399787.85 | 0.9995 | 59.30% |
|                      | 160M | 1.0000±0.0001 | 399818.00 | 0.9995 | 57.99% |
| Yates & Pearson      | 40M  | 1.0000±0.0002 | 378947.96 | 0.9474 | 100.00% |
|                      | 80M  | 1.0000±0.0002 | 384061.88 | 0.9602 | 100.00% |
|                      | 160M | 1.0000±0.0001 | 388657.51 | 0.9716 | 100.00% |
| Yates                | 40M  | 0.9703±0.0002 | 373538.20 | 0.9338 | 100.00% |
|                      | 80M  | 0.9851±0.0002 | 381142.64 | 0.9529 | 100.00% |
|                      | 160M | 0.9925±0.0001 | 387426.69 | 0.9686 | 100.00% |
| Likelihood           | 40M  | 1.0000±0.0002 | 403109.30 | 1.0078 | 0.03% |
|                      | 80M  | 1.0000±0.0002 | 400796.91 | 1.0020 | 18.61% |
|                      | 160M | 1.0000±0.0001 | 400279.85 | 1.0007 | 37.65% |
| CF Likelihood        | 40M  | 1.0000±0.0002 | 402131.56 | 1.0053 | 0.86% |
|                      | 80M  | 1.0000±0.0002 | 400290.77 | 1.0007 | 37.19% |
|                      | 160M | 1.0000±0.0001 | 400069.11 | 1.0002 | 46.85% |


## Plots for histogram fitting

### Default

<img height="250" alt="poisson_vs_gauss_test_40M_Classical" src="https://github.com/user-attachments/assets/9af979fe-00a1-4f8d-b5d0-bcf1c91f63d2" />
<img height="250" alt="poisson_vs_gauss_test_80M_Classical" src="https://github.com/user-attachments/assets/568b830e-d6fb-470e-a82a-89dff07aa877" />
<img height="250" alt="poisson_vs_gauss_test_160M_Classical" src="https://github.com/user-attachments/assets/c294b245-fe8c-48c1-a8b1-ccb7cb85c87c" />

### Yates' Correction

<img height="250" alt="poisson_vs_gauss_test_40M_Yates" src="https://github.com/user-attachments/assets/58ca0096-4c21-4109-9bbd-c01fe255fc23" />
<img height="250" alt="poisson_vs_gauss_test_80M_Yates" src="https://github.com/user-attachments/assets/88b3fb39-004a-4f93-8a17-03777948e9a4" />
<img height="250" alt="poisson_vs_gauss_test_160M_Yates" src="https://github.com/user-attachments/assets/24e6e234-326a-4495-9aad-15563f1aa721" />

### Pearson's Chisquare

<img height="250" alt="poisson_vs_gauss_test_40M_Pearson" src="https://github.com/user-attachments/assets/684c6264-1563-4888-a5c4-3647d699c10d" />
<img height="250" alt="poisson_vs_gauss_test_80M_Pearson" src="https://github.com/user-attachments/assets/62601037-c39d-430d-b2d2-5f33b072dc20" />
<img height="250" alt="poisson_vs_gauss_test_160M_Pearson" src="https://github.com/user-attachments/assets/c99a89ef-75e9-46e2-974b-13d7b4ce29b1" />

### Yates & Pearson

<img height="250" alt="poisson_vs_gauss_test_40M_Yates-Pearson" src="https://github.com/user-attachments/assets/9327fdd2-7bae-4c2f-853a-2357ada32c75" />
<img height="250" alt="poisson_vs_gauss_test_80M_Yates-Pearson" src="https://github.com/user-attachments/assets/a96e074b-0b1b-4dd0-a0d9-2066088dff69" />
<img height="250" alt="poisson_vs_gauss_test_160M_Yates-Pearson" src="https://github.com/user-attachments/assets/7f593d6b-3979-444d-9df3-7c9f8fa2ce03" />

### Nagy, Csanád, et al.

<img height="250" alt="poisson_vs_gauss_test_40M_Corrected" src="https://github.com/user-attachments/assets/4c7bbd53-4d9f-4349-981a-8f4527253b53" />
<img height="250" alt="poisson_vs_gauss_test_80M_Corrected" src="https://github.com/user-attachments/assets/cc57d52e-79ee-4a87-83cf-e184adbdabb1" />
<img height="250" alt="poisson_vs_gauss_test_160M_Corrected" src="https://github.com/user-attachments/assets/d91ca37c-81bc-4327-b9d8-646e61b322c4" />

### Likelihood

<img height="250" alt="poisson_vs_gauss_test_40M_LogLikelihood" src="https://github.com/user-attachments/assets/c3ded664-1bd1-460e-97f5-7914c4b7df8c" />
<img height="250" alt="poisson_vs_gauss_test_80M_LogLikelihood" src="https://github.com/user-attachments/assets/739f5ee6-432f-428a-8d5c-9e2bf564316b" />
<img height="250" alt="poisson_vs_gauss_test_160M_LogLikelihood" src="https://github.com/user-attachments/assets/270aba9b-2cfc-4410-92d8-2355f62041eb" />

## Plots for histogram fitting

### Default

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Classical" src="https://github.com/user-attachments/assets/613a0403-c6fe-4252-9c25-405d3aff1c22" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Classical" src="https://github.com/user-attachments/assets/3e1b415b-f98c-4e49-bb3a-c5aa089749ce" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Classical" src="https://github.com/user-attachments/assets/9f436965-5ff2-4b18-bf0c-e60b11ee89bc" />

### Yates' Correction

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Yates" src="https://github.com/user-attachments/assets/147edf58-65c0-49a5-b4ca-e7b11d493ba2" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Yates" src="https://github.com/user-attachments/assets/ff54d567-de49-4ffa-918d-af4e0ab2bc2d" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Yates" src="https://github.com/user-attachments/assets/8bfc6b3b-7769-4419-b5b9-1080d6a7ed91" />

### Pearson's Chisquare

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Pearson" src="https://github.com/user-attachments/assets/9f122f20-4de6-4bd0-b2e8-628b7e04c082" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Pearson" src="https://github.com/user-attachments/assets/4e4940ec-c3bb-46bb-ac22-1c1742d9cdfa" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Pearson" src="https://github.com/user-attachments/assets/b6fdf567-5513-4d52-9f05-a1fffd72890d" />

### Yates & Pearson

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Yates-Pearson" src="https://github.com/user-attachments/assets/146d1431-6412-44da-857d-29b5a04fd111" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Yates-Pearson" src="https://github.com/user-attachments/assets/510de90a-84df-4b96-a943-1ab0643e271b" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Yates-Pearson" src="https://github.com/user-attachments/assets/01599939-7be0-41b2-85ae-90347c802186" />

### Nagy, Csanád, et al.

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_Corrected" src="https://github.com/user-attachments/assets/7060c1c9-9e42-4b9f-894f-4cbe08eacefb" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_Corrected" src="https://github.com/user-attachments/assets/a09a0183-8790-4120-9621-86e426ade6f3" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_Corrected" src="https://github.com/user-attachments/assets/9549dfd0-5d54-4e6a-98af-0f49b724ed30" />

### Likelihood

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_LogLikelihood" src="https://github.com/user-attachments/assets/c8d16d14-6985-4d43-9386-c203c65e39f8" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_LogLikelihood" src="https://github.com/user-attachments/assets/fa2fb59b-6954-4add-b526-d7e5379fdf65" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_LogLikelihood" src="https://github.com/user-attachments/assets/5d9c258f-c925-48f5-9407-31d53020cf15" />

### CF Likelihood

<img height="250" alt="poisson_vs_gauss_test_ratio_40M_CFLikelihood" src="https://github.com/user-attachments/assets/3286b8d1-2702-4d87-9c38-a16b47a23bf6" />
<img height="250" alt="poisson_vs_gauss_test_ratio_80M_CFLikelihood" src="https://github.com/user-attachments/assets/f8e52b51-88fc-43dc-9a1f-d5782e045219" />
<img height="250" alt="poisson_vs_gauss_test_ratio_160M_CFLikelihood" src="https://github.com/user-attachments/assets/6c56e5d8-43b0-43c9-8a12-506b9c8ba821" />









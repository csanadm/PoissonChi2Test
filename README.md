# Poisson $\chi^2$ test

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
For the ratio fit:
$\chi^2_i = \frac{(A_{i}-F_{i}\cdot B_{i})^{2}}{F_i\cdot B_i\cdot(1-A_i/B_i)}$

### Yates & Pearson
One can combine the above two:<br>
$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{E_{i}}$<br>
(where the numerator is zero if $|O_{i}-E_{i}|<0.5$)<br>
For fitting a ratio, the equivalent expression is:<br>
$\chi^2_i = \frac{(|A_{i}-F_{i}\cdot B_{i}|-0.5)^{2}}{F_i\cdot B_i\cdot(1-A_i/B_i)}$

### Nagy, Csanád, et al.
In addition, one can investigate which Gaussian distribution is closest to a given Poisson distribution with parameter $\lambda$. It turns out, that this is the Gaussian with a variance shifted by 1/2, and one can include this in the chisquare definition:<br>
$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}+0.5}$
For fitting a ratio, this becomes:<br>
$\chi^2_i = \frac{(A_{i}-F_{i}\cdot B_{i})^{2}}{A_{i}\cdot(1-A_i/B_i)+0.5}$

### Log likelihood
Finally, to remedy the problems associated with assuming a Gaussian uncertainty distribution, log likelihood optimization can be utilized. In this case, one defines chisquare as:<br>
$\chi^2_i = 2\left[ E_{i} - O_{i} + O_{i} \log(O_{i}/E_{i}) \right]$
For fitting a ratio:<br>
$\chi^2_i = \left[ F_{i}\cdot B_i - A_{i} + A_{i} \log(A_{i}/(F_i\cdot B_{i})) \right]$


## Testing the $\chi^2$ options

In the attached code [**poisson_vs_gauss_test.C**](https://github.com/csanadm/PoissonChi2Test/blob/main/poisson_vs_gauss_test.C) we create a histogram with $N_{\rm bins}=400 000$ bins, and then fill it with $N_{\rm hits}$ uniformly distributed random values. In principle, if $\lambda = $N_{\rm hits}/N_{\rm bins}$, then we expect $\lambda$ hits in every bin. Furthermore, if we fit the resulting histogram with a constant, based on one of the $\chi^2$ definitions, then the fitted parameter should also be $\lambda$, within uncertainties. Finally, the $\chi^2$ should be roughly equal to $N_{\rm bins}$, and the p-value (probability, confidence level, C.L.) should have an expectation value of 50%, but certainly above 0.1% in 99.9% of the cases. We test these hypotheses below.

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





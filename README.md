# Poisson $\chi^2$ test
Testing Poisson uncertainties and $\chi^2$ definitions

## Default

$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}}$

<img height="300" alt="poisson_vs_gauss_test_40M" src="https://github.com/user-attachments/assets/40289048-bb72-42ce-aaed-72fedce94eb9" />
<img height="300" alt="poisson_vs_gauss_test_80M" src="https://github.com/user-attachments/assets/5393ecf4-4f7d-4388-b0c1-b66b160c0eaf" />
<img height="300" alt="poisson_vs_gauss_test_160M" src="https://github.com/user-attachments/assets/d8aa7ecb-5a9a-4d7a-bfc9-da2cd539efc4" />

## Yates' Correction

$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{O_{i}}$ (where the numerator is zero if $|O_{i}-E|<0.5$)

<img height="300" alt="poisson_vs_gauss_test_160M_Yates" src="https://github.com/user-attachments/assets/ce4e11a2-d848-4097-88a5-5aeee89b352f" />
<img height="300" alt="poisson_vs_gauss_test_80M_Yates" src="https://github.com/user-attachments/assets/f039fd59-e48e-4b17-b953-346764027210" />
<img height="300" alt="poisson_vs_gauss_test_40M_Yates" src="https://github.com/user-attachments/assets/1c95f79f-1135-4d23-93be-5d3785dd6d5e" />

## Pearson's Chisquare

$\chi^2_i = \frac{(O_{i}-E)^{2}}{E_{i}}$

<img height="300" alt="poisson_vs_gauss_test_80M_Pearson" src="https://github.com/user-attachments/assets/e0ca88e5-a820-4b22-84b0-bafe263cef9b" />
<img height="300" alt="poisson_vs_gauss_test_40M_Pearson" src="https://github.com/user-attachments/assets/1a1b6b22-96d8-4fbb-b13f-9b1f82b89632" />
<img height="300" alt="poisson_vs_gauss_test_160M_Pearson" src="https://github.com/user-attachments/assets/1772a22f-38a7-48d4-a5e2-2d993b8e6fe7" />

## Yates & Pearson

$\chi^2_i = \frac{(|O_{i}-E_{i}|-0.5)^{2}}{E}$ (where the numerator is zero if $|O_{i}-E|<0.5$)

<img height="300" alt="poisson_vs_gauss_test_80M_YatesMod" src="https://github.com/user-attachments/assets/72452c0a-13fc-434e-b118-cba6a1197751" />
<img height="300" alt="poisson_vs_gauss_test_40M_YatesMod" src="https://github.com/user-attachments/assets/9848681f-b7af-4555-9f25-8ff317b844f1" />
<img height="300" alt="poisson_vs_gauss_test_160M_YatesMod" src="https://github.com/user-attachments/assets/db418faa-bab3-4c01-b1d4-6891fd7fee11" />

## Nagy, Csanád, et al.

$\chi^2_i = \frac{(O_{i}-E_{i})^{2}}{O_{i}+0.5}$

<img height="300" alt="poisson_vs_gauss_test_80M_corr" src="https://github.com/user-attachments/assets/6542f7c5-65fd-4dc0-b83f-0419822dd2ba" />
<img height="300" alt="poisson_vs_gauss_test_40M_corr" src="https://github.com/user-attachments/assets/d1757ca3-4d63-427f-9cf5-d8d3fdf4921d" />
<img height="300" alt="poisson_vs_gauss_test_160M_corr" src="https://github.com/user-attachments/assets/6bff0395-5509-45a1-bdbd-7a4560b09729" />







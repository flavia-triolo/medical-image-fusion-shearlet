# Multimodal Medical Image Fusion (CT-MRI) using Shearlet Transform

Implementation and experimental evaluation of a multimodal medical image fusion framework for co-registered Computed Tomography (CT) and Magnetic Resonance Imaging (MRI), developed in MATLAB.

This work was developed as part of my Master's Degree Thesis in Mathematics at the University of Palermo (*Università degli Studi di Palermo*).

---

## Framework Overview

The fusion algorithm integrates multi-scale anisotropic harmonic analysis with edge-preserving spatial filtering:
1. **Multi-Scale Decomposition:** Discrete Shearlet Transform (DST) using a cone-adapted Parseval frame (`ShearLab2D`) to capture curvilinear contours and directional singularities.
2. **Base Layer Fusion:** Two-stage rule combining directional dispersion from local covariance matrices (GFS) and frequency-tuned saliency regularized via **Guided Image Filtering**.
3. **Detail Layer Fusion:** Absolute choose-max rule across directional high-frequency subbands.
4. **Reconstruction:** Inverse Discrete Shearlet Transform (IDST) to reconstruct the fused image.

---

## Repository Structure

* `main_fusion.m`: Core pipeline implementing the two-stage multi-scale fusion.
* `guided_filter.m`: Edge-preserving guided filtering for weight map smoothing.
* `compute_gfs_weights.m`: Local covariance matrix calculation and directional dispersion extraction via trace operator.
* `compute_petrovic_metrics.m`: Objective evaluation metrics simulating the Human Visual System (Petrovic & Xydeas edge association $Q^{AB/F}$, loss $L^{AB/F}$, artifacts $N^{AB/F}$).

---

## Dataset

The experiments were validated on 15 pairs of clinically co-registered axial brain CT and MRI images (256 × 256 pixels) obtained from the open-source repository:
* [Medical-image-fusion-dataset (GitHub)](https://github.com/wenzhezhai/Medical-image-fusion-dataset)

Images were normalized in the continuous interval $[0, 1]$ before multi-scale processing.

---

## Experimental Analysis & Results

The modular pipeline was tested iteratively across the dataset to evaluate performance, sensitivity, and robustness:

### 1. Visual Fusion Results
The shearlet-based method successfully transfers bone detail from CT and soft-tissue contrast from MRI into a single composite representation without introducing blur or contrast loss.

> *<img width="771" height="286" alt="Screenshot 2026-10-15 124304" src="https://github.com/user-attachments/assets/a668f1f3-048d-4bf0-8067-49fe0971ba25" />


---

### 2. Comparative Benchmark (Fourier vs Wavelet vs Shearlet)

The framework was evaluated against classical harmonic representation bases under identical two-stage fusion rules across all 15 clinical pairs:
* **Fourier Domain (2D FFT):** Ideal circular low-pass filter ($R_0 = 16$). Prone to spectral Gibbs ringing artifacts along anatomical boundaries.
* **Orthonormal Wavelet (2D DWT):** 4-level decomposition with Daubechies (db4). Suffering from shift-variance due to dyadic downsampling.
* **Shearlet Domain (2D DST):** Cone-adapted Parseval frame ($J=4$, $s=[1, 1, 2, 2]$). Achieves optimal sparse approximation for curvilinear singularities.

Average quantitative metrics across the dataset:
* **Edge Preservation ($Q^{AB/F}$ ↑):** **Shearlet (0.7044)** vs Fourier (0.6735) vs Wavelet (0.6279)
* **Information Loss ($L^{AB/F}$ ↓):** **Shearlet (0.2393)** vs Fourier (0.2447) vs Wavelet (0.2860)
* **Fusion Artifacts ($N^{AB/F}$ ↓):** **Shearlet (0.2056)** vs Wavelet (0.2819) vs Fourier (0.2628)

#### Visual Comparison (Case 04)
The 5-panel comparison demonstrates how the Shearlet transform successfully preserves both bone boundaries (from CT) and soft-tissue details (from MRI) while suppressing ringing and blocking artifacts:

><img width="555" height="420" alt="Screenshot 2026-10-15 125319" src="https://github.com/user-attachments/assets/7a974a42-1102-4abb-9c24-0e260b077ab3" />
 

---

### 3. Parametric Sensitivity Analysis & Weight Map Regularization

A bivariate grid search was performed over the Guided Filter spatial radius $r \in \\{3, 5, 10, 20, 30\\}$ and regularization factor $\varepsilon \in \\{10^{-6}, 10^{-5}, 10^{-4}, 10^{-3}, 10^{-1}\\}$ relative to the base layer's local variance $(\text{Var}_ {\text{median}} \approx 1.76 \times 10^{-4}, \\: \text{Var}_ {95\\%} \approx 1.50 \times 10^{-3})$.

The analysis identified three distinct operational behaviors:
1. **Under-regularized regime ($\varepsilon \le 10^{-5} \ll \text{Var}_{\text{median}}$):** The linear slope attempts to fit the abrupt step edges of the raw binary saliency map $BW_A$, producing noisy and fractured weights that introduce false edges.
2. **Transition regime ($\varepsilon \in [10^{-4}, 10^{-3}]$):** $\varepsilon$ matches the local variance of brain tissues, actively discriminating between flat areas and structural edges.
3. **Asymptotic regime ($\varepsilon \ge 10^{-2} \gg \text{Var}_{95\\%}$):** The local variance becomes negligible relative to $\varepsilon$; the filter stably converges to a smooth spatially weighted partition without step artifacts.

#### Weight Map Evolution Across Regimes
Visual evolution of the CT weight map $W_{CT}$ across the three regimes compared to the raw binary map $BW_{CT}$:

><img width="645" height="147" alt="Screenshot 2026-10-15 125640" src="https://github.com/user-attachments/assets/f3064985-1005-4722-ba24-69f09ca3792d" />


---

### 4. Robustness and Perturbation Stress-Testing
To assess stability against sensor noise, controlled artificial perturbations were added to the source images before iterative execution of the fusion pipeline:
* **Additive White Gaussian Noise (AWGN):** $\sigma^2 \in \\{0.001, 0.005, 0.01, 0.02\\}$. The uniform angular dispersion of isotropic noise across the 48 shearlet subbands enables the choose-max rule to effectively suppress noise, achieving higher PSNR and lower RMSE compared to wavelets and Fourier.
* **Bipolar Salt & Pepper Impulsive Noise:** corruption densities $d \in \\{1\\%, 2\\%, 5\\%, 10\\%\\}$.

><img width="385" height="287" alt="Screenshot 2026-10-15 125827" src="https://github.com/user-attachments/assets/938b3d98-3861-4911-b80e-f7b29ad42169" />

---

## Requirements

* MATLAB (R2020)
* Image Processing Toolbox
* [ShearLab 2D](http://www.shearlab.org/) library

---

## References

* **Algorithm Baseline:** A. Dogra and S. Kumar, *"Multi-modality medical image fusion based on guided filter and image statistics in multidirectional shearlet transform domain"*, J. Ambient Intell. Humaniz. Comput. 14.9 (2023).
* **Shearlet Framework:** G. Kutyniok, W.-Q. Lim, R. Reisenhofer, *"ShearLab 3D: Faithful digital shearlet transforms based on compactly supported shearlets"*, ACM Trans. Math. Softw. 42.1 (2016).
* **Evaluation Metrics:** V. Petrovic and C. Xydeas, *"Objective image fusion performance characterisation"*, ICCV 2 (2005).
* **Dataset:** W. Zhai, *Medical-image-fusion-dataset*, GitHub repository.

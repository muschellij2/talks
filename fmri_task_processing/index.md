# Neuroconductor Example: fMRI Task Processing
John Muschelli<br/>http://johnmuschelli.com/talks/fmri_task_processing/<br/> Johns Hopkins Bloomberg School of Public Health  
<style type="text/css">
article {
  font-size: 30pt;
}
</style>





## SPM

All of this is using the Statistical Parametric Mapping (SPM) [@penny2011statistical] software version 12:

- requies MATLAB (for this tutorial)
- All called through [`spm12r` package](https://github.com/muschellij2/spm12r)
- All code is found at https://github.com/muschellij2/talks/tree/master/fmri_task_processing


# `spm12r` Worked Example<br>Disclaimer: there is no universal fMRI pipeline<br> Many options are **specific** to this data analysis


## Data required for analysis

- One anatomical T1-weighted scan: `anat.nii.gz`
- One 4D fMRI task-related scan: `fmri.nii.gz`.  
- Information on design:
    - onsets/duration of stimuli
- Order of slices (which was first slice) 
    - or time slice measured (in ms).
- Repetition time (TR) from DICOM header/scanner/tech.

## Worked example: DICOM conversion

- DICOM to NIfTI: Convert the data using `dcm2niix` (https://github.com/rordenlab/dcm2niix).
- `dcm2niir`: https://github.com/muschellij2/dcm2niir
- `divest`: https://github.com/jonclayden/divest



## Download the data 

https://figshare.com/articles/SFO-example/5442298


```r
url = paste0("https://ndownloader.figshare.com/articles/",
             "5442298/versions/1")
# download a temporary zip file
zipfile = tempfile(fileext = ".zip")
res = httr::GET(url, write_disk(path = zipfile))

####### unzip file code (not shown) ###########

out_files = c("anat.nii.gz", "fmri.nii.gz")
```




## fMRI information


```r
fmri_filename = "fmri.nii.gz"
tr = 1.8 # seconds
hdr = neurobase::check_nifti_header(fmri_filename) # nifti header
(nslices = oro.nifti::nsli(hdr))
```

```
[1] 60
```

```r
(n_time_points = oro.nifti::ntim(hdr))
```

```
[1] 280
```

```r
hdr
```

```
NIfTI-1 format
  Type            : nifti
  Data Type       : 16 (FLOAT32)
  Bits per Pixel  : 32
  Slice Code      : 0 (Unknown)
  Intent Code     : 0 (None)
  Qform Code      : 1 (Scanner_Anat)
  Sform Code      : 2 (Aligned_Anat)
  Dimension       : 96 x 96 x 60 x 280
  Pixel Dimension : 2.25 x 2.25 x 2.55 x 1.8
  Voxel Units     : mm
  Time Units      : sec
```


# Explore the Raw Data: <br><br>http://bit.ly/neuroshiny

## Explore the Data

The first part of any preprocessing pipeline should be to use exploratory techniques to investigate  detect possible problems and artifacts.

fMRI data often contain transient spike artifacts and slow drift.

An exploratory technique such as principal components analysis (PCA) can be used.

(courtesy of Martin Lindquist)


## Types of Registration
<div style="font-size: 20pt;">

- Rigid-body registration (linear) - 6 degrees of freedom (dof)

<img src="rollpitchyaw.png" style="width: 50%; display: block; margin: auto;">
<div style="font-size: 8pt">
Image taken from [http://cnl.web.arizona.edu/imageprops.htm](http://cnl.web.arizona.edu/imageprops.htm)
</div>

- Pitch - Think of nodding ("yes")
- Yaw - Think of shaking head ("no") 
- Roll - Think of shoulder shrugging ("I don't know")
- x – left/right, y – forward/backward, z – jump up/down 

</div>

## Rigid Registration: The Math

<div style="font-size: 20pt;">

For a voxel $v$, the rigid transformation can be written as:

$$T_{\rm rigid}(v) = Rv + t$$
where $R =$
\small
$$\left[\begin{array}{ccc} \cos\beta\cos\gamma& \cos\alpha\sin\gamma + \sin\alpha\sin\beta\cos\gamma & \sin\alpha\sin\gamma - \cos\alpha\sin\beta\cos\gamma \\
-\cos\beta\sin\gamma & \cos\alpha\cos\gamma - \sin\alpha\sin\beta\sin\gamma & \sin\alpha\cos\gamma + \cos\alpha\sin\beta\sin\gamma \\
\sin\beta & -\sin\alpha\cos\beta & \cos\alpha\cos\beta \end{array}\right]$$
\normalsize

- 6 degrees of freedom
- $3$ associated with the translation vector: $t=(t_x, t_y, t_z)$
- $3$ associated with the rotation parameters: $\theta=(\alpha, \beta,\gamma)$. 

</div>


## Image Realignment: within-fMRI registration

<center>
<img src="average.png" style="width:40%; margin: auto;" alt="flow"> 
</center>

## Image Realignment 



```r
realigned = spm12_realign(filename = fmri_filename,
  time_points = seq(n_time_points),
  quality = 0.98, separation = 3,
  register_to = "mean", est_interp = "bspline4", reslice_interp = "bspline4")
# reading in the mean image
mean_img = realigned[["mean"]]
mean_nifti = readnii(mean_img)
realigned$outfiles
```

```
[1] "rfmri.nii"
```

```r
realigned$mat
```

```
[1] "fmri.mat"
```





## Image Realignment 

<center>
<img src="realign.png" style="width:40%; margin: auto;" alt="flow"> 
</center>

## Plotting the realignment parameters



<img src="index_files/figure-html/rp_plot-1.png" style="display: block; margin: auto;" />

## Slice timing correction - temporal alignment

<center>
<img src="st.png" style="width:90%; margin: auto;" alt="flow"> 
</center>

<div style="font-size: 20pt;">
(courtesy of Martin Lindquist)
</div>



## Slice timing correction - temporal alignment

- slice order: descending, dual-coil (different for ascending or interleaved)
- Need to know this from DICOM/design

```r
slice_order = c(1740, 1680, 1620, 1560, 1500, 1440, 1380, 
  1320, 1260, 1200, 1140, 1080, 1020, 960, 900, 840, 780, 
  720, 660, 600, 540, 480, 420, 360, 300, 240, 180, 120, 60, 
  0, 1740, 1680, 1620, 1560, 1500, 1440, 1380, 1320, 1260, 
  1200, 1140, 1080, 1020, 960, 900, 840, 780, 720, 660, 600, 
  540, 480, 420, 360, 300, 240, 180, 120, 60, 0)
ref_slice = 900
ta = 0 # since slice_order in ms
```


## What does this order mean?<br> <img src="slice_timing_slow.gif" style="width: 50%; display: block; margin: auto;">


## Data needed for slice timing correction

- Repetition time (from `hdr`)
- Number of time points and slices (from `hdr`)
- Slice order + need the reference slice (`ref_slice`), 
- Time between the first and the last slice within one scan (`ta`).  `ta = 0` if you give slice order in seconds/milliseconds.





## Slice timing correction - temporal alignment


```r
aimg = spm12_slice_timing(filename = realigned$outfiles,
  nslices = nslices,  
  tr = tr, slice_order = slice_order,
  time_points = seq(n_time_points),
  ta = ta, # since slice order given in ms 
  ref_slice = ref_slice, 
  prefix = "a")
print(aimg$outfile)
```




## After lice timing correction

<center>
<img src="st2.png" style="width:90%; margin: auto;" alt="flow"> 
</center>

<div style="font-size: 20pt;">
(courtesy of Martin Lindquist)
</div>


## T1 Coregistration to Mean fMRI

We then perform the coregistration of the mean image (fixed) and T1 (moving):

<center>
<img src="coreg.png" style="width:70%; margin: auto;" alt="flow"> 
</center>

## T1 Coregistration to Mean fMRI

Coregistration is estimated using `spm12_coregister_estimate`:


```r
t1_fname = "anat.nii.gz"
coreg = spm12_coregister_estimate(
  fixed = mean_img,
  moving = t1_fname, 
  cost_fun = "nmi")
coreg$outfile
```


```
[1] "anat.nii"
```

## Output file was the same: nothing happened!

- `spm12_coregister_estimate` - estimates coregistration (transforms the header)
- `spm12_coregister_reslice` - reslices the image to the same voxel dimensions (should probably be coregistered already using `estimate`)
- `spm12_coregister` - estimates and reslices

- Estimate the transformation, but do segmentation on native T1 space (better resolution)


## Anatomical MRI Segmentation 

Here we segment the image into 6 different regions, where the regions are gray matter, white matter, cerebrospinal fluid (CSF), bone, soft tissue, and the background.  


```r
seg = spm12_segment(
  filename = coreg$outfile,
  set_origin = FALSE, 
  bias_corrected = TRUE, native = TRUE,
  unmodulated = TRUE, modulated = TRUE, affine = "mni",
  sampling_distance = 1.5)
```


## Anatomical MRI Segmentation 

- `native` - native space segmentations
- `modulated` - adjusted segmentations to constrain tissue-class volumes
- `unmodulated` - unadjusted 
- `bias_corrected` - save bias-field corrected image
- `set_origin` - should AC/PC alignment be done (no because we just coregistered)

## Anatomical MRI Segmentation 

<img src="index_files/figure-html/hard_seg-1.png" style="display: block; margin: auto;" />

## Anatomical MRI Segmentation: CSF/WM/GM 

<img src="index_files/figure-html/hard_seg2-1.png" style="display: block; margin: auto;" />


## Spatial normalization to MNI

- My brain is not the same size/shape as your brain
- Want to look across subjects spatially
- Spatial normalization allows us to transform the data, stretching and scaling the data (nonlinearly) to a standard brain.
- MNI (Montreal Neurological Institute) is the most commonly used (ICBM MNI152 of some sort, http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009).  

## Spatial normalization to MNI

Affine + Non-linear transform (invertible)

<center>
<img src="nonlin.png" style="width:80%; margin: auto;" alt="flow"> 
</center>

## Spatial normalization to MNI: already done

The segmentation was done by warping the T1 to the MNI template and that transform/deformation in the segmentation output:


```r
seg$deformation
```


```
[1] "y_anat.nii"
```

## Applying spatial normalization: fMRI

We apply the deformation to the fMRI data using `spm12_normalize_write`.  


```r
bounding_box = matrix(
    c(-78, -112, -70, 
      78, 76, 85), nrow = 2, 
    byrow = TRUE)
norm = spm12_normalize_write(
  deformation = seg$deformation,
  other.files = aimg$outfile, #corrected fMRI
  bounding_box = bounding_box,
  interp = "bspline5")
```

## Applying spatial normalization: fMRI

![](index_files/figure-html/norm_show-1.png)<!-- -->

## Applying spatial normalization: Corrected T1



```r
anat_norm = spm12_normalize_write(
  deformation = seg$deformation, other.files = seg$bias_corrected,
  bounding_box = bounding_box, interp = "bspline5",
  voxel_size = c(1, 1, 1))
anat_norm$outfiles
```


```
[1] "wmanat.nii"
```

![](index_files/figure-html/anat_norm_show-1.png)<!-- -->



## Applying spatial normalization: T1, but 2x2x2


```r
anat_norm2x2x2 = spm12_normalize_write(
  deformation = seg$deformation, other.files = seg$bias_corrected,
  bounding_box = bounding_box, interp = "bspline5",
  voxel_size = c(2, 2, 2)) # note the resolution!!!
```

![](index_files/figure-html/anat_norm2x2x2_show-1.png)<!-- -->


## Spatial smoothing using a Gaussian

- Spatial smoothing should signal to noise depending on the size of activation

- Typically, the amount of smoothing is chosen a priori

- Usually global smoothing (same amount at each voxel), but can be adaptive (`adimpro` pacakge)


## Spatial smoothing using a Gaussian

- Specified using the full-width half max (FWHM) for the Gaussian smoother (not $\sigma$):  $FWHM = \sigma \sqrt{8 \log(2)}$

<img src="3dgauss.png" style="width: 50%; display: block; margin: auto;">

<div style="font-size: 20pt;">
From https://en.wikipedia.org/wiki/Gaussian_function#/media/File:Gaussian_2d.svg
</div>


## Spatial smoothing using a Gaussian


```r
smooth_norm = spm12_smooth(
  norm$outfiles[[1]], fwhm = 5, 
  prefix = "s5")
```

![](index_files/figure-html/smooth_norm_show-1.png)<!-- -->

## First level modeling: Single-subject model

In many applications, that smoothed data you will use for post-processing and analysis.  Motion correction has usually been applied above, but some realign this again.

## Conditions of the experiment (block design)

- need the onset/duration of conditions (in seconds or scans):




## First level modeling: single-subject model

- Conditions are convolved with the Hemodynamic Response Function (HRF)

<img src="spm_hrf.png" style="width: 50%; display: block; margin: auto;">

<div style="font-size: 20pt;">
https://en.wikibooks.org/wiki/SPM/Haemodynamic_Response_Function#/media/File:SPM_hemodynamic_response_function.png
</div>





## Estimate first level model

- General linear model (GLM) (not **Generalized**)
- `regressor_mat` - motion parameters and other "confounders" (not convolved with HRF) 
- `condition_list` - conditions are convolved



## Model outputs

- [Good Cheat Sheet](http://www.bobspunt.com/resources/teaching/single-subject-analysis/spmdoc/SPMdotMAT.pdf)
- beta coefficient maps of regressors and contrasts 


```r
betas = list.files(pattern = "beta.*[.]nii")
print(betas)
```

```
[1] "beta_0001.nii" "beta_0002.nii" "beta_0003.nii" "beta_0004.nii"
[5] "beta_0005.nii" "beta_0006.nii" "beta_0007.nii" "beta_0008.nii"
[9] "beta_0009.nii"
```


```r
print(first_model$spmmat)
```

```
[1] "SPM.mat"
```

- `SPM.mat` - model specification

## Contrast Manager - Creating Contrasts

- can make T-statistic of F statistic maps
- `weights` indicate which coefficients


```r
contrasts = list(
  list(name = "LeftHand", weights = c(1, rep(0, 7)),
    replicate = "none", type = "T" ),
  list(name = "RightHand", weights = c(0, 1, rep(0, 6)),
       replicate = "none", type = "T"), 
  list(name = "AllEffects",
       weights = rbind(
         c(1, rep(0, 7)),
         c(0, 1, rep(0, 6))
       ), replicate = "none", type = "F")   )
```

## Contrast Manager - Creating Contrasts


```r
contrast_res = spm12_contrast_manager(spm = first_model$spmmat,
  delete_existing = TRUE, contrast_list = contrasts)
```


```r
cons = list.files(pattern = "con.*[.]nii")
print(cons)
```

```
[1] "con_0001.nii" "con_0002.nii"
```

```r
stats = list.files(pattern = "spm(T|F).*[.]nii")
print(stats)
```

```
[1] "spmF_0003.nii" "spmT_0001.nii" "spmT_0002.nii"
```


## There is no universal fMRI pipeline

- Each step has inherent drawback and limitation (spatial resolution, artifact smoothing, etc.)
- A few different pipelines should be tested.
    - Not necessarily all combinations, but change the "knobs" a bit
- Similar to sensitivity analysis

## References

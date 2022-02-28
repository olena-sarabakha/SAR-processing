# SAR-processing
Processing codes for SAR applications.

This repository contains basic codes for some blocks of SAR processing in MATLAB language. In particular:
- matched_filter_azimuth.m is a function for matched filtering in azimuth dimension; it is used to compress the spatial chirp in the azimuth dimension.
- patch_adapt_filter.m is a function that applies the matched filter in azimuth splitting the input matrix into patches, whose dimension is determined according to the depth of focus criterion;
- resampling.m is a function that allows to resample the input matrix according to a multiplier chosen by the user;
- psf_range is a function that creates a plot of the Point Spread Function (PSF) of selected cell in the range dimension, in order to verify if the -4dB width is compliant with the expected resolution and check the symmetry of the side lobes and their level compared to the main lobe;
- psf_azimuth.m is a function that creates a plot of the Point Spread Function (PSF) of selected cell in the azimuth dimension, in order to verify if the -4dB width is compliant with the expected resolution and check the symmetry of the side lobes and their level compared to the main lobe;
- RCM_correction.m is a function that corrects the Range Cell Migration (RCM) affecting the SAR image; this effect must be corrected before focusing in azimuth;
- PD_autofocus_time.m is an autofocus function that applies the Phase Difference (PD) algorithm in the time domain to correct de-focalization due to uncertainties in the knowledge of the acquisition geometry and/or the velocity of the SAR platform; this autofocus algorithm isn't suitable in case multiple scatterers are present in the same range cell;
- PD_autofocus_freq.m is an autofocus function that applies the Phase Difference (PD) algorithm in the frequency domain to correct de-focalization due to uncertainties in the knowledge of the acquisition geometry and/or the velocity of the SAR platform; this autofocus algorithm is suitable in case multiple scatterers are presesnt in the same range cell;
- MR_autofocus.m is an autofocus function that applies the Multilook Registration algorithm to correct de-focalization due to uncertainties in the knowledge of the acquisition geometry and/or the velocity of the SAR platform;
- CO.autofocus.m is a function that applies the Contrast Optimization algorithm to correct de-focalization due to uncertainties in the knowledge of the acquisition geometry and/or the velocity of the SAR platform;
- multilook_filter_sliding.m is a function that multilooks the input image in order to reduce the speckle effect; the algorithm consists in a sliding window that averages the pixels within the window producing one averaged pixel, thus reduces the resolution of the image as a side effect; this method treats in the same way the pixels whether they belong to a homogeneous region or an edge region; the dimension of the square window is chosen by the user and determines the loss in resolution;
- frost_filter_sliding.m is a function that multilooks the input image in order to reduce the speckle effect; this method is equivalent to the multilook filter in homogeneous regions, while the regions where dominant scatterers and/or edges are present are treated in such a way that the scatterer/edge is preserved.

This repository will be periodically updated in order to improve the uploaded codes and include some new content.

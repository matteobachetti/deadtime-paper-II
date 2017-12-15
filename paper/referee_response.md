Report from Referee:

> At this point, I cannot recommend for publication "No Time for Dead Time - Use of the Fourier Amplitude Differences to Normalize Dead Time-Affected Periodograms" by Bachetti and Hupenkothen. The basic reason is that the paper as written is difficult enough to follow, contains typographical errors, misstatements and mislabelings, has poor figures, to the point that I cannot easily create my own code that reproduces the results exactly as presented in this work. An algorithmic paper like this, to be truly useful, should be sufficiently self-contained and explanatory so as to be easily reproducible.

We thank the referee for the frank report and their constructive suggestions, which allowed us to improve the paper significantly. We acknowledge that the paper was not as clear as it should have been. However. We have extensively modified the manuscript to hopefully demonstrate this. We would like to start this response with the description of a possible critical point: in order to accommodate a thorough discussion of the caveats of the method, we decided to take out the real data example from the paper (Cyg X-1). We realize that this choice might raise some concerns from the referee, as traditionally a paper might be considered incomplete without real data. However, for an algorithm paper like this one, we find that the addition of real data detracts from the main focus of the paper, which is the description of a data analysis method that works on any instrument with multiple identical detectors. The presented method will still work practically unmodified with IXPE, which will probably have a dead time only slightly shorter than NuSTAR, and of course current instruments like Astrosat. The NuSTAR detectors and dead time within NuSTAR are well enough understood to simulate very precisely, thus aside from the reassurance of seeing a real periodogram from a real BH candidate, a single example of real data periodogram would not do a better job at showing the method compared to the thousands of realizations of simulated NuSTAR power spectra we produced, that allowed us to test the efficacy of the methods in different parameter regimes.


We have also added a new discussion about several caveats of the method. Most prominantly, we discuss the effect a difference in flux between the two detectors would have, and also include the effects that a difference in count rates between the two detectors would have on the results, up to a difference higher than any reasonable sitaution found in NuSTAR data.
 
We are setting up a "Gallery" on Github showing off the FAD correction working on real data that we will link in the text, with a unique identifier on Zenodo or similar services.

Overall, we have edited the paper for clarity, expanded where necessary and completely redone the figures to improve the ease with which they can be understood.

> I realize that the authors have provided a link to a Jupyter notebook to reproduce the results and figures, and it was sufficient to highlight some of the misstatements made in this paper. However, I have some issues with the adequacy of this, as discussed below.
> In section 2.1, the discussion on creating the lightcurves is insufficient. Yes, it is reasonable to create lightcurves with an existing set of packages, but since these packages are in active development a version number for the packages should be given, or a git hash for the versions used.
> How do I know whether any version I run now or in the future will not be different?

We agree that it is good practice to refer to a specific version of the codes, so we added the version numbers used.

> More importantly "stingray.simulate.Simulator" class is *not* a description of the methods used. A few sentences describing the mechanics behind each method used is necessary. Likewise, going to the Jupyter notebook and simply seeing "from stingray" and "from hendrics" is not a description. Unless I start following a potentially long chain of imports and start reading code, "stingray" and "hendrics" are black boxes. (In general, the existing Jupyter notebook desperately needs comments in the code.) This can be rectified with short descriptions of the algorithms.

We have added some description of the algorithms, and included appropriate references where some of them were already described in detail (the simulation process is basically the same to Bachetti+2015: model periodogram -> Timmer & Koening method -> light curve -> acceptance-rejection method -> clean event lists -> dead time filter -> final event lists). We hope the new version will be easier to follow.

> This "clipped" nature of the descriptions carries throughout the paper.
> Section 2.2 the fact that the data is coming from NuSTAR is only in the section header; this is not acceptable.
> The "heasarc_pipelines" packages is only referenced as "ASCL in prep." OK, that's Astrophysics Source Code Library, but who are the authors? I tried a search at ASCL, and that package did not come up.

As we said in the opening, we eliminated the part about NuSTAR data analysis from the paper.

> Section 3, x and y are defined, without explicitly saying that you are specifically considering binned lightcurves of identical length representing the exact same times.

We have slightly expanded this explanation to make sure we state our assumptions about the nature of the light curves.

> Also in this section the expected Fourier amplitude is given as N_phot P(\nu)/4, and only *later* in the paragraph is it explicitly stated as the Leahy normalization.

We have rewritten this section in part for clarity, and have made sure to specify the normalization of the periodogram where it is first mentioned.

> (As a side note, I see figures on the Jupyter notebook with the Leahy power of white noise being plotted with values of order 100, which is not what I expect for the Leahy normalization.)

We do not fully grasp this point, and we cannot find a place in the notebook where Leahy power from _white noise_ goes to 100, in which case it might just be a misprint of the axis label. The white noise distributions we plot show precisely how the data points follow the correct chi squared distribution. Moreover, the points where the Leahy power goes to 100 in notebook 2 - if the reviewer is referring to that plot - is heavily distorted by dead time and amplified (around 300 Hz is the maximum of the deadtime positive distortion).

> Later in Section 3, A_xj is written when you mean B_xj.

We thank the reviewer for pointing out this typo.

> Section 4, Figure 3 is said when you mean Figure 2, and

We have fixed this, and have significantly changed the figures to improve clarity.

> Figure 2a is said to be Fourier Amplitudes vs. Fourier Amplitude Differences, when you actually mean the standard deviation of said quantities. (The Jupyter notebook was useful for figuring this out.)

In this case, we made a confusing mistake in both the plot and the notebook, and we thank the referee for bringing this to our attention, because it needs clarification. The only point missing in the text is that the _mean absolute values_ of the FAD and amplitudes in each frequency bin are used, and this is the reason for the apparent lack of negative values. The description with the standard deviation is somehow equivalent, with actually the same factor between the quantities. We tried both approaches in the Jupyter notebook, and we chose the plots inconsistently in the paper, but it is not what is being used in the actual algorithm (see the implementation in the file fad_correction.py in the notebooks/ directory in the repository). We have edited the text and figure to make them consistent.

> Figure 1 caption, you say "Real Fourier amplitudes". Did you mean "actual", or did you mean you're plotting only the "Real" component as opposed to the imaginary component of the amplitudes?

Real part. We have clarified this in the new version.

> The are a number of places that you simply write "spectral", and from context it is sometime obvious that you meant power spectral variability, and aren't referring to the spectrum, but with a quick read it could easily be misinterpreted as spectral energy distribution in some cases.

We changed all references to "spectral" in order to clarify that we are talking about variability only.

> In the description of applying the correction, the multiplication factor to obtain a Leahy-normalized periodogram is only applicable for a specific normalization of the FFT and whether one has applied the forward or reverse transform. I could go on.

We are not sure that we understand the point the referee makes here, but we would like to point out that the multiplication factor should work with _any_ normalization of the PDS. The way we calculate the FAD, it is done so that the white noise power spectral values in Leahy normalization flatten to ~two, but the same multiplicative factor would correct other normalizations as well: the unit is always P x multiplicative factors, in Leahy or the various rms normalizations, and the correction is just a new multiplicative factor to be included in calculating the normalization.

> Any one of these issues in isolation could be considered trivial and easily forgiven if fixed. The fact that these issues are numerous and persist throughout the entire paper greatly detracts from the clarity of the presentation.
> The Figures are for the most part very poor. The labels are too small. In some cases (Fig 3d) the ranges are poorly chosen. The order and grouping is confusing. (The left of Figure 1 is not even described in the text, the right of Fig. 1 is more naturally paired with Fig. 2a, and might have made your meaning more clear. Fig. 2b and 2c should probably separated out into their own Figure separate from Fig. 2a.) The figures are too small, and only really readable on a computer screen if zoomed in.

We changed most of the figures in order to select the most effective ones and give them more space.

> In general, there is too much explanation and description crammed into the figure caption, instead of the text. I assume this was in some attempt to obtain a "letter" length. Given rapid electronic publication these days, there often isn't a strong motivation for choosing a letter over a journal article, especially if that means sacrificing both clarity and depth.

While we understand how that perception could arise, given that our paper was otherwise not very clear, but we would like to assure the referee that this was not our intention. Our aim was for our figures to work on their own, without the need to look in the main text for essential information to understand them. We agree that the complexity of the figures meant they required too much text to describe them; our new, simpler figures should be easier to read and have shorter captions attached to them. We would like to point out that the text in the figures counts in the number of words allowed for letters, so there would be no advantage from this point of view.

> And we haven't even gotten to discussing the algorithm. The one major question I have is, what happens if the two lightcurve extractions are in fact not identical? How precisely do the two detectors have to behave identically? This is mentioned at the end of Section 5, but not adequately discussed.

In the new version, we discuss more extensively the details of this, and other caveats, in a subsection on its own. Basically, the only power spectral measurement that is heavily affected by the count rate difference in the two detectors is the single-channel periodogram. Neither the periodogram obtained by the total-intensity light curve (the sum of the two channels) or the cospectrum are very sensitive to the difference of flux in the two channels. In our simulations the parameter indicating the ratio of flux in the two channels was free to change from 0.5 to 2 (with 1 being the same flux in the two channels), and we find basically indistinguishable distributions if the fluxes in the channels are within ~30% of one another (It can be seen in the simulations in the notebooks). The simulations in figure 5 include channel ratios between 0.8 to 1.3.

> Again, this paper should be adequate within itself for me to be able to quickly script a code to reproduce the results. It is not. Once the paper is rewritten, I would be willing to try this experiment again, and see if I agree with the results that are presented.

We thank the referee for their attempt to reproduce our results. We hope that the new version of this papers simplifies this exercise significantly, and we would welcome any further suggestions.

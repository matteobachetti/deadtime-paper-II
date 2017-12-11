Report from Referee:

> At this point, I cannot recommend for publication "No Time for Dead Time - Use of the Fourier Amplitude Differences to Normalize Dead Time-Affected Periodograms" by Bachetti and Hupenkothen. The basic reason is that the paper as written is difficult enough to follow, contains typographical errors, misstatements and mislabelings, has poor figures, to the point that I cannot easily create my own code that reproduces the results exactly as presented in this work. An algorithmic paper like this, to be truly useful, should be sufficiently self-contained and explanatory so as to be easily reproducible.

We thank the referee for the frank report. We do acknowledge that the paper was not as clear as it should have been. However, we do think that a poor description should not be confused with a poor result. The extensive modifications we have done to the manuscript will hopefully demonstrate this. We would like to start this response with the description with a possible critical point: in order to accommodate a thorough discussion of the caveats of the method, we decided to take out the real data example from the paper (Cyg X-1). We realize that this choice might raise some concerns from the referee, as traditionally a paper is considered incomplete without real data. However, we find that the addition of real data does not really add anything to the main focus of the paper, which is the description of a data analysis method that works on any instrument with multiple identical detectors. Besides the reassuring fact of seeing a real periodogram from a real BH candidate, a single example of real data spectrum would not do a better job at showing the method, because it would still be a single example, subject to the usual doubts (is it a cherry-picked example?). The simulation environment we use in the paper mimicks NuSTAR sufficiently well, uses thousands of realizations of simulated NuSTAR power spectra - not just one - and the new additions in the discussion about the difference of flux between the two channels also includes effects of count rate differences that are higher than expected from any reasonable situation found in NuSTAR data (where bright sources are usually dominating the observation and are far from the detector edges). In any case, we are setting up a web page showing off the FAD correction working on real data that we will link in the text.

> I realize that the authors have provided a link to a Jupyter notebook to reproduce the results and figures, and it was sufficient to highlight some of the misstatements made in this paper. However, I have some issues with the adequacy of this, as discussed below.
> In section 2.1, the discussion on creating the lightcurves is insufficient. Yes, it is reasonable to create lightcurves with an existing set of packages, but since these packages are in active development a version number for the packages should be given, or a git hash for the versions used.
> How do I know whether any version I run now or in the future will not be different?

We do agree that it is good practice to refer to a specific version of the codes, so we added the version numbers used.

> More importantly "stingray.simulate.Simulator" class is *not* a description of the methods used. A few sentences describing the mechanics behind each method used is necessary. Likewise, going to the Jupyter notebook and simply seeing "from stingray" and "from hendrics" is not a description. Unless I start following a potentially long chain of imports and start reading code, "stingray" and "hendrics" are black boxes. (In general, the existing Jupyter notebook desperately needs comments in the code.) This can be rectified with short descriptions of the algorithms.

We have added some description of the algorithms, in some cases by referencing Bachetti+15 where some of them were already described in detail (the simulation process is basically the same: model periodogram -> Timmer & Koening method -> light curve -> acceptance-rejection method -> event lists -> dead time filter). We hope the new version will be easier to follow.

> This "clipped" nature of the descriptions carries throughout the paper.
> Section 2.2 the fact that the data is coming from NuSTAR is only in the section header; this is not acceptable.
> The "heasarc_pipelines" packages is only referenced as "ASCL in prep." OK, that's Astrophysics Source Code Library, but who are the authors? I tried a search at ASCL, and that package did not come up.

As we said in the opening, we eliminated the part about NuSTAR data analysis from the paper.

> Section 3, x and y are defined, without explicitly saying that you are specifically considering binned lightcurves of identical length representing the exact same times.
We added this information

> Also in this section the expected Fourier amplitude is given as N_phot P(\nu)/4, and only *later* in the paragraph is it explicitly stated as the Leahy normalization.

The section was basically rewritten.

> (As a side note, I see figures on the Jupyter notebook with the Leahy power of white noise being plotted with values of order 100, which is not what I expect for the Leahy normalization.)

We do not fully grasp this point, and we cannot find a place in the notebook where Leahy power from _white noise_ goes to 100, in which case it might just be a misprint of the axis label. The white noise distributions we plot show precisely how the data points follow the correct chi squared distribution. Moreover, the points where the Leahy power goes to 100 in notebook 2 - if the reviewer is referring to that plot - is heavily distorted by dead time and amplified (around 300 Hz is the maximum of the deadtime positive distortion).

> Later in Section 3, A_xj is written when you mean B_xj.

Rewritten.

> Section 4, Figure 3 is said when you mean Figure 2, and

Corrected, and figures are changed anyway.

> Figure 2a is said to be Fourier Amplitudes vs. Fourier Amplitude Differences, when you actually mean the standard deviation of said quantities. (The Jupyter notebook was useful for figuring this out.)

No, we made a confusing mistake in the plot and the notebook and this needs to be clarified. The only point missing in the text is that the _absolute values_ of the FAD and amplitudes are used, and this is the reason for the apparent lack of negative values. The description with the standard deviation is somehow equivalent, with actually the same factor between the quantities (we tried both approaches in the Jupyter notebook, and we chose the plots inconsistently), but it's not what we use in the actual algorithm (see the implementation in the file fad_correction.py in the notebooks/ directory in the repository).

> Figure 1 caption, you say "Real Fourier amplitudes". Did you mean "actual", or did you mean you're plotting only the "Real" component as opposed to the imaginary component of the amplitudes?

Real part. We tried to be clearer in the new version.

> The are a number of places that you simply write "spectral", and from context it is sometime obvious that you meant power spectral variability, and aren't referring to the spectrum, but with a quick read it could easily be misinterpreted as spectral energy distribution in some cases.

We changed all references to "spectral" in order to clarify that we are talking about variability.

> In the description of applying the correction, the multiplication factor to obtain a Leahy-normalized periodogram is only applicable for a specific normalization of the FFT and whether one has applied the forward or reverse transform. I could go on.

This point is not very clear to us, but anyway the multiplication factor should work with _any_ normalization of the PDS. The way we calculate the FAD, it is done so that the white noise power spectral values in Leahy normalization flatten to ~two, but the same multiplicative factor would correct other normalizations as well (because the unit is always P x multiplicative factors, in Leahy or the various rms normalizations, and this is just a new multiplicative factor)

> Any one of these issues in isolation could be considered trivial and easily forgiven if fixed. The fact that these issues are numerous and persist throughout the entire paper greatly detracts from the clarity of the presentation.
> The Figures are for the most part very poor. The labels are too small. In some cases (Fig 3d) the ranges are poorly chosen. The order and grouping is confusing. (The left of Figure 1 is not even described in the text, the right of Fig. 1 is more naturally paired with Fig. 2a, and might have made your meaning more clear. Fig. 2b and 2c should probably separated out into their own Figure separate from Fig. 2a.) The figures are too small, and only really readable on a computer screen if zoomed in.

We changed most of the figures in order to select the most effective ones and give them more space.

> In general, there is too much explanation and description crammed into the figure caption, instead of the text. I assume this was in some attempt to obtain a "letter" length. Given rapid electronic publication these days, there often isn't a strong motivation for choosing a letter over a journal article, especially if that means sacrificing both clarity and depth.

We respectfully disagree with this. The text in the figures counts exactly as much as the one in the main text for the ApJL rules, so there is would be no advantage whatsoever from the point of view of obtaining a letter length. We just like the figures to work on their own, without the need to look in the main text for essential information needed to understand them. The figures were quite (probably too) complicated and needed some text to describe them properly. In any case, now most of them have changed.

> And we haven't even gotten to discussing the algorithm. The one major question I have is, what happens if the two lightcurve extractions are in fact not identical? How precisely do the two detectors have to behave identically? This is mentioned at the end of Section 5, but not adequately discussed.

In the new version, we discuss more extensively the details of this, and other caveats, in a subsection on its own. Basically, the only power spectral measurement that is heavily affected by the count rate difference in the two detectors is the single-channel periodogram. Neither the periodogram obtained by the total-intensity light curve (the sum of the two channels) or the cospectrum are very sensitive to the difference of flux in the two channel. In our simulations the parameter indicating the ratio of flux in the two channels was free to change from 0.5 to 2 (with 1 being the same flux in the two channels), and we find basically indistinguishable distributions if the fluxes in the channels are within ~30% of one another (It can be seen in the simulations in the notebooks). The simulations in figure 5 include channel ratios between 0.8 to 1.3.

> I very well likely might have further comments and questions about the results. As I tried to write my own code to check the figures and results, I had some differences from the figures. This could potentially be my fault. However, after running across poorly worded description after poorly worded description, wondering how Fourier amplitudes vs. Fourier amplitude differences could be so narrowly distributed and not have negative values (i.e., Fig. 2a) until I discovered that the plots were instead standard deviation vs. standard deviation, I simply gave up.

The code to test the algorithm conveniently is (and was) contained in the fad_correction.py file in the notebooks/ directory. There, it was clear that we weren't using the standard deviation, but admittedly, we didn't write a reference to this file in the text.

> Again, this paper should be adequate within itself for me to be able to quickly script a code to reproduce the results. It is not. Once the paper is rewritten, I would be willing to try this experiment again, and see if I agree with the results that are presented.

We hope that the new version will make the referee better appreciate the algorithm.

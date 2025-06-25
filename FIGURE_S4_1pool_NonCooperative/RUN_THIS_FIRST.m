%
%  Main script to produce all figures for the model in the folder title
%  Compute Likelihood Profile and plot results
%
%  Note: No. of Ca2+ bindings sites is set by nCaSites in ComputeSetup.m 
%
%  Make sure you download the MATLAB Parallel Toolbox
%  Likelihood profile may take 1-2 hours on a 10-core CPU

Likelihood_Profile;  % Code will automatically load from COMMON folder
Figure_5B_S4B;       % Code in this folder
Figure_5C;           % Code will automatically load from COMMON folder

%
%  Main script to produce all figures for the model in the folder title
%  Compute Likelihood Profile and plot results
%
%  Note: No. of Ca2+ bindings sites is set by nCaSites in ComputeSetup.m 
%
%  Make sure you download the MATLAB Parallel Toolbox
%  Likelihood profile may take 1-2 hours on a 10-core CPU

clear;
ComputeSetup;        % Compute [Ca2+] and all other settings
Likelihood_Profile;  % Generate parameter profile (code in COMMON folder): takes some time!
Figure_5B_S4B;       % Plot the final parameter profile
Figure_5C;           % Plot model Cm vs pulse duration (code in COMMON folder)

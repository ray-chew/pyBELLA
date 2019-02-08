%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% make_plots
% Produces plots from RKLM low Mach fluid dynamics code hdf output
% for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
% dynamics"
%
% Developed by R. Klein, FU Berlin, -2019
% Modified by T. Benacchio, Politecnico di Milano, 2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Acoustic Wave
plots_acouwave('rho', 'high');
plots_acouwave('u', 'high');

plots_acouwave('rho', 'low');
plots_acouwave('u', 'low');


% Internal Waves
plots_internalwave('dY', 'NH');
plots_internalwave('dY', 'H');
plots_internalwave('dY', 'P');  

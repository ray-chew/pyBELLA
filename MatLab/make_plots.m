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

% Vortex
plots_vortex('rho')
plots_vortex('p2_n')

% Straka density current

plots_straka('dY');
plots_straka_1dcuts('dY')

% Internal Waves
plots_internalwave('dY', 'NH');
plots_internalwave('dY', 'H');
plots_internalwave('dY', 'P');  

% Baldauf-Brdar internal wave
plots_baldauf('dY');
plots_baldauf('u');
plots_baldauf('v');


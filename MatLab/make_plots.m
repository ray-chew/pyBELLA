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
plots_vortex('rho', 192)
plots_vortex('p2_n', 192)

plots_vortex_convergence('rho')
plots_vortex_convergence('p2_n')

% Straka density current

% plots_straka('dY', 400);
% plots_straka('dY', 200);
% plots_straka('dY', 100);
 plots_straka('dY', 50);
% plots_straka('dY', 25);
 plots_straka_1dcuts('dY')

% Internal Waves
plots_internalwave('dY', 'NH');
% plots_internalwave('dY', 'H');
% plots_internalwave('dY', 'H_psinc');
% plots_internalwave('dY', 'H_hyd');
plots_internalwave_cmp('H');
% plots_internalwave('dY', 'P');  
% plots_internalwave('dY', 'P_psinc');  
% plots_internalwave('dY', 'P_hyd');  
plots_internalwave_cmp('P');

% Baldauf-Brdar internal wave
plots_baldauf('dT');
plots_baldauf('u');
plots_baldauf('v');


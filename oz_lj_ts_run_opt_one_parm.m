% This program creates the Ornstein-Zernike equation example 
% with the Verlet bridge function with one free parameter, phi.  
% A free parameter phi is obtained using fminsearch.m 
%
% An interparticle interaction: the Lennard-Jones potential
% A method: a Picard iteration technique
%  
% Dr. Tsogbayar Tsednee, California State University, Northridge
% Date: Aug 26, 2017
%
function [] = oz_lj_ts_run_opt_one_parm
format long
clear;
clc;
% A free parameter phi for for the modified Verlet bridge function, which will be optimized
% Another free parameter alpha is fixed in this code
phi = 0.5; % initial guees for free parameter; you may change it

%%% finding an optimal parameters using fminsearch.m.   
[phi_opt, dpv_dpc] = fminsearch(@consistency_lj,phi);
%
[phi_opt, dpv_dpc] % phi_opt is optimal value for phi 
                   % (dpv_dpc)^(1/2) gives a criteria for a consistency
% Results 
% [phi_opt, dpv_dpc] = 0.691699218750001   0.000000000669509
%
%%%
function [dpv_dpc] = consistency_lj(phi)
%%%
Nr = 1*4096.;           % number of grid points ( = power of 2); you may change it 
L = 32.;                % length of interval ( = 8, 16, 32, etc.); you may change it 
itermax = 4000;         % max number of iteration; you may change it
tol = 10^(-12);         % tolerance of convergence; you may change it 
%
alf = 0.75000;          % damping parameter; you may change it
alpha = 1.0;            % free parameter of Verlet bridge function 
%
rho_red = 0.400;        % density in reduced units; you may change it
delta_rho_red = 0.0001; % delta\rho in numerical derivative calculation; you may change it  
T_red = 2.75;           % temperature in reduced units; you may change it
%
% Ornstein_Zernik equation solver with the Verlet bridge function
[comp_eq_hr, comp_eq_cr, eq_of_st, int_en] = oz_lj_ts(L,Nr,itermax,tol,alf,phi,alpha,rho_red,T_red);
%
[comp_eq_hr_p1, comp_eq_cr_p1, eq_of_st_p1, int_en_p1] = ...
    oz_lj_ts(L,Nr,itermax,tol,alf,phi,alpha,rho_red+delta_rho_red,T_red);
%
[comp_eq_hr_m1, comp_eq_cr_m1, eq_of_st_m1, int_en_m1] = ...
    oz_lj_ts(L,Nr,itermax,tol,alf,phi,alpha,rho_red-delta_rho_red,T_red);
%%%
% calculate dp/drho with finite difference (2-point scheme)
dp_drho = eq_of_st + rho_red*(eq_of_st_p1 - eq_of_st_m1)/(2.*delta_rho_red);
%%%

% Output ---
%    * comp_eq_hr is isothermal compressibility with h(r) 
%    * comp_eq_cr is isothermal compressibility with c(r)
%    * eq_of_st is an equation of state, (beta*p/rho)
%    * int_en is internal energy
%    * dp_drho is d(beta*p)/(d rho)
%
Output = [comp_eq_hr, comp_eq_cr, dp_drho, eq_of_st, int_en];
%
dpv_dpc = (comp_eq_hr - dp_drho).^2;
%
return




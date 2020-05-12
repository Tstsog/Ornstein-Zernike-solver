function [] = OZ_solver_one_component_plasma_hnc
%-----------------------------------------------------------------
% Solves the Orsntein-Zernike equation for one-component plasma 
% in hypernetted chain approximation
% An interparticle interaction: the Coulomb potential
% A method: a Picard iteration technique
% Uses: function lf=lsint(f)
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Reference: Kin-Chue Ng, J. Chem. Phys. 61, 2680 (1974);
% Date: February 09, 2019
% ----------------------------------------------------------------
clear;
clc;
format short
Nr = 4096; % number of grid point; you may change it
L = 128.0; % grid parameter; you may change it  
dr = L/(Nr+1); % grid parameters
itermax = 4000; tol = 10^(-14); % max number of iteration and tolerance; you may change it
%
alf = 0.750 ; % convergence acceleration parameter; you may change it
%%%
rho_red = 3/(4*pi); % number density parameter for plasma
%%%
gam = 200.0; % plasma parameter; you may change it
%%%
aa = 1.00; % parameter in Coulomb potential 
%
r = zeros(Nr,1); u_sr = zeros(Nr,1); u_lr = zeros(Nr,1); c_sr = zeros(Nr,1);
for i = 1:Nr
    r(i) = i*dr;
%    
    u_lr(i) = (gam/r(i))*erf(aa*r(i));         % long-range part of Coulomb potential
    u_sr(i) = (gam/r(i))*(1.0 - erf(aa*r(i))); % short-range part of Coulomb potential
%    
    c_sr(i) = exp( -u_sr(i) ) - 1. ;
end
%
Nk = Nr; j = (1:1:Nk)'; ii = (1:1:Nr)'; dk = pi/L; k = (j*dk);
%
hk_old = 0.;
%
ck_sr = 4.*pi*(dr^3*(Nr+1)./pi).*lsint(c_sr.*ii)./j;
uk_lr = (4.*pi*gam./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
%%%
ck = ck_sr - uk_lr;
hk = ck./(1.-rho_red.*ck); % OZ solving in k-space
%
hk_new = hk;
hk_n = alf*hk_new + (1.-alf)*hk_old;
% compute c(r) and h(r) 
c_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(ck_sr.*j);
hr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(hk_n.*j);
%%%
gam_sr = hr - c_sr; % gamma(r) = h(r) - c(r), indirect correlation function 
%%%
%Iflag = 0.;  % convergence flag, 1 indicates convergence
%
for iter = 2:itermax
%    
    iter
    %%% iteration 2
    gam_sr_old = gam_sr;
    %%%
    c_sr_old = exp(-u_sr + gam_sr_old) - gam_sr_old - 1. ;     % c(r) in hnc approximation
    ck_sr = 4.*pi.*(dr^3*(Nr+1)./pi).*lsint(c_sr_old.*ii)./j;  % c(k)
%%%
    ck = ck_sr - uk_lr;  
    hk = ck./(1.-rho_red.*ck); % OZ solving in k-space
%    
    hk_old = hk_n;
    hk_new = hk;
    hk_n = alf*hk_new + (1.-alf)*hk_old;  
%
    % compute c(r) and h(r) 
    c_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(ck_sr.*j);
    hr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(hk_n.*j);
%%% 
    gam_sr = hr - c_sr;
%%%    
    gam_sr_new = gam_sr;
%    
    if (rms(gam_sr_new - gam_sr_old) < tol);
          Iflag = 1.
        break 
    end
%
    gam_sr = gam_sr_new;
%%%
end

%%%
gr = exp( -u_sr + gam_sr); % radial distribution function
%
% average potential energy: beta*U/N = (rho/2) * int (gam/r)*h(r)*4*pi*r*r*dr
ave_pot_energy_sm = 0.;
for i = 1:Nr
    ave_pot_energy_sm = ave_pot_energy_sm + (3/2)*gam*(gr(i)-1.0)*r(i)*dr;
end
%
[gam, ave_pot_energy_sm ] % 
% Results:
% [gam, ave_pot_energy_sm ]
%   20.0000  -16.5353       vs  20  -16.53771   is from J. Chem. Phys. 61, 2680 (1974)
%   60.0000  -51.5900       vs  60  -51.59735   is from J. Chem. Phys. 61, 2680 (1974)
%  100.0000  -86.9612       vs 100  -86.97342   is from J. Chem. Phys. 61, 2680 (1974)
%  200.0000 -175.8320       vs 200 -175.85637   is from J. Chem. Phys. 61, 2680 (1974)
%%%
%plot(r, gr, '-b') % plot for pair correlation function
%axis([0. 12. 0. 2.5 ])
return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function lsint.m is taken from C.T. Kelley, B. M. Pettitt, 
%                      J. Chmp. Phys. \textbf{197}, 491 (2004)
% LSINT
% Fast sine transform with MATLAB's FFT.
%
function lf=lsint(f)
n=length(f);
ft=-fft([0,f']',2*n+2);
lf=imag(ft(2:n+1));
%
return
end
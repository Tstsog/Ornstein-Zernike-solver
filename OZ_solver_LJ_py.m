function [] = OZ_solver_LJ_py
%-----------------------------------------------------------------
% Solves the Orsntein-Zernike equation for simple liquid
% in Percus?Yevick (PY) approximation
% An interparticle interaction: the Lennard-Jones potential
% A method: a Picard iteration technique
% Uses: function lf=lsint(f)
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Reference: R. J. Baxter, JCP, 47, 4855 (1967). 
% Date: August 09, 2017
% ----------------------------------------------------------------
clear; clc;
format short
Nr = 4096; % number of grid point; you may change it
L = 128.0; % grid parameter; you may change it  
dr = L/(Nr+1); % grid parameters
itermax = 4000; tol = 10^(-14); % max number of iteration and tolerance; you may change it
%
alf = 0.550 ; % convergence acceleration parameter; you may change it
%%%
rho_red = 0.500; % The reduced density; you may change it.
T_red = 2.74;    % The reduced temperature; you may change it.
%
r = zeros(Nr,1); 
for i = 1:Nr
    r(i) = i*dr;
end
%
% Calculatite the Lennard-Jones potentail and c(r) for the first iteration
[c,U] = lj_c(r,T_red);
%%%
hk_old = 0.;
%
Nk = Nr; j = (1:1:Nk)'; ii = (1:1:Nr)'; dk = pi/L; k = (j*dk)';
%
ck = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c.*ii)./j;
hk = ck./(1.-rho_red.*ck);
%
hk_new = hk;
hk_n = alf*hk_new + (1.-alf)*hk_old;
%
cr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(ck.*j);
hr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(hk_n.*j);
G = hr - cr;
%
%Iflag = 0.;  % converfence flag, 1 indicates convergence
%
for iter = 2:itermax
%    
    iter
    %%% iter 2
    G_old = G;
    Br = PY_bridge(G_old); % Bridge function
    c_old = exp(-U + G_old + Br) - G_old - 1.; % c(r), closure
    c = (c_old);
%
    ck = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c.*ii)./j;
    hk = ck./(1.-rho_red.*ck);
%
    hk_old = hk_n;
    hk_new = hk;
    hk_n = alf*hk_new + (1.-alf)*hk_old;
%
    cr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(ck.*j);
    hr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(hk_n.*j);
    G_new = hr - cr;
%
    cor = rms(G_new - G_old);
%    
    if (cor < tol)        
        Iflag = 1.
    break
    end
    %
    G = G_new;
    %   
end


%%% Calculate isothermal compressibility,  (K^c_T)^(-1), with h(r)
sm = 0.;
for i = 1:Nr
    sm = sm + 4*pi*rho_red*hr(i)*r(i)*r(i)*dr;
end
comp_eq_hr = 1 + sm;
comp_eq_hr = comp_eq_hr^(-1.);

%%% Calculate isothermal compressibility,  K^c_T, with c(r)
sm = 0.;
for i = 1:Nr
    sm = sm + 4*pi*rho_red*cr(i)*r(i)*r(i)*dr;   
end
comp_eq_cr = 1 - sm;
%

%%% Calculation radial distribution function
gr = hr+1.; % radial distribution function
for i = 1:Nr
    if (gr(i) <= 10^(-6))
        gr(i) = 0.;
    else
        gr(i) = gr(i);
    end
end

%%% Calculate an equation of state, (beta*p/rho)
sm = 0.;
for i = 1:Nr
    sm = sm + 16*pi*rho_red*(1./T_red)*(-2.*(1./r(i))^10 + (1./r(i))^4)*gr(i)*dr;
end
eq_of_st = 1. - sm;

%%% Calculate an internal energy U ()
sm = 0.;
for i = 1:Nr
    sm = sm + (16*pi/3)*rho_red*(1./T_red)*((1./r(i))^10 - (1./r(i))^4)*gr(i)*dr;
end
int_en = 1. + sm;


%%% pressure calculation from compressibility route in PY closure: 
%%% analytical formula from JCP, 47, 4855 (1967)
press_com_py_sm1 = 0.;
for i = 1:Nr
    press_com_py_sm1  = press_com_py_sm1 + (1./2).*rho_red*cr(i)*(gr(i) - 2.) * 4.*pi*r(i)*r(i)*dr;
end
press_com_py_sm2 = 0.;
for i = 1:Nk
    press_com_py_sm2  = press_com_py_sm2 + (1./(2*pi)^(3))*( ck(i) + (1./rho_red)*log(abs(1.-rho_red.*ck(i))) ) * 4.*pi*k(i)*k(i) * dk;
end
press_com_py_over_rho = 1.0 + press_com_py_sm1 + press_com_py_sm2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
[rho_red, T_red, comp_eq_hr, comp_eq_cr, eq_of_st, int_en, press_com_py_over_rho ]
%
% 0.5000    2.7400    2.9768    2.9768    1.5215    0.2871    1.4375 

%%% plot pair correlation function
plot(r, gr, '-b') % plot for pair correlation function
axis([0. 8. 0. 2.5 ])
set(gca,'FontSize',18)
xlabel('r') % ,'fontsize',16
ylabel('g(r)' ,'Rotation', 1)

return
end

%%%
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

%%%
function Br = PY_bridge(G_old)
% PY bridge function
    Br = (log(1.+G_old) - G_old);
return 
end

%%%
function [c,U] = lj_c(r,T_red)
% Lennard-Jones potential and c(r) for the first iteration 
Nr = length(r); U = zeros(Nr,1); c = zeros(Nr,1);
for i = 1:Nr
    U(i) = (4./T_red)*( (1./r(i))^12 - (1./r(i))^6 );
    c(i) = exp(-U(i)) - 1.;
end
return
end


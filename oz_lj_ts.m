function [comp_eq_hr, comp_eq_cr, eq_of_st, int_en] = oz_lj_ts(L,Nr,itermax,tol,alf,phi,alpha,rho_red,T_red)
%function oz_lj_ts(L,Nr,itermax,tol,alf,phi,alpha,rho_red,T_red)
% This program creates the Ornstein-Zernike equation example 
% with the Verlet bridge function
% An interparticle interaction: the Lennard-Jones potential
% A method: a Picard iteration technique
%  
% Dr. Tsogbayar Tsednee, California State University, Northridge
% Date: Aug 26, 2017
% 
% A matlab code give a stable solution with follwoing values of parameters
% Nr = power of 2, (32, 64, 128, 256, 512); 
% L = 8, 16, or 32
%
% The function lsint.m is taken from 
% C.T. Kelley, B. M. Pettitt,J. Chmp. Phys. \textbf{197}, 491 (2004)
%
%
%  Usage:
%
%    oz_lj_ts(L,Nr,itermax,tol,alf,phi,alpha,rho_red,T_red)
%
%    where
%
%    * L is a length of interval ( = 8, 16, 32, etc.,)
%    * Nr is number of grid points ( = power of 2)
%    * itermax is max number of iteration
%    * tol is tolerance of convergence 
%    * alf is dumping paramater
%    * phi is free parameter of Verlet bridge function
%    * alpha is free parameter of Verlet bridge function
%    * rho_red is density in reduced units
%    * T_red is temperature in reduced units
%
% 
%  Output:
%
%    [comp_eq_hr, comp_eq_cr, eq_of_st, int_en]
%
%    where
%    * comp_eq_hr is isothermal compressibility with h(r) 
%    * comp_eq_cr is isothermal compressibility with c(r)
%    * eq_of_st is an equation of state, (beta*p/rho)
%    * int_en is internal energy
%
%
%%%
dr = L/(Nr+1);
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
Nk = Nr; j = (1:1:Nk)'; ii = (1:1:Nr)'; %dk = pi/L; %k = (j*dk)';
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
%    iter
    %%% iter 2
    G_old = G;
    Br = Verlet_bridge(r,T_red,G_old,phi,alpha); % Bridge function
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
    G = G_new;
%
    if (cor < tol)        
        Iflag = 1.
    break
    end
    
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
    if (gr(i) <= 0.000001)
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
eq_of_st = 1 - sm; 

%%% Calculate an internal energy U ()
sm = 0.;
for i = 1:Nr
    sm = sm + (16*pi/3)*rho_red*(1./T_red)*((1./r(i))^10 - (1./r(i))^4)*gr(i)*dr;
end
int_en = 1. + sm;

%[iter, phi, alpha, rho_red, T_red, comp_eq_hr, comp_eq_cr, eq_of_st, int_en ]

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
function Br = Verlet_bridge(r,T_red,G_old,phi,alpha)
% Verlet bridge function
Nr = length(r); Ua = zeros(Nr,1);
%%%
for i = 1:Nr
    if (r(i) < (2.)^(1/6.)) % 1.122462048309373
        Ua(i) = -(1./T_red);
    else
        Ua(i) = (4./T_red)*( (1./r(i))^12 - (1./r(i))^6 );
    end
end
    Br = -(0.5.*phi.*(G_old-Ua).^2)./(1.+alpha.*(G_old-Ua));
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




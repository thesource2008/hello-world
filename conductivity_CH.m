%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate conductivity of CH as a 
% function of density (gm/cm3) and Temperature 
% according to Phys. Plasmas 23, 042704 (2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = conductivity_CH(rho,T)
kb = 1.602176565E-12;                       %Since T is in eV, erg/eV
[Zbar, Z2bar, Zeff] = ionization_CH(rho,T);

ACH = 6.5095;                                  %atomic mass of CH
amu = 1.6605e-24;                           %atomic mass unit grams
ACH = ACH * amu;
ni = rho / ACH;                             %number density 1/cm3
alphaB = 5.292e-9;                          %Bohr radius in cm
r0 = (1 / alphaB) * (3 / (4 * pi * ni))^(1/3);
rs = (3 / (4 * pi * ni))^(1/3);
e = 4.8032045E-10;                          %esu
Gammai = Z2bar*e^2 / (rs * kb * T);       

ne = ni * Zbar;
me = 9.1e-28;                               %electron mass grams
hplanck = 6.626e-27;                        %Plancks constant erg S
hbar = hplanck / (2 * pi);
Tfermi = hbar^2/(2*me*kb) * (3*pi*pi*ne)^(2/3);
thetae = T / Tfermi;

gamma0 = -0.482;
gamma1a = -0.150; gamma1b = 0.275;
gamma2 = 0.193;
gamma3 = 8.364e-3;
gamma4 = -5.287e-3;
gamma5 = -3.191e-4;
gamma6 = 2.666e-5;

sigma1a = 1.00; sigma1b = 1.20;
sigma2 = -0.225;
sigma3 = -4.652e-3;
sigma4 = 3.805e-3;
sigma5 = -7.643e-5;
sigma6 = -1.391e-5;

lnGammai = log(Gammai);
lnThetae = log(thetae);

lnLambda_a = exp(gamma0 ...
            + gamma1a * (lnGammai)^1 + sigma1a*(lnThetae)^1 ...
            + gamma2 * (lnGammai)^2 + sigma2*(lnThetae)^2 ...
            + gamma3 * (lnGammai)^3 + sigma3*(lnThetae)^3 ...
            + gamma4 * (lnGammai)^4 + sigma4*(lnThetae)^4 ...
            + gamma5 * (lnGammai)^5 + sigma5*(lnThetae)^5 ...
            + gamma6 * (lnGammai)^6 + sigma6*(lnThetae)^6);
lnLambda_b = exp(gamma0 ...
            + gamma1b * (lnGammai)^1 + sigma1b*(lnThetae)^1 ...
            + gamma2 * (lnGammai)^2 + sigma2*(lnThetae)^2 ...
            + gamma3 * (lnGammai)^3 + sigma3*(lnThetae)^3 ...
            + gamma4 * (lnGammai)^4 + sigma4*(lnThetae)^4 ...
            + gamma5 * (lnGammai)^5 + sigma5*(lnThetae)^5 ...
            + gamma6 * (lnGammai)^6 + sigma6*(lnThetae)^6);
lnLambda = min(lnLambda_a, lnLambda_b);

factor1 = 20 * (2/pi)^(3/2) * kb^(7/2) * T^(5/2) / (sqrt(me) * e^4 * Zeff);
factor2 = 0.095 * (Zeff + 0.24) / (1 + 0.24 * Zeff);
factor3 = 1 / lnLambda;
kappa = factor1 * factor2 * factor3;
%To convert to SI units
kappa = kappa * 1e-9/1.16;

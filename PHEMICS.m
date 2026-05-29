function [K,w]=PHEMICS(K,fstart,B0,ne,dd_cold_ion,dd_warm_ion)
% PHEMICS : Parallel Hot-plasma EMIC Solver for Multi-ion Anisotropic Plasmas
%
% Purpose
% -------
% This function solves the linear kinetic dispersion relation for
% parallel-propagating EMIC waves in a multi-ion anisotropic plasma.
%
% Input
% -----
% K :
%     Normalized wave number:
%
%         K = k c / (omega_pe * sqrt(eps))
%
%     where eps = me/mp.
%
% fstart :
%     Initial guess of the complex normalized wave frequency:
%
%         fstart = omega / Omega_p = omega_r/Omega_p + i gamma/Omega_p
%
%     Note: fstart is dimensionless
%
% B0 :
%     Background magnetic field, T.
%
% ne :
%     Electron number density, m^-3.
%     1 cm^-3 = 1e6 m^-3.
%
% dd_cold_ion :
%     Cold ion parameters:
%
%         dd_cold_ion = [ns/ne;
%                        qs/qe;
%                        ms/mp]
%
%     Each column represents one cold ion species.
%     This input can be empty, [].
%
% dd_warm_ion :
%     Warm/hot ion parameters:
%
%         dd_warm_ion = [ns/ne;
%                        qs/qe;
%                        ms/mp;
%                        Tperp/Tpara - 1;
%                        Tpara_keV]
%
%     Each column represents one finite-temperature ion species.
%     This input cannot be empty.
%
% Output
% ------
% K :
%     Normalized wave number.
%
% w :
%     Complex normalized frequency:
%
%         w = omega/Omega_p
%
%     real(w) gives omega_r/Omega_p.
%     imag(w) gives gamma/Omega_p.
%
%
% Author: Xue Z.X.
% Code development assistance: ChatGPT (GPT-5.5 Thinking), OpenAI
% Final implementation, verification, and scientific interpretation by Xue Z.X.
% Date: 2026/05/29

%% -----------------------------
%  Input checks
% ------------------------------
if nargin < 6
    error('PHEMICS requires six inputs: K, xstart, B0, ne, dd_cold_ion, dd_warm_ion.');
end

if isempty(dd_warm_ion)
    error('dd_warm_ion cannot be empty. At least one warm/hot ion species is required.');
end

if size(dd_warm_ion,1) ~= 5
    error('dd_warm_ion must have 5 rows: [ns/ne; qs/qe; ms/mp; A; Tpara_keV].');
end

if ~isempty(dd_cold_ion) && size(dd_cold_ion,1) ~= 3
    error('dd_cold_ion must have 3 rows: [ns/ne; qs/qe; ms/mp].');
end

if dd_warm_ion(5,1) <= 0
    error('The reference warm-ion parallel temperature dd_warm_ion(5,1) must be positive.');
end

%% Physical constants
qe   = 1.602176634e-19;
me   = 9.1093837015e-31;
mp   = 1.67262192369e-27;
mu0  = 4*pi*1e-7;
eps0 = 8.8541878128e-12;

eps_m = me/mp;   % epsilon = me/mp

%% Normalized parameters
omega_pe = sqrt(ne*qe^2/(eps0*me));
Omega_e  = qe*B0/me;

p.eps = eps_m;
p.f   = omega_pe/Omega_e;

%% ions
p.cold = [];
if ~isempty(dd_cold_ion)
    for ii=1:size(dd_cold_ion,2)
        p.cold(ii).eta  = dd_cold_ion(1,ii);
        p.cold(ii).Z    = dd_cold_ion(2,ii);
        p.cold(ii).M    = dd_cold_ion(3,ii);
    end
else
    p.cold(1).eta  = 0;
    p.cold(1).Z    = 1;
    p.cold(1).M    = 1;
end

p.warm=[];
for ii=1:size(dd_warm_ion,2)
    p.warm(ii).eta     = dd_warm_ion(1,ii);
    p.warm(ii).Z       = dd_warm_ion(2,ii);
    p.warm(ii).M       = dd_warm_ion(3,ii);
    p.warm(ii).A       = dd_warm_ion(4,ii);
    p.warm(ii).T_ratio = dd_warm_ion(5,ii)/dd_warm_ion(5,1);
end

%% beta_parallel_p
Tp_J = dd_warm_ion(5,1) * 1e3 * qe;
Nhp  = p.warm(1).eta * ne;
p.beta_ref = 2*mu0*Nhp*Tp_J/B0^2;    % beta_parallel of reference hot proton
p.eta_ref  = p.warm(1).eta;    % density ratio of reference hot proton
%% Solve H-band and He-band branches
w = solveBranch(K, fstart, p);











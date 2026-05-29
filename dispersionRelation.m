function F = dispersionRelation(x, K, p)
% General parallel L-mode kinetic dispersion relation
%
% x = omega / Omega_p
% K = k c / (omega_pe * sqrt(eps))
%
% p.cold(j): pure cold ion species
% p.warm(j): finite-temperature kinetic ion species
%
% Required fields:
% p.eps
% p.f
% p.beta_ref
% p.eta_ref
% p.cold(j).eta
% p.cold(j).Z
% p.cold(j).M
% p.warm(j).eta
% p.warm(j).Z
% p.warm(j).M
% p.warm(j).A
% p.warm(j).T_ratio

eps_m = p.eps;
f     = p.f;

beta_ref = p.beta_ref;
eta_ref  = p.eta_ref;

%% -----------------------------
% Left-hand side
% ------------------------------
LHS = eps_m * (x/f)^2;

%% -----------------------------
% Electron term: cold electron
% ------------------------------
cold_e = x / (eps_m*x + 1.0);

%% -----------------------------
% Cold ion terms
% For species s:
% alpha_s = M_s / Z_s
% Omega_s / Omega_p = 1 / alpha_s
%
% cold contribution:
% Z_s * eta_s * x / (alpha_s*x - 1)
% ------------------------------
cold_ions = 0;

if isfield(p, 'cold') && ~isempty(p.cold)

    for js = 1:length(p.cold)

        eta_s = p.cold(js).eta;
        Z_s   = p.cold(js).Z;
        M_s   = p.cold(js).M;

        alpha_s = M_s / Z_s;

        cold_ions = cold_ions + ...
            Z_s * eta_s * x / (alpha_s*x - 1.0);

    end

end

%% -----------------------------
% Warm / hot kinetic ion terms
%
% zeta_s = (x - 1/alpha_s) /
%          [K * sqrt(beta_ref/eta_ref) * sqrt(Ts/Tref/Ms)]
%
% contribution:
%
% - Z_s*eta_s/alpha_s *
%   [ A_s - { (A_s+1)(1-alpha_s*x)-1 }
%     * zeta_s Z(zeta_s)/(alpha_s*x-1) ]
%
% To improve numerical stability near cyclotron resonance:
%
% zeta_s Z(zeta_s)/(alpha_s*x-1)
% = Z(zeta_s)/(alpha_s*K*C_s)
%
% where
% C_s = sqrt(beta_ref/eta_ref) * sqrt(Ts/Tref/Ms)
% ------------------------------
warm_ions = 0;

if isfield(p, 'warm') && ~isempty(p.warm)

    Cref = sqrt(beta_ref / eta_ref);

    for js = 1:length(p.warm)

        eta_s   = p.warm(js).eta;
        Z_s     = p.warm(js).Z;
        M_s     = p.warm(js).M;
        A_s     = p.warm(js).A;
        T_ratio = p.warm(js).T_ratio;

        if T_ratio <= 0
            error('Finite-temperature species must have T_ratio > 0. Put T=0 species in p.cold instead.');
        end

        alpha_s = M_s / Z_s;

        Cs = Cref * sqrt(T_ratio / M_s);

        zeta_s = (x - 1.0/alpha_s) / (K * Cs);

        Zeta_s = plasmaDispersionFunction(zeta_s);

        % Numerically stable form:
        % zeta_s*Z(zeta_s)/(alpha_s*x - 1)
        ratio_s = Zeta_s / (alpha_s * K * Cs);

        Gs = A_s - ...
            (((A_s + 1.0)*(1.0 - alpha_s*x) - 1.0) * ratio_s);

        warm_ions = warm_ions - Z_s * eta_s / alpha_s * Gs;

    end

end

%% -----------------------------
% Full dispersion relation
% ------------------------------
RHS = K^2 + cold_e + cold_ions + warm_ions;

F = RHS - LHS;

end
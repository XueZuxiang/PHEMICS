function xs = solveBranch(Kvec, x0, p, opts_in)
% SOLVEBRANCH  Track one complex EMIC dispersion branch over a K grid.
%
%   xs = SOLVEBRANCH(Kvec, x0, p)
%
%   This function solves the complex dispersion relation
%
%       F(x, K) = 0
%
%   for each normalized wave number K in Kvec, where
%
%       x = omega / Omega_p = omega_r / Omega_p + i*gamma / Omega_p.
%
%   The real part of x gives the normalized wave frequency, and the
%   imaginary part gives the normalized linear growth or damping rate.
%   A positive imaginary part corresponds to wave growth, whereas a
%   negative imaginary part corresponds to damping.
%
%
%   Inputs
%   ------
%   Kvec : array
%       Normalized wave-number grid:
%
%           K = k c / (omega_pe * sqrt(epsilon)).
%
%   x0 : complex scalar
%       Initial guess for the complex root at the first valid K point.
%       Typical choices are:
%
%           H band  : x0 = 0.5 ~ 0.8 + small imaginary part
%           He band : x0 = 0.05 ~ 0.2 + small imaginary part
%
%   p : structure
%       Plasma parameter structure used by dispersionRelation.m.
%       It usually contains electron parameters, reference beta, and the
%       cold / finite-temperature ion species definitions.
%
%   opts_in : structure
%             options for fsolve.
%
%   Output
%   ------
%   xs : complex array
%       Complex root x(K) for each K in Kvec.
%
%       real(xs) gives omega_r / Omega_p.
%       imag(xs) gives gamma / Omega_p.
%
%       If fsolve fails at a given K, or if the solution is invalid, the
%       corresponding output is set to:
%
%           NaN + i*NaN.
%
%
%   Numerical Method
%   ----------------
%   The complex equation is rewritten as two real equations:
%
%       real(F) = 0,
%       imag(F) = 0.
%
%   These two equations are solved using fsolve. The solution at the
%   previous successful K point is used as the initial guess for the next
%   K point. This continuation strategy helps track the same dispersion
%   branch across the K grid.
%
%
%   Convergence Control
%   -------------------
%   The fsolve options set the internal stopping criteria, while resTol is
%   an additional user-defined residual check:
%
%       norm(Fval) < resTol.
%
%   The initial guess y0 is updated only after a successful convergence.
%   If a K point fails, the code stores NaN and continues to the next K
%   without updating y0.
%
%
%   Required User Functions
%   -----------------------
%   This function requires:
%
%       realImagEquation.m
%       dispersionRelation.m
%       plasmaDispersionFunction.m
%
%
%   Notes
%   -----
%   1. Very small K can be difficult because |zeta| can become large.
%   2. Roots close to cold-ion cyclotron poles may be numerically delicate.
%   3. Branch tracking is sensitive to the initial guess x0.
%   4. If many NaN values appear, try using a better initial guess, a finer
%      K grid, or a bounded solver such as lsqnonlin.
%
%   Author: Xue Z.X.
%   Code development assistance: ChatGPT (GPT-5.5 Thinking), OpenAI
%   Final implementation, verification, and scientific interpretation by Xue Z.X.
%
%   See also FSOLVE, OPTIMOPTIONS.

if nargin==3
    opts = optimoptions('fsolve', ...
        'Display', 'off', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'MaxIterations', 500, ...
        'MaxFunctionEvaluations', 5000);
elseif nargin==4
    opts = opts_in;
end

xs = nan(size(Kvec)) + 1i*nan(size(Kvec));

y0 = [real(x0); imag(x0)];

resTol = 1e-7;

for ik = 1:length(Kvec)

    K = Kvec(ik);

    if ~isfinite(K) || K <= 0
        xs(ik) = nan + 1i*nan;
        warning('Invalid K = %.4f, set NaN.', K);
        continue;
    end

    fun = @(y) realImagEquation(y, K, p);

    try
        [y, Fval, exitflag] = fsolve(fun, y0, opts);
    catch
        y = [nan; nan];
        Fval = [nan; nan];
        exitflag = 0;
    end

    x = y(1) + 1i*y(2);

    if exitflag > 0 && ...
            all(isfinite(y)) && ...
            all(isfinite(Fval)) && ...
            norm(Fval) < resTol && ...
            isfinite(real(x)) && isfinite(imag(x))

        xs(ik) = x;
        y0 = y;   % only update initial guess after successful convergence

    else
        xs(ik) = nan + 1i*nan;
        warning('fsolve did not converge at K = %.4f, set NaN.', K);
    end

end
end
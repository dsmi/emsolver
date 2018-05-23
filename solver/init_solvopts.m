function sopt = init_solvopts(freq, mqs)
% sopt = init_solvopts(freq, mqs)
%
% Returns the structure with the solver options initialized with default
% values.
%

if ~exist('mqs')
    mqs = 0;
end

sopt.freq = freq;
sopt.eps = eps0;
sopt.mu = mu0;
sopt.conductivity = 5.8e7; % Copper

sopt.hf = 0;
sopt.mqs = mqs;

% Order of the quadrature used when evaluating integrals
% over triangle of the testing edge.
sopt.mqo0 = 1;
sopt.mqo1 = 1;

% Wavenumber outside the conductors
if mqs
   k = 0;
else
   k = sopt.freq * sqrt(sopt.mu * sopt.eps);
end

% Permittivity of the conductors
eps_c = (sopt.eps - j*sopt.conductivity/sopt.freq);

% Wavenumber with conductivity absorbed
k1 = sopt.freq * sqrt(sopt.mu * eps_c);

% Number of the quadrature points used when evaluating integrals
% over triangle of the source edge.
sopt.nqn0 = 8;
sopt.nqn1 = 16;

sopt.fintg_fp_0  = @(r, robs)integ_fp(k, r, robs, sopt.nqn0);
sopt.fintg_p_0   = @(r, robs)integ_p(k, r, robs, sopt.nqn0);
sopt.fintg_fxg_0 = @(r, robs)integ_fxg(k, r, robs, sopt.nqn0);

sopt.fintg_fp_1  = @(r, robs)integ_fp(k1, r, robs, sopt.nqn1);
sopt.fintg_p_1   = @(r, robs)integ_p(k1, r, robs, sopt.nqn1);
sopt.fintg_fxg_1 = @(r, robs)integ_fxg(k1, r, robs, sopt.nqn1);

function r = r_TVI(R0,distr_profile,mean_profile,alpha_profile,Tmax,dt)

% This function computes the real-time growth rate r for a single-type
% model with time-varying infectivity profile (TVI), in the cases of a
% constant duration of infectious period (a rectangular-shaped infectivity 
% profile), an exponentially-decaying infectivity profile, or a 
% gamma-shaped infectivity profile. 
% 
% Input:
%   - R0: the basic reproduction number
%   - distr_profile: the shape of the infectivity profile, which can be:
%       - 0: constant duration (rectangular), with specified mean
%       - 1: exponentially decaying, with specified mean
%       - 2: gamma-shaped, with specified mean and shape parameter
%   - mean_profile: mean of the shape omega, often called the generation
%   time
%   - alpha_profile: shape parameter of the gamma-shaped infectivity
%   profile (only used for distr_profile == 2)
%   - Tmax: largest time considered, after which the profile is trncated
%   - dt: time step used
%   - tol: a tolerance, to check that we have not truncated the infectivity
%   profile too soon.

% Reference: Pellis, Ferguson and Fraser (2010), Journal of Mathematical
% Biology

t = 0:dt:Tmax; % index tt
lt = length(t);

omega = zeros(1,lt); % generation time distribution (discretised in time)

if distr_profile == 0 % Constant duration
    dt_mean = round(2*mean_profile/dt); % Total duration is twise the mean, then approximated to be a multiple of the time step
    omega(1:dt_mean) = 1/dt; % Total area under omega integrates to 1
elseif distr_profile == 1 % Exponential shape
    omega =  exp(-t/mean_profile) / mean_profile;
elseif distr_profile == 2
    lambda = 1/mean_profile;
    omega = lambda^alpha_profile * t.^(alpha_profile-1) .* exp(-lambda*t) / gamma(alpha_profile);
else
    error('Aborting: ill-specified infectivity profile')
end
if mean_profile > Tmax
    error('Aborting: Tmax is too small')
end

fun = @(x) ( R0 * ( dt * trapz(omega.*exp(-x*t)) ) - 1 );

r0 = fun(0);
rInf = fun(100);
if ( sign(r0) == sign(rInf) ) % Ensure there is a zero, before I start searching for one
    r = NaN;
else
    r = fzero( fun, [ 0 100 ] );
end


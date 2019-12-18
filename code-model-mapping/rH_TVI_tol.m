function r = rH_TVI_tol(AvCasesGen,Rg,distr_profile,mean_profile,alpha_profile,Tmax,dt,tol)

% This function computes the approximate real-time growth rate r for a 
% model with households and a time-varying infectivity profile (TVI) with
% a constant duration of infectious period (a rectangular-shaped
% infectivity profile), an exponentially-decaying infectivity profile, or a
% gamma-shaped infectivity profile. To be precise, whether the computation 
% of the real-time growth rate is exact or not depends on how the average
% number of cases in each generation, AvCasesGen, has been pre-computed. A
% tolerance parameter checks the infectivity profile is not truncted too
% early.
% 
% Input:
%   - AvCasesGen: the average number of cases in each generation of a
%   within-household epidemic
%   - Rg: the average number of global infectious contacts a singla case
%   makes (the R0 for global contacts, i.e. outisde the household)
%   - distr_profile: the shape of the infectivity profile, which can be:
%       - 0: constant duration (rectangular), with specified mean
%       - 1: exponentially decaying, with specified mean
%       - 2: gamma-shaped, with specified mean and shape parameter
%   - mean_profile: mean of the shape omega, often called the generation
%   time
%   - alpha_profile: shape parameter of the gamma-shaped infectivity
%   profile (only used for distr_profile == 2)
%   - Tmax: largest tiem considered, after which the profile is trncated
%   - dt: time step used
%   - tol: a tolerance, to check that we have not truncated the infectivity
%   profile too soon.

% Reference: Pellis, Ferguson and Fraser (2010), Journal of Mathematical
% Biology

t = 0:dt:Tmax; % index tt
lt = length(t); 
nh = length(AvCasesGen); % number of generations

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

% figure(1);
% clf;
% plot(t,omega,'r');

% omega
% integral = dt*trapz(omega)
% last = omega(lt)

if ( omega(lt) > tol ) || ( 1 - dt*trapz(omega) > tol )
    warning(['Careful! Tmax may be too early! omega(',num2str(Tmax),') = ',num2str(omega(lt))]);
end

afs = sum(AvCasesGen) - 1; % Average final size (afs), excluding the initial case

% I calculate 2 approximate values: the best approximation is not really a 
% true lower bound, as you can get values lower than that, but it usually 
% is, while the upper bound is a real upper bound but is often very high.
fun_lower = @(x) ( Rg * sum( AvCasesGen .* ( dt * trapz(omega.*exp(-x*t)) ).^(1:nh) ) - 1 );
fun_upper = @(x) ( Rg * ( ( dt * trapz(omega.*exp(-x*t)) ) + afs * ( dt * trapz(omega.*exp(-x*t)) )^2 ) - 1 );

r_lower0 = fun_lower(0);
r_lowerInf = fun_lower(100);
if ( sign(r_lower0) == sign(r_lowerInf) ) % Ensure there is a zero, before I start searching for one
    r_lower = NaN;
else
    r_lower = fzero( fun_lower, [ 0 100 ] );
end
r_upper0 = fun_upper(0);
r_upperInf = fun_upper(100);
if (sign(r_upper0) == sign(r_upperInf) )
    r_upper = NaN;
else
    r_upper = fzero( fun_upper, [ 0 100 ] );
end

r = [ r_lower, r_upper ]; % Return the bracket, rather than a single value


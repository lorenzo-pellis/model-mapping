function Rg = find_Rg_from_r_1typeH_varHsize_Sellke(PI_single,desired_r,Alphag,Gamg,Rh_H,Alphah,Gamh,eta,den_str,nsim)

% This function tries to find the value for the overall global infectivity,
% Rg, so that the real-time growth rate r computed from this value of Rg
% matches a desired, pre-specified, value of r, and it is used in the model
% mapping paper to compare models by matching r, instead of R0. Because
% of the assumption of a time-varying, non-random, Gamma-shaped infectivity 
% profile, there is no analytical expression for r, so Monte Carlo
% simulations are run instead, using the Sellke construction.
% This particular file is for a single-type households model. There is a
% different one for a 2-type model.
% There are variants of these files that use C to compute the Sellke
% epidemics.
%
% Input: all the household-related parameters:
%   - PI_single: the size-biased household size distribution (single-type)
%   - Alphag: shape parameter of the global infectivity profile
%   - Gamg: scale parameter of the global infectivity profile
%   - Rh_H: the mean number of infectious contact an infectives
%   makes with all other household members, throughout his/her entire
%   infectious period (the within-household R0)
%   - Alphag: shape parameter of the within-household infectivity profile
%   - Gamg: scale parameter of the within-household infectivity profile
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default in the model mapping paper, but can be changed in general
%   - den_str: a string that specifies the form of frequency dependence:
%       - 'n': transmission has the form beta/n^eta
%       - 'n-1': transmission has the form beta/(n-1)^eta
%   - nsim: the number of Monte Carlo simulations used
% 
% Output:
%   - Rg: the average number cases a single case makes through global 
%   contacts only, during the entire infectious period, in a fully 
%   susceptible population, i.e. the basic reproduction number for global 
%   contacts only
% 
% Note 1: because Rg is found using fzero, and each call requires finding r
% using fzero, I should not run Monte Carlo epidemics in the inner loop
% (stochastic variability might prevent convergence, unless the number of
% simulations i extremely high - and computational efficiency would be an
% issue). Therefore, all Monte Carlo epidemics are run at the beginning and
% stored in a large array "inftimesall".
% 
% Note 2: The code for the Monte Carlo epidemic gives both the final size
% (z) and the times of infection of each infective (inftimes), so
% when I compute r I should use, for each run, only the times of
% individuals that get infected, so that infectives' indices run from 1 to 
% z. I have found it significantly faster to use
% NaN for those not infected (and for indices out of household composition
% bounds) and use logical indices to exclude those when adding the
% contribution of all infectives together. If I do so, I do not need to
% store how many get infected in each epidemic.
% 
% Note 3: When I use logical indices as described in Note 2, if I run
% only 1 epidemic (e.g. for debugging purposes), the indices get all out of
% place because Matlab skips singleton dimensions. Therefore, I
% counterintuitively need to permute all dimensions so that the index of
% current simulation is the last dimension in the big array "inftimesall"
% 
% Reference: Supplementary Methods of
% Pellis, L. et al (2020), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 09-11-2019 

optsRg = optimset('TolX',1e-5);
optsr = []; % optimset('display','iter');
% % For debugging purposes:
% optsRg = optimset('TolX',1e-5,'Display','iter');
% optsr = optimset('TolX',1e-10);

nmax = length(PI_single); % Largest household size

% zall = zeros(maxHsize,nsim); % Final sizes for each Sellke epidemic for each household size
inftimesall = NaN(nmax,nmax,nsim); % Time of infections of all cases in all Sellke epidemics
% Each element of inftimesall is: 
% (household size, # infective, # simulations)

% First, run all Monte Carlo epidemics. Each batch of nsim epidemics need
% to be run for each household size.
for in = 1:nmax
    if strcmp(den_str,'n-1')
        if in == 1
            lambda = 0; % It doesn't matter, as within-household transmission is a household of size 1 is irrelevant
        else
            lambda = Rh_H/(in-1)^eta;
        end
    else
        lambda = Rh_H/in^eta;
    end
    for isim = 1:nsim
        [z, inftimes] = multitypeSellke_Gamma_1init(1, in, lambda, Alphah, Gamh, 1);
        % The function to run the Sellke construction spits
        % out also the final size. However, I do not use it
        % (see Note 2 above), so I do not need to store it
        inftimesall(in,1:z,isim) = inftimes(1,1:z);
        % Note that:
        % 1) Only the infection times of those that got
        % infected are copied, to preserve the NaNs in
        % the initialised "inftimesall" to discard via
        % logical indexing when the sum is made
        % 2) I need to permute the dimensions so that
        % isim is the last one (or I have troubles when
        % I run only 1 simulation as singleton
        % dimensions are skipped).
    end
    disp(['Household size ',num2str(in),' of ',num2str(nmax),' done!']);
end

% Prepare some initial guesses
rguess = desired_r; % I need a sensible starting value for r: given I aim for desired_r, I just go for that
Rgguess = ( ( Gamg + desired_r ) / Gamg )^Alphag; % A starting guess for Rg

% Prepare and run fzero
Rgfun = @(x) ( r_1typeH_varHsize_Sellke_precomputed(PI_single,x,Alphag,Gamg,inftimesall,nsim,rguess,optsr) - desired_r );
Rg = fzeroinf(Rgfun,0,optsRg,2,[],Rgguess);
% fzeroinf tips:
% 1) For debugging reasons, I pass optsRg
% 2) fzeroinf starts from the guess and then expands the search interval by
% a factor. Here I use 2, but using 10 (default) makes little difference.

disp(['Rg = ',num2str(Rg)]);

end

function rhat = r_1typeH_varHsize_Sellke_precomputed(PI_single,Rg,Alphag,Gamg,inftimesall,nsim,rguess,optsr)

% This function computes the real-time growth rate for a households models 
% with time-varying, non-random, Gamma-shaped infectivity profile. 
% Therefore it requires Monte Carlo simulations of all within-household 
% epidemics. Rather than running them, they have been precomputed and the 
% times of infection of all cases are passed as input.
% Other inputs as described in the calling function.
% 
% Output:
%   - rhat: an unbiased estimate of the real-time growth rate, more precise
%   the larger nsim is

% A guess for r is not needed, as it is provided (giving I aim for desired_r, that's what I go for)
% rguess = 0; % Just any initial guess: makes little difference
rfun = @(s) ( spectral_radius_discountedNGM(s,PI_single,Rg,Alphag,Gamg,inftimesall,nsim) - 1 );
rhat = fzeroinf(rfun,-min(min(Gamg)),optsr,[],[],rguess);
% fzeroinf tips:
% 1) For debugging reasons, I pass optsr
% 2) fzeroinf starts from the guess and then expands the search interval by
% a factor. Here I use 2, but using 10 (default) makes little difference.
% 
% For debugging purposes:
% disp(rhat)
end


function y = spectral_radius_discountedNGM(x,PI_single,Rg,Alphag,Gamg,inftimesall,nsim)

% This function computes the Laplace transform of the infectivity profile.
% When x = 0, this should give Rg. The parameters are as
% defined in the calling function.
nmax = length(PI_single);
temp = 0;
for in = 1:nmax
    temp = temp + PI_single(in) * sum( exp( -x * inftimesall(in,~isnan(inftimesall(in,:,:))) ) );
end 
y = Rg * (Gamg / (x+Gamg))^Alphag * temp / nsim;

end




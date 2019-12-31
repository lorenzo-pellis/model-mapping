function Rg = find_Rg_from_r_2typeH_varHsize_Sellke(PI,desired_r,NGM_G_shape,Alphag,Gamg,Rh,shapeHvector,eta,Alphah,Gamh,nsim)

% This function tries to find the value for the overall global infectivity,
% Rg, so that the real-time growth rate r computed from this value of Rg
% matches a desired, pre-specified, value of r, and it is used in the model
% mapping paper to compare models by matching r, instead of R0. Because
% of the assumption of a time-varying, non-random, Gamma-shaped infectivity 
% profile, there is no analytical expression for r, so Monte Carlo
% simulations are run instead, using the Sellke construction.
% This particular file is for a 2-type households model. There is a
% different one for a 1-type model.
% There are variants of these files that use C to compute the Sellke
% epidemics.
%
% Input: all the household-related parameters:
%   - PI: a (n+1)x(n+1)x2 double matrix. The first layer has elements (i,j)
%   giving the probability that the household of a randomly selected adult
%   has (i-1) adults and (j-1) children. The second layer for the household
%   of a randomly selected child. Note: row 1 in layer 1 is 0 and column 1
%   in layer 2 is 0.
%   - NGM_G_shape: the next-generation matrix of global contacts up to a
%   multiplicative factor Rg, which is what we want to find
%   - Alphag: 2x2 matrix of the shape parameters of the global infectivity 
%   profiles from either type to either type
%   - Gamg: 2x2 matrix of the scale parameters of the global infectivity
%   profiles from either type to either type
%   - Rh: the mean number of infectious contact an infectives
%   makes with all other household members, throughout his/her entire
%   infectious period (the within-household R0)
%   - shapeHvector = [ h_ratio, thetaH, psiH, phiH ], containing relevant 
%   household parameters:
%       - h_ratio: c_a/c_c, where c_a and c_c are, respectively, the number
%       of contacts an adult and a child make per unit of time
%       - thetaH: household assortativity, i.e. fraction of contacts a
%       child made with other children. This is always set to NaN, as
%       random mixing is assumed in the household, so that the
%       assortativity is household-composition-dependent (it takes the
%       value of the fraction of (other) children in the household) and
%       therefore needs to be computed for each household composition
%       - psiH: relative susceptibility of children in the household (same
%       as global, by assumption that it is a property of the individual)
%       - phiH: relative infectivity of children in the household (same
%       as global, by assumption that it is a property of the individual)
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default and never changed
%   - Alphah: 2x2 matrix of the shape parameters of the within-household
%   infectivity profiles from either type to either type
%   - Gamh: 2x2 matrix of the scale parameters of the within-household
%   infectivity profiles from either type to either type
%   - nsim: the number of Monte Carlo simulations used
% 
% Output:
%   - Rg: the global infectivity factor, which then multiplies NGM_G_shape
%   to give the next-generation matrix of global contacts
% 
% Note 1: because Rg is found using fzero, and each call requires finding r
% using fzero, I should not run Monte Carlo epidemics in the inner loop
% (stochastic variability might prevent convergence, unless the number of
% simulations i extremely high - and computational efficiency would be an
% issue). Therefore, all Monte Carlo epidemics are run at the beginning and
% stored in a large array "inftimesall".
% 
% Note 2: The code for the Monte Carlo epidemic gives both stratified
% final sizes (z) and times of infection of each infective (inftimes), so
% when I compute r I should use, for each run, only the times of
% individuals that get infected, so adult indices run from 1 to z(1) and
% child indices from 1 to z(2). I have found it significantly faster to use
% NaN for those not infected (and for indices out of household composition
% bounds) and use logical indices to exclude those when adding the
% contribution of all infectives together. If I do so, I do not need to
% store how many get infected in each epidemic.
% 
% Note 3: When I use logical indices as described in note 2, if I run
% only 1 epidemic (e.g. for debugging purposes), the indices get all out of
% place because Matlab skips singleton dimensions. Therefore, I
% counterintuitively need to permute all dimensions so that the index of
% current simulation is the last dimension in the big array "inftimesall"
% 
% Reference: Supplementary Methods of
% Pellis, L. et al (2019), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 06/11/2019 

% For debugging purposes:
optsRg = optimset('TolX',1e-5);
% optsRg = optimset('TolX',1e-5,'Display','iter');
optsr = []; % optimset('display','iter');
% optsr = optimset('TolX',1e-10);


h_ratio = shapeHvector(1);
thetaH = shapeHvector(2);
psiH = shapeHvector(3);
phiH = shapeHvector(4);

[ max_inA,max_inC,~ ] = size(PI); % It needs to read all three dimensions, but I know the third is 2.
maxHsize = max(max_inA,max_inC) - 1; % Largest possible number of infectives of one single type

Lh1to1 = zeros(2); % This is the matrix of average 1-to-1 infection pressure 
% (i.e. infectious rates integrated over time, referred to in the paper as \Lambda_h)
% zall = zeros(2,nsim,2,max_inA,max_inC); % Stratified final sizes storing all Sellke epidemics (not needed - see Note 2 above)
inftimesall = NaN(max_inA,max_inC,2,2,maxHsize,nsim); % Time of infections of all cases in all Sellke epidemics
% Each element of inftimesall is: (# adults, # children, type of primary case,
% type of infective, # infective of that type, # simulations)

% First, run all Monte Carlo epidemics. Each batch of nsim epidemics need
% to be run for each household composition and each type of initial
% infective. For each household composition I first compute the 1-to-1
% total infectivity in the matrix Lh1to1
for inA = 1:max_inA
    for inC = 1:max_inC
        if ( PI(inA,inC,1) ~= 0 || PI(inA,inC,2) ~= 0 )
            nA = inA - 1;
            nC = inC - 1;
            nH = [nA,nC];
            assert(sum(nH) <= maxHsize);
            if isnan(shapeHvector(2)) % I can't use thetaH here, because otherwise I change its value at the first iteration and then it's not NaN any longer
                thetaH = ( nC - 1 ) / ( ( nC - 1 ) + nA * h_ratio ); % If shapeHvector(2) is NaN, it means I want random mixing in the households
            end
            if nA == 0
                if nC > 1
                    Lh1to1(2,2) = Rh * psiH * thetaH * phiH / ( nC - 1 )^eta;
                end
            elseif nA == 1
                if nC == 1
                    Lh1to1(1,2) = Rh * ( 1 - thetaH ) * phiH / nA^eta;
                    Lh1to1(2,1) = Rh * psiH * ( 1 - thetaH ) * ( nC/nA ) / nC^eta;
                elseif nC > 1
                    Lh1to1(1,2) = Rh * ( 1 - thetaH ) * phiH / nA^eta;
                    Lh1to1(2,1) = Rh * psiH * ( 1 - thetaH ) * ( nC/nA ) / nC^eta;
                    Lh1to1(2,2) = Rh * psiH * thetaH * phiH / ( nC - 1 )^eta;
                end
            else
                if nC == 0
                    % If there are no children, then thetaH is not well
                    % defined, but I don't care, since it is
                    % multiplied by nC = 0. So I'll leave as it is.
                    Lh1to1(1,1) = Rh * ( h_ratio - ( 1 - thetaH ) * ( nC/nA ) ) / ( nA - 1 )^eta;
                elseif nC == 1
                    Lh1to1(1,1) = Rh * ( h_ratio - ( 1 - thetaH ) * ( nC/nA ) ) / ( nA - 1 )^eta;
                    Lh1to1(1,2) = Rh * ( 1 - thetaH ) * phiH / nA^eta;
                    Lh1to1(2,1) = Rh * psiH * ( 1 - thetaH ) * ( nC/nA ) / nC^eta;
                else
                    Lh1to1(1,1) = Rh * ( h_ratio - ( 1 - thetaH ) * ( nC/nA ) ) / ( nA - 1 )^eta;
                    Lh1to1(1,2) = Rh * ( 1 - thetaH ) * phiH / nA^eta;
                    Lh1to1(2,1) = Rh * psiH * ( 1 - thetaH ) * ( nC/nA ) / nC^eta;
                    Lh1to1(2,2) = Rh * psiH * thetaH * phiH / ( nC - 1 )^eta;
                end
            end
            for tp = 1:2
                if nH(tp) > 0
                    for isim = 1:nsim
                        [z, inftimes] = multitypeSellke_Gamma_1init(2, nH', Lh1to1, Alphah, Gamh, tp);
                        % The function to run the Sellke construction spits
                        % out also the final size. However, I do not use it
                        % (see Note 2 above), so I do not need to store it
                        for b = 1:2
                            inftimesall(inA,inC,tp,b,1:z(b),isim) = permute(inftimes(b,1:z(b)),[ 4 5 1 2 3 ]);
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
                    end
                end
            end
            disp(['Household composition (',num2str(nA),',',num2str(nC),') of (',num2str(max_inA-1),',',num2str(max_inC-1),') done!']);
        end
    end
end
% Prepare some initial guesses
rguess = 2*desired_r; % I need a sensible starting value for r: given I aim for desired_r, I just go for that
Rgguess = ( ( max(max(Gamg)) + desired_r ) / max(max(Gamg)) )^(min(min(Alphag))) / max(eig(NGM_G_shape)) ; % A starting guess for Rg

% Prepare and run fzero
Rgfun = @(x) ( r_2typeH_varHsize_Sellke_precomputed(PI,x*NGM_G_shape,Alphag,Gamg,inftimesall,nsim,rguess,optsr) - desired_r );
Rg = fzeroinf(Rgfun,0,optsRg,2,[],Rgguess);
% fzeroinf tips:
% 1) For debugging reasons, I pass optsRg
% 2) fzeroinf starts from the guess and then expands the search interval by
% a factor. Here I use 2, but using 10 (default) makes little difference.

disp(['Rg = ',num2str(Rg)]);

end

function rhat = r_2typeH_varHsize_Sellke_precomputed(PI,NGM_G,Alphag,Gamg,inftimesall,nsim,rguess,optsr)

% This function computes the real-time growth rate for a multitype 
% households models with time-varying, non-random, Gamma-shaped infectivity 
% profile. Therefore it requires Monte Carlo simulations of all 
% within-household epidemics. Rather than running them, they have been
% precomputed and the times of infection of all cases are passed as input.
% Other inputs as described in the calling function.
% 
% Output:
%   - rhat: an unbiased estimate of the real-time growth rate, more precise
%   the larger nsim is

% A guess for r is not needed, as it is provided (giving I aim for desired_r, that's what I go for)
% rguess = 0; % Just any initial guess: makes little difference
rfun = @(s) ( spectral_radius_discountedNGM(s,PI,NGM_G,Alphag,Gamg,inftimesall,nsim) - 1 );
rhat = fzeroinf(rfun,-min(min(Gamg)),optsr,[],[],rguess);
% fzeroinf tips:
% 1) For debugging reasons, I pass optsr
% 2) fzeroinf starts from the guess and then expands the search interval by
% a factor. Here I use 2, but using 10 (default) makes little difference.
% 
% For debugging purposes:
% disp(rhat)
end


function y = spectral_radius_discountedNGM(x,PI,NGM_G,Alphag,Gamg,inftimesall,nsim)

% This function computes the matrix of Laplace transforms of the
% infectivity profiles from either type to eather type and its dominant
% eigenvalue. When x = 0, this should give the NGM. The parameters are as
% defined in the calling function.
[ max_inA,max_inC,ntypes ] = size(PI); % ntypes = 2
Phi = zeros(ntypes);
for a = 1:ntypes % Type of infectee
    for tp = 1:ntypes % Type of the primary case
        temp = zeros(ntypes,1);
        for b = 1:ntypes % Type of infector
            for inA = 1:max_inA
                for inC = 1:max_inC
                    temp(b) = temp(b) + PI(inA,inC,tp) * sum( exp( -x * inftimesall(inA,inC,tp,b,~isnan(inftimesall(inA,inC,tp,b,:,:))) ) );
                end
            end
        end
        Phi(a,tp) = ( NGM_G(a,:).*(Gamg(a,:)./(x+Gamg(a,:))).^Alphag(a,:) ) * temp / nsim;
    end
end
y = max(eig(Phi));

end




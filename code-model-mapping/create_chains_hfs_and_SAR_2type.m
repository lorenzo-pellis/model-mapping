function [ av_gen_chains, hfs_by_comp_by_ic, SAR_by_comp_by_ic ] = create_chains_hfs_and_SAR_2type(PI,Rh,shapeHvector,eta,bincoeffmat)

% The main purpose of this function is to create av_gen_chains (the other
% outputs are less important, but created here for convenience), which is
% key for the computation of R0 for model AH: in the paper, \mu^p_{t,i}
% denotes element av_gen_chain(t,i,p), where t is the index of the type
% (age class), i is the index of the generation in the within-household
% epidemic and p is the index of the type of the primary case.
%
% Input: all the household-related parameters:
%   - PI: a (n+1)x(n+1)x2 double matrix. The first layer has elements (i,j)
%   giving the probability that the household of a randomly selected adult
%   has (i-1) adults and (j-1) children. The second layer for the household
%   of a randomly selected child. Note: row 1 in layer 1 is 0 and column 1
%   in layer 2 is 0.
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
%   - bincoeffmat: a precomputed matrix of binomial coefficients
% 
% Output:
%   - av_gen_chains: this is a 2xnx2 array giving the average number of
%   cases of each type (adult or children, index 1), in each generation of
%   the within-household epidemic (index 2) started with a primary case of
%   each type (index 3). The average is both with respect to the stochastic
%   dynamics of the epidemic within a household of specific composition,
%   and also across all household compositions, in proportions given by the
%   probability that the household of a randomly selected adult and child
%   has a specific composition (array PI in input).
%   - hfs_by_comp_by_ic: nxnx2 array, with element (i,j,k) giving the
%   average final size of an epidemic within a household with (i-1) adults
%   and (j-1) children started with an adult (k=1) or child (k=2) primary
%   case (hfs = household final size). Note that this cannot be derived
%   directly by av_gen_chains, because the latter is already averaged over
%   the household composition distribution.
%   - SAR_by_comp_by_ic: nxnx2 array, similar to hfs_by_comp_by_ic, but
%   reporting the secondary attack rate, i.e. the fraction of initial
%   susceptibles in the household that end up infected.
% 
% Reference: Supplementary Methods, Section 1.2.4 of
% Pellis, L. et al (2020), Nature Communications
%
% Note: this function relies on the function "create_avRFchain_6dim". Old
% versions relied on a recursive method I had developed, but that turned
% out to be less efficient.
%
% Author: Lorenzo Pellis
% Last update: 31-10-2019 

h_ratio = shapeHvector(1);
thetaH = shapeHvector(2);
psiH = shapeHvector(3);
phiH = shapeHvector(4);

[ max_inA,max_inC,~ ] = size(PI); % It needs to read all three dimensions, but I know the third is 2.

% Create the chains of infectives in each generation
% max_gen = ( max_inA - 1 ) + ( max_inC - 1 ); % This is the general case (you could have all cases of one type and then all case fo the other) 
max_gen = max_inA - 1; % This is the real max number of generations in practice, because households are not larger than 8
av_gen_chains = zeros(2,max_gen,2); % Each layer is the average chain when you start with a single case of each type (first a, second c)
av_gen_chains_temp = zeros(2,max_gen,2);
Lh1to1 = zeros(2); % This is the matrix of average 1-to-1 infection pressure 
% (i.e. infectious rates integrated over time, referred to in the paper as \Lambda_h)
A = eye(2);
hfs_by_comp_by_ic = zeros(max_inA,max_inC,2);
SAR_by_comp_by_ic = zeros(max_inA,max_inC,2);

for inA = 1:max_inA
    for inC = 1:max_inC
        if ( PI(inA,inC,1) ~= 0 || PI(inA,inC,2) ~= 0 )
            nA = inA - 1;
            nC = inC - 1;
            nH = [nA,nC];
            Ngen = sum(nH);
            if Ngen <= max_gen
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
                for tprim = 1:2
                    if all(A(tprim,:)<=nH)
                        av_chain_6dim = create_avRFchain_6dim(nH,A(tprim,:),Lh1to1,bincoeffmat);
                        av_gen_chains_temp(:,1:Ngen,tprim) = av_gen_chains_temp(:,1:Ngen,tprim) + PI(inA,inC,tprim) * av_chain_6dim;
                        hfs_by_comp_by_ic(inA,inC,tprim) = sum(sum(av_chain_6dim));
                        if (nA+nC) > 1
                            SAR_by_comp_by_ic(inA,inC,tprim) = ( hfs_by_comp_by_ic(inA,inC,tprim) - 1 ) / ( nA + nC - 1 );
                        end
                    end
                end
            end
        end
    end
end
av_gen_chains = av_gen_chains_temp;


function afs_compbiased = afs_Addy_2type_compbiased(H,PI,Rh,shapeHvector,pesc,eta)

% This function computes the average final size for each type, expressed as
% fraction of the total numbers of each type, of a within-household
% epidemic where individuals also escape global infection from outside, and
% allowing a distribution of household compositions. The hard part is to
% compute the average final outcome for each household composition, which
% is done in the function "afs_Addy_2type_const". What I do here is only
% preparing the right infection rates for each household composition and
% average over the correct (i.e. composition-biased) household composition
% distribution. The method is restricted to 2 types (adults and children) 
% only, and assumes constant total individual infectivity (constant 
% duration of infectious period, or a non-random time-since-infection 
% infectivity profile).
% 
% Input:
%   - H: an (n+1)x(n+1) matrix with element (i,j) giving the probability 
%   that a randomly selected household has (i-1) adults and (j-1) children,
%   where n is the largest household size.
%   - PI: a (n+1)x(n+1)x2 double matrix. The first layer has elements (i,j)
%   giving the probability that the household of a randomly selected adult
%   has (i-1) adults and (j-1) children. The second layer for the household
%   of a randomly selected child. Note: row 1 in layer 1 is 0 and column 1
%   in layer 2 is 0 (can be computed from H, but I've already done it 
%   before, so I might as well pass it directly).
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
%   - pesc = [q_a,q_c]: a vector with the probabilities of susceptible
%   adults and children, respectively, escaping global infection from 
%   outside the household.
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default and never changed
% 
% Reference: Supplementary Methods, Section 1.2.7 of
% Pellis, L. et al (2020), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 29-12-2019 


h_ratio = shapeHvector(1);
thetaH = shapeHvector(2);
psiH = shapeHvector(3);
phiH = shapeHvector(4);

afs_temp = zeros(2,1); % column vector
[ max_inA,max_inC ] = size(H);
for inA = 1:max_inA % index for the number of adults (= nA + 1)
    for inC = 1:max_inC % index for the number of children (= nC + 1)
        if H(inA,inC) ~= 0 % Do this only if there are households with this composition (the bottom-right half of H is all 0, as there is a cut-off on the household size)
            nA = inA - 1;
            nC = inC - 1;
            h = [nA,nC];
            Lh1to1 = zeros(2); % This is the matrix of average 1-to-1 infection pressure 
            % (i.e. infectious rates integrated over time, referred to in the paper as \Lambda_h)
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
            afs_temp = afs_temp + reshape(PI(inA,inC,:),[2,1]) .* afs_Addy_2type_const(h,Lh1to1,pesc);
        end
    end
end
afs_compbiased = afs_temp;

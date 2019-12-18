function R0 = R0H_varHsize_check(PI,Rh,Rg)

% This function computes R0 for a (single-type) model with households, but
% it's mostly for cross-checking the more efficient function ("R0H") that I
% use in the main code, where the chains of average cases in each
% generation are first created and then R0 is computed directly from them.
% Here instead I create the chains inside the function (in 2 or 3 ways), to
% check they are correct and which form is most efficient.
% 
% Input:
%   - PI: size-biased household size distribution, i.e. PI(n) gives the 
%   probability that a randomly selected individual in the population 
%   belongs to a household of size n
%   - Rh: the mean number of infectious contact an infectives
%   makes with all other household members, throughout his/her entire
%   infectious period (the within-household R0)
%   - Rg is the average number cases a single case makes through global 
%   contacts only, during the entire infectious period, in a fully 
%   susceptible population, i.e. the basic reproduction number for global 
%   contacts only.
% 
% Note: R0 could be computed directly from the average, which is what I do
% in the main code (function "R0H")

max_gen = length(PI);

mus = create_chain_1type_check(PI,Rh,1,'n-1'); % Here eta is 1, but I don't want to generalise it
M = zeros(max_gen);
M(1,1:max_gen) = Rg*mus;
for i = 2:max_gen
    M(i,i-1) = 1;
end

R0 = max(eig(M));


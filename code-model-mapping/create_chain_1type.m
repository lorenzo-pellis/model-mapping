function av_gen_chain = create_chain_1type(PI,Rh,eta,str_den)

% This function creates the average number of cases in each generation of
% an epidemic in a small population (which is used to then compute R0 for
% the households model). The average is over the stochastic dynamics and 
% the household size distribution. This code assumes a pure Reed-Frost 
% model, as the core function called ("create_avRFchain_3dim") assumes a
% pure Reed-Frost model.
% 
% Input:
%   - PI: size-biased household size distribution, i.e. PI(n) gives the 
%   probability that a randomly selected individual in the population 
%   belongs to a household of size n
%   - Rh: the mean number of infectious contact an infectives
%   makes with all other household members, throughout his/her entire
%   infectious period (the within-household R0)
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default in the model mapping paper, but can be changed in general
%   - str_den: a string that specifies the form of frequency dependence:
%       - 'n': transmission has the form beta/n^eta
%       - 'n-1': transmission has the form beta/(n-1)^eta
% 
% Output: av_gen_chain is row vector of length N = largest household size,
% with element av_gen_chain(n) giving the average number of cases in
% generation n starting from 1 initial case (so av_gen_chain(1) = 1 by
% definition). The average is over the stochastic dynamics and the
% size-biased household size distribution.
% 
% Author: Lorenzo Pellis
% Last update: 05-06-2019 

max_gen = length(PI);
av_gen_chain_temp = zeros(1,max_gen);
for n = 1:max_gen
    if PI(n) ~= 0
        if strcmp(str_den,'n-1')
            lambda = Rh/(n-1)^eta; % 1-to-1 escaping probability
        else
            lambda = Rh/n^eta;
        end
        av_gen_chain_temp(1:n) = av_gen_chain_temp(1:n) + PI(n) * create_avRFchain_3dim(n,1,lambda);
    end
end
av_gen_chain = av_gen_chain_temp;

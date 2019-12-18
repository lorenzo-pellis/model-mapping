function TP = get_pAA_from_Rh(distr,Rh)

% This function obtains the value of the adult-to-adult transmission
% probability from the value overall infectivity parameter within the
% household Rh, i.e. the mean number of infectious contact an infectives
% makes with all other household members, throughout his/her entire
% infectious period (the "R0" within the household, just not all infectious
% contacts will result in an infection)
% 
% Reference: Supplementary Methods, Section 1.2.3, and Supplementary Figure 2 of
% Pellis, L. et al (2019), Nature Communications
% 
% Note: in the paper, Rh is called \beta_h
%
% Author: Lorenzo Pellis
% Last update: 31/10/2019 

max_size = length(distr);
distr_2up = [ 0, distr(2:max_size) / sum( distr(2:max_size) ) ]; % This is the household size distribution, conditional on at least 2 people in the household

temp = 0;
for n = 2:max_size
    lambda1to1 = Rh / ( n - 1 ); % This is the 1-to-1 infection pressure  
    % (i.e. person-to-person infection rate integrated over time)
    temp = temp + distr_2up(n) * exp( -lambda1to1 );
end
TP = 1 - temp;
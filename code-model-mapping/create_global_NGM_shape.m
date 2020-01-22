function S = create_global_NGM_shape(F_ratio,psi,phi,theta,g_ratio)

% This function creates the shape of the next generation matrix for global
% infectious contacts. It's not yet a NGM (needs multiplying by a scaling 
% factor).
% 
% Input:
%   - F_ratio = F_c/F_a: ratio of fraction of children divided by fraction
%   of adults in the population
%   - psi: relative susceptibility of children versus adults
%   - phi: relative infectivity of children versus adults
%   - theta: assorativity of children, defined as the fraction the contacts
%   of a child that are made with other children
%   - g_ratio = c_a/c_c: ratio between the number of contacts an adults 
%   makes, divided by those a child makes, per unit of time (notice the
%   ratio is the opposite from F_ratio)
% 
% Reference: matrix K in Supplementary Methods, Section 1.1.5 of
% Pellis, L. et al (2020), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 12-05-2019 

S = [       g_ratio - ( 1 - theta ) * F_ratio            ( 1 - theta ) * phi            
                psi * ( 1 - theta ) * F_ratio              psi * theta * phi    ];

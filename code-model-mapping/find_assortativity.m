function theta = find_assortativity(v,F_ratio,psi,phi,g_ratio,pAA,failsymb,testplot)

% This function searches for the correct value of the assortativity for 
% model A, such that the dominant eigenvector of the next generation matrix
% of model A is the same as the one "given" by model AH. This is restricted
% to a 2-type model (adults and children), as relevant for the model
% mapping paper below.
% 
% Input:
%   - v = [v_a,v_c]: vector with fractions of adults and children in each
%   generation of model AH (the "truth"), which this code attempts to match
%   - F_ratio = F_c/F_a: ratio of fraction of children divided by fraction
%   of adults in the population
%   - psi: relative susceptibility of children versus adults
%   - phi: relative infectivity of children versus adults
%   - g_ratio = c_a/c_c: ratio between the number of contacts an adults 
%   makes, divided by those a child makes, per unit of time (notice the
%   ratio is the opposite from F_ratio)
%   - pAA: infectious adult to susceptible adult probability of infection,
%   a measure of the within-household infectivity
%   - failsymb: a specified output for when no assortativity can be found
%   - testplot: a Boolean input, to request a plot (to increase
%   understanding)
%
% Output: the new "mapped" value of the assortativity for the closest model
% A to model AH
% 
% Note: this file contains 2 functions, and some of the variables span both
% of them (so they are highlighted in a different colour)
% 
% Reference: Supplementary Methods, Sections 1.3 of
% Pellis, L et al (2019), Nature Communications
%
% Author: Lorenzo Pellis
% Last update: 26/05/2019 

if testplot

    % First I explore and plot the range of the assortativity
    ass = 0:0.1:1;
    lt = length(ass);
    w = zeros(2,lt);

    for x = 1:lt
        w(:,x) = get_dominant_eigenvector(create_global_NGM_shape(F_ratio,psi,phi,ass(x),g_ratio));
    end
    
    figure(98)
    clf;
    hold on;
    set(gca,'ColorOrder',[ 0 0 1; 1 0 0 ]); % Set the first line in blue (adults) and the second in red (children)
    plot(ass,w,'Linewidth',2);
    plot([0,1],[v(1),v(1)],'b-.','Linewidth',2);
    plot([0,1],[v(2),v(2)],'r-.','Linewidth',2);
    % title(['Eigenvectors v^{AH} and v^A(\theta^A) for p_{AA} = ',num2str(pAA),' and \psi = ',num2str(psiG)],'Fontsize',16);
    title(['Proportions of adults and children (p_{AA} = ',num2str(pAA),' and \psi = ',num2str(psi),')'],'Fontsize',16);
    xlabel('Assortativity of children - A model (\theta^A)','Fontsize',12);
    ylabel('Fraction of incidence','Fontsize',12);
    legend('v^A_a','v^A_c','v^{AH}_a','v^{AH}_c','Location','NorthEast');
    hold off;
    drawnow();
end    

% Find the assortativity
funvc = @(a) (v(2) - getvc(a)); % Vector components sum to 1, so it's sufficient to match one
vc0 = funvc(0);
vc1 = funvc(1);
if sign(vc0) == sign(vc1) % Test if there is a potential solution
    theta = failsymb; 
else
    theta = fzero(funvc,[0,1]);
end


function vc = getvc(a)
    u = get_dominant_eigenvector(create_global_NGM_shape(F_ratio,psi,phi,a,g_ratio));
    vc = u(2);
end

end

function av_gen_6dim = create_avRFchain_6dim(n,a,Lh1to1,bincoeffmat)

% This function computes the average number of cases of each type in each
% generation of a Reed-Frost model, in a specific household with n=[n1,n2] 
% cases of each type starting with a=[a1,a2] initial infectives of each 
% type (2 types only).
% 
% Inputs:
%   - n = [n1,n2] gives the numbers of each type
%   - a = [a1,a2] is the vector of initial cases (in general I'll have only
%   [1,0] or [0,1]
%   - Lh1to1 is a matrix of average 1-to-1 infection pressure (i.e. 
%   infectious rates integrated over time, referred to in the model mapping 
%   paper as \Lambda_h), so a next generation matrix for the infectious 
%   contacts would be:
%        (n1-1) * Lh1to1(1,1)       n1  * Lh1to1(1,2)
%           n2  * Lh1to1(2,1)    (n2-1) * Lh1to1(2,2)
%   - Rg is a 2-by-2 matrix representing how many individuals of each type 
%   are infected by an individual of each type in the community: they'll be
%   primary cases in new household.
%
% Output:
%   - av_gen_6dim: a 2xMaxk matrix, where Maxk is the largest possible 
%   generation index (the name 6dim comes from how the calculation is done,
%   not to describe what this array is).
% 
% Reference: Supplementary Methods, Section 1.2.4 of
% Pellis, L. et al (2020), Nature Communications
%
% Note: this function uses a 6-dimensional array binP_6dim, where element
% binP_6dim(ia1,ia2,im1,im2,is1,is2) represents the probability that m1 and
% m2 susceptibles of types 1 and 2 (here adults and children), out of the
% initial s1 and s2, escape infection from a1 and a2 infectious cases of
% types 1 and 2. Given some of these numbers can be 0, I adopt the
% convention of using letter i to denote the associated index, which is
% always 1 more than the actual number (so, ia1 is the index corresponding
% to a1=ia1-1 infectious cases of type 1).
% 
% Methodology reference: Appendix A of Pellis, Ball & Trapman (2012), 
% Mathematical Biosciences (only the single-type model, so only needs a 
% 3-dimensional array)
%
% Author: Lorenzo Pellis
% Last update: 31-10-2019 


% Although the method in Pellis, Ball & Trapman (2012), Mathematical
% Biosciences is general, in this context we assume simple Reed-Frost 
% model, because the model mapping paper assumes a time-since-infection 
% model, so that all infectives have the same non-random total infectivity.

% For the Reed-Frost model, we have:
Q = exp(-Lh1to1); % matrix of the 1-to-1 escaping probabilities (it's not the exponential of a matrix, it's elementwise)

in1 = n(1)+1;
in2 = n(2)+1;
Maxk = sum(n); % Largest generation index if we start counting from 1
binP_6dim = zeros(in1,in2,in1,in2,in1,in2);

% Fill in the matrix of binomial escaping probabilities
% Older versions were computing each term, but this is costly, so I now use
% a pre-computed matrix of binomial coefficients
for ia1 = 1:in1
    for ia2 = 1:in2
        q1 = Q(1,1)^(ia1-1) * Q(1,2)^(ia2-1);
        q2 = Q(2,1)^(ia1-1) * Q(2,2)^(ia2-1);
        for is1 = 1:(in1-ia1+1)
            for is2 = 1:(in2-ia2+1)
                for im1 = 1:is1
                    for im2 = 1:is2
%                         binP_6dim_Ball(ia1,ia2,im1,im2,is1,is2) = nchoosek(is1-1,im1-1) * q1^(im1-1) * (1-q1)^(is1-im1) * nchoosek(is2-1,im2-1) * q2^(im2-1) * (1-q2)^(is2-im2);
                        binP_6dim(ia1,ia2,im1,im2,is1,is2) = bincoeffmat(is1,im1) * q1^(im1-1) * (1-q1)^(is1-im1) * bincoeffmat(is2,im2) * q2^(im2-1) * (1-q2)^(is2-im2);
                    end
                end
            end
        end
    end
end

% Now work out the average number of cases mu_6dim, where
% mu_6dim(t,k,ia1,ia2,is1,is2) gives the average number of cases of type t
% in generation k of an epidemic started with a1=ia1-1 and a2=ia2-1 initial
% infectives of type 1 and 2, and s1=is1-1 and s2=is2-1 initial
% susceptibles of types 1 and 2.

% First, set up the initial conditions for the recursive relation
mu_6dim = zeros(2,Maxk,in1,in2,in1,in2);
for ia1 = 1:in1
    for ia2 = 1:in2
        for is1 = 1:(in1-ia1+1)
            for is2 = 1:(in2-ia2+1)
                mu_6dim(1,1,ia1,ia2,is1,is2) = ia1-1;
                mu_6dim(2,1,ia1,ia2,is1,is2) = ia2-1;
            end
        end
    end
end
% All other initial conditions are 0, which are automatically set by
% initialising mu_6dim

% Now write the recursive relation
for ik = 2:Maxk
    for ia1 = 1:in1
        for ia2 = 1:in2
            for is1 = 1:(in1-ia1+1)
                for is2 = 1:(in2-ia2+1)
                    temp1 = 0;
                    temp2 = 0;
                    for ii = 1:is1
                        for ij = 1:is2
                            temp1 = temp1 + mu_6dim(1,ik-1,ii,ij,is1-ii+1,is2-ij+1) * binP_6dim(ia1,ia2,is1-ii+1,is2-ij+1,is1,is2);
                            temp2 = temp2 + mu_6dim(2,ik-1,ii,ij,is1-ii+1,is2-ij+1) * binP_6dim(ia1,ia2,is1-ii+1,is2-ij+1,is1,is2);
                        end
                    end
                    mu_6dim(1,ik,ia1,ia2,is1,is2) = temp1;
                    mu_6dim(2,ik,ia1,ia2,is1,is2) = temp2;
                end
            end
        end
    end
end

% Finally, I'm only interested in the chains with a=[a1,a2] initial
% infectives out of a population of composition n=[n1,n2]
av_gen_6dim = mu_6dim(:,:,a(1)+1,a(2)+1,in1-a(1),in2-a(2));

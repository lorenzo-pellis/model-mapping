function RFgen = create_avRFchain_3dim(n,a,lambda)

% This function creates the average number of cases in each generation of
% an epidemic in a small population of size n with a initial infectives. 
% This is used to then compute R0 for the households model. This code 
% assumes a pure Reed-Frost model.
% 
% Input:
%   - n: the size of the group
%   - a: the number of initial infectives
%   - lambda: the 1-to-1 infection pressure, i.e. the total infectivity
%   from a single infective to a single susceptible.
% 
% Output: av_gen_chain is row vector of length N = largest household size,
% with element av_gen_chain(n) giving the average number of cases in
% generation n starting from 1 initial case (so av_gen_chain(1) = 1 by
% definition). The average is over the stochastic dynamics and the
% size-biased household size distribution.
% 
% Author: Lorenzo Pellis
% Last update: 05/06/2019 

global bincoeffmat;

in = n+1;
binP_3dim = zeros(in,in,n);
for ib = 1:in
    q = exp(-(ib-1)*lambda);
    for is = 1:(in-ib+1)
        for im = 1:is
            binP_3dim(ib,im,is) = bincoeffmat(is,im) * q^(im-1) * (1-q)^(is-im);
        end
    end
end

mu_3dim = zeros(n,in,in);
for ib = 1:in
    for is = 1:(in-ib+1)
        mu_3dim(1,ib,is) = ib-1;
    end
end

for ik = 2:n
    for ib = 2:in
        for is = 2:(in-ib+1)
            temp = 0;
            for ij = 2:(is-ik+2)
                temp = temp + mu_3dim(ik-1,ij,is-ij+1) * binP_3dim(ib,is-ij+1,is);
            end
            mu_3dim(ik,ib,is) = temp;
        end
    end
end

RFgen = mu_3dim(:,a+1,in-a)';        
e = 0;

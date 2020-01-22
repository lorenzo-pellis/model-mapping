function RFgen = matrix_avRFchain_den(TOT,i0,R0,distribution,iota,alpha,den_str)

% This function computes the average number of cases in each generation of
% an epidemic in a small group by computing probabilities for each chain
% using an efficient matrix multiplication technique. It uses the algorithm
% in Picard and Lefevre (1990) for the "collective" Reed-Frost model, which
% is even more general than the randomised Reed-Frost model.
% 
% Inputs:
% TOT = group size (including initial infective)
% i0 = number of initial infectives
% R0 = mean of the 1-to-all total infectivity. 
% distribution = distribution of the total infectivity
%   0 = constant, leading to the simple Reed-Frost model (e.g. constant duration of infection)
%   1 = exponentially distributed (e.g. Markovian case of constant recovery rate)
%   2 = Gamma distributed
% iota = mean of the distribution: usually irrelevant as it can be absorbed in R0
% alpha = shape parameter of distribution, if Gamma
% den_str = a string (either 'n-1' or 'n'), telling if the 1-to-1
% infectivity is R0/TOT or R0/(TOT-1).
% 
% Construct a matrix P of transition probabilities from each state (s,a) to each state (m,i)
% I use the function map_h_state to map pair of (susc,inf) into a single
% index. The matrix P is here called GMP, and its construction is built in,
% instead of being relegated to another function. The method is the
% following:
% 
% Pesc = probability vector of 1,2,..,s susceptibles escaping out of the
% initial s. They fill columns (see the example matrix), except for the
% first element at the top (all infected), which is found as 1-sum(Pesc).
% q(k) = mgf(k*lambda,distribution,iota,alpha)
% Q(k) = mgf(k*lambda,distribution,iota,alpha)^a, where a is the number of
% initial infectives
% 
% Current generation (s,a) <--> h
% Next generation    (m,i) <--> k
% 
% Example with 3 people, a represents the number of infectives existing at one time:
    % a = 1 --> prob 0 cases out of 0   1 out of 1      2 out of 2  --> in this line nobody escapes, i.e. k = 0
    % n = 3             0               0 out of 1      1 out of 2  --> k=1
    %                   0                   0           0 out of 2  --> k=2

    % a = 2 -->     0 out of 0      1 out of 1      nothing
    % n = 2             0           0 out of 1      nothing
    %               nothing         nothing         nothing

    % a = 3 -->     0 out of 0      nothing         nothing
    % n = 1         nothing         nothing         nothing
    %               nothing         nothing         nothing
% 
% Last update: 12-05-2019 

if i0==round(i0)
else
    error('Problem: non-integer initial number of infectives');
end

max = (TOT+2)*(TOT+1)/2-1;
GMP = zeros(max);

if strcmp(den_str,'n-1')
    av_Phi = R0/(TOT-1); % Phi is the (random) total one-to-one infectivity, av_Phi is its average.
else
    av_Phi = R0/TOT;
end
lambda = av_Phi/iota; % If there is constant infectivity, then lambda is the one-to-one infection rate
s0 = TOT-i0; % initial number of susceptibles 
q = zeros(1,TOT-1); % vector of probabilities of an entire set of s susceptbles escaping infection from a single infective. 
% Note that the I use the index s = 1,...,Ms (0 susc case is not counted),
% instead of is = 1,...,ls=Ms+1, so q(1) is the prob of 1 susc escaping
% It goes up to TOT-1 
for s = 1:TOT-1
    q(s)=mgf(s*lambda,distribution,iota,alpha);
end
a = 0; % This is a particular case, because I don't solve the triangular system and I exclude the case where they are all susceptibles
for is = 1:TOT % Exclude s = TOT
    s = is-1;
    k = map_gen_2to1(TOT,s,a);
    GMP(k,k) = 1; % There is no transition in this case, so h=k and P=1
end
for a = 1:TOT % Now, for every a we redefine the maximum number of susceptibles
    Ms = TOT-a;
    Q = zeros(1,Ms);
    for ss = 1:Ms
        Q(ss)=q(ss)^a;
    end
    for is = 1:(TOT-a+1) % Now I consider all values of s from 0 to TOT-a
        s = is-1;
        k = map_gen_2to1(TOT,s,a);
        %%%%%%% Here's the triangular system of Picard and Lefevre %%%%%%%
        A = zeros(s);
        b = zeros(s,1);
        for l = 1:s
            for j = l:s
                A(l,j)=prod(j:-1:(j-l+1));
            end
            b(l,1)=prod(s:-1:(s-l+1)) * Q(l);
        end
        Pesc = A\b; % Column vector of escaping probabilities
        h = map_gen_2to1(TOT,0,s); % Particular case where nobody escapes...
        GMP(h,k) = 1-sum(Pesc);
        for m = 1:s
            i = s-m;
            h = map_gen_2to1(TOT,m,i);
            GMP(h,k) = Pesc(m);
        end
    end
end

av_gen = zeros(1,s0+1);
av_gen(1) = i0;
Pvec = zeros(max,1);
k = map_gen_2to1(TOT,s0,i0);
Pvec(k)=1;
for g = 2:s0+1
    Pvec = GMP * Pvec;
    temp = 0;
    for k = 1:max
        pair = map_gen_1to2(TOT,k);
        temp = temp + pair(2)*Pvec(k); % pair(2) = i;
    end
    av_gen(g) = temp;
end
RFgen = av_gen;
            
        

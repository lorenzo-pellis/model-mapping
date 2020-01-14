function [z, inftimes]=multitypeSellke_Gamma_1init(m, n, Lh1to1, Alpha, Gam, tpc)

% This code runs a Monte Carlo epidemic in a small group with individuals 
% of m different types, using the Sellke construction (see references 
% below) and gamma-shaped infectivity profiles between each pair of types. 
% The epidemic is initiated by a single individual of a specified type.
% 
% Input:
%   - n = column vector [n(1);n(2);...;n(m)] of the number of individuals 
%   of each type
%   - Lh1to1: a matrix containing the 1-to-1 total infectivities, i.e.
%   element (i,j) gives the total infectivity a specified infective of type
%   j exerts on a specified individual of type i
%   - Alpha = matrix of alpha parameters of the per-type Gamma TVI
%   - Gam = matrix of the scale parameters of the per-type Gamma TVI
%   - tpc = type of the primary case
%   [ - In future versions: i0 = column vector of number of initial cases 
%   of each type ]
% 
% Output:
%   - z = column vector of stratified final size (including initial cases)
%   - inftimes = m x n matrix, with the times of infection of all type-a
%   individuals in row a
% 
% Reference: Supplementary Methods, Section 2.3.2 of
% Pellis, L et al (2019), Nature Communications
% 
% Methodological references: 
% Sellke (1983), Journal of Applied Probability
% House, Ross & Sirl (2013), Proceedings of the Royal Society A
% book by Andersson & Britton (2000)
%
% Author: Lorenzo Pellis
% Last update: 31/10/2019 

% rng(7); % Sets the seed for the random number generator, to aid reproducibility
% However, for the paper, the seed is set in the main code
% "Model_Mapping_code_match_r", so this should stay commented out

global opts1;
guessintsize = max(max(Alpha./Gam)); % A fixed amount to identify an initial guess for an fzero call
tinf = 0; % Current time, which will progress as we follow the epidemic

% opts = optimset('Display','iter','TolX',1e-5);
% opts = optimset('TolX',1e-3); % Using TolX = 1e-3 gives an error of 1e-5 or 1e-6 compared to To;X = 1e-10, but it was ~20s instead of ~22s 

% m = length(n); % Number of types
N = max(n); % Size of largest group
inftimes=-ones(m,N); % Initialise matrix of infections times
inftimes(tpc,1)=0; % The initial infective is infected at time 0
Qm = Inf(m,N); % Matrix containing the resilience thresholds
z=zeros(m,1); % Vector of number of each type currently infected - will end up giving the final sizes 
Qm(tpc,1) = 0; % Initial case has already been infected
z(tpc) = 1; % Add initial case case to final size counter
s = n;
s(tpc) = s(tpc)-1; % Current number of susceptibles of each type
tbracket = -ones(m,2); % Interval bracketing the time of infection of the weakest susceptible of each type
nextQ = NaN(m,1); % Resilience threshold of weakest remaining susceptible (next to be infected) of each type
Qrand = -log(rand(1,sum(s))); % Draw all random threshold in one go for efficiency
Qcount = 0; % Counter to run through the thresholds
for im = 1:m
    if s(im) == 0
        nextQ(im) = Inf;
    else
        Qm(im,(z(im)+1):n(im)) = sort(Qrand((Qcount+1):(Qcount+s(im)))); % Distribute resilience thresholds to individuals of each type 
        Qcount = Qcount + s(im); % Move on the counter to use new thresholds
        tbracket(im,:)=[0 Inf]; % Initially, the times of infection can be anything
        nextQ(im) = Qm(im,z(im)+1); % The smallest threshold for those of type im
    end
end
infpres = Lh1to1(:,tpc); % Initial infection pressure acting on a susceptible of each type comes from the initial infective
nextpotinfs = infpres > nextQ; % Qm(sub2ind(size(Qm),1:m,[z+1]'))';
% Avoid doing costly gamma computations if you know somebody won't get
% infected even if challenged with the maximum total infectivity (t->Inf)
while any(nextpotinfs) && all(z<=n) % If at least one gets infected & ???
    nextpotinftimes = Inf(m,1); % Initialise potential times of infection for the weakest individual of each type
    for im = 1:m
        if nextpotinfs(im) % Look for time of infection only if the weakest individual of this type can get infected at some point
            % We look for the time at which the total pressure from all
            % infectors matches the resilience threshold of the weakest of type im
            funr = @(s) ( cuminfpress(s,m,Lh1to1(im,:),inftimes,Alpha(im,:),Gam(im,:)) - Qm(im,z(im)+1) );
            % The time of the next infection can only be larger than the
            % time where we have currently arrived (gammaguess(im,1), which
            % is also set to tinf for every im):
            if isinf(tbracket(im,2))
                % If no time for an infectious contact has been computed 
                % before for the weakest susceptible of type im,
                % tbracket(im,2) = Inf, so I use fzeromin with minimum = tinf
                % (i.e. current time) and initial guess = tinf + guessintsize
                % (i.e. a fixed amount, ideally roughly equal to the generation time)
                
                % To test/debug, uncomment next line:
%                 disp(['Search in [ ',num2str(tbracket(im,1)),', Inf ], with initial guess ',num2str(tinf+guessintsize)]);
                nextpotinftimes(im) = fzeromin(funr,tbracket(im,1),opts1,1,[],tinf+guessintsize);
            else
                % If we have already computed a time for a potential 
                % infectious contact for the weakest suscpetible of type  
                % im, tbracket(im,2) is finite and equal to this time. It 
                % means this guy was supposed to be infected at this time,  
                % but someone else has been infected before. This guy will 
                % then certainly be infected, but we know it will be sooner 
                % than tbracket(im,2), because the additional contribution 
                % to the cumulative infection pressure of that someone else
                % can only erode the resilience of our guy faster.
                % Therefore, tbracket(im,2) is an upper bound, and I should
                % search in tbracket(im,:) = [ tinf, tbracket(im,2) ] using
                % fzero. HOWEVER, numerical errors might creep in, so I
                % still use fzeromin and use tbracket(im,2) as a starting 
                % guess, which fzeromin uses as the starting upper bound. 

                % To test/debug, uncomment next line:
%                 disp(['Search in [ ',num2str(tbracket(im,1)),', ',num2str(tbracket(im,2)),' ]']);
                % nextpotinftimes(im) = fzero(funr,tbracket(im,:),opts1);
                nextpotinftimes(im) = fzeromin(funr,tbracket(im,1),opts1,1,[],tbracket(im,2));
            end
            tbracket(im,2) = nextpotinftimes(im); % At each iteration, a 
            % set of potential times of infection has been computed for the
            % weakest susceptible of each type. The system will be updated
            % to the time of the first infection time. Then all the others 
            % need to be recalculated, but can only be smaller than the 
            % ones just calculated, which I keep as an upper bound for the
            % future time calculations.
        end
    end
    [tinf,index] = min(nextpotinftimes); % First event gives time and identity of the newly infected case
    z(index) = z(index) + 1; % Increase final size count in that group
    inftimes(index,z(index)) = tinf; % Store the new time of infection
    if ( z(index) == n(index) )
        nextQ(index) = Inf;
    else
        nextQ(index) = Qm(index,z(index)+1); % Update weakest individual to the next in that group
    end
    infpres = infpres + Lh1to1(:,index); % Update total infection pressure vector
    nextpotinfs = infpres > nextQ; % Recompute the logical vector of which type will see a new infection with the current infection pressure
    tbracket(:,1) = tinf; % Increase lowest bound of future times of infection to current time
    tbracket(index,2) = Inf; % Clear upper bound of next infection time, but only for the type of the individual just infected
end


    function y = cuminfpress(t,m,L,inftimes,alpha,gam)
    % This subroutine computes the cumulative infection pressure experienced by
    % a specific individual at time t, due to all individuals of each type,
    % accounting for their respective infectious times, for a gamma-shaped
    % infectivity profile.
    % 
    % Input:
    %   - t = the time at which the alleged susceptiple is challenged
    %   - m = number of types
    %   - L = row vector of 1-to-1 total infection pressures: this is the row
    %   of the matrix Lh1to1 corresponding to the type of the susceptible being
    %   challenged, so that the contribution of all types of infectors can be
    %   added up
    %   - inftimes = mxN matrix, where m is the number of types and N is the
    %   largest number of individuals of the same type: this will contain
    %   information only until time t
    %   - alpha = shape parameter of the gamma pdf
    %   - gam = scale parameter of the gamma pdf
    % 
    % Output:
    %   - y = the total infection pressure accumulated up to time t

    y = 0;
    for im = 1:m
        y = y + L(im) * sum(gammainc(gam(im)*(t-inftimes(im,inftimes(im,:)>=0)),alpha(im)));
    end



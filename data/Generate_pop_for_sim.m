% This script creates the household distribution files that I use for the 
% individual-based stochastic simulation written in C when running the
% model mapping procedue of Pellis, L. et al (2020), Nature Communications.
% 
% Input data:
%   - XX_H_structure_ModelMapping.txt: A text file with the table of 
%       household structure composition of the country considered (XX = GB 
%       for Great Britain, SL for Sierra Leone, SA for South Africa). The 
%       structure of this file mimics the form of Supplementary Tables 1, 4
%       and 6 of the paper. The file for GB is created in Excel directly
%       from data. Those for SL and SA by running the "ReadDHShhdata.m"
%       code.
% 
% Outputs:
%   - XX_H_structure_ModelMapping_sim_e5.txt: A text file similar to the
%       input one, to check how different the population structure used in
%       the simulation is from the one used in the actual mapping procedure
%   - H_GB_5.dat: A text file with the structure of schools (which are 
%       not used at all, but the simulation requires them - due to past
%       choices, in the simulation "H" stands for workplaces/schools and 
%       "W" for households). The file contains: the total number of
%       children, the total number of schools and the largest size (100),
%       and then a list with an index for each school and the sublist of
%       each of the children in the school.
%   - W_GB_5.dat: A text file with the structure of households used in the 
%       simulation. The file contains: the total number of adults and
%       children, the total number of households and the largest number of 
%       adults and children in a household, and finally a list with an 
%       index for each household, its size, and the sublist of adults 
%       (smaller indices) and children (largest indices) in each household,
%       with "-1" used to fill in the remaining gaps to the household size.
% 
% Note: No seed for the random number generator was set, so it's impossible
% to reproduce exactly the population structure used in the paper. However,
% the files used for the simulation are provided and any newly generated
% file is very similar.
% 
% Update: 21-01-2020

clearvars;
country = 'GB'; % Great Britain
% country = 'SL'; % Sierra-Leone
% country = 'SA'; % South-Africa

% Folder and workspace options
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    wrk_path = [base_dir,'\data\'];
    code_path = [base_dir,'\code-model-mapping\'];
else
    wrk_path = [base_dir,'/data/'];
    code_path = [base_dir,'/code-model-mapping/'];
end
cd(wrk_path);

% Load the file with the information about the household size composition
input_distr = [country,'_H_structure_ModelMapping.txt'];
H = load(input_distr);

tot = 100000; % Total population size: DO NOT change this!

[ max_inA,max_inC ] = size(H);
max_nA = max_inA - 1;
max_nC = max_inC - 1;
iH = 1;
nH = sum(sum(H>0)); 
Hcomp = NaN(nH,2);
HcompDistr = NaN(nH,1);
for inA = 1:max_inA
    for inC = 1:max_inC
        nA = inA-1;
        nC = inC-1;
        if H(inA,inC) > 0
            Hcomp(iH,:) = [ nA, nC ];
            HcompDistr(iH) = H(inA,inC);
            iH = iH+1;
        end
    end
end
assert(iH==nH+1);

% First normalise the distribution (in case it doesn't add up to 1) and
% then construct a cumulative distribution from it:
HcompCumDistr = distr_to_cumdistr(HcompDistr/sum(HcompDistr));
% Check whether the cumulative distribution arrives exactly at 1 or not:
err1 = 1 - sum(HcompDistr);
err2 = 1 - HcompCumDistr(end);
if ( err2 > 0 )
    warning([ 'The cumulative distribution is slightly defective and I had to add ', num2str(err2), ' to reach 1' ]);
elseif ( err2 < 0 )
    warning([ 'The cumulative distribution is slightly in excess and I removed ', num2str(-err2), ' to bring the largest value down to 1' ]);
end
HcompCumDistr(end) = 1; % If for any reason the distribution does not reach exactly 1, force it to.

N = 0;
NA = 0;
NC = 0;
i = 0;
indices = [];
while N < tot
    i = i + 1;
    j = draw_random_deviates_from_empirical_cum_distribution(HcompCumDistr);
    indices(i) = j;
    na = Hcomp(j,1);
    nc = Hcomp(j,2);
    hsize = sum(Hcomp(j,:));
    N = N + hsize;
    NA = NA + na;
    NC = NC + nc;
end
if N > tot
    N = N - hsize;
    NA = NA - na;
    NC = NC - nc;
    R = tot - N;
    j = draw_random_deviates_from_empirical_cum_distribution(HcompCumDistr);
    na = Hcomp(j,1);
    nc = Hcomp(j,2);
    hsize = sum(Hcomp(j,:));
    while hsize ~= R
        j = draw_random_deviates_from_empirical_cum_distribution(HcompCumDistr);
        na = Hcomp(j,1);
        nc = Hcomp(j,2);
        hsize = sum(Hcomp(j,:));
    end
    indices(i) = j;
    N = N + hsize;
    NA = NA + na;
    NC = NC + nc;
end
NH = i; % Number of households
assert(N==tot);

Cindices = randperm(NC) - 1;
filename = ['W_',country,'_',num2str(round(log10(tot))),'.dat'];
fileID = fopen(filename,'w');
fprintf(fileID,'%d\t%d\t%d\t%d\t%d\r\n',NA,NC,NH,max_nA,max_nC); % Replace "\r\n" with "\r" in Windows
currA = NC;
currCind = 1;
newH = zeros(max_inA,max_inC);
for i = 1:NH
    na = Hcomp(indices(i),1);
    nc = Hcomp(indices(i),2);
    newH(na+1,nc+1) = newH(na+1,nc+1) + 1;
    fprintf(fileID,'%d\t%d\t%d',i-1,na,nc);
    for j = 1:na
        fprintf(fileID,'\t%d',currA);
        currA = currA + 1;
    end
    for j = (na+1):max_nA;
        fprintf(fileID,'\t-1');
    end
    for j = 1:nc
        fprintf(fileID,'\t%d',Cindices(currCind));
        currCind = currCind + 1;
    end
    for j = (nc+1):max_nC;
        fprintf(fileID,'\t-1');
    end
    fprintf(fileID,'\r\n');
end
fclose(fileID);

newH = newH/NH;
newfilename = [country,'_H_structure_ModelMapping_sim_e',num2str(round(log10(tot))),'.txt'];
dlmwrite(newfilename,newH,'delimiter','\t');

filename = ['H_',country,'_',num2str(round(log10(tot))),'.dat'];
fileID = fopen(filename,'w');
Ssize = 100;
NS = ceil(NC/Ssize);
fprintf(fileID,'%d\t%d\t100\r\n',NC,NS);
curr = 0;
for i = 1:(NS-1)
    fprintf(fileID,'%d\t%d',i-1,Ssize);
    for j = 1:Ssize
        fprintf(fileID,'\t%d',curr);
        curr = curr + 1;
    end
    fprintf(fileID,'\r\n');
end
% Last "school"
lefttoplace = NC - curr;
leftovers = NS * Ssize - NC;
fprintf(fileID,'%d\t%d',NS-1,lefttoplace);
for j = 1:lefttoplace
    fprintf(fileID,'\t%d',curr);
    curr = curr + 1;
end
for j = 1:leftovers
    fprintf(fileID,'\t-1');
end
fprintf(fileID,'\r\n');
fclose(fileID);

% Compare summary statistics of input distirbution and distribution for simulation
cd(code_path);

F1 = create_population_composition(H);
H_single = create_1type_distr(H);
PI_single = create_1type_size_biased_distr(H_single);
muH1 = sum(H_single.*(1:length(H_single)));
sb_muH1 = sum(PI_single.*(1:length(PI_single)));

F2 = create_population_composition(newH);
H_single = create_1type_distr(newH);
PI_single = create_1type_size_biased_distr(H_single);
muH2 = sum(H_single.*(1:length(H_single)));
sb_muH2 = sum(PI_single.*(1:length(PI_single)));

% Print out some summary statistics of the population structure now
% generated ("New"), and compared with the results used for in the model 
% mapping paper ("Old"):
disp( ' ' )
disp( '[ F_a, F_c ]' )
disp( [ 'Old: ', num2str( F1 ) ] )
disp( [ 'New: ', num2str( F2 ) ] )
disp( ' ' )
disp( 'muH:' )
disp( [ 'Old: ', num2str( muH1 ) ] )
disp( [ 'New: ', num2str( muH2 ) ] )
disp( ' ' )
disp( 'sb_muH:' )
disp( [ 'Old: ', num2str( sb_muH1 ) ] )
disp( [ 'New: ', num2str( sb_muH2 ) ] )

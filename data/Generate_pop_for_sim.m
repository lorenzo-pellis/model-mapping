% This script creates the household distribution file that I use for the IB
% stochastic simulation in C. 

clearvars;
% country = 'GB';
% country = 'SL'; % Sierra-Leone
country = 'SA'; % South-Africa

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

input_distr = [country,'_H_structure_ModelMapping.txt'];
H = load(input_distr);

tot = 100000;

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

HcompCumDistr = distr_to_cumdistr(normalise_distribution(HcompDistr));

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
fprintf(fileID,'%d\t%d\t%d\t%d\t%d\r\n',NA,NC,NH,max_nA,max_nC); % \r for Windows
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
muH1 = sum(H_single.*[1:length(H_single)]);
sb_muH1 = sum(PI_single.*[1:length(PI_single)]);

F2 = create_population_composition(newH);
H_single = create_1type_distr(newH);
PI_single = create_1type_size_biased_distr(H_single);
muH2 = sum(H_single.*[1:length(H_single)]);
sb_muH2 = sum(PI_single.*[1:length(PI_single)]);

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

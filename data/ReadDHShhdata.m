% This script reads the DHS household data and creates the table describing
% the household size composition of SL and SA, which look like
% Supplementary Tables 4 and 6 of:
% Pellis, L. et al (2020), Nature Communications
% 
% Note: the DHS household data is not provided and needs to be obtained
% directly from www.dhsprogram.com (free registration). The files required
% are the "Household member recode" files:
%   - SLPR51FL.SAV for Sierra Leone, available at
% https://www.dhsprogram.com/data/dataset/Sierra-Leone_Standard-DHS_2008.cfm
%   - ZAPR31FL.SAV for South Africa, available at
% https://www.dhsprogram.com/data/dataset/South-Africa_Standard-DHS_1998.cfm
% Both files need to be converted using the code in file "DHSconvert.R" to
% .CSV files and then opened in Excel and saved as .xlsx files.
% 
% Input data:
%   - SLPR51FL.xlsx: An Excel file obtained as described above (i.e. not
%       provided here) for Sierra Leone.
%   - ZAPR31FL.xlsx: An Excel file obtained as described above (i.e. not
%       provided here) for South Africa.
% 
% Outputs:
%   - XX_H_structure_ModelMapping.txt: A text file with the table of 
%       household structure composition of the country considered (XX = SL 
%       for Sierra Leone, SA for South Africa). The structure of this file 
%       mimics the form of Supplementary Tables 4 and 6 of the paper. 
% 
% Update: 22-01-2020

% Switches:
useweights = true; % if true, I use weights as I should do according to DHS explanations if false, I discard weights.
discard_largehh = true; % if true, I just chop the distribution and renormalise; if false, I absorb hh larger than max into max
handlenan = 'C'; % What to do with NaN values? Exclude the house (E), put 
% case as child (C), put case as adult (A), give a random age (R).

oldestchild = 18; % Largest age of a child (included)
largesthh = 17; % Largest household size (included)

% country = 'GB'; % Great Britain
country = 'SL'; % Sierra Leone
% country = 'SA'; % South Africa

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

if strcmp(country,'SL');
    data = xlsread('SLPR51FL.xlsx','A2:E41986');
elseif strcmp(country,'SA')
    data = xlsread('ZAPR31FL.xlsx','A2:E52907');
else
    error('GB not working here...');
end
[ ld, nc ] = size(data);
ih = 0; % hh index
ch = data(1,2); % current household id in the data
cc = data(1,1); % current cluster (just to check there is no case where id remains the same but cluster changes)
ic = 1; % current line index
nanflag = false;
indhnew1 = find(data(:,1)~=cc,1);
indhnew2 = find(data(:,2)~=ch,1);
indhnew = min([indhnew1,indhnew2]);
Hcell = {}; % create an empty cell array
while ~isempty(indhnew)
    assert(data(indhnew-1,1)==cc);
    vec = data(ic:(indhnew-1),5);
    if any(isnan(vec))
        nanflag = true;
        switch handlenan
            case 'A'
                vec(isnan(vec)) = oldestchild+1;
            case 'C'
                vec(isnan(vec)) = oldestchild;
            case 'E'
                ic = indhnew;
                ch = data(ic,2);
                cc = data(ic,1);
                indhnew1 = find(data(ic:end,1)~=cc,1);
                indhnew2 = find(data(ic:end,2)~=ch,1);
                indhnew = ic - 1 + min([indhnew1,indhnew2]);
                break;
            case 'R'
                vec(isnan(vec)) = randi(100);
            otherwise
                error('Invalid option about how to handle NaNs');
        end
    end
    ih = ih+1;
    Hcell{ih,1} = ih;
    Hcell{ih,2} = data(ic,3); % the weight
    Hcell{ih,3} = data(ic,4); % the hh size
    Hcell{ih,4} = vec; % the age of the members
    assert((indhnew-ic)==length(Hcell{ih,4}));
    Hcell{ih,5} = nanflag;
    nanflag = false;
    Hcell{ih,6} = sum(Hcell{ih,4}<=oldestchild);
    Hcell{ih,7} = sum(Hcell{ih,4}> oldestchild);
    assert((Hcell{ih,6}+Hcell{ih,7})==Hcell{ih,3});
    
    ic = indhnew;
    ch = data(ic,2);
    cc = data(ic,1);
    indhnew1 = find(data(ic:end,1)~=cc,1);
    indhnew2 = find(data(ic:end,2)~=ch,1);
    indhnew = ic - 1 + min([indhnew1,indhnew2]);
    cc = data(ic,1);
end
nh = ih;

%%

Hmax = max([Hcell{:,3}]);
if discard_largehh
    if largesthh < Hmax
        Hmat = zeros(largesthh+1);

        if useweights
            for ih = 1:nh
                if Hcell{ih,3} <= largesthh
                    pos = [Hcell{ih,7},Hcell{ih,6}];
                    Hmat(pos(1)+1,pos(2)+1) = Hmat(pos(1)+1,pos(2)+1)+Hcell{ih,2};
                end
            end
            Hmat = Hmat / sum(sum(Hmat));
        else
            for ih = 1:nh
                if Hcell{ih,3} <= largesthh
                    pos = [Hcell{ih,7},Hcell{ih,6}];
                    Hmat(pos(1)+1,pos(2)+1) = Hmat(pos(1)+1,pos(2)+1)+1;
                end
            end
            Hmat = Hmat / sum(sum(Hmat));
        end
    end
else
    Hmat = zeros(Hmax+1);

    if useweights
        for ih = 1:nh
            pos = [Hcell{ih,7},Hcell{ih,6}];
            Hmat(pos(1)+1,pos(2)+1) = Hmat(pos(1)+1,pos(2)+1)+Hcell{ih,2};
        end
        Hmat = Hmat / sum(sum(Hmat));
    else
        for ih = 1:nh
            pos = [Hcell{ih,7},Hcell{ih,6}];
            Hmat(pos(1)+1,pos(2)+1) = Hmat(pos(1)+1,pos(2)+1)+1;
        end
        Hmat = Hmat / sum(sum(Hmat));
    end
end
if strcmp(country,'SL');
    dlmwrite('SL_H_structure_ModelMapping.txt',Hmat,'\t');
elseif strcmp(country,'SA')
    dlmwrite('SA_H_structure_ModelMapping.txt',Hmat,'\t');
else
    error('GB not working here');
end

cd(code_path);

F = create_population_composition(Hmat);
F_ratio = F(2)/F(1); % Ratio N_C / N_A. Careful: numerator and denominator are swapped with respect to g_ratio and h_ratio...
g_ratio = 1; thetaG = NaN;
thetaG_random_mix = 1 / ( 1 + g_ratio / F_ratio ); % This is the assortativity that should give random mixing: theta = N_C g_C / ( N_A g_A + N_C g_C )
if isnan(thetaG)
    thetaG = thetaG_random_mix;
end
PI = create_2type_size_biased_distr(Hmat);
H_single = create_1type_distr(Hmat);
PI_single = create_1type_size_biased_distr(H_single);
PI_type = create_type_biased_distr(PI);
muH = sum(H_single.*(1:length(H_single)));
sb_muH = sum(PI_single.*(1:length(PI_single)));
mu_type = sum(PI_type.*[1:length(H_single);1:length(H_single)],2);
muH1 = F*mu_type;
disp('Matrix for the household composition:')
disp(Hmat)
disp(['Fraction of adults and children:     ',num2str(F)])
disp(' ')
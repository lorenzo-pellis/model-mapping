% This script is the main code to run the model mapping procedure to
% compare model with age and households (AH) to models with age (A),
% households (H) and random mixing (U). Please refer to publication:
% Pellis, L. et al (2019), Nature Communications
%
% Date: 15-11-2019
% 
% Fixed parameters (uncomment applicable values in the code below):
%   - country: this selects the household composition structure (how adults
%           and children are distributed in households of different size)
%   - population: this selects the global assortativity (thetaG) and the
%           ratio of the number of contacts of children VS adults (g_ratio,
%           in the paper denoted by gamma). Global random mixing is
%           achieved for thetaG = fraction of children in the population,
%           so it's set as NaN to start with and computed based on the
%           population structure
%   - R0: the basic reproduction number
%   - phiG: the relative infectivity of children VS adults 

clearvars;

% Switches
Activate_checks = 1; % If 0, all cross check mechanisms are de-activated, to make everything a bit faster
Activate_C_codes = 0; % If 0, C codes call are de-activated
Activate_workspace_saving = 1; % If 0, the workspace is not saved automatically
Activate_continue = 1; % This allows continuing from the last run, if anything makes the code crash or you need to stop it
Activate_deletefile = 1; recycle('off'); % If 1 the code eletes the output file of the stochastic simulation after use
warning('off','MATLAB:nearlySingularMatrix');

% country = 'GB'; % Great Britain
% country = 'SL'; % Sierra Leone
country = 'SA'; % South Africa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time-independent parameters

pop = '2ran'; thetaG = NaN; g_ratio = 1; % Random (which for GB is thetaG = 0.2273)
% pop = 'm4r'; thetaG = 0.4; g_ratio = 1; 
% pop = 'm4UK'; thetaG = 0.4; g_ratio = 0.75; 
% pop = 'm5r'; thetaG = 0.5; g_ratio = 1; 
% pop = 'm5UK'; thetaG = 0.5; g_ratio = 0.75; 
% pop = 'UK'; thetaG = 0.58; g_ratio = 0.75; % UK
% pop = 'ass'; thetaG = 0.7; g_ratio = 0.75; % More extreme than UK

R0 = 1.5;
phiG = 1; % relative global infectivity of children versus adults

pAA_min = 0;
if strcmp(country,'SL')
    pAA_max = 0.475;
    dpAA = 0.025;
elseif strcmp(country,'SA')
    pAA_max = 0.63;
    dpAA = 0.033;
else
    pAA_max = 0.95;
    dpAA = 0.05;
end
pAA_vec = pAA_min:dpAA:pAA_max;
l1 = length(pAA_vec);
psiG_min = 1;
psiG_max = 4;
dpsiG = 0.2;
psiG_vec = psiG_min:dpsiG:psiG_max;
l2 = length(psiG_vec);

psiH_vec = psiG_vec; % assume the same susceptibility in households as in the community. Here it is a vector
phiH = phiG;
% thetaG = NaN; % Assortativity of children in the community. If NaN then it's random, i.e. I compute it after I computed F
thetaH = NaN; % Assortativity of children in the household. If NaN then it's random, i.e. it depends on the household size
% g_ratio = 1;%0.75; % Ratio of the contact rate gA / gC
h_ratio = g_ratio;%0.75; % Ratio of the contact rate hA / hC
c_ratio = g_ratio; % For the A model, I use the overall ratio of contacts
n_init_inf = 50;
mapfail = NaN; % This is the value to use if the mapping procedure fails in finding a suitable assortativity

% Real-time-related parameters

Tg = 2.85; % only to calculate r
% Tg corresponds to mean_profile, because we ust TVI model, i.e. distr_profile = 2
alpha_profile = 9; % only to calculate r
lambda_profile = alpha_profile / Tg;
Tmax = 50;
dt = 0.1;
tol = 0.01;

% Folder and workspace options
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-model-mapping\'];
    temp_path = [base_dir,'\code-model-mapping\temp\'];
    wrksp_path = [base_dir,'\output-workspaces\',country,'\'];
    runcommand = 'MRCModelMapping_Win.exe';
    tool_path = [base_dir,'\tools\'];
    check_path = [code_path,'\check-codes\'];
else
    code_path = [base_dir,'/code-model-mapping/'];
    temp_path = [base_dir,'/code-model-mapping/temp/'];
    wrksp_path = [base_dir,'/output-workspaces/',country,'/'];
    runcommand = './ModelMapping_Mac';
    tool_path = [base_dir,'/tools/'];
    check_path = [code_path,'/check-codes/'];
end
% Test immediately if both folders exist and you got the operating system right...
cd(wrksp_path);
cd(code_path);
if Activate_checks
    addpath(check_path);
end
if Activate_C_codes
    addpath(genpath(tool_path))
end

workspace_emergency_name = 'workspace_emergency_save';

% Population parameters

% H = probability distribution of households with each possible structure:
% H(i,j) has (i-1) adults and (j-1) children
if strcmp(country,'SL');
    input_distr = 'SL_H_structure_ModelMapping.txt';
elseif strcmp(country,'SA')
    input_distr = 'SA_H_structure_ModelMapping.txt';
else
    input_distr = 'GB_H_structure_ModelMapping.txt';
end
H = load(input_distr);

% Don't change the following
den_str = 'n-1'; % string containing the denominator I want to use to scale the 1-to-1 within-household infectivities
eta = 1; % exponent of the denominator in the freqency dependent transmission in households (n-1)^eta

%%%%%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saved variables
Rh_vec = zeros(1,l1); % Transmission probability from adult to adult in the household of a randomly chosen adult
Rg = zeros(l1,l2); % For each run, the 1-to-all total infectivities from j to i are LambdaG(i,j) = Rg * shape(i,j)
R0_G = zeros(l1,l2);
R0_AH = zeros(l1,l2);
v_AH = zeros(l1,l2,2);
vH_AH = zeros(l1,l2,2);
hfs_AH = zeros(l1,l2);
zAH = zeros(l1,l2);
zAHsim = zeros(l1,l2);
zAHsim_check = zeros(l1,l2);
tAHsim = zeros(l1,l2);
piAHsim = zeros(l1,l2);
afs_each_typeAH = zeros(l1,l2,2);
afs_each_typeAH_sim_check = zeros(l1,l2,2);
r_AH_best = zeros(l1,l2);
r_AH_high = zeros(l1,l2);
theta_A = zeros(l1,l2); % Assortativity of model A
Rg_A = zeros(l1,l2); % Global transmission of model A
zAt = zeros(l1,l2);
afs_each_typeAt = zeros(l1,l2,2);
r_A = zeros(l1,l2);
R0_A = zeros(l1,l2);
zA = zeros(l1,l2);
zAsim = zeros(l1,l2);
tAsim = zeros(l1,l2);
piAsim = zeros(l1,l2);
afs_each_typeA = zeros(l1,l2,2);
Rh_H = zeros(l1,l2);
hfs_H = zeros(l1,l2);
Rg_H = zeros(l1,l2);
r_H_best = zeros(l1,l2);
r_H_high = zeros(l1,l2);
R0_H = zeros(l1,l2);
zH = zeros(l1,l2);
zHsim = zeros(l1,l2);
zHsim_check = zeros(l1,l2);
tHsim = zeros(l1,l2);
piHsim = zeros(l1,l2);
NGM_G_stored = zeros(l1,l2,2,2);
NGM_At_stored = zeros(l1,l2,2,2);
NGM_A_stored = zeros(l1,l2,2,2);
hfs_by_ic = zeros(l1,l2,2);
av_obs_h_size = zeros(1,l2);
SAR_app = zeros(l1,l2);
SAR = zeros(l1,l2);
SAR_by_ic = zeros(l1,l2,2);
hfs_check = zeros(l1,l2);
inVSout = zeros(l1,l2,2);
inVSout_prop = zeros(l1,l2,2);
inVSout_A = zeros(l1,l2,2);
inVSout_A_prop = zeros(l1,l2,2);
diffA_H = zeros(l1,l2);
diffAH_H = zeros(l1,l2);
diffAH_A = zeros(l1,l2);
rel_diffA_H = zeros(l1,l2);
rel_diffAH_H = zeros(l1,l2);
rel_diffAH_A = zeros(l1,l2);
eflag = zeros(l1,l2); % Error flag

% Recycled variables
shapeG = zeros(2);
shapeH = zeros(2);
NGM_G_shape = zeros(2);
NGM_G = zeros(2);
NGM_At = zeros(2);
shape_A = zeros(2);
NGM_shape_A = zeros(2);
R1_temp = 0;
R2 = 0;
R3_temp = 0;
r_AH_temp = [];
r_H_temp = [];

% Population parameters initialisation:
F = [0,0]; % fractions of adults and children
% H (loaded before) = probability distribution of households with each possible structure: H(i,j) has i-1 adults and j-1 children
PI = []; % double matrix (one layer per type) for the size biased distribution of household sizes. PI_{a,c}^{(t)} = prob that the house of a randomly selected type t has a adults and c children
H_single = []; % vector of household size distribution. H_n = prob that a randomly selected household has size n
PI_single = []; % vector of size-biased household size distribution. PI_n = prob that the house of a randomly selected individual has size n
PI_type = []; % This gives the size distribution of the household of a randomly selected adult or child (it has 2 rows). PI_n^{t} = prob that the house of a r.s. type t has size n
muH = 0; % average size of a r.s. household
sb_muH = 0; % average size of the house of a r.s. individual
mu_type = 0; % average size of the house of a r.s. adult or child
muV = 0; % This is the average size of a household during the exp phase (it's not a property of the population, 
% it requires the dominant eigenvector [adults,children] as they occur during the exp growing phase)

%%%%%%%%%%%%%% Population structure %%%%%%%%%%%%%%%%%%%%%

F = create_population_composition(H);
F_ratio = F(2)/F(1); % Ratio N_c / N_a. Careful: numerator and denominator are swapped with respect to g_ratio and h_ratio...
thetaG_random_mix = 1 / ( 1 + g_ratio / F_ratio ); % This is the assortativity that should give random mixing: theta = N_C g_C / ( N_A g_A + N_C g_C )
if isnan(thetaG) % If I want random mixing, I initially select thetaG = NaN, and compute it here
    thetaG = thetaG_random_mix;
end
PI = create_2type_size_biased_distr(H);
H_single = create_1type_distr(H);
PI_single = create_1type_size_biased_distr(H_single);
PI_type = create_type_biased_distr(PI);
muH = sum(H_single.*[1:length(H_single)]);
sb_muH = sum(PI_single.*[1:length(PI_single)]);
mu_type = sum(PI_type.*[1:length(H_single);1:length(H_single)],2);
maxHsize = length(H_single);

% Pre-compilation of a table of binomial coefficients.
% I need them for future calculations, but they are expensive to compute,
% so I will do it now once and store such values in a table. Such table is
% a global variable, so I can access it from all functions that need it.
global bincoeffmat;
[ max_inA,max_inC ] = size(H);
bincoeffmat = NaN(size(H));
for inA = 1:max_inA
    for inC = 1:max_inC
        nA = inA-1;
        nC = inC-1;
        if nA >= nC
            bincoeffmat(inA,inC) = nchoosek(nA,nC);
        end
    end
end

% % In case it is needed, I can create a unique list of all household
% % compositions. However I do not use it in this version of the code
% iH = 1;
% nH = sum(sum(H>0)); 
% Hcomp = NaN(nH,2);
% HcompDistr = NaN(nH,1);
% HcompDistrBias = NaN(nH,2);
% for inA = 1:max_inA
%     for inC = 1:max_inC
%         nA = inA-1;
%         nC = inC-1;
%         if H(inA,inC) > 0
%             Hcomp(iH,:) = [ nA, nC ];
%             HcompDistr(iH) = H(inA,inC);
%             HcompDistrBias(iH,:) = Hcomp(iH,:) * HcompDistr(iH);
%             iH = iH+1;
%         end
%     end
% end
% assert(iH==nH+1);
% HcompDistrBias(:,1) = HcompDistrBias(:,1) / sum( HcompDistrBias(:,1) );
% HcompDistrBias(:,2) = HcompDistrBias(:,2) / sum( HcompDistrBias(:,2) );


%%%%%%%%%%%%%% Population structure used in the simulation (pop = 1e5)
% When I run my individual-based stochastic simulations, I have to generate
% a synthetic population of households (the same one is used for all
% simulations). Because it's randomly generated, stochastic variations
% might lead to slightly different summary description of the population
% from those that I derived from data and use in this script. Here I
% calculate the same quantities as above, but for the population used in
% the simulation to check potential discrepancies.
if strcmp(country,'SL')
    Hsim = load('SL_H_structure_ModelMapping_sim_e5.txt');
elseif strcmp(country,'SA')
    Hsim = load('SA_H_structure_ModelMapping_sim_e5.txt');
else
    Hsim = load('GB_H_structure_ModelMapping_sim_e5.txt');
end
PIsim = create_2type_size_biased_distr(Hsim);
H_single_sim = create_1type_distr(Hsim);
PI_single_sim = create_1type_size_biased_distr(H_single_sim);
PI_type_sim = create_type_biased_distr(PIsim);
muH_sim = sum(H_single_sim.*[1:length(H_single_sim)]);
sb_muH_sim = sum(PI_single_sim.*[1:length(PI_single_sim)]);
mu_type_sim = sum(PI_type_sim,2); % Just to check

% The output of the mapping procedure (this script) is saved in a workspace
% with the following name:
workspace_name = [ 'R0', num2str(R0*10,'%02d'),...
    '_pAA', num2str(ceil(pAA_min*100),'%02d'),'_',num2str(ceil(dpAA*100),'%02d'),'_',num2str(ceil(pAA_max*100),'%02d'),...
    '_psi', num2str(psiG_min*10,'%02d'), '_', num2str(dpsiG*10,'%02d'),'_', num2str(psiG_max*10,'%02d'),...
    '_phi', num2str(phiG*10,'%02d'), '_theta',num2str(round(thetaG*10),'%02d'),...
    '_gammaG', num2str(g_ratio*100,'%03d'), '_H', num2str(h_ratio*100,'%03d') ];
if Activate_C_codes % It takes much longer if I also compute peak incidence and time to the peak from simulations. If I do, the file name shows it
    workspace_name = [ workspace_name, '_100sim_e5_init', num2str(n_init_inf,'%03d') ];
end
workspace_name_temp = [ workspace_name, '_temp' ];
if Activate_continue % Trick to avoid recalculationg everything if my code has crashed half-way through
    if ~exist([temp_path,workspace_name_temp,'.mat'],'file')
        for countdown = 5:-1:1
            disp(['No previous run attempt with these parameters was found. Start the run from scratch in ',num2str(countdown),' seconds...']);
            pause(1);
        end
        Activate_continue = false;
        flag_entry_loop = false; % Indicate I want to skip all loops already calculated
        i1_temp = l1;
        i2_temp = l2;
        TotTime = 0;
    else
        load([temp_path,workspace_name_temp]); % Load workspace saved when the code failed
        for countdown = 5:-1:1
            disp(['Previous run attempt found: continue from loop ',num2str(i2),'-',num2str(i1),' of ',num2str(l2),'-',num2str(l1),' in ',num2str(countdown),' seconds...']);
            pause(1);
        end
        Activate_continue = true;
        flag_entry_loop = true; % Indicate I want to skip all loops already calculated
        i1_temp = i1;
        i2_temp = i2;
    end
else
    flag_entry_loop = false;
    i1_temp = l1;
    i2_temp = l2;
    TotTime = 0;
end

tStart = tic;
SingleRun = tic;

for i2 = 15:l2 % External loop is for psi (so figures in the paper are computed in rows)
    if Activate_continue && ( i2 < i2_temp )
    else
        psiG = psiG_vec(i2);
        psiH = psiH_vec(i2);

        % Construct the shape of the global NGM. Then you only need to multiply
        % it by a coefficient
        NGM_G_shape = create_global_NGM_shape(F_ratio,psiG,phiG,thetaG,g_ratio);

        % Note: I wanted to allow for assortativity in the household, but it
        % turns out the shape depends on the household composition...
        % So I use thetaH = NaN for random mixing within each household 
        % (computed differently for each household composition).
        % To pass all the information I need for when I have to compute the
        % assortativity of each household composition, I store all relevant
        % parameters in a vector:
        shapeHvector = [ h_ratio, thetaH, psiH, phiH ];

        for i1 = 10:l1 % Internal loop is for amount of within-household adult-to-adult transmission probability
            if flag_entry_loop && Activate_continue && ( i1 <= i1_temp )
            else
                disp(['Loop ',num2str(i2),'-',num2str(i1),' of ',num2str(l2),'-',num2str(l1)]);
                pAA = pAA_vec(i1);

                %%%%%%%%%%%%%%%%%% Find the right Rh for this pAA %%%%%%%%%%%%%%%
                if pAA == 0
                    Rh = 0;
                else
                    funpAA = @(x) ( get_pAA_from_Rh(H_single,x) - pAA );
                    if Activate_checks
                        Rh1 = funpAA(0);
                        Rh2 = funpAA(100);
                        if ( sign(Rh1) == sign(Rh2) )
                            Rh = NaN;
                            eflag(i1,i2) = eflag(i1,i2) + 1;
                            error('Cannot get the required Rh. Abort everything!')
                        else
                            Rh = fzero( funpAA, [0 100] );
                        end
                    else
                        Rh = fzeromin( funpAA, 0, [], [], [], 1 );
                    end
                end
                Rh_vec(i1) = Rh;

                % Just to check Rh has been found correctly (I had issues with using the wrong distribution, so Rh and pAA were not 1-to-1 related)
                if Activate_checks
                    TP = get_pAA_from_Rh(H_single,Rh);
                    if ( TP - pAA > 1e-6 )
                        warning([ 'Careful! Checked value of pAA failed. [ pAA, TP ] = [ ',num2str(pAA),', ',num2str(TP),' ].'] );
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% AH model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                disp('Start AH model');
                % Construct the chains of average number of cases in each
                % generation (it's a 2 x Hsize x 2 tensor). See inside of
                % function for details. The second and third outputs are
                % not very important, but it's cheaper to compute them 
                % all in one go. av_gen_chains is key for computing R0
                [ av_gen_chains_2type, hfs_by_comp_by_ic, SAR_by_comp_by_ic ] = create_chains_hfs_and_SAR_2type(PI,Rh,shapeHvector,eta,bincoeffmat); % Here is where I decide how the infectivity scales with the household size

                if Activate_checks
                    % A consistency check: the sum of the average chains of
                    % infectives in each generation should give the same
                    % result as the average final size of a
                    % within-household epidemic, which I can compute using
                    % Addy, Longini & Haber (1991), Biometrics
                    check1 = zeros(2);
                    check2 = zeros(2);
                    for bb = 1:2 % First start with 1 initial adult, then 1 initial child
                        init_cases = zeros(1,2);
                        init_cases(bb) = 1;
                        check1(:,bb) = afsn_Addy_2type_distr_initinf_check(PI(:,:,bb),Rh,shapeHvector,init_cases,[1;1],eta);
                        for aa = 1:2
                            check2(aa,bb) = sum( av_gen_chains_2type(aa,:,bb) );
                        end
                    end
                    check1 = check1 + eye(2); % Still to add the initial infectives
                    diff12 = abs(check1-check2);
                    if any( diff12 >= 1e-6 )
                        warning('Possible errors in computing the chains of infectives: the total doesn''t match!');
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                end

                %%% Find the Rg needed to achieve the selected R0.
                % shape = diag( F.^(-1) ) * NGM_shape; % This is the key relationship between shape and NGM_shape
                funAH = @(x) ( R0AH(av_gen_chains_2type,x*NGM_G_shape) - R0 );
                if Activate_checks
                    x1 = funAH(0);
                    x2 = funAH(100);
                    if ( sign(x1) == sign(x2) )
                        Rg(i1,i2) = NaN;
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                        error('Cannot get the required R0. Abort everything!')
                    else
                        Rg(i1,i2) = fzero( funAH, [0 100] );
                    end
                else
                    Rg(i1,i2) = fzeromin( funAH, 0, [], [], [], 1 );
                end
                %%% We have found Rg from R0
                %%% Find the Rg needed to achieve the selected r
%                 Rg(i1,i2) = find_Rg_from_r_2typeH_varHsize_Sellke(PI,desired_r,NGM_G_shape,alpha_profile*ones(2),lambda_profile*ones(2),Rh,shapeHvector,eta,alpha_profile*ones(2),lambda_profile*ones(2),nsimr);
%                 disp(['Computation of Rg from r for loop ',num2str(i2),'-',num2str(i1),' of ',num2str(l2),'-',num2str(l1)]);
                %%% We have found Rg from r
                
                NGM_G = Rg(i1,i2) * NGM_G_shape;
                NGM_G_stored(i1,i2,:,:) = NGM_G;
                % Now that we have found the right Rg, recalculate R0 to
                % confirm it coincides with what we want:
                R0_AH(i1,i2) = R0AH(av_gen_chains_2type,Rg(i1,i2)*NGM_G_shape);

                % Contribution of within-household VS between-household transmission
                R0_G(i1,i2) = max(eig(NGM_G)); % This is the R0 only due to global infectious contacts
                inVSout(i1,i2,1) = R0_AH(i1,i2) - R0_G(i1,i2);
                inVSout(i1,i2,2) = R0_G(i1,i2);
                inVSout_prop(i1,i2,:) = inVSout(i1,i2,:) / R0_AH(i1,i2); % Proportion of transmission inside the household compared to outside

                % Calculate the final size for AH, zAH(i1,i2);
                tic; % Time this, because it's a potentially expensive process
                tempAH = afs_2typeH_varHsize(H,PI,F,NGM_G,Rh,shapeHvector,eta); % This is still quite inefficient, but I'm not sure I can improve it
                toc;
                if Activate_checks
                    tic;
                    tempAHcheck = afs_2typeH_varHsize_check(H,PI,F,muH,NGM_G,Rh,shapeHvector,eta); % This is still quite inefficient, but I'm not sure I can improve it
                    toc;
                    if sum(tempAH-tempAHcheck)>1e-10
                        warning('Careful, the two methods to compute the final size zAH are giving different answers...');
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                end
                afs_each_typeAH(i1,i2,:) = tempAH;
                zAH(i1,i2) = F * tempAH;
                if ( sum(zAH(i1,i2)) < 0.0001 )
                    zAH(i1,i2) = 0;
                end

                if Activate_checks % Check the discrepancy in final size if I use Hsim, rather than H for the population composition
                    tempAHsim_check = afs_2typeH_varHsize_check(Hsim,PI,F,muH,NGM_G,Rh,shapeHvector,eta); % This is still quite inefficient, but I'm not sure I can improve it
                    afs_each_typeAH_sim_check(i1,i2,:) = tempAHsim_check;
                    zAHsim_check(i1,i2) = F * tempAHsim_check;
                    if ( sum(zAHsim_check(i1,i2)) < 0.0001 )
                        zAHsim_check(i1,i2) = 0;
                    end
                end

                %%%%%%%%%%%%%%%%%%%% computing v_AH %%%%%%%%%%%%%%%%%%%%%%%

                [ R1, v_inc ] = R0AH_and_vAH(av_gen_chains_2type,NGM_G);
                assert(abs(R0-R1)<1e-10); % Check the computed R0 is what it should be
                if Activate_checks
                    [ R2, v2_inc, NGM_At ] = R0AH_vAH_and_NGM_check(av_gen_chains_2type,NGM_G);
                    if ( R0 - R2 >= 1e-6 ) % bypass numerical errors
                        warning('Problem: R0 and R2 should be the same!')
                        disp([ R0, R2 ])
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                    if any( abs( v_inc - v2_inc ) >= 1e-6 )
                        warning('Problem in computing v_AH: the two approaches don''t match!')
                        disp([ v_inc, v2_inc ])
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                    NGM_At_stored(i1,i2,:,:) = NGM_At; % Store the "true NGM", but I'm not sure it's right
                    % In the saved workspaces for the paper, I had actually
                    % stored the "true NGM", but I'm not convinced it's 
                    % right. In any case, this output is not used.
                    [ R3, v3_inc ] = get_dominant_eigenpair( NGM_At ); % However, NGM_At has the correct R0 and dominant eigenvector v_inc
                    if ( R0 - R3 >= 1e-6 ) % bypass numerical errors
                        warning('Problem: R0 and R3 should be the same!')
                        disp([ R0, R3 ])
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                    if any( abs( v_inc - v3_inc ) >= 1e-6 )
                        warning('Problem in computing v_AH: the true NGM approach seems to fail!')
                        disp([ v_inc, v3_inc ])
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                end
                v_AH(i1,i2,1) = v_inc(1);
                v_AH(i1,i2,2) = v_inc(2);
                % Proportions of adults/children among the household primary cases:
                vH_inc = NGM_G * v_inc / sum( NGM_G * v_inc ); 
                vH_AH(i1,i2,1) = vH_inc(1);
                vH_AH(i1,i2,2) = vH_inc(2);
                % Average generation chains we expect to see during the exponentially growing phase:
                av_H_gen = sum( av_gen_chains_2type(:,:,1) * vH_inc(1) + av_gen_chains_2type(:,:,2) * vH_inc(2) );
                % Average final size in a within-household epidemic 
                % (household final size, or "hfs") computed during the
                % exponentially growing phase (to be compared with 
                % av_obs_h_size (see below) = 3.0677 in random mixing and 
                % psi=phi=1) during the exp growing phase (no multiple 
                % reintorductions) and is needed to compute Rh_H later on):
                hfs_AH(i1,i2) = sum( av_H_gen ); 
                hfs_by_ic(i1,i2,1) = sum( sum( av_gen_chains_2type(:,:,1) ) ); % hsf by initial case, to be compared with mu^a = 2.73
                hfs_by_ic(i1,i2,2) = sum( sum( av_gen_chains_2type(:,:,2) ) ); % hsf by initial case, to be compared with mu^c = 4.20
                if Activate_checks
                    hfs_check(i1,i2) = hfs_by_ic(i1,i2,1) * vH_inc(1) + hfs_by_ic(i1,i2,2) * vH_inc(2);
                    if ( abs( hfs_AH(i1,i2) - hfs_check(i1,i2) ) > 10e-7 )
                        warning('Careful! hfs_AH and hfs_check values are different')
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                end
                % Secondary attack rates (SAR) are the same as household
                % final sizes (hfs), but expressed in terms of proportions.
                % However, the proportion needs ot be calculated for each
                % household composition and then averaged across the
                % household composition distribution (hence why I computed
                % the SAR before, inside the same function for the chains)
                SAR_by_ic(i1,i2,1) = sum( sum( SAR_by_comp_by_ic(:,:,1) .* PI(:,:,1) ) );
                SAR_by_ic(i1,i2,2) = sum( sum( SAR_by_comp_by_ic(:,:,2) .* PI(:,:,2) ) );
                SAR(i1,i2) = SAR_by_ic(i1,i2,1) * vH_inc(1) + SAR_by_ic(i1,i2,2) * vH_inc(2);
                % This is the household size distribution of a randomly
                % chosen primary case during the exponential phase (recall
                % that PI_type is the size distribution of the household a 
                % randomly selected adult/child):
                PI_V = vH_inc' * PI_type;
                % Observed average household size of an infected household
                % during the exponentially growing phase:
                av_obs_h_size(i2) = sum( PI_V .* (1:length(PI_V)) ); 
                % Prepare the initial infectives for the simulation, for
                % model AH and model A:
                n_init_inf_A_AH = round( vH_inc(1) * n_init_inf );
                n_init_inf_C_AH = n_init_inf - n_init_inf_A_AH;
                n_init_inf_A_A = round( v_inc(1) * n_init_inf );
                n_init_inf_C_A = n_init_inf - n_init_inf_A_A;

                %%%%%%%%%% A little extra stuff, not really relevant for official results
                % Check the approsimate value of r. Because the generation 
                % time distribution is the same for adults and children, 
                % then I can ignore the fact that I have 2 types
                r_AH_temp = rH_TVI_tol(av_H_gen,R0_G(i1,i2),2,Tg,alpha_profile,Tmax,dt,tol); % This returns two values hopefully bracketing r
                r_AH_best(i1,i2) = r_AH_temp(1); % Best estimate (typically a lower bound, but not always)
                r_AH_high(i1,i2) = r_AH_temp(2); % Higher bound

                if Activate_checks
                    % Calculate the final size for the A model using the 
                    % true NGM. First, not sure the NGM_At is computed
                    % correctly. Second, even if correct, there is no
                    % reason the predicted final size should be the same!
                    tempAt = afs_multitype_NGM(F,NGM_At);
                    afs_each_typeAt(i1,i2,:) = tempAt;
                    zAt(i1,i2) = F * tempAt;
                    
                    % Proportion on transmission within the household 
                    % compared to outside, among adults only - not sure
                    % this is correct, as it relies on NGM_At, but in the
                    % workspaces of the paper I have saved this. However,
                    % it's not used for any results.
                    inVSout_A(i1,i2,1) = NGM_At(1,1) - NGM_G(1,1);
                    inVSout_A(i1,i2,2) = NGM_G(1,1);
                    inVSout_A_prop(:,:,1) = inVSout_A(:,:,1) ./ sum(inVSout_A,3);       
                    inVSout_A_prop(:,:,2) = inVSout_A(:,:,2) ./ sum(inVSout_A,3);       
                end

                %%%%%%%%%%%% Call C code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if Activate_C_codes % The stochastic simulations are in C
                    disp('Run AH simulations');
                    simtime = tic; % Check how long they take
                    cmdlAH = [ runcommand,' param.txt ',country,' ', num2str(Rg(i1,i2),'%.3f'), ' ', num2str(Rh,'%.3f'), ' ', num2str(psiG,'%.1f'), ' ', ...
                         num2str(phiG,'%.1f'), ' ', num2str(thetaG,'%.3f'),' ', num2str(g_ratio,'%.2f'),' ', num2str(h_ratio,'%.2f'),' ',...
                         num2str(n_init_inf_A_AH,'%d'),' ', num2str(n_init_inf_C_AH,'%d') ];
                    [status result] = system( cmdlAH ); % cmdlAH is the command line to run the executable with the right arguments
                    % status is 0 if the command is executed correctly
                    % result contains what the C code spits out in the standard output (not used here)

                    output_file_name_AH = [ country,'_Rg', num2str(Rg(i1,i2),'%.3f'), '_Rw', num2str(Rh,'%.3f'), '_sigma', num2str(psiG,'%.1f'), ...
                        '_rho', num2str(phiG,'%.1f'), '_ass', num2str(thetaG,'%.3f'), '_gammaG', num2str(g_ratio,'%.2f'),'_H', num2str(h_ratio,'%.2f'),...
                        '__averages.dat' ]; % The simulation creates a data file (.dat, but it's just a text file) with this name
                    [labels,HOW_MANY,data] = readColData( output_file_name_AH, 16, 0, 1 ); % Read the Excel file with a function written by someone else.
                    zAHsim(i1,i2) = data(3)/100; % Average final size from the simulation (just to cross-check the analytical result). Simulation gives percentage, now turned in a fraction
                    tAHsim(i1,i2) = data(7); % Time to the peak
                    piAHsim(i1,i2) = data(5); % Peak incidence, expressed in percentages
                    if Activate_checks && ( abs( zAH(i1,i2) - zAHsim(i1,i2) ) > 0.01 ) % Check analytical and simulated final size are not too different from each other
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                        warning('Careful! Analytical and simulated final size does not coincide for model AH');
                    end
                    tsims = toc(simtime);
                    disp(['Done in ',datestr(tsims/86400,'SS,FFF'),' seconds']);
                end

        %%%%%%%%%%%%%%%%%%%%%%%%% A model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                disp('Start A model');
                % We assume psi and phi are correct both in the household and in
                % the community, but we need to estimate theta_A. This can be an
                % issue when things are perfectly "symmetric" or when the matrix is
                % reducible. Check we are not in a strange case:
                if Activate_checks
                    theta_A(i1,i2) = find_assortativity(v_inc,F_ratio,psiG,phiG,g_ratio,pAA,mapfail,true); % Request a plot
                else
                    if ( psiG == 1 && phiG == 1 ) % In this case a valid assortativity can never be found.
                        theta_A(i1,i2) = thetaG; % Return the specified globala ssortativity for model AH
                    else
                        theta_A(i1,i2) = find_assortativity(v_inc,F_ratio,psiG,phiG,g_ratio,pAA,mapfail,false);
                    end
                end                

                if ( theta_A(i1,i2) == mapfail ) || ( isnan(mapfail) && isnan(theta_A(i1,i2)) ) % If no suitable assortativity can be found, keep track of this
                    Rg_A(i1,i2) = mapfail;
                    NGM_A = mapfail;
                    NGM_A_stored(i1,i2,:,:) = mapfail;
                    r_A(i1,i2) = mapfail;
                    R0_A(i1,i2) = mapfail;
                    afs_each_typeA(i1,i2,:) = mapfail;
                    zA(i1,i2) = mapfail;
                    zAsim(i1,i2) = mapfail;
                    tAsim(i1,i2) = mapfail;
                    piAsim(i1,i2) = mapfail;
                else
                    % Construct the shape of the NGM matrix for model A
                    NGM_shape_A = create_global_NGM_shape(F_ratio,psiG,phiG,theta_A(i1,i2),g_ratio);

                    %%% Calculate Rg_A, the total global infectivity for model A, from R0
                    funA = @(b) ( max(eig(b*NGM_shape_A)) - R0 );
                    if Activate_checks
                        b1 = funA(0);
                        b2 = funA(100); % It's unlikely Rg_A > 100, but this bit is only for cross-checking when debugging
                        if ( sign( b1 ) == sign( b2 ) )
                            Rg_A(i1,i2) = NaN;
                            eflag(i1,i2) = eflag(i1,i2) + 1;
                            error('Cannot get the required R0 for model A. Abort mapping!');
                        else
                            Rg_A(i1,i2) = fzero( funA, [0,100] );
                        end
                    else % When not debugging, use fzeromin, which searches a bit better and tries much larger values than 100, if needed
                        Rg_A(i1,i2) = fzeromin( funA, 0, [], [], [], 1 );
                    end
                    %%% We have found Rg_A from R0                    
                    %%% Calculate Rg_A, the total global infectivity for model A, from r
%                     Rg_A(i1,i2) = ( ( 1 + desired_r / lambda_profile )^alpha_profile ) / max(eig(NGM_shape_A));
                    %%% We have found Rg_A from r                    

                    % Finally obtain the NGM for model A
                    NGM_A = Rg_A(i1,i2) * NGM_shape_A;
                    NGM_A_stored(i1,i2,:,:) = NGM_A;
                    if NGM_A(1,1) < 0 % I found a value for the assortativity, but it leads to negative adult-to-adult reproduction number
                        r_A(i1,i2) = mapfail;
                        R0_A(i1,i2) = mapfail;
                        afs_each_typeA(i1,i2,:) = mapfail;
                        zA(i1,i2) = mapfail;
                        zAsim(i1,i2) = mapfail;
                        tAsim(i1,i2) = mapfail;
                        piAsim(i1,i2) = mapfail;
                    else
                        % Calculate the real-time growth rate for model A
                        r_A(i1,i2) = r_TVI(R0,2,Tg,alpha_profile,Tmax,dt);
                        R0_A(i1,i2) = max(eig(NGM_A)); % This should be the same as R0

                        % Calculate the final size for model A using the NGM of model A
                        tempA = afs_multitype_NGM(F,NGM_A);
                        afs_each_typeA(i1,i2,:) = tempA;
                        zA(i1,i2) = F * tempA;

                        %%%%%%%%%%%% Call C code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        if Activate_C_codes % The stochastic simulations are in C
                            disp('Run A simulations');
                            simtime = tic; % Check how long they take
                            cmdlA = [ runcommand,' param.txt ',country,' ', num2str(Rg_A(i1,i2),'%.3f'), ' ', num2str(0,'%.3f'), ' ', num2str(psiG,'%.1f'), ' ', ...
                                 num2str(phiG,'%.1f'), ' ', num2str(theta_A(i1,i2),'%.3f'),' ', num2str(c_ratio,'%.2f'),' ', num2str(1,'%.2f'),' ',...
                                 num2str(n_init_inf_A_A,'%d'),' ', num2str(n_init_inf_C_A,'%d')  ];
                            [status result] = system( cmdlA ); % cmdlA is the command line to run the executable with the right arguments
                            % status is 0 if the command is executed correctly
                            % result contains what the C code spits out in the standard output (not used here)

                            output_file_name_A = [ country,'_Rg', num2str(Rg_A(i1,i2),'%.3f'), '_Rw', num2str(0,'%.3f'), '_sigma', num2str(psiG,'%.1f'), ...
                                '_rho', num2str(phiG,'%.1f'), '_ass', num2str(theta_A(i1,i2),'%.3f'), '_gammaG', num2str(c_ratio,'%.2f'),'_H', num2str(1,'%.2f'),...
                                '__averages.dat' ]; % The simulation creates a data file (.dat, but it's just a text file) with this name
                            [labels,HOW_MANY,data] = readColData( output_file_name_A, 16, 0, 1 ); % Read the Excel file with a function written by someone else.
                            zAsim(i1,i2) = data(3)/100; % Average final size from the simulation (just to cross-check the analytical result). Simulation gives percentage, now turned in a fraction
                            tAsim(i1,i2) = data(7); % Time to the peak
                            piAsim(i1,i2) = data(5); % Peak incidence, expressed in percentages
                            if Activate_checks && ( abs( zA(i1,i2) - zAsim(i1,i2) ) > 0.01 ) % Check analytical and simulated final size are not too different from each other
                                warning('Careful! Analytical and simulated final size does not coincide for model A');
                                eflag(i1,i2) = eflag(i1,i2) + 1;
                            end
                            tsims = toc(simtime);
                            disp(['Done in ',datestr(tsims/86400,'SS,FFF'),' seconds']);
                        end
                    end
                end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% H model %%%%%%%%%%%%%%%%%%%%%%%%%

                disp('Start H model');
                % Calculate Rh_H such that the average final size of a household
                % epidemic during the exp growing phase is the same as hfs_AH
                % But we need to use the household size distribution we observe 
                % during the exp growing phase: PI_V, not PI
                fun_hfs = @(q) ( hfs_AH(i1,i2) - afsn_Ball_sizebiased_den(PI_V,q,0,1,1,eta,den_str) );
                if Rh == 0
                    Rh_H(i1,i2) = 0;
                else
                    if Activate_checks
                        q1 = fun_hfs(0);
                        q2 = fun_hfs(10);
                        if ( sign(q1) == sign(q2) ) % Check if I can find a solution before searching
                            Rh_H(i1,i2) = NaN;
                            eflag(i1,i2) = eflag(i1,i2) + 1;
                            error('Cannot get the required final size in households during the exponentially growing phase. Abort mapping!')
                        else
                            Rh_H(i1,i2) = fzero( fun_hfs, [0 10] );
                        end
                    else
                        Rh_H(i1,i2) = fzeromin( fun_hfs, 0, [], [], [], 1 );
                    end
                end
                if Activate_checks
                    hfs_H(i1,i2) = afsn_Ball_sizebiased_den(PI_V,Rh_H(i1,i2),0,1,1,eta,den_str);
                    if abs( hfs_AH(i1,i2) - hfs_H(i1,i2) ) > 1e-6
                        warning('Careful! The average household epidemic size in models AH and H is different!')
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                    hfs_H_check = afsn_Addy_sizebiased_den(PI_V,1,Rh_H(i1,i2),0,1,1,eta,den_str);
                    if abs( hfs_H_check - hfs_H(i1,i2) ) > 1e-6
                        warning( 'Careful! The average household epidemic size of model H computed in two different ways do not coincide!' )
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end                
                    len = length(PI_V);
                    PI_V_2up = zeros(1,len);
                    PI_V_2up(2:len) = PI_V(2:len) / sum(PI_V(2:len));
                    hfs_by_size = zeros(1,len);
                    for n = 1:len
                        hfs_by_size(n) = afsn_Ball(n,1,Rh_H(i1,i2)/(n-1),0,1,1);
                    end
                    hfs_H_check2 = PI_V * hfs_by_size';
                    if abs( hfs_H_check2 - hfs_H(i1,i2) ) > 1e-6
                        warning( 'Careful! The average household epidemic size of model H computed in two different ways do not coincide!' )
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end                
                end
                
                % Given Rh_H, compute the average number of cases in each generation
                av_gen_chain_1type = create_chain_1type(PI_single,Rh_H(i1,i2),1,'n-1');

                %%% Find Rg_H for model H from R0
                funR0H = @(r) ( R0H(av_gen_chain_1type,r) - R0 ); % r = Rg
                if Activate_checks
                    r1 = funR0H(0);
                    r2 = funR0H(10);
                    if ( sign(r1) == sign(r2) )
                        Rg_H(i1,i2) = NaN;
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                        error('Cannot get the required R0 for model H. Abort mapping!')
                    else
                        Rg_H(i1,i2) = fzero( funR0H, [0 10] );
                    end
                    % Try computing R0H in different ways
                    funR0H_check = @(r) ( R0H_varHsize_check(PI_single,Rh_H(i1,i2),r) - R0 ); % Aiming for r = Rg
                    Rg_H_check = fzero( funR0H_check, 1 );
                    if abs( Rg_H(i1,i2) - Rg_H_check ) > 1e-10
                        warning( 'Careful! The values of R0 for model H computed in two different ways do not coincide!' )
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end                
                else
                    Rg_H(i1,i2) = fzeromin( funR0H, 0, [], [], [], 1 ); % Slightly more robust function if I don't know a good interval bracketing the solution
                end
                %%% We have found Rg_H from R0                                    
                %%% Find Rg_H for model H from r
%                 Rg_H(i1,i2) = find_Rg_from_r_1typeH_varHsize_Sellke(PI_single,desired_r,alpha_profile,lambda_profile,Rh_H(i1,i2),alpha_profile,lambda_profile,eta,den_str,nsimr);
                %%% We have found Rg_H from r                    
                
                % Calculate the final size for model H
                zH(i1,i2) = afsH_varHsize_Addy_den(PI_single,Rg_H(i1,i2),Rh_H(i1,i2),0,1,1,1,'n-1');

                if Activate_checks
                    if isnan(zH(i1,i2))
                        error( 'Could not compute the final size for model H correctly! Abort mapping' );
                    else
                        zH_check = afsH_varHsize_Britton_den(PI_single,Rg_H(i1,i2),Rh_H(i1,i2),0,1,1,1,'n-1');
                        if abs( zH(i1,i2) - zH_check ) > 1e-6
                            warning( 'Careful! The values of final size for model H computed in two different ways do not coincide!' )
                            eflag(i1,i2) = eflag(i1,i2) + 1;
                        end
                    end
                    % Now try to compute the final size using the household
                    % size distribution used in the simulation, to check
                    % the difference:
                    zHsim_check(i1,i2) = afsH_varHsize_Addy_den(PI_single_sim,Rg_H(i1,i2),Rh_H(i1,i2),0,1,1,1,'n-1');
                end

                % Calculate the real-time growth rate for model H
                r_H_temp = rH_TVI_tol(av_H_gen,Rg_H(i1,i2),2,Tg,alpha_profile,Tmax,dt,tol);
                r_H_best(i1,i2) = r_H_temp(1);
                r_H_high(i1,i2) = r_H_temp(2);

                %%%%%%%%%%%% Call C code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if Activate_C_codes % The stochastic simulations are in C
                    disp('Run H simulations');
                    simtime = tic; % Check how long they take
                    cmdlH = [ runcommand,' param.txt ',country,' ', num2str(Rg_H(i1,i2),'%.3f'), ' ', num2str(Rh_H(i1,i2),'%.3f'), ' ', num2str(1,'%.1f'), ' ', ...
                         num2str(1,'%.1f'), ' ', num2str(thetaG_random_mix,'%.3f'),' ', num2str(1,'%.2f'),' ', num2str(1,'%.2f'),' ',...
                         num2str(n_init_inf_A_AH,'%d'),' ', num2str(n_init_inf_C_AH,'%d') ];
                    [status result] = system( cmdlH ); % cmdlH is the command line to run the executable with the right arguments
                    % status is 0 if the command is executed correctly
                    % result contains what the C code spits out in the standard output (not used here)

                    output_file_name_H = [ country,'_Rg', num2str(Rg_H(i1,i2),'%.3f'), '_Rw', num2str(Rh_H(i1,i2),'%.3f'), '_sigma', num2str(1,'%.1f'), ...
                        '_rho', num2str(1,'%.1f'), '_ass', num2str(thetaG_random_mix,'%.3f'), '_gammaG', num2str(1,'%.2f'),'_H', num2str(1,'%.2f'),...
                        '__averages.dat' ]; % The simulation creates a data file (.dat, but it's just a text file) with this name
                    [labels,HOW_MANY,data] = readColData( output_file_name_H, 16, 0, 1 ); % Read the Excel file with a function written by someone else.
                    zHsim(i1,i2) = data(3)/100; % Average final size from the simulation (just to cross-check the analytical result). Simulation gives percentage, now turned in a fraction
                    tHsim(i1,i2) = data(7); % Time to the peak
                    piHsim(i1,i2) = data(5); % Peak incidence, expressed in percentages
                    if Activate_checks && ( abs( zH(i1,i2) - zHsim(i1,i2) ) > 0.01 ) % Check analytical and simulated final size are not too different from each other
                        disp('Careful! Analytical and simulated final size does not coincide for the H model');
                        eflag(i1,i2) = eflag(i1,i2) + 1;
                    end
                    tsims = toc(simtime);
                    disp(['Done in ',datestr(tsims/86400,'SS,FFF'),' seconds']);
                end
                if Activate_C_codes && ( eflag(i1,i2) > 0 )
                    beep;
                end

        %%%%%%%%%%%%%%%%%%%% Dysplay output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp(['R0:      ',num2str(R0)]);
%                 disp(['r:      ',num2str(desired_r)]);
                disp(['phi:     ',num2str(phiG)]);
                disp(['psi:     ',num2str(psiG)]);
                disp(['theta:   ',num2str(thetaG)]);
                disp(['theta_A: ',num2str(theta_A(i1,i2))]);
                disp(['pAA:     ',num2str(pAA)]);
                disp(['Rh:      ',num2str(Rh)]);
                disp(['Rh_H:    ',num2str(Rh_H(i1,i2))]);
                disp(['Rg:      ',num2str(Rg(i1,i2))]);
                disp(['Rg_A:    ',num2str(Rg_A(i1,i2))]);
                disp(['Rg_H:    ',num2str(Rg_H(i1,i2))]);
                disp(['R0_AH:   ',num2str(R0_AH(i1,i2))]);
                disp(['R0_A:    ',num2str(R0_A(i1,i2))]);
                disp(['R0_H:    ',num2str(R0_H(i1,i2))]);
                disp(['zAH:     ',num2str(zAH(i1,i2))]);
                disp(['zAHsim:  ',num2str(zAHsim(i1,i2))]);
                disp(['zA:      ',num2str(zA(i1,i2))]);
                disp(['zAsim:   ',num2str(zAsim(i1,i2))]);
                disp(['zH:      ',num2str(zH(i1,i2))]);
                disp(['zHsim:   ',num2str(zHsim(i1,i2))]);
                disp(['tAHsim:  ',num2str(tAHsim(i1,i2))]);
                disp(['tAsim:   ',num2str(tAsim(i1,i2))]);
                disp(['tHsim:   ',num2str(tHsim(i1,i2))]);
                disp(['piAHsim: ',num2str(piAHsim(i1,i2))]);
                disp(['piAsim:  ',num2str(piAsim(i1,i2))]);
                disp(['piHsim:  ',num2str(piHsim(i1,i2))]);
                disp(['hfs_AH:  ',num2str(hfs_AH(i1,i2))]);
                disp(['SAR:     ',num2str(SAR(i1,i2))]);
                disp(['hfs if the first case is an adult:  ',num2str(hfs_by_ic(i1,i2,1))]);
                disp(['hfs if the first case is a child:   ',num2str(hfs_by_ic(i1,i2,2))]);
                disp(['SAR if the first case is an adult:  ',num2str(SAR_by_ic(i1,i2,1))]);
                disp(['SAR if the first case is a child:   ',num2str(SAR_by_ic(i1,i2,2))]);
                disp(['Contribution of household VS global overall transmission:                ',num2str(inVSout(i1,i2,1)),' ',num2str(inVSout(i1,i2,2))]);
                disp(['Proportional contribution of household VS global overall transmission:   ',num2str(inVSout_prop(i1,i2,1)),' ',num2str(inVSout_prop(i1,i2,2))]);
                disp(['Contribution of household VS global A-to-A transmission:                 ',num2str(inVSout_A(i1,i2,1)),' ',num2str(inVSout_A(i1,i2,2))]);
                disp(['Proportional contribution of household VS global A-to-A transmission:    ',num2str(inVSout_A_prop(i1,i2,1)),' ',num2str(inVSout_A_prop(i1,i2,2))]);
                disp(['Activate checks: ',num2str(Activate_checks),' --- Error flag: ',num2str(eflag(i1,i2))])
                disp(' ');
                if ( theta_A(i1,i2) == mapfail ) || ( isnan(mapfail) && isnan(theta_A(i1,i2)) ) || ( NGM_A(1,1) < 0 )
                    if ( NGM_A(1,1) < 0 )
                        disp('Assortativity found but leading to negative AA rate. Map to A failed!');
                    else
                        disp('Map to A failed!')
                    end
                end
                disp(['Run ',num2str(i2),'-',num2str(i1),' of ',num2str(l2),'-',num2str(l1),' completed!']);
                disp(' ');

                tElapsed = toc(tStart);
                tStart = tic;
                TotTime = TotTime + tElapsed;
                save([temp_path,workspace_name_temp]); % Save in the workspace all the values that have been computed till now, in case something goes wrong
                if Activate_C_codes && Activate_deletefile % Clear up the outputs of the individual-based stochastic simulation
                    delete( output_file_name_AH );
                    delete( output_file_name_A );
                    delete( output_file_name_H );
                end
            end
        end
        flag_entry_loop = false;
    end
end

OneRunTime = toc(SingleRun);
disp(['Total time needed (if in a single run): ',datestr(OneRunTime/86400, 'dd:HH:MM:SS,FFF')]);
disp(['Total time needed (estimated by adding separate "continue" runs): ',datestr(TotTime/86400, 'dd:HH:MM:SS,FFF')]);

save(workspace_emergency_name); % Just in case something goes wrong with names and paths, do an emergency save
cd(wrksp_path);
if Activate_workspace_saving
    save(workspace_name); % Save the proper workspace
end
cd(code_path);
if Activate_checks
    rmpath(check_path);
end
if Activate_C_codes
    rmpath(genpath(tool_path))
end

%% Make a cute plot of the final size appear when you're done!
[X,Y] = meshgrid(pAA_vec,psiG_vec);

afs_limits = [ min([min(min(zAH)),min(min(zA)),min(min(zAt)),min(min(zH))]), max([max(max(zAH)),max(max(zA)),max(max(zAt)),max(max(zH))]) ];

figure(10)
clf;
surf(X,Y,zAH');
% zlim([0,1]);
zlim(afs_limits);
title('afs AH model','Fontsize',16);
xlabel('pAA','Fontsize',12);
ylabel('\psi','Fontsize',12);
zlabel('zAH','Fontsize',12);


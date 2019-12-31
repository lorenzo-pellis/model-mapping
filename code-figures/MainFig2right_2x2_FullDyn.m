% This is the code to create the right part of Figure 2 of the main text of 
% Pellis, L et al (2019), Nature Communications
% 
% It relies on a hand-made function called "mycomplexfig"
% It also relies on files containing the epidemic dynamics, in the
% subfolder "simulation-dynamics". If such files are not available, this 
% code runs the stochastic epidemics from a slightly different executable 
% file than the one used in the main model mapping code (i.e. with "_dyn" 
% at the end of the name)
% 
% Update: 31-12-2019

close all; % close all figures
clearvars;
Run_simulation_anyway = 0; % If true, stochastic simulations for the 4 models and 4 parameter combinations are re-run 
Activate_save_fig = 1; % If true, figures are saved
Activate_plot_from_new_simulations = 0; 
% If 0, I make plots from pre-computed and saved simulations (folder code-figures\simulation-dynamics\saved-simulation-dynamics\)
% If 1, I make plots from newly computed workspaces (folder code-figures\simulation-dynamics\)

country = 'GB';
R0 = 2;
phiG = 1;
test_pAA = [ 0.15, 0.75 ];
test_psi = [ 2.4, 1.2 ];
first_subletter = 'e';

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-figures\'];
    fig_path = [base_dir,'\output-figures\main\'];
    tool_path = [base_dir,'\tools\'];
    sim_path = [base_dir,'\code-figures\simulation-dynamics\'];
    runcommand = 'MRCModelMapping_Win_dyn.exe';
    wrksp_path = [base_dir,'\saved-workspaces\GB\'];
    if Activate_plot_from_new_simulations
        saved_sim_path = sim_path;
        ncols = 17; % The simulations used to generate the figure in the main paper 
        % were old and with a slightly different input (18 cols instead of 17)
    else
        saved_sim_path = [sim_path,'\saved_simulation_dynamics\'];
        ncols = 18; % The simulations used to generate the figure in the main paper 
        % were old and with a slightly different input (18 cols instead of 17)
    end
else
    code_path = [base_dir,'/code-figures/'];
    fig_path = [base_dir,'/output-figures/main/'];
    tool_path = [base_dir,'/tools/'];
    sim_path = [base_dir,'/code-figures/simulation-dynamics'];
    runcommand = './ModelMapping_Mac_dyn';
    wrksp_path = [base_dir,'/saved-workspaces/GB/'];
    if Activate_plot_from_new_simulations
        saved_sim_path = sim_path;
        ncols = 17; % The simulations used to generate the figure in the main paper 
        % were old and with a slightly different input (18 cols instead of 17)
    else
        saved_sim_path = [sim_path,'/saved_simulation_dynamics/'];
        ncols = 18; % The simulations used to generate the figure in the main paper 
        % were old and with a slightly different input (18 cols instead of 17)
    end
end
addpath(genpath(tool_path));
% Pre-load baseline workspace to get basic variable values
wrks_name = 'R020_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050';
load([wrksp_path,wrks_name]); 
cd(sim_path); % Work in the directory where the codes for figure simulations are

totpop = 100000;
D.X = (10:50)';
D.Y = cell(2,2);
subtitle_list = cell(2,2);
M = zeros(41,8);

for i = 1:2
    for j = 1:2
        ia = find( abs(pAA_vec - test_pAA(j))< 1e-10, 1 );
        ib = find( abs(psiG_vec - test_psi(i))< 1e-10, 1 );
        
        % Model AH
        fname = [ country, '_Rg',num2str(Rg(ia,ib),'%.3f'), '_Rw', num2str(Rh_vec(ia),'%.3f'), '_sigma', num2str(psiG_vec(ib),'%.1f'), ...
            '_rho', num2str(phiG,'%.1f'), '_ass', num2str(thetaG,'%.3f'), '_gammaG', num2str(g_ratio,'%.2f'),'_H', num2str(h_ratio,'%.2f'),...
            '__sync_large.dat' ];
        cmdlAH = [ runcommand,' param.txt ',country,' ', num2str(Rg(ia,ib),'%.3f'), ' ', num2str(Rh_vec(ia),'%.3f'), ' ', num2str(psiG_vec(ib),'%.1f'), ' ', ...
             num2str(phiG,'%.1f'), ' ', num2str(thetaG,'%.3f'),' ', num2str(g_ratio,'%.2f'),' ', num2str(h_ratio,'%.2f'),' ',...
             num2str(n_init_inf_A_AH,'%d'),' ', num2str(n_init_inf_C_AH,'%d') ];
        if Run_simulation_anyway
            if exist(fname,'file')
                disp('File for dynamics of model AH already exists, but run_code option is active, so I''m running simulations again...');
            else
                disp('File for dynamics of model AH does not exist: running simulations to generate it...');
            end
            [status result] = system( cmdlAH ); % cmdlAH is the command line to run the executable with the right arguments
        else                
            if exist(fname,'file')
                disp('File for dynamics of model AH already exists: no need to run simulations...');
            else
                disp('File for dynamics of model AH does not exist: running simulations to generate it...');
                [status result] = system( cmdlAH ); % cmdlAH is the command line to run the executable with the right arguments
            end
        end
        cd(saved_sim_path);
        [labels,times,data] = readColData( fname, ncols, 0, 0 );
        cd(sim_path);
        inc = data(:,4) * 100 / totpop;
        cuminc = data(:,5) * 100 / totpop;
        
        M(:,1) = inc(11:51);
        M(:,5) = cuminc(11:51);
        
        % Model A
        fname = [ country, '_Rg', num2str(Rg_A(ia,ib),'%.3f'), '_Rw0.000_sigma', num2str(psiG_vec(ib),'%.1f'), ...
            '_rho', num2str(phiG,'%.1f'), '_ass', num2str(theta_A(ia,ib),'%.3f'), '_gammaG', num2str(c_ratio,'%.2f'),'_H', num2str(1,'%.2f'),...
            '__sync_large.dat' ];
        cmdlA = [ runcommand,' param.txt ',country,' ', num2str(Rg_A(ia,ib),'%.3f'), ' ', num2str(0,'%.3f'), ' ', num2str(psiG_vec(ib),'%.1f'), ' ', ...
             num2str(phiG,'%.1f'), ' ', num2str(theta_A(ia,ib),'%.3f'),' ', num2str(g_ratio,'%.2f'),' ', num2str(h_ratio,'%.2f'),' ',...
             num2str(n_init_inf_A_A,'%d'),' ', num2str(n_init_inf_C_A,'%d') ];
        if Run_simulation_anyway
            if exist(fname,'file')
                disp('File for dynamics of model A already exists, but run_code option is active, so I''m running simulations again...');
            else
                disp('File for dynamics of model A does not exist: running simulations to generate it...');
            end
            [status result] = system( cmdlA ); % cmdlA is the command line to run the executable with the right arguments
        else                
            if exist(fname,'file')
                disp('File for dynamics of model A already exists: no need to run simulations...');
            else
                disp('File for dynamics of model A does not exist: running simulations to generate it...');
                [status result] = system( cmdlA ); % cmdlA is the command line to run the executable with the right arguments
            end
        end
        cd(saved_sim_path);
        [labels,times,data] = readColData( fname, ncols, 0, 0 );
        cd(sim_path);
        inc = data(:,4) * 100 / totpop;
        cuminc = data(:,5) * 100 / totpop;
        
        M(:,2) = inc(11:51);
        M(:,6) = cuminc(11:51);
        
        % Model H
        fname = [ country, '_Rg', num2str(Rg_H(ia,ib),'%.3f'), '_Rw', num2str(Rh_H(ia,ib),'%.3f'),...
            '_sigma1.0_rho', num2str(phiG,'%.1f'), '_ass', num2str(thetaG,'%.3f'), '_gammaG', num2str(g_ratio,'%.2f'),'_H', num2str(h_ratio,'%.2f'),...
            '__sync_large.dat' ];
        cmdlH = [ runcommand,' param.txt ',country,' ', num2str(Rg_H(ia,ib),'%.3f'), ' ', num2str(Rh_H(ia,ib),'%.3f'), ' ', num2str(1,'%.1f'), ' ', ...
             num2str(1,'%.1f'), ' ', num2str(thetaG_random_mix,'%.3f'),' ', num2str(1,'%.2f'),' ', num2str(1,'%.2f'),' ',...
             num2str(n_init_inf_A_AH,'%d'),' ', num2str(n_init_inf_C_AH,'%d') ];
        if Run_simulation_anyway
            if exist(fname,'file')
                disp('File for dynamics of model H already exists, but run_code option is active, so I''m running simulations again...');
            else
                disp('File for dynamics of model H does not exist: running simulations to generate it...');
            end
            [status result] = system( cmdlH ); % cmdlH is the command line to run the executable with the right arguments
        else                
            if exist(fname,'file')
                disp('File for dynamics of model H already exists: no need to run simulations...');
            else
                disp('File for dynamics of model H does not exist: running simulations to generate it...');
                [status result] = system( cmdlH ); % cmdlH is the command line to run the executable with the right arguments
            end
        end
        cd(saved_sim_path);
        [labels,times,data] = readColData( fname, ncols, 0, 0 );
        cd(sim_path);
        inc = data(:,4) * 100 / totpop;
        cuminc = data(:,5) * 100 / totpop;
        
        M(:,3) = inc(11:51);
        M(:,7) = cuminc(11:51);
        
        % Model U
        fname = [ country, '_Rg', num2str(R0,'%.3f'), '_Rw0.000_sigma1.0',...
            '_rho', num2str(phiG,'%.1f'), '_ass', num2str(thetaG,'%.3f'), '_gammaG', num2str(g_ratio,'%.2f'),'_H', num2str(h_ratio,'%.2f'),...
            '__sync_large.dat' ];
        cmdlU = [ runcommand,' param.txt ',country,' ', num2str(R0,'%.3f'), ' ', num2str(0,'%.3f'), ' ', num2str(1,'%.1f'), ' ', ...
             num2str(1,'%.1f'), ' ', num2str(thetaG_random_mix,'%.3f'),' ', num2str(1,'%.2f'),' ', num2str(1,'%.2f'),' ',...
             num2str(n_init_inf_A_AH,'%d'),' ', num2str(n_init_inf_C_AH,'%d') ];
        if Run_simulation_anyway
            if exist(fname,'file')
                disp('File for dynamics of model U already exists, but run_code option is active, so I''m running simulations again...');
            else
                disp('File for dynamics of model U does not exist: running simulations to generate it...');
            end
            [status result] = system( cmdlU ); % cmdlU is the command line to run the executable with the right arguments
        else                
            if exist(fname,'file')
                disp('File for dynamics of model U already exists: no need to run simulations...');
            else
                disp('File for dynamics of model U does not exist: running simulations to generate it...');
                [status result] = system( cmdlU ); % cmdlU is the command line to run the executable with the right arguments
            end
        end
        cd(saved_sim_path);
        [labels,times,data] = readColData( fname, ncols, 0, 0 );
        cd(sim_path);
        inc = data(:,4) * 100 / totpop;
        cuminc = data(:,5) * 100 / totpop;
        
        M(:,4) = inc(11:51);
        M(:,8) = cuminc(11:51);
        
        % Store output
        D.Y{i,j} = M;
        subtitle_list{i,j} = ['p_{aa} = ',num2str(pAA_vec(ia),'%1.2f'),', \psi = ',num2str(psiG_vec(ib))];
    end
end
        
%% Figure
cd(code_path)
% To constuct the figure, I pass 3 structures:
%   - D = Data
%   - T = Text (title(s), axe labels, legend text, etc.)
%   - L = Layout information for the figure spatial structure

% %%%%%% Data
D.colorlist = [ 1 0 0; 0 1 1; 1 1 0; 0 0 1 ]; % Order of lines: AH, A, H, U, but plots are A, AH, U, H
% U = 0 0 1; A = 0 1 1; E = 0 1 0; H = 1 1 0; AH = 1 0 0
D.linelist = {'-','--'};

% %%%%%% Text
T.labelx = 'Time (days)';
T.labely = 'Daily incidence (% of population)';
T.labely2 = 'Cumulative incidence (thousands)';

T.subtitles = subtitle_list;
T.first_subletter = first_subletter;

% %%%%%% Layout
L.fig_width_cm = 8; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7.5; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 2;
L.ncols = 2;
% Title
L.xlabel_height_cm = 0.7;
L.ylabel_width_cm = 0.7;
L.ylabel_xpos = 5/8;
L.y2axis = true;
L.y2label_width_cm = 0.7;
% Legend
% The legend is created outside, by hand, in PowerPoint
% Subtitles
L.subtitle_height_cm = 0.5;

T.subtitle_font_size = 6;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'dyn', D, L, T );

if Activate_save_fig
    cd(fig_path);
    export_fig Main_Figure2right -pdf -nocrop -transparent;
    cd(code_path);
end
rmpath(genpath(tool_path));

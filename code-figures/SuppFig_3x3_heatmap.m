% This is the code to create many figures in the supplementary text of 
% Pellis, L et al (2019), Nature Communications
% 
% It generates a 3x3 plot with heatmaps. There is no title, but subtitles 
% are generated for each row and column.
% It relies on a hand-made function called "mycomplexfig.m"
% 
% Update: 13/10/2019

close all; % close all figures
clearvars; % clear all variables
Activate_save_fig = 0; % If true, figures are saved
Activate_C_codes = 1;

Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0
% Choose here the variable to plot and the figure letter (top left)
which_output = 'vc'; figletter = 'A';
% which_output = 'vhc'; figletter = 'B';
% which_output = 'SAR'; figletter = 'C';
% which_output = 'Th'; figletter = 'D';
% which_output = 'ass'; figletter = 'A';

country = 'GB'; % Great Britain
% country = 'SL'; % Sierra Leone
% country = 'SA'; % South Africa

R0vals = [ 1.5 2 4 ]; R0name = 'R0basic';
% R0vals = [ 1.1 1.3 1.5 ]; R0name = 'R0down';
% R0vals = [ 1.7 2 2.3 ]; R0name = 'R0middle';
% R0vals = [ 2.7 3.2 4 ]; R0name = 'R0up';
% popfig = '2ran'; thetaGval = NaN; gamval = 1; % Random
% popfig = 'm4r'; thetaGval = 0.4; gamval = 1; 
% popfig = 'm4UK'; thetaGval = 0.4; gamval = 0.75; 
% popfig = 'm5r'; thetaGval = 0.5; gamval = 1; 
% popfig = 'm5UK'; thetaGval = 0.5; gamval = 0.75; 
popfig = 'UK'; thetaGval = 0.58; gamval = 0.75; % UK
% popfig = 'ass'; thetaGval = 0.7; gamval = 0.75; % More extreme than UK
phivals = [ 1 1.5 2 ];

rvals = [ 0.14552 0.25282 0.52588 ]; % These are the values of r corresponding to the values of R0 = 1.5, 2 and 4
% nsimval = 100;


% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-figures\'];
    fig_path = [base_dir,'\output-figures\supp\',country,'\',popfig,'\'];
    tool_path = [base_dir,'\tools\'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'\output-workspaces\',country,'\'];
    else
        wrksp_path = [base_dir,'\saved-workspaces\',country,'\'];
    end
    if use_match_r
        wrksp_path = [wrksp_path,'match-r\'];
    end
else
    code_path = [base_dir,'/code-figures/'];
    fig_path = [base_dir,'/output-figures/supp/',country,'/',popfig,'/'];
    tool_path = [base_dir,'/tools/'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'/output-workspaces/',country,'/'];
    else
        wrksp_path = [base_dir,'/saved-workspaces/',country,'/'];
    end
    if use_match_r
        wrksp_path = [wrksp_path,'match-r/'];
    end
end
if strcmp(country,'SL') % If country is Sierra Leone
    if isnan(thetaGval) % If mixing is random
        if ~use_match_r
            preloadedwrksp = 'R020_pAA00_03_48_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050';
        else
            preloadedwrksp = 'r025282_pAA00_03_48_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
        end
    else
        if ~use_match_r
            preloadedwrksp = 'R020_pAA00_03_48_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050';
        else
            preloadedwrksp = 'r025282_pAA00_03_48_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050_rfixed_nsimr00100';
        end
    end
elseif strcmp(country,'SA') % If country is Aouth Africa
    if isnan(thetaGval) % If mixing is random
        if ~use_match_r
            preloadedwrksp = 'R020_pAA00_04_63_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050';
        else
            preloadedwrksp = 'r025282_pAA00_04_63_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
        end
    else
        if ~use_match_r
            preloadedwrksp = 'R020_pAA00_04_63_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050';
        else
            preloadedwrksp = 'r025282_pAA00_04_63_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050_rfixed_nsimr00100';
        end
    end
else % If country is Great Britain
    if isnan(thetaGval) % If mixing is random
        if ~use_match_r
            preloadedwrksp = 'R020_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050';
        else
            preloadedwrksp = 'r025282_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
        end
    else
        if ~use_match_r
            preloadedwrksp = 'R020_pAA00_05_95_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050';
        else
            preloadedwrksp = 'r025282_pAA00_05_95_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050_rfixed_nsimr00100';
        end
    end
end
load([wrksp_path,preloadedwrksp]);
cd(code_path); % Work in the directory where the codes for figures are
x_vec = pAA_vec;
y_vec = psiG_vec;

if strcmp(which_output,'ass')
    y_vec(1) = []; % Remove \psiG = 1 when working with assortativity
end
[X,Y] = meshgrid(x_vec,y_vec);
lx = length(x_vec);
ly = length(y_vec);
lR0 = length(R0vals);
lphi = length(phivals);
if isnan(thetaGval)
    thetaGval = thetaG; 
end

% Build figure name
fname = ['3x3_',which_output,'_',country,'_',popfig,'_',R0name];
if use_match_r
    fname = [fname,'_match_r'];
end
    
wrks_list = cell(lR0,lphi);
out_list = cell(lR0,lphi);
% outtitle_list = cell(lR0,lphi);
coltitle_list = cell(lR0);
rowtitle_list = cell(lphi);
out_min = NaN(lR0,lphi);
out_max = NaN(lR0,lphi);

for ip = 1:lphi
    for iR0 = 1:lR0
        % Build workspace name
        if use_match_r
            wname = [ 'r', num2str(round(rvals(iR0)*100000),'%06d') ];
        else
            wname = [ 'R0', num2str(R0vals(iR0)*10,'%02d') ];
        end
        wname = [ wname,'_pAA', num2str(ceil(pAA_min*100),'%02d'),'_',num2str(ceil(dpAA*100),'%02d'),'_',num2str(ceil(pAA_max*100),'%02d'),...
            '_psi', num2str(psiG_min*10,'%02d'), '_', num2str(dpsiG*10,'%02d'),'_', num2str(psiG_max*10,'%02d'),...
            '_phi', num2str(phivals(ip)*10,'%02d'), '_theta',num2str(round(thetaGval*10),'%02d'),...
            '_gammaG', num2str(gamval*100,'%03d'), '_H', num2str(gamval*100,'%03d') ];
        if Activate_C_codes
            wname = [ wname, '_100sim_e5_init', num2str(n_init_inf,'%03d') ];
        end
        if use_match_r
            wname = [ wname, '_rfixed_nsimr',num2str(nsimval,'%05d') ];
        end
        wrks_list{ip,iR0} = wname;
        
        try_name = strcat(wrksp_path,wrks_list{ip,iR0});
        load(try_name);

        if strcmp(which_output,'vc')
            out_list{ip,iR0} = 100*v_AH(:,:,2)';
        elseif strcmp(which_output,'vhc')
            out_list{ip,iR0} = 100*vH_AH(:,:,2)';
        elseif strcmp(which_output,'SAR')
            out_list{ip,iR0} = 100*SAR';
        elseif strcmp(which_output,'Th')
            out_list{ip,iR0} = 100*inVSout_prop(:,:,1)';
        elseif strcmp(which_output,'ass')
            out_list{ip,iR0} = theta_A(:,2:end)';
        else
            disp('Incorrect output');
            beep;
        end
%         outtitle_list{iR0,ip} = ['R_0 = ',num2str(R0,'%1.1f'),', \phi = ',num2str(phiG,'%1.1f')];
%         outtitle_list{ip,iR0} = ['R_0 = ',num2str(R0),', \phi = ',num2str(phiG)];
        coltitle_list{iR0} = ['R_0 = ',num2str(R0)];
        rowtitle_list{ip} = ['\phi = ',num2str(phiG)];
        out_min(ip,iR0) = min(min(out_list{ip,iR0}));
        out_max(ip,iR0) = max(max(out_list{ip,iR0}));
    end
end

%% Figure

outmin_min = min(min(out_min)); outmax_max = max(max(out_max));
if strcmp(which_output,'vc')
    title = 'Percentage of children in incidence';
elseif strcmp(which_output,'vhc')
    title = 'Percentage of children among first household cases';
elseif strcmp(which_output,'SAR')
    title = 'Household secondary attack rate (%)';
elseif strcmp(which_output,'Th')
    title = 'Percentage of transmission in household';
elseif strcmp(which_output,'ass')
    title = ['Mapped assortativity \theta^A of A model (input: \theta^{AH} = ',num2str(thetaG,'%.2f'),')'];
else
    title = 'Not sure what you want me to plot!';
    beep;
end


% %%%%%% Data
D.X = X;
D.Y = Y;
D.Z = out_list;
if strcmp(which_output,'ass')
    D.Zover = out_list;
    if strcmp(popfig,'2ran')
        D.Zcont = thetaG;
    else
        D.Zcont = [];
    end
    D.Clim = [ 0, 2/3 ]; % Manual input, to plot ass for 2ran and UK side by side. Same extreme used for SA and SL for consistency
else
    D.Clim = [ outmin_min, outmax_max ];
end
D.Yticks = 0:0.5:4;

if strcmp(country,'SL')
    D.Xticks = [ 0 0.1 0.2 0.3 0.4 ];
elseif strcmp(country,'SA')
    D.Xticks = [ 0 0.15 0.3 0.45 0.6 ];
else
    D.Xticks = [ 0 0.2 0.4 0.6 0.8 ];
end

% %%%%%% Text
% Annoyingly, spaces in text are weird and system-dependent (differently in
% Windows and Mac), so the first lines show how it should work, while
% the following lines offer a solution for how the result should look, 
% (which works for me in Windows)
% T.labelx = 'Adult-to-adult within-household transmission probability (p_{aa})';
% T.labely = 'Relative susceptibility of children (\psi)';
T.labelx = 'Adult-to-adult within-household transmission probability{(p_{aa})}';
T.labely = 'Relative susceptibility of children versus adults { (\psi)}';

% T.subtitles = outtitle_list;
T.coltitles = coltitle_list;
T.rowtitles = rowtitle_list;
T.figletter = figletter;

% %%%%%% Layout
L.fig_width_cm = 7; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 3;
L.ncols = 3;
% Title
L.title_height_cm = 0;
L.xlabel_height_cm = 0.8;
L.ylabel_width_cm = 0.8;
T.labely_ycorr = -0.01;
% Legend
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0.7;
% Subtitles
% L.subtitle_height_cm = 0;
L.coltitle_height_cm = 0.5;
L.rowtitle_width_cm = 0.4;
L.colorbar = true;
% T.figletter_ycorr = 0.005;

T.axes_font_size = 6;
T.subtitle_font_size = 6;
T.subaxes_font_size = 5;
T.legend_font_size = 5;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'heatmap', D, L, T );

if Activate_save_fig
    addpath(genpath(tool_path));
    cd(fig_path);
    expcmd = ['export_fig ',fname,' -pdf -nocrop -transparent'];
    eval(expcmd);
    cd(code_path);
    rmpath(genpath(tool_path));
end

% This is the code to create Figure 1 of the main text of 
% Pellis, L. et al (2020), Nature Communications
% 
% It relies on a hand-made function called "mysubplot_outputs"
% 
% Update: 22-01-2020

close all; % close all figures
clearvars; % clear all variables
Activate_save_fig = 1; % If true, figures are saved
Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
fig_base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
% Names of folders are preceded by "fig_" because there is a risk that
% loading a workspace might override the names used in this scrip
if ispc
    fig_code_path = [fig_base_dir,'\code-figures\'];
    fig_fig_path = [fig_base_dir,'\output-figures\main\'];
    fig_tool_path = [fig_base_dir,'\tools\'];
    if Activate_plot_from_new_workspaces
        fig_wrksp_path = [fig_base_dir,'\output-workspaces\GB\'];
    else
        fig_wrksp_path = [fig_base_dir,'\saved-workspaces\GB\'];
    end
else
    fig_code_path = [fig_base_dir,'/code-figures/'];
    fig_fig_path = [fig_base_dir,'/output-figures/main/'];
    fig_tool_path = [fig_base_dir,'/tools/'];
    if Activate_plot_from_new_workspaces
        fig_wrksp_path = [fig_base_dir,'/output-workspaces/GB/'];
    else
        fig_wrksp_path = [fig_base_dir,'/saved-workspaces/GB/'];
    end
end
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
% Loaded workspaces generated from other computers, with a different folder
% structure might have saved anonymous functions with paths not recognised 
% when this script is run, so the line above turns off that warning.
wrks_name = 'R020_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050';
% wrks_name = 'R020_pAA00_05_95_log2psi-200_025_200_phi15_theta02_gammaG100_H100_100sim_e5_init050';
% wrks_name = 'R040_pAA00_05_95_psi10_02_40_phi20_theta07_gammaG075_H075_100sim_e5_init050';
load(strcat(fig_wrksp_path,wrks_name));
cd(fig_code_path); % Work in the directory where the codes for figures are

% Annoyingly, spaces in text are weird and system-dependent (differently in
% Windows and Mac), so the first lines show how it should work, while
% the following lines offer a solution for how the result should look, 
% (which works for me in Windows)
% % Correct spacing (R2019a on Windows):
labelx = 'Relative susceptibility of children versus adults (\psi)'; y_vec = psiG_vec;
labely = '    Adult-to-adult within-household transmission probability (p_{aa})'; x_vec = pAA_vec;
% Fudged spacing (R2016a on Mac):
% labelx = 'Relative susceptibility of children versus adults { (\psi)}'; y_vec = psiG_vec;
% labely = '    Adult-to-adult within-household transmission probability{(p_{aa})}'; x_vec = pAA_vec;

[X,Y] = meshgrid(x_vec,y_vec);

zAH_perc = 100*zAH;
zA_perc = 100*zA;
zH_perc = 100*zH;
zlimits = [ min([min(min(zAH_perc)),min(min(zA_perc)),min(min(zH_perc))]), max([max(max(zAH_perc)),max(max(zA_perc)),max(max(zH))]) ];
pilimits = [ min([min(min(piAHsim)),min(min(piAsim)),min(min(piHsim))]), max([max(max(piAHsim)),max(max(piAsim)),max(max(piHsim))]) ];
tlimits = [ min([min(min(tAHsim)),min(min(tAsim)),min(min(tHsim))]), max([max(max(tAHsim)),max(max(tAsim)),max(max(tHsim))]) ];

SARlimits = [ min(min(SAR)), max(max(SAR)) ];

rho = v_AH(:,:,1) ./ v_AH(:,:,2);
rholimits = [ min(min(rho)), max(max(rho)) ];

inVSout_prop_limits = [ min(min(inVSout_prop(:,:,1))), max(max(inVSout_prop(:,:,1))) ];

iVSo = inVSout_prop(:,:,1) ./ inVSout_prop(:,:,2);
iVSolimits = [ min(min(iVSo)), max(max(iVSo)) ];

figure_list(:,:,1) = zAH_perc';
model_list{1,:} = 'model  AH';
climits_list(1,:) = zlimits;
clabel_list{1,:} = ['    Final size   '; '(% of population)'];

figure_list(:,:,2) = piAHsim';
model_list{2,:} = 'model  AH';
climits_list(2,:) = pilimits;
clabel_list{2,:} = ['Daily peak incidence';'  (% of population) '];

figure_list(:,:,3) = tAHsim';
model_list{3,:} = 'model  AH';
climits_list(3,:) = tlimits;
clabel_list{3,:} = ['      Time to peak     ';'(mean generation times)'];

figure_list(:,:,4) = zA_perc';
model_list{4,:} = 'model  A';
climits_list(4,:) = zlimits;
clabel_list{4,:} = 'Final size (% of population)';

figure_list(:,:,5) = piAsim';
model_list{5,:} = 'model  A';
climits_list(5,:) = pilimits;
clabel_list{5,:} = 'Daily peak incidence (% of population)';

figure_list(:,:,6) = tAsim';
model_list{6,:} = 'model  A';
climits_list(6,:) = tlimits;
clabel_list{6,:} = 'Time to peak (units: generation times)';

figure_list(:,:,7) = zH_perc';
model_list{7,:} = 'model  H';
climits_list(7,:) = zlimits;
clabel_list{7,:} = 'Final size (% of population)';

figure_list(:,:,8) = piHsim';
model_list{8,:} = 'model  H';
climits_list(8,:) = pilimits;
clabel_list{8,:} = 'Daily peak incidence (% of population)';

figure_list(:,:,9) = tHsim';
model_list{9,:} = 'model  H';
climits_list(9,:) = tlimits;
clabel_list{9,:} = 'Time to peak (units: generation times)';

figure_handles = mysubplot_outputs( 22, 25, 'Comparison between models'' outputs',...
    labely, labelx, model_list, X, Y, figure_list, climits_list, clabel_list );

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    export_fig Main_Figure1 -pdf -nocrop -transparent;
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end

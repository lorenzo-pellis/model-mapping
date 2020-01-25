% This is the code to create Figure 3 of the supplementary text of 
% Pellis, L et al (2020), Nature Communications
% 
% Update: 11-10-2019

close all; % close all figures
clearvars;
Activate_save_fig = 1; % If true, figures are saved
Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
fig_base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
% Names of folders are preceded by "fig_" because there is a risk that
% loading a workspace might override the names used in this scrip
if ispc
    fig_code_path = [fig_base_dir,'\code-figures\'];
    fig_fig_path = [fig_base_dir,'\output-figures\supp\'];
    fig_tool_path = [fig_base_dir,'\tools\'];
    if Activate_plot_from_new_workspaces
        fig_wrksp_path = [fig_base_dir,'\output-workspaces\assortativity\'];
    else
        fig_wrksp_path = [fig_base_dir,'\saved-workspaces\assortativity\'];
    end
else
    fig_code_path = [fig_base_dir,'/code-figures/'];
    fig_fig_path = [fig_base_dir,'/output-figures/supp/'];
    fig_tool_path = [fig_base_dir,'/tools/'];
    if Activate_plot_from_new_workspaces
        fig_wrksp_path = [fig_base_dir,'/output-workspaces/assortativity/'];
    else
        fig_wrksp_path = [fig_base_dir,'/saved-workspaces/assortativity/'];
    end
end
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
% Loaded workspaces generated from other computers, with a different folder
% structure might have saved anonymous functions with paths not recognised 
% when this script is run, so the line above turns off that warning.
wrksp_name = 'GB_R020_pAA50_50_50_psiGcustom_phi10_theta02_gammaG100_H100_plot_v';
% wrksp_name = 'GB_R020_pAA20_30_80_psi06_04_18_phi10_theta02_gammaG100_H100_plot_v';
% wrksp_name = 'SL_R020_pAA20_30_80_psi04_02_18_phi10_theta04_gammaG100_H100_plot_v';
cd(fig_code_path)

fig_name = 'plot_v';
load([fig_wrksp_path,wrksp_name]);
tot_plots = l2;
switch tot_plots
    case 1
        nrows = 1;
        ncols = 1;
        width = 20;
        height = 20;
    case 2
        nrows = 1;
        ncols = 2;
        width = 20;
        height = 20;
    case 3
        nrows = 2;
        ncols = 2;
        width = 20;
        height = 20;
    case 4
        nrows = 2;
        ncols = 2;
        width = 20;
        height = 20;
    case 5
        nrows = 2;
        ncols = 3;
        width = 20;
        height = 20;
    case 6
        nrows = 2;
        ncols = 3;
        width = 20;
        height = 20;
    case 7
        nrows = 3;
        ncols = 3;
        width = 20;
        height = 20;
    case 8
        nrows = 4;
        ncols = 2;
        width = 28;
        height = 15;
    case 9
        nrows = 3;
        ncols = 3;
        width = 20;
        height = 20;
    case 10
        nrows = 3;
        ncols = 4;
        width = 20;
        height = 20;
    case 11
        nrows = 3;
        ncols = 4;
        width = 20;
        height = 20;
    case 12
        nrows = 3;
        ncols = 4;
        width = 26;
        height = 20;
    otherwise
        warning('Too many subplots: not plotting anything...');
end

for i1 = 1:l1
    for i2 = 1:l2
        p = i2;%( i1 - 1 ) * l1 + i2;
        vAH_list(:,p) = all_vAH(:,i1,i2);
        v_list(:,:,p) = all_v(:,:,i1,i2);
        vtitle_list(p,:) = ['\psi = ',num2str(psiG_vec(i2),'%1.2f')];
    end
    ass_handles = mysubplot_find_v( width, height, '', ...
        'Assortativity of children in model A   (\theta^A)',...
        'Fractions adults (blue) and children (red)',...
        vtitle_list, ass, v_list, vAH_list );
end

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    export_fig plot_v -pdf -nocrop -transparent;
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end

% This is the code to create Figure 3 of the supplementary text of 
% Pellis, L et al (2019), Nature Communications
% 
% Update: 11/10/2019

close all; % close all figures
clearvars;
Activate_save_fig = 1; % If true, figures are saved
Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-figures\'];
    fig_path = [base_dir,'\output-figures\supp\'];
    tool_path = [base_dir,'\tools\'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'\output-workspaces\assortativity\'];
    else
        wrksp_path = [base_dir,'\saved-workspaces\assortativity\'];
    end
else
    code_path = [base_dir,'/code-figures/'];
    fig_path = [base_dir,'/output-figures/supp/'];
    tool_path = [base_dir,'/tools/'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'/output-workspaces/assortativity/'];
    else
        wrksp_path = [base_dir,'/saved-workspaces/assortativity/'];
    end
end
wrksp_name = 'GB_R020_pAA50_50_50_psiGcustom_phi10_theta02_gammaG100_H100_plot_v';
% wrksp_name = 'GB_R020_pAA20_30_80_psi06_04_18_phi10_theta02_gammaG100_H100_plot_v';
% wrksp_name = 'SL_R020_pAA20_30_80_psi04_02_18_phi10_theta04_gammaG100_H100_plot_v';
cd(code_path)

fig_name = 'plot_v';
load([wrksp_path,wrksp_name]);
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
    addpath(genpath(tool_path));
    cd(fig_path);
    export_fig plot_v -pdf -nocrop -transparent;
    cd(code_path);
    rmpath(genpath(tool_path));
end

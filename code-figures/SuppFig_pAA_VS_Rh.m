% This is the code to create Figure 1 of the supplementary text of 
% Pellis, L. et al (2020), Nature Communications
% 
% Update: 10-10-2019

close all;
clearvars;
Activate_save_fig = 1; % If true, figures are saved

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
fig_base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
% Names of folders are preceded by "fig_" because there is a risk that
% loading a workspace might override the names used in this scrip
if ispc
    fig_code_path = [fig_base_dir,'\code-figures\'];
    fig_temp_path = [fig_base_dir,'\code-model-mapping\'];
    fig_fig_path = [fig_base_dir,'\output-figures\supp\'];
    fig_tool_path = [fig_base_dir,'\tools\'];
else
    fig_code_path = [fig_base_dir,'/code-figures/'];
    fig_temp_path = [fig_base_dir,'/code-model-mapping/'];
    fig_fig_path = [fig_base_dir,'/output-figures/supp/'];
    fig_tool_path = [fig_base_dir,'/tools/'];
end

cd(fig_temp_path); % Population structure and some of the relevant codes are in this folder
input_distr = 'GB_H_structure_ModelMapping.txt';
H = load(input_distr);
H_single = create_1type_distr(H);

Rh_min = 0;
Rh_max = 10;
dRh = 0.01;
Rh = Rh_min:dRh:Rh_max;
lRh = length(Rh);
pAA = zeros(1,lRh);

for iRh = 1:lRh
    pAA(iRh) = get_pAA_from_Rh(H_single,Rh(iRh));
end

figure(1)
clf;
hold on;
plot(Rh,pAA,'b','Linewidth',2)
for n = 2:8
    plot(Rh,1-exp(-Rh/(n-1)),'b-.','Linewidth',1)
end
% title('Adult-to-adult within-household transmission probability','Fontsize',16);
set(gca,'fontsize',14);
xlabel('\beta_h','Fontsize',18);
ylabel('p_{aa}','Fontsize',18);

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    export_fig pAA_VS_Rh -pdf -nocrop -transparent;
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end


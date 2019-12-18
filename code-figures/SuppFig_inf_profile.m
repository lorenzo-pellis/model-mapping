% This is the code to create Figure 1 of the supplementary text of 
% Pellis, L et al (2019), Nature Communications
% 
% Update: 10/10/2019

close all;
clearvars;
Activate_save_fig = 1; % If true, figures are saved

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-figures\'];
    fig_path = [base_dir,'\output-figures\supp\'];
    tool_path = [base_dir,'\tools\'];
else
    code_path = [base_dir,'/code-figures/'];
    fig_path = [base_dir,'/output-figures/supp/'];
    tool_path = [base_dir,'/tools/'];
end
cd(code_path)

supp_min = 0;
supp_max = 6;
xmax = 12;
dx = 0.01;
x = supp_min:dx:supp_max;
x0 = supp_max:dx:xmax;
lx = length(x);
c = 0.65;
alpha = 5;
lambda = 3;

y = gampdf(x,alpha,1/lambda); % First argument is the shape, then the scale (careful how Matlab defines the scale parameter)
y(lx) = 0; % Artificially put the last point to 0

figure;
clf;
h = axes('Xtick',[],'Xticklabel',[],'Ytick',[],'Yticklabel',[]);
xlim([0 13]);
hold on;
ax = annotation('arrow',[ 0.13 0.13+0.775 ], [ 0.11 0.11 ]);
ax.LineWidth = 1;
ay = annotation('arrow',[ 0.13 0.13 ], [ 0.11 0.11+0.815 ]);
ay.LineWidth = 1;
plot(x,y,'r','Linewidth',1)
plot(x0,zeros(length(x0)),'r','Linewidth',1)
plot([ 0 supp_max ], [ c c ],'k','Linewidth',1);
plot([ supp_max supp_max ], [ c 0 ],'k','Linewidth',1);
plot([ supp_max xmax ], [ c c ],'k--','Linewidth',1);
plot([ xmax xmax ], [ c 0 ],'k--','Linewidth',1);
text(-0.8,c,'$c$','Fontsize',24','interpreter','latex'); 
text(supp_max-0.3,-0.05,'$T$','Fontsize',24,'interpreter','latex'); 
text(xmax-0.3,-0.05,'$T''$','Fontsize',24,'interpreter','latex'); 
hold off;

if Activate_save_fig
    addpath(genpath(tool_path));
    cd(fig_path);
    export_fig infprof -pdf -nocrop -transparent;
    cd(code_path);
    rmpath(genpath(tool_path));
end

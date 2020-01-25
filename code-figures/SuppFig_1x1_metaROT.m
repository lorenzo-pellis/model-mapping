% This is the code to create Figure 25 in the supplementary text of 
% Pellis, L et al (2020), Nature Communications
% 
% It generates a 1x1 graph with a summary of the rule of thumb for each 
% studied population and each studies threshold.
% 
% It relies on a hand-made function called "mycomplexfig"
% It also relies on pre-computed workspaces where information about the 
% rule of thumb is stored. Such workspaces are contained in the subfolder
% called "rule-of-thumb". If they are not available, one needs to first run
% "Analyse_data_for_RuleOfThumb.m"
% 
% Update: 14-10-2019

close all; % close all figures
clearvars;
Activate_save_fig = 1; % If true, figures are saved

usez = true; % If true, I include the final size (z) in the overall acceptance region plot, otherwise not
usepi = true; % If true, I include the peak incidence (pi) in the overall acceptance region plot, otherwise not
uset = true;  % If true, I include the time to the peak (t) in the overall acceptance region plot, otherwise not
use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0

clist = {'GB','SA','SL'};
elist = [ 0.01, 0.05, 0.1 ];

current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
fig_base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
% Names of folders are preceded by "fig_" because there is a risk that
% loading a workspace might override the names used in this scrip
if ispc
    fig_code_path = [fig_base_dir,'\code-figures\'];
    fig_fig_path = [fig_base_dir,'\output-figures\supp\'];
    fig_tool_path = [fig_base_dir,'\tools\'];
else
    fig_code_path = [fig_base_dir,'/code-figures/'];
    fig_fig_path = [fig_base_dir,'/output-figures/supp/'];
    fig_tool_path = [fig_base_dir,'/tools/'];
end
cd(fig_code_path); % Work in the directory where the code for figures are

% Build figure name
fname = ['1x1_metaROT_'];
if usez
    fname = [fname,'z'];
end
if usepi
    fname = [fname,'pi'];
end
if uset
    fname = [fname,'t'];
end
if use_match_r
    fname = [fname,'_match_r'];
end


for e = 1:3
    tolval = elist(e);
    for c = 1:3
        country = clist{c};
        
        if tolval == 0.01
            if strcmp(country,'GB')
                myROTcoeff = 0.092;
                myROTinter = -0.048;
            elseif strcmp(country,'SA')
                myROTcoeff = 0.074;%0.072;
                myROTinter = -0.051;%-0.048;
            elseif strcmp(country,'SL')
                myROTcoeff = 0.068;%0.07;
                myROTinter = -0.052;%-0.062;
            else
                myROTcoeff = 0;
                myROTinter = 0;
            end
        elseif tolval == 0.05
            if strcmp(country,'GB')
                myROTcoeff = 0.24;%0.236;
                myROTinter = -0.1;%-0.094;
            elseif strcmp(country,'SA')
                myROTcoeff = 0.225;%0.228;
                myROTinter = -0.16;%-0.166;
            elseif strcmp(country,'SL')
                myROTcoeff = 0.22;
                myROTinter = -0.18;%-0.188;
            else
                myROTcoeff = 0;
                myROTinter = 0;
            end
        elseif tolval == 0.1
            if strcmp(country,'GB')
                myROTcoeff = 0.368;
                myROTinter = -0.14;
            elseif strcmp(country,'SA')
                myROTcoeff = 0.347;%0.348;
                myROTinter = -0.23;%-0.22;
            elseif strcmp(country,'SL')
                myROTcoeff = 0.34;
                myROTinter = -0.26;
            else
                myROTcoeff = 0;
                myROTinter = 0;
            end
        else
            myROTcoeff = 0;
            myROTinter = 0;
        end
        myROT = [ myROTinter, myROTcoeff];

        datax(e,c) = 100*myROTinter;
        datay(e,c) = 100*myROTcoeff';
    end
end

%% Figure

% %%%%%% Data
D.X = datax;
D.Y = datay;
D.colorlist = {'g','r','b'};
D.limx = [-28,-3];
% D.limy = [0,100];

% %%%%%% Text
T.figletter = 'D';
T.title = 'Rule of thumb extensions';
T.labelx = '\beta_0';
T.labely = '\beta_1';

% %%%%%% Layout
L.fig_width_cm = 7; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 portrait ~ 25cm
L.screen_scale = 2.5;
% Panels
L.nrows = 1;
L.ncols = 1;
% Title
L.title_height_cm = 0.6;
L.title_ypos = 2/3;
L.xlabel_height_cm = 1;
L.ylabel_width_cm = 0.8;
% Legend
L.legend = true;
T.legentries = {'Great Britain','South Africa','Sierra Leone',' 1% accuracy',' 5% accuracy','10% accuracy'};
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0;
L.colorbar = false;
% Subtitles
L.subtitle_height_cm = 0;
T.title_font_size = 8;
T.axes_font_size = 8;
T.subaxes_font_size = 6;


L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'metaROT', D, L, T );

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    expcmd = ['export_fig ',fname,' -pdf -nocrop -transparent'];
    eval(expcmd);
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end
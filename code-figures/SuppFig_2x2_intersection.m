% This is the code to create many figures in the supplementary text of 
% Pellis, L et al (2020), Nature Communications
% 
% It generates a 2x2 graph with the simplest-model-acceptance-region plot 
% for each of the 3 outputs and then the overall simplest-model-acceptance-
% region plot, where the simplest model is selected for which all 3 outputs
% are within the selected tolerance. There is no title, but subtitles 
% are generated for each row and column.
% It relies on a hand-made function called "mycomplexfig.m"
% 
% Update: 16-01-2020

close all; % close all figures
clearvars; % clear all variables
Activate_save_fig = 1; % If true, figures are saved
Activate_C_codes = 1;

Activate_plot_from_new_workspaces = 1; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0

country = 'GB'; % Great Britain
% country = 'SL'; % Sierra Leone
% country = 'SA'; % South Africa

popfig = '2ran'; thetaGval = NaN; gamval = 1; % Random
% popfig = 'm4r'; thetaGval = 0.4; gamval = 1; 
% popfig = 'm4UK'; thetaGval = 0.4; gamval = 0.75; 
% popfig = 'm5r'; thetaGval = 0.5; gamval = 1; 
% popfig = 'm5UK'; thetaGval = 0.5; gamval = 0.75; 
% popfig = 'UK'; thetaGval = 0.58; gamval = 0.75; % UK
% popfig = 'ass'; thetaGval = 0.7; gamval = 0.75; % More extreme than UK

%%% 5% threshold - variable phi
phi_in_title = 1;
tolval = 0.05;
% R0val = 1.5; phiGfig = 1; first_subletter = 'a';
R0val = 2; phiGfig = 1; first_subletter = 'e';
% R0val = 4; phiGfig = 1; first_subletter = 'i';
% R0val = 1.5; phiGfig = 2; first_subletter = 'm';
% R0val = 2; phiGfig = 2; first_subletter = 'q';
% R0val = 4; phiGfig = 2; first_subletter = 'u';

%%% phi = 1 - variable thresholds
% phi_in_title = 0;
% phiGfig = 1;
% R0val = 1.5; tolval = 0.01; first_subletter = 'a';
% % R0val = 2; tolval = 0.01; first_subletter = 'e';
% % R0val = 4; tolval = 0.01; first_subletter = 'i';
% % R0val = 1.5; tolval = 0.1; first_subletter = 'm';
% % R0val = 2; tolval = 0.1; first_subletter = 'q';
% % R0val = 4; tolval = 0.1; first_subletter = 'u';

rval = 0.25282; % [ 0.14552 0.25282 0.52588 ] are the values corresponding to R0 = 1.5, 2 and 4;
% nsimval = 100;
if phi_in_title
    title = [ 'Simplest model acceptance region: R_0 = ',num2str(R0val),', \phi = ',num2str(phiGfig) ];
else
    title = [ 'Simplest model acceptance region: R_0 = ',num2str(R0val),', \epsilon = ',num2str(100*tolval),'%' ];
end

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
fig_base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
% Names of folders are preceded by "fig_" because there is a risk that
% loading a workspace might override the names used in this scrip
if ispc
    fig_code_path = [fig_base_dir,'\code-figures\'];
    fig_fig_path = [fig_base_dir,'\output-figures\supp\',country,'\',popfig,'\'];
    fig_tool_path = [fig_base_dir,'\tools\'];
    if Activate_plot_from_new_workspaces
        fig_wrksp_path = [fig_base_dir,'\output-workspaces\',country,'\'];
    else
        fig_wrksp_path = [fig_base_dir,'\saved-workspaces\',country,'\'];
    end
    if use_match_r
        fig_wrksp_path = [fig_wrksp_path,'match-r\'];
    end
else
    fig_code_path = [fig_base_dir,'/code-figures/'];
    fig_fig_path = [fig_base_dir,'/output-figures/supp/',country,'/',popfig,'/'];
    fig_tool_path = [fig_base_dir,'/tools/'];
    if Activate_plot_from_new_workspaces
        fig_wrksp_path = [fig_base_dir,'/output-workspaces/',country,'/'];
    else
        fig_wrksp_path = [fig_base_dir,'/saved-workspaces/',country,'/'];
    end
    if use_match_r
        fig_wrksp_path = [fig_wrksp_path,'match-r/'];
    end
end
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
% Loaded workspaces generated from other computers, with a different folder
% structure might have saved anonymous functions with paths not recognised 
% when this script is run, so the line above turns off that warning.
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
load([fig_wrksp_path,preloadedwrksp]);
cd(fig_code_path); % Work in the directory where the codes for figures are
x_vec = pAA_vec;
y_vec = psiG_vec;
[X,Y] = meshgrid(x_vec,y_vec);
lx = length(x_vec);
ly = length(y_vec);
if isnan(thetaGval)
    thetaGval = thetaG; 
end


% Workspace name
if use_match_r
    wname = [ 'r', num2str(round(r*100000),'%06d') ];
else
    wname = [ 'R0', num2str(R0val*10,'%02d') ];
end
wname = [ wname,'_pAA', num2str(pAA_min*100,'%02d'),'_',num2str(round(dpAA*100),'%02d'),'_',num2str(round(pAA_max*100),'%02d'),...
    '_psi', num2str(psiG_min*10,'%02d'), '_', num2str(dpsiG*10,'%02d'),'_', num2str(psiG_max*10,'%02d'),...
    '_phi', num2str(phiGfig*10,'%02d'), '_theta',num2str(round(thetaGval*10),'%02d'),...
    '_gammaG', num2str(gamval*100,'%03d'), '_H', num2str(gamval*100,'%03d') ];
if Activate_C_codes
    wname = [ wname, '_100sim_e5_init', num2str(n_init_inf,'%03d') ];
end
if use_match_r
    wname = [ wname, '_rfixed_nsimr',num2str(nsimval,'%05d') ];
end
load([fig_wrksp_path,wname]);

% Build figure name
if phi_in_title
    fname = ['2x2_inters_',country,'_',popfig,'_R0', num2str(R0val*10,'%02d'),'_phi', num2str(phiGfig*10,'%02d'),'_fixed_eps', num2str(tolval*100,'%02d')];
else
    fname = ['2x2_inters_',country,'_',popfig,'_R0', num2str(R0val*10,'%02d'),'_eps', num2str(tolval*100,'%02d'),'_fixed_phi', num2str(phiGfig*10,'%02d')];
end
if use_match_r
    fname = [fname,'_match_r'];
end

ztol = tolval;
sztol = tolval;
pitol = tolval;
ttol = tolval;
zabstol = ztol;
zreltol = ztol;
szabstol = sztol;
szreltol = sztol;
piabstol = pitol;
pireltol = pitol;
tabstol = ttol;
treltol = ttol;

szAH_A = afs_each_typeAH(:,:,1);
szAH_C = afs_each_typeAH(:,:,2);
szA_A = afs_each_typeA(:,:,1);
szA_C = afs_each_typeA(:,:,2);

% Final size differences
dz_AH_H = zAH - zH; % dz = difference in final size
dz_AH_A = zAH - zA;
dz_AH_At = zAH - zAt;
dz_A_At = zA - zAt;
zplot_end = max([max(max(abs(dz_AH_H))),max(max(abs(dz_AH_A))),max(max(abs(dz_AH_At))),max(max(abs(dz_A_At)))]);

rel_dz_AH_H = ( zAH - zH ) ./ ( ( zAH + zH ) / 2 );
rel_dz_AH_A = ( zAH - zA ) ./ ( ( zAH + zA ) / 2 );
rel_dz_AH_U = ( zAH - zH(1,1) ) ./ ( ( zAH + zH(1,1) ) / 2 );
rel_dz_AH_At = ( zAH - zAt ) ./ ( ( zAH + zAt ) / 2 );
rel_dz_A_At = ( zA - zAt ) ./ ( ( zA + zAt ) / 2 );
zplot_rel_end = max([max(max(abs(rel_dz_AH_H))),max(max(abs(rel_dz_AH_A))),max(max(abs(rel_dz_AH_At))),max(max(abs(rel_dz_A_At)))]);
zplot_rel_end_perc = ceil( zplot_rel_end / zreltol ) * zreltol;

zrelregA = abs(rel_dz_AH_A) < zreltol;
zrelregH = abs(rel_dz_AH_H) < zreltol;
zrelregU = abs(rel_dz_AH_U) < zreltol;
zrelreg = -1 * ones(lx,ly);
zrelreg_check = -1 * ones(lx,ly);
for ix = 1:lx
    for iy = 1:ly
        if zrelregU(ix,iy)
            zrelreg(ix,iy) = 0;
        elseif zrelregA(ix,iy)
            zrelreg(ix,iy) = 1;
        elseif zrelregH(ix,iy)
            zrelreg(ix,iy) = 3;
        else
            zrelreg(ix,iy) = 4;
        end

        if ( xor(zrelregA(ix,iy),zrelregH(ix,iy)) ) && (~zrelregU(ix,iy))
            zrelreg_check(ix,iy) = 0;
        elseif ( zrelregA(ix,iy) && zrelregH(ix,iy) && (~zrelregU(ix,iy)) )
            zrelreg_check(ix,iy) = 1;
            zrelreg(ix,iy) = 2;
        end
    end
end

if Activate_C_codes
    % Peak incidence differences
    dpi_AH_H = piAHsim - piHsim; % dpi = difference in peak incidence
    dpi_AH_A = piAHsim - piAsim;
    piplot_end = max([max(max(abs(dpi_AH_H))),max(max(abs(dpi_AH_A)))]);

    rel_dpi_AH_H = ( piAHsim - piHsim ) ./ ( ( piAHsim + piHsim ) / 2 );
    rel_dpi_AH_A = ( piAHsim - piAsim ) ./ ( ( piAHsim + piAsim ) / 2 );
    rel_dpi_AH_U = ( piAHsim - piHsim(1,1) ) ./ ( ( piAHsim + piHsim(1,1) ) / 2 );
    piplot_rel_end = max([max(max(abs(rel_dpi_AH_H))),max(max(abs(rel_dpi_AH_A)))]);
    piplot_rel_end_perc = ceil( piplot_rel_end / pireltol ) * pireltol;

    pirelregA = abs(rel_dpi_AH_A) < pireltol;
    pirelregH = abs(rel_dpi_AH_H) < pireltol;
    pirelregU = abs(rel_dpi_AH_U) < pireltol;
    pirelreg = -1 * ones(lx,ly);
    pirelreg_check = -1 * ones(lx,ly);
    for ix = 1:lx
        for iy = 1:ly
            if pirelregU(ix,iy)
                pirelreg(ix,iy) = 0;
            elseif pirelregA(ix,iy)
                pirelreg(ix,iy) = 1;
            elseif pirelregH(ix,iy)
                pirelreg(ix,iy) = 3;
            else
                pirelreg(ix,iy) = 4;
            end

            if ( xor(pirelregA(ix,iy),pirelregH(ix,iy)) ) && (~pirelregU(ix,iy))
                pirelreg_check(ix,iy) = 0;
            elseif ( pirelregA(ix,iy) && pirelregH(ix,iy) && (~pirelregU(ix,iy)) )
                pirelreg_check(ix,iy) = 1;
                pirelreg(ix,iy) = 2;
            end
        end
    end

    % Time of peak differences
    dt_AH_H = tAHsim - tHsim; % dt = difference in peak time
    dt_AH_A = tAHsim - tAsim;
    tplot_end = max([max(max(abs(dt_AH_H))),max(max(abs(dt_AH_A)))]);

    rel_dt_AH_H = ( tAHsim - tHsim ) ./ ( ( tAHsim + tHsim ) / 2 );
    rel_dt_AH_A = ( tAHsim - tAsim ) ./ ( ( tAHsim + tAsim ) / 2 );
    rel_dt_AH_U = ( tAHsim - tHsim(1,1) ) ./ ( ( tAHsim + tHsim(1,1) ) / 2 );
    tplot_rel_end = max([max(max(abs(rel_dt_AH_H))),max(max(abs(rel_dt_AH_A)))]);
    tplot_rel_end_perc = ceil( tplot_rel_end / treltol ) * treltol;

    trelregA = abs(rel_dt_AH_A) < treltol;
    trelregH = abs(rel_dt_AH_H) < treltol;
    trelregU = abs(rel_dt_AH_U) < treltol;
    trelreg = -1 * ones(lx,ly);
    trelreg_check = -1 * ones(lx,ly);
    for ix = 1:lx
        for iy = 1:ly
            if trelregU(ix,iy)
                trelreg(ix,iy) = 0;
            elseif trelregA(ix,iy)
                trelreg(ix,iy) = 1;
            elseif trelregH(ix,iy)
                trelreg(ix,iy) = 3;
            else
                trelreg(ix,iy) = 4;
            end

            if ( xor(trelregA(ix,iy),trelregH(ix,iy)) ) && (~trelregU(ix,iy))
                trelreg_check(ix,iy) = 0;
            elseif ( trelregA(ix,iy) && trelregH(ix,iy) && (~trelregU(ix,iy)) )
                trelreg_check(ix,iy) = 1;
                trelreg(ix,iy) = 2;
            end
        end
    end

    % Overall rejection regions (intersection)
    relregA = all(cat(3,zrelregA,pirelregA,trelregA),3);
    relregH = all(cat(3,zrelregH,pirelregH,trelregH),3);
    relregU = all(cat(3,zrelregU,pirelregU,trelregU),3);

else
    relregA = zrelregA;
    relregH = zrelregH;
    relregU = zrelregU;
    pirelreg = zeros(size(zrelreg));
    trelreg = zeros(size(zrelreg));
end

relreg = -1 * ones(lx,ly);
relreg_check = -1 * ones(lx,ly);
for ix = 1:lx
    for iy = 1:ly
        if relregU(ix,iy)
            relreg(ix,iy) = 0;
        elseif relregA(ix,iy)
            relreg(ix,iy) = 1;
        elseif relregH(ix,iy)
            relreg(ix,iy) = 3;
        else
            relreg(ix,iy) = 4;
        end

        if ( xor(relregA(ix,iy),relregH(ix,iy)) ) && (~relregU(ix,iy))
            relreg_check(ix,iy) = 0;
        elseif ( relregA(ix,iy) && relregH(ix,iy) && (~relregU(ix,iy)) )
            relreg_check(ix,iy) = 1;
            relreg(ix,iy) = 2;
        end
    end
end
zrelreg(isnan(zA)) = -1 * zrelreg(isnan(zA));
pirelreg(isnan(zA)) = -1 * pirelreg(isnan(zA));
trelreg(isnan(zA)) = -1 * trelreg(isnan(zA));
relreg(isnan(zA)) = -1 * relreg(isnan(zA));

output_list{1,1} = zrelreg';
subtitle_list{1,1} = 'Final size';

output_list{1,2} = pirelreg';
subtitle_list{1,2} = 'Daily peak incidence';

output_list{2,1} = trelreg';
subtitle_list{2,1} = 'Time to peak';

output_list{2,2} = relreg';
subtitle_list{2,2} = 'Overall (intersection)';


%% Figure

% %%%%%% Data
D.X = x_vec;
D.Y = y_vec;
D.Z = output_list;
D.Clim = [ -4.5, 4.5 ];
D.Climleg = [ -0.5, 4.5 ];
wm = 0.7; % white mix shade
D.colormap = [ 1 wm wm; 1 1 wm; wm 1 wm; wm 1 1; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0 ];

% %%%%%% Text
% Annoyingly, spaces in text are weird and system-dependent (differently in
% Windows and Mac), so the first lines show how it should work, while
% the following lines offer a solution for how the result should look, 
% (which works for me in Windows)
% Correct spacing (R2019a on Windows):
T.labelx = 'Adult-to-adult within-household transmission probability (p_{aa})';
T.labely = 'Relative susceptibility of children (\psi)';
% % Fudged spacing (R2016a on Mac):
% T.labelx = 'Adult-to-adult within-household transmission probability{(p_{aa})}';
% T.labely = 'Relative susceptibility of children versus adults { (\psi)}';
T.clabels = { 'Unstructured', 'Age', 'Either', 'Households', 'Both' };

T.subtitles = subtitle_list;
T.first_subletter = first_subletter;

% %%%%%% Layout
L.fig_width_cm = 7; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 landscape ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 2;
L.ncols = 2;
% Title
L.title_height_cm = 0;
L.xlabel_height_cm = 0.8;
L.ylabel_width_cm = 0.8;
T.labely_ycorr = -0.02;
% Legend
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0.7;
% Subtitles
L.subtitle_height_cm = 0.5;
L.colorbar = true;
L.hseparator_cm = 0;

T.axes_font_size = 6;
T.subtitle_font_size = 6;
T.subaxes_font_size = 5;
T.legend_font_size = 5;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'OAR_grid', D, L, T );

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    expcmd = ['export_fig ',fname,' -pdf -nocrop -transparent'];
    eval(expcmd);
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end

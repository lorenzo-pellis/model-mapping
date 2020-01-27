% This is the code to create many figures in the supplementary text of 
% Pellis, L. et al (2020), Nature Communications
% 
% It generates a 2x2 plot with heatmaps for a single output and for 3
% models in position top-left, top-right, and bottom-left. Then ont he
% bottom-right the simplest-model-acceptance-region for that specific
% output is plotted.
% It relies on a hand-made function called "mycomplexfig_mix.m"
% 
% Update: 13-10-2019

close all; % close all figures
clearvars;
Activate_save_fig = 1; % If true, figures are saved
Activate_C_codes = 1;

Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0
impose_same_clim = true; % If true, for easier comparison between different 
% parameter values (e.g. R0), a fixed limit for the color scale is used.

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

R0val = 1.5; first_subletter = 'a'; % a or m
% R0val = 2; first_subletter = 'e'; % e or q
% R0val = 4; first_subletter = 'i'; % i or u
% Choose between z, pi, t, sz
output = 'z'; sameclim = [ 40 100 ];
% output = 'pi'; sameclim = [ 1.5 17.5 ];
% output = 't'; sameclim = [ 4.8 15 ];
% output = 'sz'; sameclim = [ 30 100 ];

rval = 0.25282; % [ 0.14552 0.25282 0.52588 ] are the values corresponding to R0 = 1.5, 2 and 4;
% nsimval = 100;
phiGfig = 1; % The value of the relative infectivity of children VS adults 
tolval = 0.05; % The value of the tolerance used in the simplest-model-acceptance-region plot


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
            preloadedwrksp = 'r025282_pAA00_03_48_psi10_02_40_phi10_theta04_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
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
addpath(genpath(fig_tool_path));
cd(fig_code_path); % Work in the directory where the codes for figures are
x_vec = pAA_vec;
y_vec = psiG_vec;
[X,Y] = meshgrid(x_vec,y_vec);
lx = length(x_vec);
ly = length(y_vec);
if isnan(thetaGval)
    thetaGval = thetaG; 
end

% Build figure name
fname = ['2x2_mix_',output,'_',country,'_',popfig,'_R0', num2str(R0val*10,'%02d'),'_phi', num2str(phiGfig*10,'%02d')];
if use_match_r
    fname = [fname,'_match_r'];
end

% Workspace name
if use_match_r
    wname = [ 'r', num2str(round(rval*100000),'%06d') ];
else
    wname = [ 'R0', num2str(R0val*10,'%02d') ];
end
wname = [ wname,'_pAA', num2str(pAA_min*100,'%02d'),'_',num2str(dpAA*100,'%02d'),'_',num2str(pAA_max*100,'%02d'),...
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

% Stratified final size differences
dszAH_A_A = szAH_A - szA_A;
dszAH_A_C = szAH_C - szA_C;
szplot_end = max([max(max(abs(dszAH_A_A))),max(max(abs(dszAH_A_C)))]);

rel_dsz_AH_A_A = ( szAH_A - szA_A ) ./ ( ( szAH_A + szA_A ) / 2 );
rel_dsz_AH_A_C = ( szAH_C - szA_C ) ./ ( ( szAH_C + szA_C ) / 2 );
szplot_rel_end = max([max(max(abs(rel_dsz_AH_A_A))),max(max(abs(rel_dsz_AH_A_C)))]);
szplot_rel_end_perc = ceil( szplot_rel_end / szreltol ) * szreltol;

zAH_perc = 100*zAH;
zA_perc = 100*zA;
zH_perc = 100*zH;
zlimits = [ min([min(min(zAH_perc)),min(min(zA_perc)),min(min(zH_perc))]), max([max(max(zAH_perc)),max(max(zA_perc)),max(max(zH_perc))]) ];
pilimits = [ min([min(min(piAHsim)),min(min(piAsim)),min(min(piHsim))]), max([max(max(piAHsim)),max(max(piAsim)),max(max(piHsim))]) ];
tlimits = [ min([min(min(tAHsim)),min(min(tAsim)),min(min(tHsim))]), max([max(max(tAHsim)),max(max(tAsim)),max(max(tHsim))]) ];
szAH_A_perc = 100*szAH_A;
szAH_C_perc = 100*szAH_C;
szA_A_perc = 100*szA_A;
szA_C_perc = 100*szA_C;
szlimits = [ min([min(min(szAH_A_perc)),min(min(szA_A_perc)),min(min(szAH_C_perc)),min(min(szA_C_perc))]), max([max(max(szAH_A_perc)),max(max(szA_A_perc)),max(max(szAH_C_perc)),max(max(szA_C_perc))]) ];

if strcmp(output,'z')
    output_list{1,1} = zAH_perc';
    subtitle_list{1,1} = 'Model AH';
    % clabel_list{1,1} = 'Final size (% of population)';

    output_list{1,2} = zA_perc';
    subtitle_list{1,2} = 'Model A';

    output_list{2,1} = zH_perc';
    subtitle_list{2,1} = 'Model H';

    output_list{2,2} = zrelreg';
%     subtitle_list{2,2} = 'Simplest model acceptance regions';
    subtitle_list{2,2} = 'Acceptance regions';

    title = [ 'Average epidemic final size (%): R_0 = ',num2str(R0val),', \phi = ',num2str(phiGfig) ];
    limc = zlimits;
elseif strcmp(output,'pi')
    output_list{1,1} = piAHsim';
    subtitle_list{1,1} = 'Model AH';
    % clabel_list{1,1} = 'Final size (% of population)';

    output_list{1,2} = piAsim';
    subtitle_list{1,2} = 'Model A';

    output_list{2,1} = piHsim';
    subtitle_list{2,1} = 'Model H';

    output_list{2,2} = pirelreg';
    subtitle_list{2,2} = 'Acceptance regions';

    title = ['Average daily peak incidence (%): R_0 = ',num2str(R0val),', \phi = ',num2str(phiGfig) ];
    limc = pilimits;
elseif strcmp(output,'t')
    output_list{1,1} = tAHsim';
    subtitle_list{1,1} = 'Model AH';
    % clabel_list{1,1} = 'Final size (% of population)';

    output_list{1,2} = tAsim';
    subtitle_list{1,2} = 'Model A';

    output_list{2,1} = tHsim';
    subtitle_list{2,1} = 'Model H';

    output_list{2,2} = trelreg';
    subtitle_list{2,2} = 'Acceptance regions';

    title = [ 'Average time to peak incidence (generations): R_0 = ',num2str(R0val),', \phi = ',num2str(phiGfig) ];
    limc = tlimits;
elseif strcmp(output,'sz')
    output_list{1,1} = szAH_A_perc';
    subtitle_list{1,1} = 'Adults - model AH';
    % clabel_list{1,1} = 'Final size (% of population)';

    output_list{1,2} = szA_A_perc';
    subtitle_list{1,2} = 'Adults - model A';

    output_list{2,1} = szAH_C_perc';
    subtitle_list{2,1} = 'Children - model AH';

    output_list{2,2} = szA_C_perc';
    subtitle_list{2,2} = 'Children - model A';

    title = [ 'Age-stratified final size (% of age class): R_0 = ',num2str(R0val),', \phi = ',num2str(phiGfig) ];
    limc = szlimits;
else
    disp( 'Problems: no valid output specified' );
end


%% Figure

% %%%%%% Data
D.X = x_vec;
D.Y = y_vec;
D.Z = output_list;
if impose_same_clim
    D.Clim = sameclim;
else
    D.Clim = limc;
end

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
T.clabels = { 'U', 'A', 'A or H', 'H', 'A and H' };

% if strcmp(output,'sz')
%     T.title_xcorr = -0.01;
% end
% T.title_xcorr = 0.035;
% T.labelx_xcorr = 0.05;
% T.labely_ycorr = 0.01;

T.subtitles = subtitle_list;
T.first_subletter = first_subletter;

% %%%%%% Layout
L.fig_width_cm = 7; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 2;
L.ncols = 2;
% Title
L.title_height_cm = 0;
L.xlabel_height_cm = 0.8;
L.ylabel_width_cm = 0.7;
% T.labely_ycorr = -0.03;
% Legend
L.vrlegend_width_cm = 0.2;
% L.hblegend_height_cm = 0.7;
% Subtitles
L.subtitle_height_cm = 0.5;
L.subtitle_ypos = 0.55; 
% Sublegend
L.subcolorbar = true;
L.htsublegend_height_cm = 0.5;
% L.subcolorbar_bottom_limit_in_box = 1/4;
L.hseparator_cm = 0;

%%%%% Text size
T.axes_font_size = 6;
T.subtitle_font_size = 6;
T.subaxes_font_size = 5;
T.legend_font_size = 5;
% if strcmp(output,'t')
%     T.title_xcorr = -0.03;
% end

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

if strcmp(output,'sz')
    subplot_handles = mycomplexfig( 'heatmap', D, L, T );
else
    subplot_handles = mycomplexfig_mix( D, L, T );
end

if Activate_save_fig
    cd(fig_fig_path);
    expcmd = ['export_fig ',fname,' -pdf -nocrop -transparent'];
    eval(expcmd);
    cd(fig_code_path);
end
rmpath(genpath(fig_tool_path));

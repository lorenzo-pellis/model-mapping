% This is the code to create the left part of Figure 2 of the main text of 
% Pellis, L et al (2020), Nature Communications
% 
% It relies on a hand-made function called "mycomplexfig"
% 
% Update: 09-10-2019

close all; % close all figures
clearvars; % clear all variables
Activate_save_fig = 1; % If true, figures are saved
Activate_C_codes = true;
Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0
use_marks = true; xm = [ 0.15, 0.75 ]; ym = [ 3, 1.2 ]; firstm = 'e';

tolval = 0.05;
first_subletter = 'a';

% Path stuff
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
    if use_match_r
        fig_wrksp_path = [fig_wrksp_path,'match-r\'];
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
    if use_match_r
        fig_wrksp_path = [fig_wrksp_path,'match-r/'];
    end
end
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
% Loaded workspaces generated from other computers, with a different folder
% structure might have saved anonymous functions with paths not recognised 
% when this script is run, so the line above turns off that warning.
if use_match_r
    wrks_name = 'r025282_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
else
    wrks_name = 'R020_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050';
end
load([fig_wrksp_path,wrks_name]); 
cd(fig_code_path); % Work in the directory where the codes for figures are

x_vec = pAA_vec; lx = length(x_vec);
y_vec = psiG_vec; ly = length(y_vec);

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

relreg = -1 * ones(lx,ly); % Initialise to a value not used
relreg_check = -1 * ones(lx,ly); % Initialise to a value not used
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

output_list{1,1} = zrelreg';
subtitle_list{1,1} = 'Final size';

output_list{1,2} = pirelreg';
subtitle_list{1,2} = 'Peak daily incidence';

output_list{2,1} = trelreg';
subtitle_list{2,1} = 'Time to peak';

output_list{2,2} = relreg';
subtitle_list{2,2} = 'Overall (intersection)';


%% Figure
% To constuct the figure, I pass 3 structures:
%   - D = Data
%   - T = Text (title(s), axe labels, legend text, etc.)
%   - L = Layout information for the figure spatial structure

[X,Y] = meshgrid(x_vec,y_vec);

% %%%%%% Data
D.X = x_vec;
D.Y = y_vec;
D.Z = output_list;
D.Clim = [ -0.5, 4.5 ];
D.colormap = [ 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0 ];

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
L.fig_height_cm = 7.5; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 2;
L.ncols = 2;
% Title
L.xlabel_height_cm = 0.7;
L.ylabel_width_cm = 0.8;
T.labely_ycorr = -0.02;
% Legend
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0.7;
% Subtitles
L.subtitle_height_cm = 0.5;
L.colorbar = true;

T.title_font_size = 7.5;
T.axes_font_size = 6;
T.subtitle_font_size = 6;
T.subaxes_font_size = 4;
T.legend_font_size = 5;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'OAR_grid', D, L, T );

if use_marks
    sph = subplot_handles(2,2);
    im = 0;
    orange = [ 1 0.6 0 ];
    corr = [ 0, 0.01; 0 0; 0 0.03; 0 0 ];
    for bm = 1:2
        for am = 1:2
            im = im+1;
            lm = char(firstm+im-1); % letter
            axes( sph );
            bw = 1.5 * ( x_vec(2) - x_vec(1) ); % width
            bh = 1.5 * ( y_vec(2) - y_vec(1) ); % height
            rectangle( 'position',[ xm(am) - bw / 2, ym(bm) - bh / 2, bw, bh ],...
                'EdgeColor', orange,...
                'Linewidth',1.5,'FaceColor','w' );
            ht = text( xm(am) + corr(im,1), ym(bm) - bh/3  + corr(im,2), lm, 'Color', orange, 'EdgeColor', 'none',...
                'FontSize',L.screen_scale * 5,'Fontweight','bold',...
                'HorizontalAlignment','center','VerticalAlignment','baseline' );
        end
    end
end

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    export_fig Main_Figure2left -pdf -nocrop -transparent;
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end

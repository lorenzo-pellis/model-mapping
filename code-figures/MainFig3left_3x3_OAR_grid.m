% This is the code to create the left part of Figure 3 of the main text of 
% Pellis, L et al (2020), Nature Communications
% 
% It relies on a hand-made function called "mycomplexfig"
% 
% Update: 09-10-2019

close all; % close all figures
clearvars;
Activate_C_codes = 1;
Activate_save_fig = 1; % If true, figures are saved
Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

usez = true; % If true, I include the final size (z) in the overall acceptance region plot, otherwise not
usepi = true; % If true, I include the peak incidence (pi) in the overall acceptance region plot, otherwise not
uset = true;  % If true, I include the time to the peak (t) in the overall acceptance region plot, otherwise not
useSAR = true; % If true, the SAR contours are overlaid on the plot
use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0
use_marks = 1;

country = 'GB';

R0vals = [ 1.5 2 4 ]; R0name = 'R0basic';
rvals = [ 0.14552 0.25282 0.52588 ]; % In case I want to match r
pop = '2ran'; thetaGval = 0.2273; gamval = 1; % Random
% pop = 'UK'; thetaGval = 0.58; gamval = 0.75; % UK
phivals = [ 1 1.5 2 ];
tolval = 0.05;
figletter = 'A';

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
        fig_wrksp_path = [fig_wrksp_path,'Match-r\'];
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
        fig_wrksp_path = [fig_wrksp_path,'Match-r/'];
    end
end
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
% Loaded workspaces generated from other computers, with a different folder
% structure might have saved anonymous functions with paths not recognised 
% when this script is run, so the line above turns off that warning.
% 
% Start by loading one workspace, to check paths are correct and to load
% values of common variables:
if use_match_r
    wrks_name = 'r025282_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
else
    wrks_name = 'R020_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050';
end
load([fig_wrksp_path,wrks_name]); 
cd(fig_code_path); % Work in the directory where the codes for figures are

x_vec = pAA_vec;
y_vec = psiG_vec;
if use_marks
    mark(1).plot = 1;
    mark(1).pos = [ 0.1, 2; 0.2, 2; 0.3, 2; 0.4, 2 ];
    mark(1).n = 4;
    mark(1).symb = '1';
    mark(1).corr = [ 0 0 ];

    mark(2).plot = 6;
    if strcmp(pop,'2ran')
        mark(2).pos = [ 0.3, 1.05 ];
    else % UK
        mark(2).pos = [ 0.3, 1.05 ];
    end
    mark(2).n = 1;
    mark(2).symb = '2';
    mark(2).corr = [ 0 0 ];

    mark(3).plot = 7;
    mark(3).pos = [ 0.3, 3.8; 0.35, 3.2; 0.4, 2.8; 0.45, 2.4; 0.5, 2; 0.55, 1.8; 0.6, 1.6; 0.65, 1.4; 0.7, 1.2 ];
    mark(3).n = 9;
    mark(3).symb = '3';
    mark(3).corr = [ 0 0 ];
end

[X,Y] = meshgrid(x_vec,y_vec);
lx = length(x_vec);
ly = length(y_vec);
lR0 = length(R0vals);
lphi = length(phivals);

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

for ip = 1:lphi
    for iR0 = 1:lR0
        % Build workspace name
        if use_match_r
            wrks_name = [ 'r', num2str(round(rvals(iR0)*100000),'%06d') ];
        else
            wrks_name = [ 'R0', num2str(R0vals(iR0)*10,'%02d') ];
        end
        wrks_name = [ wrks_name,'_pAA', num2str(round(pAA_min*100),'%02d'),'_',num2str(round(dpAA*100),'%02d'),'_',num2str(round(pAA_max*100),'%02d'),...
            '_psi', num2str(psiG_min*10,'%02d'), '_', num2str(dpsiG*10,'%02d'),'_', num2str(psiG_max*10,'%02d'),...
            '_phi', num2str(phivals(ip)*10,'%02d'), '_theta',num2str(round(thetaGval*10),'%02d'),...
            '_gammaG', num2str(gamval*100,'%03d'), '_H', num2str(gamval*100,'%03d') ];
        if Activate_C_codes
            wrks_name = [ wrks_name, '_100sim_e5_init', num2str(n_init_inf,'%03d') ];
        end
        if use_match_r
            wrks_name = [ wrks_name, '_rfixed_nsimr',num2str(nsimval,'%05d') ];
        end
        wrks_list{ip,iR0} = wrks_name;

        try_name = strcat(fig_wrksp_path,wrks_list{ip,iR0});
        load(try_name);

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

            if ( usez && usepi && uset )        
                relregA = all(cat(3,zrelregA,pirelregA,trelregA),3);
                relregH = all(cat(3,zrelregH,pirelregH,trelregH),3);
                relregU = all(cat(3,zrelregU,pirelregU,trelregU),3);
            elseif ( usez && usepi && ~uset )
                relregA = all(cat(3,zrelregA,pirelregA),3);
                relregH = all(cat(3,zrelregH,pirelregH),3);
                relregU = all(cat(3,zrelregU,pirelregU),3);
            elseif ( usez && ~usepi && uset )
                relregA = all(cat(3,zrelregA,trelregA),3);
                relregH = all(cat(3,zrelregH,trelregH),3);
                relregU = all(cat(3,zrelregU,trelregU),3);
            elseif ( ~usez && usepi && uset )
                relregA = all(cat(3,trelregA,pirelregA),3);
                relregH = all(cat(3,trelregH,pirelregH),3);
                relregU = all(cat(3,trelregU,pirelregU),3);
            elseif ( usez && ~usepi && ~uset )
                relregA = zrelregA;
                relregH = zrelregH;
                relregU = zrelregU;
            elseif ( ~usez && usepi && ~uset )
                relregA = pirelregA;
                relregH = pirelregH;
                relregU = pirelregU;
            elseif ( ~usez && ~usepi && uset )
                relregA = trelregA;
                relregH = trelregH;
                relregU = trelregU;
            else
                disp('Did you really mean not to use any output whatsoever for the overall acceptance region? Really?');
            end
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

        zrelreg(isnan(zA)) = -1 * zrelreg(isnan(zA));
        pirelreg(isnan(zA)) = -1 * pirelreg(isnan(zA));
        trelreg(isnan(zA)) = -1 * trelreg(isnan(zA));
        relreg(isnan(zA)) = -1 * relreg(isnan(zA));
        
        %%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        relreg_list{ip,iR0} = relreg';
        if use_match_r
            reltitle_list{ip,iR0} = ['r = ', num2str(desired_r,'%1.3f'),', \phi = ',num2str(phiG,'%1.1f')];
        else
            reltitle_list{ip,iR0} = ['R_0 = ',num2str(R0,'%1.1f'),', \phi = ',num2str(phiG,'%1.1f')];
        end
        if useSAR
            overlaid_list{ip,iR0} = SAR'*100;
        else
            overlaid_list{ip,iR0} = inVSout_prop(:,:,1)'*100;
        end
    end
end

%% Figure
% To constuct the figure, I pass 3 structures:
%   - D = Data
%   - T = Text (title(s), axe labels, legend text, etc.)
%   - L = Layout information for the figure spatial structure

% %%%%%% Data
D.X = x_vec;
D.Y = y_vec;
D.Z = relreg_list;
D.Clim = [ -4.5, 4.5 ];
D.Climleg = [ -0.5, 4.5 ];
wm = 0.7; % white mix shade, to highlight where the map fails
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

T.subtitles = reltitle_list;
T.figletter = figletter;

% %%%%%% Layout
L.fig_width_cm = 7; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 8; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 3;
L.ncols = 3;
% Title
L.xlabel_height_cm = 0.7;
L.ylabel_width_cm = 0.7;
T.labely_ycorr = -0.02;
% Legend
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0.7;
L.colorbar = true;
% Subtitles
L.subtitle_height_cm = 0.4;
T.subaxes_font_size = 4;
T.legend_font_size = 5;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'OAR_grid', D, L, T );

if use_marks
    sph = subplot_handles(2,2);
    orange = [ 1 0.6 0 ];
    bw = 1.5 * ( x_vec(2) - x_vec(1) ); % width
    bh = 1.5 * ( y_vec(2) - y_vec(1) ); % height
    for im = 1:3
        axes( subplot_handles( mark(im).plot ) );
        for jm = 1:mark(im).n
            xm = mark(im).pos(jm,1);
            ym = mark(im).pos(jm,2);
            rectangle( 'position',[ xm - bw / 2, ym - bh / 2, bw, bh ],...
                'EdgeColor', orange,...
                'Linewidth',1.5,'FaceColor','w' );
            ht = text( xm + mark(im).corr(1), ym - bh/3  + mark(im).corr(2), mark(im).symb, 'Color', orange, 'EdgeColor', 'none',...
                'FontSize',L.screen_scale * 4,'Fontweight','bold',...
                'HorizontalAlignment','center','VerticalAlignment','baseline' );
        end
    end
end

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    export_fig Main_Figure3left -pdf -nocrop -transparent;
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end


% This is the code to create many figures in the supplementary text of 
% Pellis, L et al (2019), Nature Communications
% 
% It generates a 3x3 graph with the overall simplest-model-acceptance-
% region plots for 3 values of R0 and 3 values of the relative infectivity
% of children VS adults (phi). There is no title, but subtitles are 
% generated for each row and column.
% It relies on a hand-made function called "mycomplexfig.m"
% 
% Update: 31-12-2019

close all; % close all figures
clearvars;
Activate_save_fig = 1; % If true, figures are saved
Activate_C_codes = 1;

Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

usez = true; % If true, I include the final size (z) in the overall acceptance region plot, otherwise not
usepi = true; % If true, I include the peak incidence (pi) in the overall acceptance region plot, otherwise not
uset = true;  % If true, I include the time to the peak (t) in the overall acceptance region plot, otherwise not
useSAR = true; % If true, the SAR contours are overlaid on the plot
use_match_r = 0; % If false, I match R0; if true, I use the r correspondent to the desired R0
use_marks = 0; % If true, I add a numbered mark in selected locations
use_intermediate = 0; % Leave this as default choice. If you want intermediate contact patterns, choose them below
use_log2psi = 0; % If true, I let log2psi run from -2 to 2 with steps of 0.2.

country = 'GB'; % Great Britain
% country = 'SL'; % Sierra Leone
% country = 'SA'; % South Africa

R0vals = [ 1.5 2 4 ]; R0name = 'R0basic';
% R0vals = [ 1.1 1.3 1.5 ]; R0name = 'R0down';
% R0vals = [ 1.7 2 2.3 ]; R0name = 'R0middle';
% R0vals = [ 2.7 3.2 4 ]; R0name = 'R0up';
popfig = '2ran'; thetaGval = NaN; gamval = 1; % Random
% popfig = 'm4r'; thetaGval = 0.4; gamval = 1; tolval = 0.05; figletter = 'A'; use_intermediate = true;
% popfig = 'm4UK'; thetaGval = 0.4; gamval = 0.75; tolval = 0.05; figletter = 'C'; use_intermediate = true;
% popfig = 'm5r'; thetaGval = 0.5; gamval = 1; tolval = 0.05; figletter = 'B'; use_intermediate = true;
% popfig = 'm5UK'; thetaGval = 0.5; gamval = 0.75; tolval = 0.05; figletter = 'D'; use_intermediate = true;
% popfig = 'UK'; thetaGval = 0.58; gamval = 0.75; % UK
% popfig = 'ass'; thetaGval = 0.7; gamval = 0.75; % UK
if ~use_log2psi
    phivals = [ 1 1.5 2 ];
else
    phivals = [ 0.5 1 2 ];
end
rvals = [ 0.14552 0.25282 0.52588 ]; % These are the values of r corresponding to the values of R0 = 1.5, 2 and 4
nsimval = 100;
if strcmp(country,'SL')
    if strcmp(popfig,'2ran')
        tolval = 0.05; figletter = 'A';
    else % UK
        tolval = 0.05; figletter = 'B';
    end
elseif strcmp(country,'SA')
    if strcmp(popfig,'2ran')
        tolval = 0.05; figletter = 'A';
    else % UK
        tolval = 0.05; figletter = 'B';
    end
else % GB        
    if ~use_intermediate
        tolval = 0.01; figletter = 'A'; % A or D
%         tolval = 0.05; figletter = 'B'; % B or E
%         tolval = 0.1; figletter = 'C'; % C or F
    end
end

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-figures\'];
    fig_path = [base_dir,'\output-figures\supp\',country,'\',popfig,'\'];
    tool_path = [base_dir,'\tools\'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'\output-workspaces\',country,'\'];
    else
        wrksp_path = [base_dir,'\saved-workspaces\',country,'\'];
    end
    if use_match_r
        wrksp_path = [wrksp_path,'match-r\'];
    end
else
    code_path = [base_dir,'/code-figures/'];
    fig_path = [base_dir,'/output-figures/supp/',country,'/',popfig,'/'];
    tool_path = [base_dir,'/tools/'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'/output-workspaces/',country,'/'];
    else
        wrksp_path = [base_dir,'/saved-workspaces/',country,'/'];
    end
    if use_match_r
        wrksp_path = [wrksp_path,'match-r/'];
    end
end
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
if use_log2psi
    if isnan(thetaGval) % If mixing is random
        preloadedwrksp = 'R020_pAA00_05_95_log2psi-200_025_200_phi10_theta02_gammaG100_H100_100sim_e5_init050';
    else
        preloadedwrksp = 'R020_pAA00_05_95_log2psi-200_025_200_phi10_theta06_gammaG075_H075_100sim_e5_init050';
    end
end
load([wrksp_path,preloadedwrksp]);
cd(code_path); % Work in the directory where the codes for figures are
x_vec = pAA_vec;
if use_log2psi
    y_vec = logpsiG_vec;
else
    y_vec = psiG_vec;
end
[X,Y] = meshgrid(x_vec,y_vec);
lx = length(x_vec);
ly = length(y_vec);
lR0 = length(R0vals);
lphi = length(phivals);
if isnan(thetaGval)
    thetaGval = thetaG; 
end


% if strcmp(country,'SL')
%     if isnan(thetaGval)
%         if ~use_match_r
%             preloadedwrksp = 'R020_pAA00_03_48_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050';
%         else
%             preloadedwrksp = 'r025282_pAA00_03_48_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
%         end
%     else
%         if ~use_match_r
%             preloadedwrksp = 'R020_pAA00_03_48_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050';
%         else
%             preloadedwrksp = 'r025282_pAA00_03_48_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050_rfixed_nsimr00100';
%         end
%     end
% elseif strcmp(country,'SA')
%     if isnan(thetaGval)
%         if ~use_match_r
%             preloadedwrksp = 'R020_pAA00_04_63_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050';
%         else
%             preloadedwrksp = 'r025282_pAA00_04_63_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
%         end
%     else
%         if ~use_match_r
%             preloadedwrksp = 'R020_pAA00_04_63_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050';
%         else
%             preloadedwrksp = 'r025282_pAA00_04_63_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050_rfixed_nsimr00100';
%         end
%     end
% else
%     if isnan(thetaGval)
%         if ~use_match_r
%             preloadedwrksp = 'R020_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050';
%         else
%             preloadedwrksp = 'r025282_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050_rfixed_nsimr00100';
%         end
%     else
%         if ~use_match_r
%             preloadedwrksp = 'R020_pAA00_05_95_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050';
%         else
%             preloadedwrksp = 'r025282_pAA00_05_95_psi10_02_40_phi10_theta06_gammaG075_H075_100sim_e5_init050_rfixed_nsimr00100';
%         end
%     end
% end
% 
% % Workspace folder
% if ispc
%     wrks_path = ['C:\Users\Lorenzo\Desktop\Model Mapping\Workspaces\',country,'\'];
%     if use_match_r
%         wrks_path = [wrks_path,'Match r\'];
%     end
% else
%     wrks_path = ['/Volumes/Macintosh data/lpellis/Desktop/3. Current projects/Model Mapping/Workspaces/',country,'/'];
%     if use_match_r
%         wrks_path = [wrks_path,'Match r/'];
%     end
% end
% load([wrks_path,preloadedwrksp]);
% x_vec = pAA_vec;
% y_vec = psiG_vec;



if use_marks
    mark(1).plot = 1;
    mark(1).pos = [ 0.15, 2];
    mark(1).n = 1;
    mark(1).symb = '1';
    mark(1).corr = [ 0 0 ];

    mark(2).plot = 2;
    if strcmp(popfig,'2ran')
        mark(2).pos = [ 0.45, 1.4 ];
    else % UK
        mark(2).pos = [ 0.45, 1.4 ];
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

% Build figure name
if ~use_log2psi
    fname = ['3x3_OAR_ROT_'];
else
    fname = ['3x3_OAR_ROT_log2psi_'];
end
if use_marks
    fname = [ fname, 'marks_' ];
end
if useSAR
    fname = [fname,'overSAR'];
else
    fname = [fname,'overTh'];
end
fname = [ fname,'_tol',num2str(round(tolval*100),'%02d'),'_'];
if usez
    fname = [fname,'z'];
end
if usepi
    fname = [fname,'pi'];
end
if uset
    fname = [fname,'t'];
end
fname = [fname,'_',country,'_',popfig,'_',R0name];
if use_match_r
    fname = [fname,'_match_r'];
end

for ip = 1:lphi
    for iR0 = 1:lR0
        % Build workspace name
        if use_match_r
            wname = [ 'r', num2str(round(rvals(iR0)*100000),'%06d') ];
        else
            wname = [ 'R0', num2str(R0vals(iR0)*10,'%02d') ];
        end
        if ~use_log2psi
            wname = [ wname,'_pAA', num2str(ceil(pAA_min*100),'%02d'),'_',num2str(ceil(dpAA*100),'%02d'),'_',num2str(ceil(pAA_max*100),'%02d'),...
                '_psi', num2str(psiG_min*10,'%02d'), '_', num2str(dpsiG*10,'%02d'),'_', num2str(psiG_max*10,'%02d'),...
                '_phi', num2str(phivals(ip)*10,'%02d'), '_theta',num2str(round(thetaGval*10),'%02d'),...
                '_gammaG', num2str(gamval*100,'%03d'), '_H', num2str(gamval*100,'%03d') ];
        else
            wname = [ wname,'_pAA', num2str(round(pAA_min*100),'%02d'),'_',num2str(round(dpAA*100),'%02d'),'_',num2str(round(pAA_max*100),'%02d'),...
                '_log2psi', num2str(logpsiG_min*100,'%-3d'), '_', num2str(logdpsiG*100,'%03d'),'_', num2str(logpsiG_max*100,'%03d'),...
                '_phi', num2str(phivals(ip)*10,'%02d'), '_theta',num2str(round(thetaGval*10),'%02d'),...
                '_gammaG', num2str(gamval*100,'%03d'), '_H', num2str(gamval*100,'%03d') ];
        end
        if Activate_C_codes
            wname = [ wname, '_100sim_e5_init', num2str(n_init_inf,'%03d') ];
        end
        if use_match_r
            wname = [ wname, '_rfixed_nsimr',num2str(nsimval,'%05d') ];
        end
        wrksp_list{ip,iR0} = wname;

        try_name = strcat(wrksp_path,wrksp_list{ip,iR0});
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

        if ~use_log2psi
            zrelreg(isnan(zA)) = -1 * zrelreg(isnan(zA));
            pirelreg(isnan(zA)) = -1 * pirelreg(isnan(zA));
            trelreg(isnan(zA)) = -1 * trelreg(isnan(zA));
            relreg(isnan(zA)) = -1 * relreg(isnan(zA));
        else
            zrelreg(isnan(zA)) = -1 * zrelreg(isnan(zA)) - 1;
            pirelreg(isnan(zA)) = -1 * pirelreg(isnan(zA)) - 1;
            trelreg(isnan(zA)) = -1 * trelreg(isnan(zA)) - 1;
            relreg(isnan(zA)) = -1 * relreg(isnan(zA)) - 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         zrelreg_list{iR0,ip} = zrelreg';
%         zreltitle_list{iR0,ip} = ['R_0 = ',num2str(R0,'%1.1f'),', \phi = ',num2str(phiG,'%1.1f')];
%         pirelreg_list{iR0,ip} = pirelreg';
%         pireltitle_list{iR0,ip} = ['R_0 = ',num2str(R0,'%1.1f'),', \phi = ',num2str(phiG,'%1.1f')];
%         trelreg_list{iR0,ip} = trelreg';
%         treltitle_list{iR0,ip} = ['R_0 = ',num2str(R0,'%1.1f'),', \phi = ',num2str(phiG,'%1.1f')];
        relreg_list{ip,iR0} = relreg';
        if use_match_r
%             reltitle_list{ip,iR0} = ['r = ', num2str(desired_r,'%1.3f'),', \phi = ',num2str(phiG,'%1.1f')];
            coltitle_list{iR0} = ['r = ',num2str(desired_r,'%1.3f')];
        else
%             reltitle_list{ip,iR0} = ['R_0 = ',num2str(R0,'%1.1f'),', \phi = ',num2str(phiG,'%1.1f')];
            coltitle_list{iR0} = ['R_0 = ',num2str(R0)];
        end
        rowtitle_list{ip} = ['\phi = ',num2str(phiG)];
        relregA_list{ip,iR0} = relregA';
        if useSAR
            overlaid_list{ip,iR0} = SAR'*100;
        else
            overlaid_list{ip,iR0} = inVSout_prop(:,:,1)'*100;
        end
    end
end

%% Figure

% %%%%%% Data
D.X = x_vec;
D.Y = y_vec;
D.Z = relreg_list;
D.nZover = 3;
D.Zover = { overlaid_list, overlaid_list, relregA_list };
D.Zcont = { [10 40 70 ], [ 20 30 50 60 80 90 ], 0.05 };
D.Zlinewidth = { 0.3, 0.3, 0.6 };
D.Zshowtext = { 'on', 'off', 'off' };
D.Climleg = [ -0.5, 4.5 ];
wm = 0.7; % white mix shade
if ~use_log2psi
    D.Clim = [ -4.5, 4.5 ];
    D.colormap = [ 1 wm wm; 1 1 wm; wm 1 wm; wm 1 1; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0 ];
    D.Yticks = 1:0.5:4;
else
    D.Clim = [ -5.5, 4.5 ];
    D.colormap = [ 1 wm wm; 1 1 wm; wm 1 wm; wm 1 1; wm wm 1; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0 ];
    D.Yticks = -2:2;
end

if strcmp(country,'SL')
    D.Xticks = [ 0 0.1 0.2 0.3 0.4 ];
elseif strcmp(country,'SA')
    D.Xticks = [ 0 0.15 0.3 0.45 0.6 ];
else
    D.Xticks = [ 0 0.2 0.4 0.6 0.8 ];
end

% %%%%%% Text
% Annoyingly, spaces in text are weird and system-dependent (differently in
% Windows and Mac), so the first lines show how it should work, while
% the following lines offer a solution for how the result should look, 
% (which works for me in Windows)
% T.labelx = 'Adult-to-adult within-household transmission probability (p_{aa})';
% T.labely = 'Relative susceptibility of children (\psi)';
T.labelx = 'Adult-to-adult within-household transmission probability{(p_{aa})}';
if ~use_log2psi
    T.labely = 'Relative susceptibility of children versus adults { (\psi)}';
else
    T.labely = 'Log_2 relative susceptibility of children versus adults{ (\psi)}   ';
end
T.clabels = { 'Unstructured', 'Age', 'Either', 'Households', 'Both' };

% T.subtitles = reltitle_list;
T.coltitles = coltitle_list;
T.rowtitles = rowtitle_list;
T.figletter = figletter;

% %%%%%% Layout
L.fig_width_cm = 7; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 3;
L.ncols = 3;
% Title
L.title_height_cm = 0;
L.xlabel_height_cm = 0.8;
L.ylabel_width_cm = 0.8;
T.labely_ycorr = -0.01;
% Legend
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0.7;
% Subtitles
% L.subtitle_height_cm = 0;
L.coltitle_height_cm = 0.5;
L.rowtitle_width_cm = 0.4;
L.colorbar = true;

T.axes_font_size = 6;
T.subtitle_font_size = 6;
T.subaxes_font_size = 5;
T.legend_font_size = 5;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'OAR_over', D, L, T );

if use_marks
    sph = subplot_handles(2,2);
    orange = [ 1 0.75 0 ];
    bw = 1.6 * ( x_vec(2) - x_vec(1) ); % width
    bh = 1.6 * ( y_vec(2) - y_vec(1) ); % height
    for im = 1:3
        axes( subplot_handles( mark(im).plot ) );
        for jm = 1:mark(im).n
            xm = mark(im).pos(jm,1);
            ym = mark(im).pos(jm,2);
            rectangle( 'position',[ xm - bw / 2, ym - bh / 2, bw, bh ],...
                'EdgeColor', orange,...
                'Linewidth',1.5,'FaceColor','w' );
            ht = text( xm + mark(im).corr(1), ym - bh/3  + mark(im).corr(2), mark(im).symb, 'Color', orange, 'EdgeColor', 'none',...
                'FontSize',L.screen_scale * 3.4,'Fontweight','bold',...
                'HorizontalAlignment','center','VerticalAlignment','baseline' );
        end
    end
end

if Activate_save_fig
    addpath(genpath(tool_path));
    cd(fig_path);
    expcmd = ['export_fig ',fname,' -pdf -nocrop -transparent'];
    eval(expcmd);
    cd(code_path);
    rmpath(genpath(tool_path));
end

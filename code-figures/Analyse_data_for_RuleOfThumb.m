% This code analyses the overall simplest model acceptance regions and the
% values of the SAR along the line separating the regions in parameter
% space where households should be included or are not necessary for
% accurate predictions. This is then used to generate the Rule Of Thumb
% depicted in Figure 4 in the Main text of
% Pellis, L. et al (2019), Nature Communications
% 
% Update: 31-12-2019

close all; % close all figures
clearvars;
Activate_C_codes = 1; % If true, the ROT includes using the peak incidence and time to the peak
Activate_fig_check = 0; % If true, the code plots a figure for every loop, to check the code is doing the right thing
Activate_clean_output = 1; % If true, only the main contour line is followed, discarding small closed contour lines
% that appear in plots with significant noise, e.g. due to stochastic variability
Activate_plot_from_new_workspaces = 0; 
% If 0, I make plots from pre-computed and saved workspaces (folder saved-workspaces)
% If 1, I make plots from newly computed workspaces (folder output-workspaces)

usez = true; % If true, I include the final size (z) in the overall acceptance region plot, otherwise not
usepi = true; % If true, I include the peak incidence (pi) in the overall acceptance region plot, otherwise not
uset = true;  % If true, I include the time to the peak (t) in the overall acceptance region plot, otherwise not
use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0

country = 'GB'; % Great Britain
% country = 'SL'; % Sierra-Leone
% country = 'SA'; % South-Africa

psirange = [ 1 4 ];
R0vals = [ 1.1 1.3 1.5 1.7 2 2.3 2.7 3.2 4 ]; R0name = 'R0all';
% R0vals = [ 1.3 1.5 1.7 2 2.3 2.7 3.2 ]; R0name = 'R0all-2';
% R0vals = [ 1.5 1.7 2 2.3 2.7 ]; R0name = 'R0all-4';
% R0vals = [ 1.5 2 4 ]; R0name = 'R0basic';
% R0vals = [ 1.1 1.3 1.5 ]; R0name = 'R0down';
% R0vals = [ 1.7 2 2.3 ]; R0name = 'R0middle';
% R0vals = [ 2.7 3.2 4 ]; R0name = 'R0up';
% R0vals = [ 1.3 1.5 1.7 ]; R0name = 'R0best3';
% R0vals = [ 1.3 1.5 1.7 2 ]; R0name = 'R0best4';
% R0vals = [ 1.3 1.5 1.7 2 2.3 ]; R0name = 'R0best5';
popROT = '2ran'; thetaGval = NaN; gamval = 1; % Random
% popROT = 'm4r'; thetaGval = 0.4; gamval = 1; 
% popROT = 'm4UK'; thetaGval = 0.4; gamval = 0.75; 
% popROT = 'm5r'; thetaGval = 0.5; gamval = 1; 
% popROT = 'm5UK'; thetaGval = 0.5; gamval = 0.75; 
% popROT = 'UK'; thetaGval = 0.58; gamval = 0.75; % UK
% popROT = 'ass'; thetaGval = 0.7; gamval = 0.75; % highly assortative, only for SL
phivals = [ 1 1.5 2 ];
tolval = 0.05; % Possible values: 0.01, 0.05 and 0.1

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-figures\'];
    ROT_path = [base_dir,'\output-rule-of-thumb\'];
    tool_path = [base_dir,'\tools\'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'\output-workspaces\',country,'\'];
    else
        wrksp_path = [base_dir,'\saved-workspaces\',country,'\'];
    end
else
    code_path = [base_dir,'/code-model-mapping/'];
    ROT_path = [base_dir,'/output-rule-of-thumb/'];
    tool_path = [base_dir,'/tools/'];
    if Activate_plot_from_new_workspaces
        wrksp_path = [base_dir,'/output-workspaces/',country,'/'];
    else
        wrksp_path = [base_dir,'/saved-workspaces/',country,'/'];
    end
end
% Start by loading the workspace with random mixing (which depends on the 
% population), to check paths are correct, store the right assortativity 
% for random mixing and to load values of common variables:
if strcmp(country,'GB')
    wrks_name = 'R020_pAA00_05_95_psi10_02_40_phi10_theta02_gammaG100_H100_100sim_e5_init050';
elseif strcmp(country,'SA')
    wrks_name = 'R020_pAA00_04_63_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050';
elseif strcmp(country,'SL')
    wrks_name = 'R020_pAA00_03_48_psi10_02_40_phi10_theta05_gammaG100_H100_100sim_e5_init050';
else
    error('No valid country specified!');
end
load([wrksp_path,wrks_name]); 
if isnan(thetaGval)
    thetaGval = thetaG;
end
addpath(tool_path)
cd(code_path); % Work in the directory where the codes for figures are

labelx = 'p_{AA}'; x_vec = pAA_vec;
labely = '\psi'; y_vec = psiG_vec;

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
% Pre-compilation of a table of binomial coefficients.
% I need them for future calculations, but they are expensive to compute,
% so I will do it now once and store such values in a table. Such table is
% aglobal variable, so I can access it from all functions that need it.
global bincoeffmat;
if ( ~exist('bincoeffmat','var') || all( size(bincoeffmat) ~= size(H) ) )
    % If for some reason it's not available, construct the matrix of binomial coefficients
    [ max_inA,max_inC ] = size(H);
    bincoeffmat = NaN(size(H));
    for inA = 1:max_inA
        for inC = 1:max_inC
            nA = inA-1;
            nC = inC-1;
            if inA >= inC
                bincoeffmat(inA,inC) = nchoosek(nA,nC);
            end
        end
    end
end

% Build ROT file name (output of the analysis)
ROT_name = [ 'ROT_tol',num2str(round(tolval*100),'%02d'),'_'];
if usez
    ROT_name = [ROT_name,'z'];
end
if usepi
    ROT_name = [ROT_name,'pi'];
end
if uset
    ROT_name = [ROT_name,'t'];
end
ROT_name = [ROT_name,'_',country,'_',popROT,'_',R0name];
ROT_name = [ROT_name,'_psirange',num2str(round(psirange(1)*10),'%02d'),'_',num2str(round(psirange(2)*10),'%02d')];
if Activate_clean_output
    ROT_name = [ROT_name,'_clean'];
else
    ROT_name = [ROT_name,'_dirty'];
end

nrowsROT = 5;
dotsROT = NaN(nrowsROT,0);
% dotsROT is the main variable where I store the output. Each columns
% contains [ R0, phi, pAA, psi, SAR ]^T, where pAA and psi come from the
% countour function, SAR is calculated in each of these values, and then I
% loop through each value of phi (for fixed R0), and then through each R0.
progressbar(0,0,0);
for iR0 = 1:lR0
    R0 = R0vals(iR0);
    for ip = 1:lphi
        % Load workspace
        wrks_name = [ 'R0', num2str(R0vals(iR0)*10,'%02d') ];
        wrks_name = [ wrks_name,'_pAA', num2str(pAA_min*100,'%02d'),'_',num2str(ceil(dpAA*100),'%02d'),'_',num2str(ceil(pAA_max*100),'%02d'),...
            '_psi', num2str(psiG_min*10,'%02d'), '_', num2str(dpsiG*10,'%02d'),'_', num2str(psiG_max*10,'%02d'),...
            '_phi', num2str(phivals(ip)*10,'%02d'), '_theta',num2str(round(thetaGval*10),'%02d'),...
            '_gammaG', num2str(gamval*100,'%03d'), '_H', num2str(gamval*100,'%03d') ];
        if Activate_C_codes
            wrks_name = [ wrks_name, '_100sim_e5_init', num2str(n_init_inf,'%03d') ];
        end
        load( [wrksp_path, wrks_name] );
         
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

            if ( usez && usepi && uset )        
                maxreldiffAH_A = max(max(abs(rel_dz_AH_A),abs(rel_dpi_AH_A)),abs(rel_dt_AH_A));
            elseif ( usez && usepi && ~uset )
                maxreldiffAH_A = max(abs(rel_dz_AH_A),abs(rel_dpi_AH_A));
            elseif ( usez && ~usepi && uset )
                maxreldiffAH_A = max(abs(rel_dz_AH_A),abs(rel_dt_AH_A));
            elseif ( ~usez && usepi && uset )
                maxreldiffAH_A = max(abs(rel_dpi_AH_A),abs(rel_dt_AH_A));
            elseif ( usez && ~usepi && ~uset )
                maxreldiffAH_A = abs(rel_dz_AH_A);
            elseif ( ~usez && usepi && ~uset )
                maxreldiffAH_A = abs(rel_dpi_AH_A);
            elseif ( ~usez && ~usepi && uset )
                maxreldiffAH_A = abs(rel_dt_AH_A);
            else
                error('Did you really mean not to use any output whatsoever for the plot? Really?!?')
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
                error('Did you really mean not to use any output whatsoever for the plot? Really?!?')
            end
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
                
        if Activate_fig_check
            figure(10);
            clf; 
            cm = colormap( [ 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0 ] );
            imagesc( x_vec, y_vec, relreg' )
            set(gca,'YDir','normal')
            caxis( [ 0 4 ] );
            set( gca, 'TickDir', 'out', 'TickLength', [ 0.02, 0.025 ] );
            hold on
            contour( x_vec, y_vec, maxreldiffAH_A'*100 , [ tolval, tolval ]*100, 'ShowText', 'on', 'Color', 'k', 'Linewidth', 2, 'LabelSpacing', 100 );
            contour( x_vec, y_vec, SAR'*100, [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ]*100, 'Color', 'k', 'Linewidth', 1 );
            xlabel(labelx);
            ylabel(labely);
        end
        Craw = contourc( x_vec, y_vec, maxreldiffAH_A' , [ tolval, tolval ] );
        
        currind = 1;
        ind1 = [];
        while currind <= size(Craw,2)
            assert( Craw(1,currind) == tolval );
            ind1 = [ ind1, currind ];
            currind = currind + Craw(2,currind) + 1;
        end
        if Activate_clean_output
            [ m, ind2 ] = max( Craw(2,ind1) );
            C = Craw(:,(ind1(ind2)+1):(ind1(ind2)+Craw(2,ind1(ind2))));
        else
            Craw(:,ind1) = [];
            C = Craw;
        end
        C(:,C(2,:)<psirange(1)) = [];
        C(:,C(2,:)>psirange(2)) = [];
        
        %%%%%%%%%% Run the code to compute the exact SAR along the line
        R0 = R0vals(iR0);
        phiG = phivals(ip);

        lC = size(C,2);
        MROT = NaN(nrowsROT,lC);
        MROT(1,:) = R0;
        MROT(2,:) = phiG;
        MROT(3:4,:) = C;
        
        for iC = 1:lC            

            pAA = C(1,iC);
            psiG = C(2,iC);
            psiH = psiG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% Same bits as main code in Model_Mapping_code.m %%%%%%%%%%%%%%%

            % Construct the shape of the global NGM. Then you only need to multiply
            % it by a coefficient
            NGM_G_shape = create_global_NGM_shape(F_ratio,psiG,phiG,thetaG,g_ratio);

            % Note: I wanted to allow for assortativity in the household, but it
            % turns out the shape depends on the household composition...
            % So I use thetaH = NaN for random mixing within each household 
            % (computed differently for each household composition).
            % To pass all the information I need for when I have to compute the
            % assortativity of each household composition, I store all relevant
            % parameters in a vector:
            shapeHvector = [ h_ratio, thetaH, psiH, phiH ];

            %%%%%%%%%%%%%%%%%% Find the right Rh for this pAA %%%%%%%%%%%%%%%
            if pAA == 0
                Rh = 0;
            else
                funpAA = @(x) ( get_pAA_from_Rh(H_single,x) - pAA );
                Rh = fzeromin( funpAA, 0, [], [], [], 1 );
            end

            % Construct the chains of average number of cases in each
            % generation (it's a 2 x Hsize x 2 tensor). See inside of
            % function for details. The second and third outputs are
            % not very important, but it's cheaper to compute them 
            % all in one go. av_gen_chains is key for computing R0
            [ av_gen_chains_2type, hfs_by_comp_by_ic, SAR_by_comp_by_ic ] = create_chains_hfs_and_SAR_2type(PI,Rh,shapeHvector,eta,bincoeffmat); % Here is where I decide how the infectivity scales with the household size

            % Find the Rg needed to achieve the selected R0.
            % shape = diag( F.^(-1) ) * NGM_shape; % This is the key relationship between shape and NGM_shape
            funAH = @(x) ( R0AH(av_gen_chains_2type,x*NGM_G_shape) - R0 );
            Rg = fzeromin( funAH, 0, [], [], [], 1 );

%             Lh = Phrates1to1_byHcomp_MRC(Hcomp,Rh,shapeHvector,eta); % 1to1 rates in households, by household composition
            disp(['Loop ',num2str(iC),' of ',num2str(lC)]);
            NGM_G = Rg * NGM_G_shape;
            R0_AH = R0AH(av_gen_chains_2type,Rg*NGM_G_shape);

            %%%%%%%%%%%%%%%%%%%% computing v_AH %%%%%%%%%%%%%%%%%%%%%%%

            [ R1, v_inc ] = R0AH_and_vAH(av_gen_chains_2type,NGM_G);
            % Proportions of adults/children among the household primary cases:
            vH_inc = NGM_G * v_inc / sum( NGM_G * v_inc ); 
            % Average generation chains we expect to see during the exponentially growing phase:
%             av_H_gen = sum( av_gen_chains_2type(:,:,1) * vH_inc(1) + av_gen_chains_2type(:,:,2) * vH_inc(2) );
            SAR_ROT_temp = [ sum( sum( SAR_by_comp_by_ic(:,:,1) .* PI(:,:,1) ) ), sum( sum( SAR_by_comp_by_ic(:,:,2) .* PI(:,:,2) ) ) ];
            SAR_ROT = SAR_ROT_temp * vH_inc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            MROT(5,iC) = SAR_ROT;
            
            frac3 = iC/lC;
            frac2 = (frac3-1+ip)/lphi;
            frac1 = (frac2-1+iR0)/lR0;
            
            progressbar(frac1,frac2,frac3);
        end
        
        dotsROT = [ dotsROT, MROT ];
        
    end
end
rmpath(tool_path)

cd(ROT_path);
save(ROT_name);
cd(code_path); % Return to the coding directory



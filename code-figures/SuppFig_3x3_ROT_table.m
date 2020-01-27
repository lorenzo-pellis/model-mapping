% This is the code to create Figure 24 in the supplementary text of 
% Pellis, L. et al (2020), Nature Communications
% 
% It generates a 3x3 graph with, in each panel, a table with mulitple
% values of R0 (each row for a different starting value and each column for
% a different ending value), to identify the intervals in R0 values over
% which to do a linear regression. The numbers in each cell show the
% coefficients (top: regression coefficient; bottom: intercept) and the
% colour in each cell the adjusted R^2. Each panel refers to a different
% value of the relative infectivity \phi (row) and range of relative 
% susceptibility \psi (column) of children VS adult. There is no title, but
% subtitles are generated for each row and column.
% It relies on a hand-made function called "mycomplexfig.m"
% 
% Update: 26-01-2020

close all; % close all figures
clearvars;
Activate_save_fig = 1; % If true, figures are saved
Activate_C_codes = 1;
Activate_clean_output = 1;

Activate_plot_from_new_ROT_analysis = 0; 
% If 0, I make plots from pre-computed and saved workspaces for the rule-of-thumb analysis (folder saved-rule-of-thumb)
% If 1, I make plots from newly computed workspaces for the rule-of-thumb analysis (folder output-rule-of-thumb)

usez = true; % If true, I include the final size (z) in the overall acceptance region plot, otherwise not
usepi = true; % If true, I include the peak incidence (pi) in the overall acceptance region plot, otherwise not
uset = true;  % If true, I include the time to the peak (t) in the overall acceptance region plot, otherwise not
useSAR = true;
use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0

country = 'GB'; % Great Britain
% country = 'SL'; % Sierra Leone
% country = 'SA'; % South Africa

R0vals = [ 1.1 1.3 1.5 1.7 2 2.3 2.7 3.2 4 ]; R0name = 'R0all';
% R0vals = [ 1.5 2 4 ]; R0name = 'R0basic';
% R0vals = [ 1.1 1.3 1.5 ]; R0name = 'R0down';
% R0vals = [ 1.7 2 2.3 ]; R0name = 'R0middle';
% R0vals = [ 2.7 3.2 4 ]; R0name = 'R0up';
popfig = '2ran'; thetaGval = NaN; gamval = 1; % Random
% popfig = 'm4r'; thetaGval = 0.4; gamval = 1; 
% popfig = 'm4UK'; thetaGval = 0.4; gamval = 0.75; 
% popfig = 'm5r'; thetaGval = 0.5; gamval = 1; 
% popfig = 'm5UK'; thetaGval = 0.5; gamval = 0.75; 
% popfig = 'UK'; thetaGval = 0.58; gamval = 0.75; % UK
% popfig = 'ass'; thetaGval = 0.7; gamval = 0.75; % UK
phivals = [ 1 1.5 2 ];
% endphi = phivals; % I explore phi = 1 to 1, 1.5 or 2
psirange = [ 1 4 ]; % Don't modify this
endpsi = [ 2 3 4 ]; % I explore psi = 1-2, 1-3 or 1-4
tolval = 0.05;
figletter = 'A';

toprow = R0vals(3:end);
leftcol = R0vals(1:(end-2));

rvals = [ 0.14552 0.25282 0.52588 ]; % These are the values of r corresponding to the values of R0 = 1.5, 2 and 4
% nsimval = 100;

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
    if Activate_plot_from_new_ROT_analysis
        fig_ROT_path = [fig_base_dir,'\output-rule-of-thumb\'];
    else
        fig_ROT_path = [fig_base_dir,'\saved-rule-of-thumb\'];
    end
else
    fig_code_path = [fig_base_dir,'/code-figures/'];
    fig_fig_path = [fig_base_dir,'/output-figures/supp/',country,'/',popfig,'/'];
    fig_tool_path = [fig_base_dir,'/tools/'];
    if Activate_plot_from_new_ROT_analysis
        fig_ROT_path = [fig_base_dir,'/output-rule-of-thumb/'];
    else
        fig_ROT_path = [fig_base_dir,'/saved-rule-of-thumb/'];
    end
end
cd(fig_code_path); % Work in the directory where the code for figures are
addpath(genpath(fig_tool_path));
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
% Loaded workspaces generated from other computers, with a different folder
% structure might have saved anonymous functions with paths not recognised 
% when this script is run, so the line above turns off that warning.

% Workspace name
wname = [ 'ROT_tol',num2str(round(tolval*100),'%02d'),'_'];
if usez
    wname = [wname,'z'];
end
if usepi
    wname = [wname,'pi'];
end
if uset
    wname = [wname,'t'];
end
wname = [wname,'_',country,'_',popfig,'_',R0name];
wname = [wname,'_psirange',num2str(round(psirange(1)*10),'%02d'),'_',num2str(round(psirange(2)*10),'%02d')];
if Activate_clean_output
    wname = [wname,'_clean'];
else
    wname = [wname,'_dirty'];
end
load([fig_ROT_path,wname]);

% Build figure name
fname = ['3x3_ROT_table_'];
if useSAR
    fname = [fname,'overSAR'];
else
    fname = [fname,'overTh'];
end
fname = [ fname,'_',country,'_tol',num2str(round(tolval*100),'%02d'),'_'];
if usez
    fname = [fname,'z'];
end
if usepi
    fname = [fname,'pi'];
end
if uset
    fname = [fname,'t'];
end
fname = [fname,'_',popfig,'_',R0name];
if use_match_r
    fname = [fname,'_match_r'];
end

lx = length(toprow);
ly = length(leftcol);
lendpsi = length(endpsi);
lendphi = length(phivals);

progressbar(0,0,0,0)
minMMcolor = Inf;
maxMMcolor = -Inf;
for i = 1:lendphi
    for j = 1:lendpsi
        % dotsROT is the main variable storing the output. Each columns
        % contains [ R0, phi, pAA, psi, SAR ]^T, where pAA and psi come 
        % from the countour function, SAR is calculated in each of these 
        % values, and then I loop through each value of phi (for fixed R0),
        % and then through each R0.

        % Select only specific values of phi and ranges of values of psi
        cM = ( dotsROT(2,:) == phivals(i) ) & ( endpsi(j)-1 <= dotsROT(4,:) ) & ( dotsROT(4,:) <= endpsi(j) );
        M{i,j} = dotsROT(:,cM); % Pick only the columns of dotsROT I'm interested in
        
        for iy = 1:ly
            for ix = 1:lx
                cMM = ( leftcol(iy) <= M{i,j}(1,:) ) & ( M{i,j}(1,:) <= toprow(ix) );
                MM{i,j} = M{i,j}(:,cMM); % Pick only the columns corresponding to values of R0 inside the interval of interest
                if ix < iy
                    MMvals1{i,j}(iy,ix) = 0;
                    MMvals2{i,j}(iy,ix) = 0;
                    MMtext{i,j}{iy,ix} = '';
                    MMcolor{i,j}(iy,ix) = NaN;
                else
                    aa = MM{i,j}(1,:)'; % R0
                    bb = MM{i,j}(5,:)'; % SAR (in fraction of population)
                    mdl = fitlm(aa,bb);
                    valuem = mdl.Coefficients.Estimate(1);
                    valueq = mdl.Coefficients.Estimate(2);
                    MMvals1{i,j}(iy,ix) = valuem;
                    MMvals2{i,j}(iy,ix) = valueq;
                    MMtext{i,j}{iy,ix} = [ ' ',num2str(valueq,'%1.3f');num2str(valuem,'%+1.3f') ];
                    R2color = mdl.Rsquared.Adjusted;
                    MMcolor{i,j}(iy,ix) = R2color;
                    if R2color > maxMMcolor
                        maxMMcolor = R2color;
                    end
                    if R2color < minMMcolor
                        minMMcolor = R2color;
                    end                    
                end
                frac4 = ix/lx;
                frac3 = ((iy-1) + frac4) / ly;
                frac2 = ((j-1) + frac3) / lendpsi;
                frac1 = ((i-1) + frac2) / lendphi;
                
                progressbar(frac1,frac2,frac3,frac4);
            end
        end
        [ maxcolval, maxrowind ] = max( MMcolor{i,j} );
        [ maxval, maxcolind ] = max( maxcolval );
        MMmaxind{i,j} = [ maxrowind(maxcolind), maxcolind ];
        MMmaxval{i,j} = [ MMvals1{i,j}(maxrowind(maxcolind),maxcolind), MMvals2{i,j}(maxrowind(maxcolind),maxcolind) ];
        MMmaxR2{i,j} = maxval;
        [ mincolval, mincolind ] = min( MMcolor{i,j} );
        [ minval, minrowind ] = min( mincolval );
        MMminind{i,j} = [ minrowind, mincolind(minrowind) ];
        MMminval{i,j} = [ MMvals1{i,j}(minrowind,mincolind(minrowind)), MMvals2{i,j}(minrowind,mincolind(minrowind)) ];
        MMminval{i,j} = minval;
    end
end

%% Figure
     
% %%%%%% Data
D.X = toprow;
D.Y = leftcol;
% D.Z = { MMvals1, MMvals2, MMcolor };
D.Z = { MMtext, MMcolor };
D.Clim = [ floor(minMMcolor*10)/10, 1];
D.Climleg = D.Clim;
D.colormap = 'hot';

% %%%%%% Text
T.figletter = figletter;
T.title = ['\psi = 1 - 2                      \psi = 2 - 3                      \psi = 3 - 4'];
T.labely = ['    \phi = 2                      \phi = 1.5                      \phi = 1'];
T.labelx = ['Adjusted R^2 (corresponding to linear regression coefficient in each cell)'];

% T.title_xcorr = 0.035;
% T.labelx_xcorr = 0.05;
T.labely_ycorr = -0.03;

% %%%%%% Layout
L.fig_width_cm = 7; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 3;
L.ncols = 3;
% Title
T.title_font_size = 6;
L.title_height_cm = 0.5;
L.xlabel_height_cm = 0.4;
L.xlabel_ypos = -1;
L.ylabel_width_cm = 0.5;
% Legend
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0.7;
L.colorbar = true;
L.colorbar_bottom_limit_in_box = 1;
% T.cbar_string = 'Adjusted R^2';
% T.cbar_string_font_size = 6;

% Table
T.table_incells_font_size = 2.1;
T.table_outcells_font_size = 2.5;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'ROT_table', D, L, T );

if Activate_save_fig
    cd(fig_fig_path);
    expcmd = ['export_fig ',fname,' -pdf -nocrop -transparent'];
    eval(expcmd);
    cd(fig_code_path);
end
rmpath(genpath(fig_tool_path));

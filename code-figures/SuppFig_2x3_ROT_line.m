% This is the code to create Figure 25 in the supplementary text of 
% Pellis, L et al (2019), Nature Communications
% 
% It generates a 2x3 graph with, in each panel, the linear regression
% to work out the rule-of-thumb. Each column for a different country and
% each row for random or assortative mixing. There is no title, but
% subtitles are generated for each row and column.
% 
% It relies on a hand-made function called "mycomplexfig"
% It also relies on pre-computed workspaces where information about the 
% rule of thumb is stored. Such workspaces are contained in the subfolder
% called "rule-of-thumb". If they are not available, one needs to first run
% "Analyse_data_for_RuleOfThumb.m"
% 
% Update: 14/10/2019

close all; % close all figures
clearvars;
Activate_clean_output = 1;
Activate_save_fig = 0; % If true, figures are saved
Activate_plot_from_new_ROT_analysis = 0; 
% If 0, I make plots from pre-computed and saved workspaces for the rule-of-thumb analysis (folder saved-rule-of-thumb)
% If 1, I make plots from newly computed workspaces for the rule-of-thumb analysis (folder output-rule-of-thumb)

usez = true; % If true, I include the final size (z) in the overall acceptance region plot, otherwise not
usepi = true; % If true, I include the peak incidence (pi) in the overall acceptance region plot, otherwise not
uset = true;  % If true, I include the time to the peak (t) in the overall acceptance region plot, otherwise not
useSAR = true;
use_match_r = false; % If false, I match R0; if true, I use the r correspondent to the desired R0

% country = 'GB';
% % country = 'SA'; % South-Africa
% % country = 'SL'; % Sierra-Leone

R0vals = [ 1.1 1.3 1.5 1.7 2 2.3 2.7 3.2 4 ]; R0name = 'R0all';
% R0vals = [ 1.5 2 4 ]; R0name = 'R0basic';
% R0vals = [ 1.1 1.3 1.5 ]; R0name = 'R0down';
% R0vals = [ 1.7 2 2.3 ]; R0name = 'R0middle';
% R0vals = [ 2.7 3.2 4 ]; R0name = 'R0up';
% R0vals = [ 1.3 1.5 1.7 ]; R0name = 'R0best3';
% R0vals = [ 1.3 1.5 1.7 2 ]; R0name = 'R0best4';
% R0vals = [ 1.3 1.5 1.7 2 2.3 ]; R0name = 'R0best5';
% R0vals = [ 1.3 1.5 1.7 2 2.3 2.7 ]; R0name = 'R0best6';
% R0vals = [ 1.1 1.3 1.5 1.7 2 2.3 2.7 ]; R0name = 'R0best7';
% R0vals = [ 1.1 1.3 1.5 1.7 2 2.3 2.7 3.2 ]; R0name = 'R0best8';
% % pop = '2ran'; thetaGval = NaN; gamval = 1; % Random
% % pop = 'm4r'; thetaGval = 0.4; gamval = 1; 
% % pop = 'm4UK'; thetaGval = 0.4; gamval = 0.75; 
% % pop = 'm5r'; thetaGval = 0.5; gamval = 1; 
% % pop = 'm5UK'; thetaGval = 0.5; gamval = 0.75; 
% pop = 'UK'; thetaGval = 0.58; gamval = 0.75; % UK
% % pop = 'ass'; thetaGval = 0.7; gamval = 0.75; % UK
psirange = [ 1 4 ]; % Don't change this
maxprod = 3; % I consider more reliable only phi*psi<maxprod
% tolval = 0.01; minmaxR0 = [ 1.3 4 ]; figletter = 'A';
tolval = 0.05; minmaxR0 = [ 1.3 3.2 ]; figletter = 'B';
% tolval = 0.1; minmaxR0 = [ 1.3 2.5 ]; figletter = 'C';
rvals = [ 0.14552 0.25282 0.52588 ]; % These are the values of r corresponding to the values of R0 = 1.5, 2 and 4
% nsimval = 100;

clist = {'GB','SA','SL'};
plist = {'2ran','UK'};
tlist = [ NaN, 0.58 ];
glist = [ 1, 0.75 ];

% Path stuff
current_dir = cd;
eval('cd ..'); % Move to the folder 1 level up, which is assumed to be the "base" folder
base_dir = cd; % This is assumed to be the self-contained folder with all relevant files and subfolders
if ispc
    code_path = [base_dir,'\code-figures\'];
    fig_path = [base_dir,'\output-figures\supp\'];
    tool_path = [base_dir,'\tools\'];
    if Activate_plot_from_new_ROT_analysis
        ROT_path = [base_dir,'\output-rule-of-thumb\'];
    else
        ROT_path = [base_dir,'\saved-rule-of-thumb\'];
    end
else
    code_path = [base_dir,'/code-figures/'];
    fig_path = [base_dir,'/output-figures/supp/'];
    tool_path = [base_dir,'/tools/'];
    if Activate_plot_from_new_ROT_analysis
        ROT_path = [base_dir,'/output-rule-of-thumb/'];
    else
        ROT_path = [base_dir,'/saved-rule-of-thumb/'];
    end
end
cd(code_path); % Work in the directory where the code for figures are

for s = 1:3
    country = clist{s};
    for p = 1:2
        pop = plist{p};
        thetaGval = tlist(p);
        gamval = glist(p);

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
                myROTcoeff = 0.24;
                myROTinter = -0.1;
%                 myROTcoeff = 0.235;%0.236; % Current best for metaROT; % best for fitting (between assortative and random)
%                 myROTinter = -0.09;%-0.094;
            elseif strcmp(country,'SA')
                myROTcoeff = 0.225;
                myROTinter = -0.16;%-0.166;
%                 myROTcoeff = 0.226;%0.228;
%                 myROTinter = -0.165;
            elseif strcmp(country,'SL')
                myROTcoeff = 0.22;
                myROTinter = -0.18;%-0.188;
%                 myROTcoeff = 0.223;%0.222;
%                 myROTinter = -0.19;%-0.188;
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

        toprow = R0vals(3:end);
        leftcol = R0vals(1:(end-2));

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
        wname = [wname,'_',country,'_',pop,'_',R0name];
        wname = [wname,'_psirange',num2str(round(psirange(1)*10),'%02d'),'_',num2str(round(psirange(2)*10),'%02d')];
        if Activate_clean_output
            wname = [wname,'_clean'];
        else
            wname = [wname,'_dirty'];
        end
        load([ROT_path,wname]);

        whichnearall = ( dotsROT(2,:) .* dotsROT(4,:) <= maxprod );
        whichfarall = ~whichnearall;
        whichnearfit = ( minmaxR0(1) <= dotsROT(1,:) ) & ( dotsROT(1,:) <= minmaxR0(2) ) & whichnearall;
        whichfarfit = ( minmaxR0(1) <= dotsROT(1,:) ) & ( dotsROT(1,:) <= minmaxR0(2) ) & whichfarall;

        a = dotsROT(1,whichfarfit)';
        b = dotsROT(nrowsROT,whichfarfit)';
        % Fitting best line
        mdl = fitlm(a,b);
        valuesg = mdl.Coefficients.Estimate;
        % plot([0,4],[valuesg(1),valuesg(1)+valuesg(2)*4],':','color',0.6*ones(1,3),'Linewidth',1);

        a = dotsROT(1,whichnearfit)';
        b = dotsROT(nrowsROT,whichnearfit)';
        % Fitting best line
        mdl = fitlm(a,b);
        valuesk = mdl.Coefficients.Estimate;
        % plot([0,4],[valuesk(1),valuesk(1)+valuesk(2)*4],'k:','Linewidth',1);

        % plot([0,4],[myROT(1),myROT(1)+myROT(2)*4],'k','Linewidth',2);

        datax{p,s} = dotsROT(1,:)';
        datay{p,s} = 100*dotsROT(nrowsROT,:)';
        dataz{p,s} = whichnearall;
        datag{p,s} = 100*valuesg;
        datak{p,s} = 100*valuesk;
        datam{p,s} = 100*myROT;
        strcountry = {'UK','South Africa','Sierra Leone'};
        strmix = {'Random','Assortative'};

    end
end

% Build figure name
fname = ['2x3_ROT_line_'];
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
fname = [fname,'_',R0name];
if use_match_r
    fname = [fname,'_match_r'];
end

%% Figure

% %%%%%% Data
D.X = datax;
D.Y = datay;
D.Z = dataz;
D.nZover = 3;
D.Zover = { datag, datak, datam };
D.limx = [0,4.1];
D.limy = [0,100];

% %%%%%% Text
T.labelx = 'R_0';
T.labely = ['SAR (%) at ',num2str(round(tolval*100)),'% relative error threshold'];

% T.title_xcorr = 0.035;
% T.labelx_xcorr = 0.05;
% T.labely_ycorr = 0.01;

T.coltitles = {['Great Britain',char(10),'\beta_1 = ',num2str(datam{1,1}(2),'%2.1f'),', \beta_0 = ',num2str(datam{1,1}(1),'%2.1f'),sprintf('\\fontsize{25} ')],...
    ['South Africa',char(10),'\beta_1 = ',num2str(datam{1,2}(2),'%2.1f'),', \beta_0 = ',num2str(datam{1,2}(1),'%2.1f'),sprintf('\\fontsize{25} ')],...
    ['Sierra Leone',char(10),'\beta_1 = ',num2str(datam{1,3}(2),'%2.1f'),', \beta_0 = ',num2str(datam{1,3}(1),'%2.1f'),sprintf('\\fontsize{25} ')]};
T.rowtitles = {'Random','Assortative'};
T.figletter = figletter;

% %%%%%% Layout
L.fig_width_cm = 12; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 2;
L.ncols = 3;
% Title
L.title_height_cm = 0;
L.xlabel_height_cm = 1;
L.ylabel_width_cm = 1;
% Legend
if tolval == 0.01
    L.legend = true;
    T.legentries = {'Data (unlikely)','Data (more likely)','Regression (unlikely)','Regression (more likely)','Empirical rule of thumb'};
else
    L.legend = false;
end
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0;
L.colorbar = false;
% Subtitles
L.subtitle_height_cm = 0;
T.title_font_size = 8;
T.legend_font_size = 4;


L.coltitle_height_cm = 0.8;
L.coltitle_ypos = 1/2;
L.rowtitle_width_cm = 0.4;

T.axes_font_size = 8;
T.subtitle_font_size = 8;
T.subaxes_font_size = 6;
T.legend_font_size = 6;

L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'ROT_line', D, L, T );

if Activate_save_fig
    addpath(genpath(tool_path));
    cd(fig_path);
    expcmd = ['export_fig ',fname,' -pdf -nocrop -transparent'];
    eval(expcmd);
    cd(code_path);
    rmpath(genpath(tool_path));
end
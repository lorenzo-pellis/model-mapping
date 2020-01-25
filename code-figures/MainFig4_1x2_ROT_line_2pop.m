% This is the code to create Figure 4 of the main text of 
% Pellis, L et al (2020), Nature Communications
% 
% It relies on a hand-made function called "mycomplexfig"
% It also relies on pre-computed workspaces where information about the 
% rule of thumb is stored. Such workspaces are contained in the subfolder
% called "rule-of-thumb". If they are not available, one needs to first run
% "Analyse_data_for_RuleOfThumb.m"
% 
% Update: 25-01-2020

close all; % close all figures
clearvars;
Activate_clean_output = 1; % Clean output (leave as it is: details explained in Supplementary Discussion, Section 2.4)
Activate_save_fig = 1; % If true, figures are saved
Activate_plot_from_new_ROT_analysis = 0; 
% If 0, I make plots from pre-computed and saved workspaces for the rule-of-thumb analysis (folder saved-rule-of-thumb)
% If 1, I make plots from newly computed workspaces for the rule-of-thumb analysis (folder output-rule-of-thumb)

usez = true; % If true, I include the final size (z) in the overall acceptance region plot, otherwise not
usepi = true; % If true, I include the peak incidence (pi) in the overall acceptance region plot, otherwise not
uset = true;  % If true, I include the time to the peak (t) in the overall acceptance region plot, otherwise not

clist = {'GB','SL'}; % List of countries appearing in the plot
plist = {'UK'}; % List of population mixing patterns
R0vals = [ 1.1 1.3 1.5 1.7 2 2.3 2.7 3.2 4 ]; R0name = 'R0all';
psirange = [ 1 4 ]; % Range of values of psi considered
minmaxR0 = [ 1.3 3.2 ]; % This gives the smallest and largest values of R0 
% included in the linear regression (i.e. esclude extremes R0 = 1.1 and 4)
tolval = 0.05;

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
    if Activate_plot_from_new_ROT_analysis
        fig_ROT_path = [fig_base_dir,'\output-rule-of-thumb\'];
    else
        fig_ROT_path = [fig_base_dir,'\saved-rule-of-thumb\'];
    end
else
    fig_code_path = [fig_base_dir,'/code-figures/'];
    fig_fig_path = [fig_base_dir,'/output-figures/main/'];
    fig_tool_path = [fig_base_dir,'/tools/'];
    if Activate_plot_from_new_ROT_analysis
        fig_ROT_path = [fig_base_dir,'/output-rule-of-thumb/'];
    else
        fig_ROT_path = [fig_base_dir,'/saved-rule-of-thumb/'];
    end
end
cd(fig_ROT_path); % Work in the directory where the analysis for the ROT appears
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
% Loaded workspaces generated from other computers, with a different folder
% structure might have saved anonymous functions with paths not recognised 
% when this script is run, so the line above turns off that warning.

for s = 1:length(clist)
    country = clist{s};
    for p = 1:length(plist);
        pop = plist{p};

        % Values of slope (myROTcoeff) and intercept (myROTinter) are found
        % by trial and error, to eyeball a decent fit which works well for
        % both random and assortative mixing.
        if tolval == 0.01
%             if strcmp(country,'GB')
%                 myROTcoeff = 0.092;
%                 myROTinter = -0.048;
%             elseif strcmp(country,'SA')
%                 myROTcoeff = 0.074;%0.072;
%                 myROTinter = -0.051;%-0.048;
%             elseif strcmp(country,'SL')
%                 myROTcoeff = 0.068;%0.07;
%                 myROTinter = -0.052;%-0.062;
%             else
%                 myROTcoeff = 0;
%                 myROTinter = 0;
%             end
        elseif tolval == 0.05
            if strcmp(country,'GB')
                myROTcoeff = 0.24;
                myROTinter = -0.1;
%             elseif strcmp(country,'SA')
%                 myROTcoeff = 0.225;
%                 myROTinter = -0.16;
            elseif strcmp(country,'SL')
                myROTcoeff = 0.22;
                myROTinter = -0.18;
            else
                myROTcoeff = 0;
                myROTinter = 0;
            end
        elseif tolval == 0.1
%             if strcmp(country,'GB')
%                 myROTcoeff = 0.368;
%                 myROTinter = -0.14;
%             elseif strcmp(country,'SA')
%                 myROTcoeff = 0.347;%0.348;
%                 myROTinter = -0.23;%-0.22;
%             elseif strcmp(country,'SL')
%                 myROTcoeff = 0.34;
%                 myROTinter = -0.26;
%             else
%                 myROTcoeff = 0;
%                 myROTinter = 0;
%             end
        else
            myROTcoeff = 0;
            myROTinter = 0;
        end
        myROT = [ myROTinter, myROTcoeff];

        % Workspace name
        wrks_name = [ 'ROT_tol',num2str(round(tolval*100),'%02d'),'_'];
        if usez
            wrks_name = [wrks_name,'z'];
        end
        if usepi
            wrks_name = [wrks_name,'pi'];
        end
        if uset
            wrks_name = [wrks_name,'t'];
        end
        wrks_name = [wrks_name,'_',country,'_',pop,'_',R0name];
        wrks_name = [wrks_name,'_psirange',num2str(round(psirange(1)*10),'%02d'),'_',num2str(round(psirange(2)*10),'%02d')];
        if Activate_clean_output
            wrks_name = [wrks_name,'_clean'];
        else
            wrks_name = [wrks_name,'_dirty'];
        end
        load([fig_ROT_path,wrks_name]);

        % "near" refers to points from values of phi and psi that are more
        % 'reasonable', where by reasonable I mean: phi*psi <= 3
        whichnearall = ( dotsROT(2,:) .* dotsROT(4,:) <= 3 ); % indices of 'reasonable' values
        whichfarall = ~whichnearall; % other ('unreasonable') values
        whichnearfit = ( minmaxR0(1) <= dotsROT(1,:) ) & ( dotsROT(1,:) <= minmaxR0(2) ) & whichnearall;
        whichfarfit = ( minmaxR0(1) <= dotsROT(1,:) ) & ( dotsROT(1,:) <= minmaxR0(2) ) & whichfarall;

        % Fit to 'unreasonable' points ("near")
        a = dotsROT(1,whichfarfit)';
        b = dotsROT(nrowsROT,whichfarfit)';
        % Fitting best line
        mdl = fitlm(a,b);
        valuesg = mdl.Coefficients.Estimate;

        % Fir to 'reasonable' points ("far")
        a = dotsROT(1,whichnearfit)';
        b = dotsROT(nrowsROT,whichnearfit)';
        % Fitting best line
        mdl = fitlm(a,b);
        valuesk = mdl.Coefficients.Estimate;

        datax{p,s} = dotsROT(1,:)';
        datay{p,s} = 100*dotsROT(nrowsROT,:)';
        dataz{p,s} = whichnearall;
        datag{p,s} = 100*valuesg;
        datak{p,s} = 100*valuesk;
        datam{p,s} = 100*myROT;
    end
end


%% Figure
cd(fig_code_path)
% To constuct the figure, I pass 3 structures:
%   - D = Data
%   - T = Text (title(s), axe labels, legend text, etc.)
%   - L = Layout information for the figure spatial structure

% %%%%%% Data
D.X = datax;
D.Y = datay;
D.Z = dataz;
D.nZover = 3;
D.Zover = { datag, datak, datam };
D.limx = [0,4.1];
D.limy = [0,100];

% %%%%%% Text
T.title = ['A)  Great Britain: \beta_1 = ',num2str(round(datam{1,1}(2)),'%d'),', \beta_0 = ',num2str(round(datam{1,1}(1)),'%d'),...
    '   B)  Sierra Leone: \beta_1 = ',num2str(round(datam{1,2}(2)),'%d'),', \beta_0 = ',num2str(round(datam{1,2}(1)),'%d'),'      '];
T.title = ['A)  Great Britain: \eta_1 = ',num2str(round(datam{1,1}(2)),'%d'),', \eta_0 = ',num2str(round(datam{1,1}(1)),'%d'),...
    '   B)  Sierra Leone: \eta_1 = ',num2str(round(datam{1,2}(2)),'%d'),', \eta_0 = ',num2str(round(datam{1,2}(1)),'%d'),'      '];
T.labelx = 'R_0';
T.labely = ['SAR (%) at ',num2str(round(tolval*100)),'% relative error threshold'];

% %%%%%% Layout
L.fig_width_cm = 9; % Text width of A4 portrait ~ 15cm
L.fig_height_cm = 7; % Text width of A4 portrait ~ 25cm
L.screen_scale = 3;
% Panels
L.nrows = 1;
L.ncols = 2;
% Title
L.title_height_cm = 0.5;
L.xlabel_height_cm = 1;
L.ylabel_width_cm = 0.8;
% Legend
L.legend = true;
T.legentries = {'Data (unlikely)','Data (more likely)','Regression (unlikely)','Regression (more likely)','Empirical rule of thumb'};
L.vrlegend_width_cm = 0.2;
L.hblegend_height_cm = 0;
L.colorbar = false;
% Subtitles
L.subtitle_height_cm = 0;


T.title_font_size = 7.5;
T.axes_font_size = 6;
T.subaxes_font_size = 6;


L = set_mycomplexfig_layout( L );
T = set_mycomplexfig_text( T, L );
D = set_mycomplexfig_data( D );

subplot_handles = mycomplexfig( 'ROT_line', D, L, T );

if Activate_save_fig
    addpath(genpath(fig_tool_path));
    cd(fig_fig_path);
    export_fig Main_Figure4 -pdf -nocrop -transparent;
    cd(fig_code_path);
    rmpath(genpath(fig_tool_path));
end


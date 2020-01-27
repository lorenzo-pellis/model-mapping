function L = set_mycomplexfig_layout( L )
% 
% This is a hard-made functions written in support of "mycomplexfig", a
% function I use to make plots that have multiple panels showing variants
% of the same output. I mostly use it to create the figures used in the
% model mapping paper:
% 
% Reference: Pellis, L. et al (2020), Nature Communications
% 
% Update: 22-06-2019

%%%%%% Layout
% Dimensions
if ~isfield(L,'fig_width_cm')
    L.fig_width_cm = 15; % Text width of A4 portrait ~ 15cm
end
if ~isfield(L,'fig_height_cm')
    L.fig_height_cm = 15; % Text width of A4 portrait ~ 25cm
end
if ~isfield(L,'screen_scale')
    L.screen_scale = 2;
end
L.fig_width = L.screen_scale * L.fig_width_cm;
L.fig_height = L.screen_scale * L.fig_height_cm;

% Panels
if ~isfield(L,'fig_width_cm')
    L.nrows = 1;
end
if ~isfield(L,'fig_width_cm')
    L.ncols = 1;
end

% Title
if ~isfield(L,'title_height_cm')
    L.title_height_cm = 0;
end
if ~isfield(L,'title_ypos')
    L.title_ypos = 1/2; % Title y position in the rectangle devoted to it
end
L.title_height = L.screen_scale * L.title_height_cm / L.fig_height;

% Labels
if ~isfield(L,'xlabel_height_cm')
    L.xlabel_height_cm = 0;
end
if ~isfield(L,'ylabel_width_cm')
    L.ylabel_width_cm = 0;
end
if ~isfield(L,'xlabel_ypos')
    L.xlabel_ypos = 0.3;
end
if ~isfield(L,'ylabel_xpos')
    L.ylabel_xpos = 3/4;
end
L.xlabel_height = L.screen_scale * L.xlabel_height_cm / L.fig_height;
L.ylabel_width = L.screen_scale * L.ylabel_width_cm / L.fig_width;

if ~isfield(L,'y2axis')
    L.y2axis = false;
end
if ~isfield(L,'y2label_width_cm')
    L.y2label_width_cm = 0;
end
if ~isfield(L,'y2label_xpos')
    L.y2label_xpos = 3/8;
end
L.y2label_width = L.screen_scale * L.y2label_width_cm / L.fig_width;

% Legend
if ~isfield(L,'vllegend_width_cm')
    L.vllegend_width_cm = 0;
end
if ~isfield(L,'vrlegend_width_cm')
    L.vrlegend_width_cm = 0;
end
if ~isfield(L,'htlegend_height_cm')
    L.htlegend_height_cm = 0;
end
if ~isfield(L,'hblegend_height_cm')
    L.hblegend_height_cm = 0;
end
L.vllegend_width = L.screen_scale * L.vllegend_width_cm / L.fig_width;
L.vrlegend_width = L.screen_scale * L.vrlegend_width_cm / L.fig_width;
L.htlegend_height = L.screen_scale * L.htlegend_height_cm / L.fig_height;
L.hblegend_height = L.screen_scale * L.hblegend_height_cm / L.fig_height;

if ~isfield(L,'colorbar')
    L.colorbar = false;
end
if ~isfield(L,'legend')
    L.legend = false;
end
if ~isfield(L,'colorbar_text_out')
    L.colorbar_text_out = true;
end
if ~isfield(L,'colorbar_fraction_of_box')
    L.colorbar_fraction_of_box = 1/3;
end
if ~isfield(L,'colorbar_bottom_limit_in_box')
    L.colorbar_bottom_limit_in_box = 1 - L.colorbar_fraction_of_box;
end


% % Margins (to accommodate the axes values only left and bottom
% if ~isfield(L,'lmargin_cm')
%     L.lmargin_cm = 0.4;
% end
% if ~isfield(L,'bmargin_cm')
%     L.bmargin_cm = 0.2;
% end
% L.lmargin = L.screen_scale * L.lmargin_cm / L.fig_width;
% L.bmargin = L.screen_scale * L.bmargin_cm / L.fig_width;

% Subtitles
if ~isfield(L,'subtitle_height_cm')
    L.subtitle_height_cm = 0;
end
if ~isfield(L,'subtitle_ypos')
    L.subtitle_ypos = 1/2; 
end
L.subtitle_height = L.screen_scale * L.subtitle_height_cm / L.fig_height;

% Sublabels
if ~isfield(L,'xsublabel_height_cm')
    L.xsublabel_height_cm = 0;
end
if ~isfield(L,'ysublabel_width_cm')
    L.ysublabel_width_cm = 0;
end
if ~isfield(L,'xsublabel_ypos')
    L.xsublabel_ypos = 1/2;
end
if ~isfield(L,'ysublabel_xpos')
    L.ysublabel_xpos = 1/2;
end
L.xsublabel_height = L.screen_scale * L.xsublabel_height_cm / L.fig_height;
L.ysublabel_width = L.screen_scale * L.ysublabel_width_cm / L.fig_width;

% Sublegend
if ~isfield(L,'vlsublegend_width_cm')
    L.vlsublegend_width_cm = 0;
end
if ~isfield(L,'vrsublegend_width_cm')
    L.vrsublegend_width_cm = 0;
end
if ~isfield(L,'htsublegend_height_cm')
    L.htsublegend_height_cm = 0;
end
if ~isfield(L,'hbsublegend_height_cm')
    L.hbsublegend_height_cm = 0;
end
L.vlsublegend_width = L.screen_scale * L.vlsublegend_width_cm / L.fig_width;
L.vrsublegend_width = L.screen_scale * L.vrsublegend_width_cm / L.fig_width;
L.htsublegend_height = L.screen_scale * L.htsublegend_height_cm / L.fig_height;
L.hbsublegend_height = L.screen_scale * L.hbsublegend_height_cm / L.fig_height;

if ~isfield(L,'subcolorbar')
    L.subcolorbar = false;
end
if ~isfield(L,'subcolorbar_text_out')
    L.subcolorbar_text_out = true;
end
if ~isfield(L,'subcolorbar_fraction_of_box')
    L.subcolorbar_fraction_of_box = 1/2;
end
if ~isfield(L,'subcolorbar_bottom_limit_in_box')
    L.subcolorbar_bottom_limit_in_box = 1/5;
end

% Column title
if ~isfield(L,'coltitle')
    L.coltitle = [];
end
if ~isfield(L,'coltitle_height_cm')
    L.coltitle_height_cm = 0;
end
if ~isfield(L,'coltitle_ypos')
    L.coltitle_ypos = 0.4; 
end
L.coltitle_height = L.screen_scale * L.coltitle_height_cm / L.fig_height;

% Row title
if ~isfield(L,'rowtitle')
    L.rowtitle = [];
end
if ~isfield(L,'rowtitle_width_cm')
    L.rowtitle_width_cm = 0;
end
if ~isfield(L,'rowtitle_xpos')
    L.rowtitle_xpos = 2/3; 
end
L.rowtitle_width = L.screen_scale * L.rowtitle_width_cm / L.fig_width;


% Subplot separator
if ~isfield(L,'vseparator_cm')
    L.vseparator_cm = 0.25;
end
if ~isfield(L,'hseparator_cm')
    L.hseparator_cm = 0.25;
end
L.vseparator = L.screen_scale * L.vseparator_cm / L.fig_height;
L.hseparator = L.screen_scale * L.hseparator_cm / L.fig_height;

% Subplot margins
if ~isfield(L,'lsubmargin_cm')
    L.lsubmargin_cm = 0;
end
if ~isfield(L,'bsubmargin_cm')
    L.bsubmargin_cm = 0;
end
if ~isfield(L,'rsubmargin_cm')
    L.rsubmargin_cm = 0;
end
if ~isfield(L,'tsubmargin_cm')
    L.tsubmargin_cm = 0;
end
L.lsubmargin = L.screen_scale * L.lsubmargin_cm / L.fig_width;
L.bsubmargin = L.screen_scale * L.bsubmargin_cm / L.fig_width;
L.rsubmargin = L.screen_scale * L.rsubmargin_cm / L.fig_height;
L.tsubmargin = L.screen_scale * L.tsubmargin_cm / L.fig_height;


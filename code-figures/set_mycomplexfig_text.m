function T = set_mycomplexfig_text( T, L )
% 
% This is a hard-made functions written in support of "mycomplexfig", a
% function I use to make plots that have multiple panels showing variants
% of the same output. I mostly use it to create the figures used in the
% model mapping paper:
% 
% Reference: Pellis et al (2019), Nature Communications
% 
% Update: 22/06/2019

scalingfactor = L.screen_scale;

%%%%%% Text
if ~isfield(T,'figletter')
    T.figletter = [];
end
if ~isfield(T,'title')
    T.title = [];
end
if ~isfield(T,'labelx')
    T.labelx = [];
end
if ~isfield(T,'labely')
    T.labely = [];
end
if ~isfield(T,'labely2')
    T.labely2 = [];
end
if ~isfield(T,'first_subletter')
    T.first_subletter = [];
end
if ~isfield(T,'subtitles')
    T.subtitles = cell(L.nrows,L.ncols);
end
if ~isfield(T,'coltitles')
    T.coltitles = cell(L.ncols);
end
if ~isfield(T,'rowtitles')
    T.rowtitles = cell(L.nrows);
end
if ~isfield(T,'clabels')
    T.clabels = [];
end
if ~isfield(T,'cbar_string')
    T.cbar_string = [];
end
if ~isfield(T,'legentries')
    T.legentries = [];
end

%%%%%% Position

% Title
if ~isfield(T,'title_xcorr')
    T.title_xcorr = 0;
end
if ~isfield(T,'title_ycorr')
    T.title_ycorr = 0;
end

% Labels
if ~isfield(T,'labelx_xcorr')
    T.labelx_xcorr = 0;
end
if ~isfield(T,'labelx_ycorr')
    T.labelx_ycorr = 0;
end
if ~isfield(T,'labely_xcorr')
    T.labely_xcorr = 0;
end
if ~isfield(T,'labely_ycorr')
    T.labely_ycorr = 0;
end
if ~isfield(T,'labely2_xcorr')
    T.labely2_xcorr = 0;
end
if ~isfield(T,'labely2_ycorr')
    T.labely2_ycorr = 0;
end

% Letter
if ~isfield(T,'subletter_xcorr')
    T.subletter_xcorr = 0;
end
if ~isfield(T,'subletter_ycorr')
    T.subletter_ycorr = 0;
end

% Subletter
if ~isfield(T,'figletter_xcorr')
    T.figletter_xcorr = 0;
end
if ~isfield(T,'figletter_ycorr')
    T.figletter_ycorr = 0;
end

% Subtitles
if ~isfield(T,'subtitles_xcorr')
    T.subtitles_xcorr = 0;
end
if ~isfield(T,'subtitles_ycorr')
    T.subtitles_ycorr = 0;
end

% Coltitles
if ~isfield(T,'coltitles_xcorr')
    T.coltitles_xcorr = 0;
end
if ~isfield(T,'coltitles_ycorr')
    T.coltitles_ycorr = 0;
end

% Rowtitles
if ~isfield(T,'rowtitles_xcorr')
    T.rowtitles_xcorr = 0;
end
if ~isfield(T,'rowtitles_ycorr')
    T.rowtitles_ycorr = 0;
end


%%%%% Text size
if ~isfield(T,'title_font_size')
    T.title_font_size = 7.5;
end
if ~isfield(T,'axes_font_size')
    T.axes_font_size = 6;
end
if ~isfield(T,'subtitle_font_size')
    T.subtitle_font_size = 5;
end
if ~isfield(T,'subaxes_font_size')
    T.subaxes_font_size = 4;
end
if ~isfield(T,'legend_font_size')
    T.legend_font_size = 5;
end
if ~isfield(T,'cbar_string_font_size')
    T.cbar_string_font_size = 4;
end
if ~isfield(T,'table_incells_font_size')
    T.table_incells_font_size = 3;
end
if ~isfield(T,'table_outcells_font_size')
    T.table_outcells_font_size = 3;
end
T.title_font_size_sc = scalingfactor * T.title_font_size;
T.axes_font_size_sc = scalingfactor * T.axes_font_size;
T.subtitle_font_size_sc = scalingfactor * T.subtitle_font_size;
T.subaxes_font_size_sc = scalingfactor * T.subaxes_font_size;
T.legend_font_size_sc = scalingfactor * T.legend_font_size;
T.cbar_string_font_size_sc = scalingfactor * T.cbar_string_font_size;
T.table_incells_font_size_sc = scalingfactor * T.table_incells_font_size;
T.table_outcells_font_size_sc = scalingfactor * T.table_outcells_font_size;




function subplot_handles = mysubplot_outputs( fig_width, fig_height, title, xaxis, yaxis, model_list, XData, YData, ZData_list, CLim_list, clabel_list )
% 
% This function creats multiple subplots in a grid but with only 1 big
% x-axis and 1 big y-axis. There are multiple versions of this function,
% adapted to every single type of figure I'm interested in producing
% 
% Input:
%   - fig_width: width of figure in cm
%   - fig_height: height of figure in cm
%   - title: big title for the entire figure (if there is one)
%   - xaxis: label for the x-axes (one large one at the bottom)
%   - yaxis: label for the y-axes (one large one at the left)
%   - model_list: label for the type of model used
%   - XData: values on the x-axes
%   - YData: values on the y-axes
%   - ZData_list: plotted values of the heatmaps
%   - CLim_list: limits of colour axes
%   - clabel_list: label for the colour axes
% 
% In particular, this code has been created when working on the paper:
% Pellis et al (2019), Nature Communications
% 
% Update: 11/10/2019

set(0, 'Units', 'Centimeters' );

cov = 1; % Coverage = shrinking factor to decide how much of the actual drawing of the main plot covers the area dedicated to it
subcov = 1; % subCoverage = shrinking factor to decide how much of the actual drawing of each subplot covers the subarea dedicated to it
xcorr = 1; % Horizontal correction of the position of the x-axis label
[ n_subplots, ~ ] = size( model_list );
[ ~, ~, check ] = size( ZData_list );
if n_subplots ~= check
    error('Error, number of subtitles different from number of plots');
end

%%%%%%%%%%%%%%% Figure attributes
% The figure is made with a lot of rectangles ("placeholders") where to put
% all figure elements. However, then they are adapted manually, so all is a
% bit hacky...
title_height = 1.5; % put 0 if you don't want a title
title_font_size = 18;
title_ypos = 1/3; % Title y position in the rectangle devoted to it

xlabel_height = 2;
ylabel_width = 1.5;
axes_font_size = 18;
xlabel_ypos = 1/3;
ylabel_xpos = 2/3;

vllegend_width = 0; % space for vertical left legend
vrlegend_width = 1; % space for vertical right legend
htlegend_height = 2; % space for horizontal top legend
hblegend_height = 0; % space for horizontal bottom legend

%%%%%%%%%%%%%%% Subplots attributes
subtitle_height = 0.8;
subtitle_font_size = 18;
subtitle_ypos = 1/4; % Title y position in the rectangle devoted to it
subfig_letters_size = 14; % Letters identifying the subplot

xsublabel_height = 0;
ysublabel_width = 0.5;
subaxes_font_size = 12;
subaxes_values_font_size = 12;
xsublabel_ypos = 0;
ysublabel_xpos = 0;

vlsublegend_width = 0;
vrsublegend_width = 0.3;
htsublegend_height = 0;
hbsublegend_height = 0;


screen = get(0,'ScreenSize');
screen_width = screen(3);
screen_height = screen(4);

xLeft = ( screen_width - fig_width ) - 0.3; yTop = ( screen_height - fig_height ) - 2.1;

switch n_subplots
    case 1
        nrows = 1;
        ncols = 1;
    case 2
        nrows = 1;
        ncols = 2;
    case 3
        nrows = 1;
        ncols = 3;
    case 4
        nrows = 2;
        ncols = 2;
    case 5
        nrows = 2;
        ncols = 3;
    case 6
        nrows = 2;
        ncols = 3;
    case 7
        nrows = 3;
        ncols = 3;
    case 8
        nrows = 3;
        ncols = 3;
    case 9
        nrows = 3;
        ncols = 3;
    otherwise
        warning('Too many subplots: not plotting anything...');
end

figure
clf;
colormap jet
set(gcf,'Units','centimeters')
set(gcf,'Position',[ xLeft yTop fig_width fig_height ])

% outax = outer axes handle, for the axes along the border of the actual figure
outax = axes( 'position', [ 0 0 1 1 ], 'visible', 'off', 'Color', 'none' );
% Area = part of the figure where the suplots are drawn (excluding big
% title, big legends and big axes labels)
xArea = ylabel_width + vllegend_width;
Area_width = cov * ( fig_width - xArea - vrlegend_width );
yArea = xlabel_height + hblegend_height;
Area_height = cov * ( fig_height - yArea - title_height - htlegend_height );

% Put big axes labels, big title and big legend(s)
text( ( vllegend_width + ylabel_width * ( 1 - ylabel_xpos ) ) / fig_width, ( yArea + Area_height / 2 ) / fig_height, yaxis, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', axes_font_size );
text( ( xArea + vlsublegend_width + ysublabel_width - xcorr + Area_width / 2 ) / fig_width, ( hblegend_height + xlabel_height * xlabel_ypos ) / fig_height, xaxis, 'HorizontalAlignment', 'center', 'FontSize', axes_font_size);
% text( ( xArea + vlsublegend_width + ysublabel_width - xcorr + Area_width / 2 ) / fig_width, ( yArea + Area_height + htlegend_height + title_height * title_ypos ) / fig_height, title, 'HorizontalAlignment', 'center', 'FontSize', title_font_size);

% subplot = area dedicated to each subplot (including title, axes labels
% and legends)
% subArea = the are dedicated to each drawn part of the subplot (once you 
% remove the subplot title, axes labels and legends)
nsubplots = nrows * ncols;
subplot_width = Area_width / ncols; %Aw / ncols * fig_width;
subplot_height = Area_height / nrows; %Ah / nrows * fig_height;
spw = subplot_width / fig_width;
sph = subplot_height / fig_height;
subArea_width = subcov * ( subplot_width - ysublabel_width - vlsublegend_width - vrsublegend_width );
subArea_height = subcov * ( subplot_height - subtitle_height - xsublabel_height - hbsublegend_height - htsublegend_height );
sAw = subArea_width / fig_width;
sAh = subArea_height / fig_height;

n = 0;
xsubplot = zeros(1,nsubplots);
ysubplot = zeros(1,nsubplots);
sahandle = zeros(1,nsubplots); % subaxes handle
subplot_handles = zeros(1,nsubplots); % subplot handle
for i = 1:nrows
    for j = 1:ncols
        n = n+1;
        axes( outax );
        xsubplot = xArea + ( j - 1 ) * subplot_width; % / fig_width;
        ysubplot = yArea + Area_height - i * subplot_height; % / fig_width;
        
        xsp = xsubplot / fig_width;
        ysp = ysubplot / fig_height;
        
        xsubArea = xsubplot + ysublabel_width + vlsublegend_width;
        ysubArea = ysubplot + xsublabel_height + hbsublegend_height;
        xsubfigureID = xsubplot;
        xsubtitle = xsubplot + subplot_width / 2;
        ysubtitle = ysubplot + hbsublegend_height + xsublabel_height + subArea_height + htsublegend_height + subtitle_height * subtitle_ypos;
        xsubxlabel = xsubplot + vlsublegend_width + ysublabel_width + subArea_width / 2;
        ysubxlabel = ysubplot + hbsublegend_height + xsublabel_height * xsublabel_ypos;
        xsubylabel = xsubplot + vlsublegend_width + ysublabel_width * ysublabel_xpos;
        ysubylabel = ysubplot + hbsublegend_height + xsublabel_height + subArea_height / 2;
        
        % text( xsubylabel / fig_width, ysubylabel / fig_height, 'y subaxis', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', subaxes_font_size );
        % text( xsubxlabel / fig_width, ysubxlabel / fig_height, 'x subaxis', 'HorizontalAlignment', 'center', 'FontSize', subaxes_font_size);
        if j == 1 % figure letter in line 1 needs to be adjusted
            text( xsubfigureID / fig_width - 0.02, ysubtitle / fig_height + 0.01, [char(n+96),')'], 'HorizontalAlignment', 'left', 'FontSize', subfig_letters_size, 'Fontweight', 'bold');
        else
            text( xsubfigureID / fig_width, ysubtitle / fig_height + 0.01, [char(n+96),')'], 'HorizontalAlignment', 'left', 'FontSize', subfig_letters_size, 'Fontweight', 'bold');
        end
        % text( xsubtitle / fig_width, ysubtitle / fig_height, subtitles_list(n,:), 'HorizontalAlignment', 'center', 'FontSize', subtitle_font_size);
              
        xsA = xsubArea / fig_width;
        ysA = ysubArea / fig_height;
        sahandle(n) = axes( 'position', [ xsp, ysp, spw, sph ], 'visible', 'off' );
        subplot_handles(n) = axes( 'position', [ xsA, ysA, sAw, sAh ], 'FontSize', subaxes_values_font_size );
        [cont,cont] = contourf(XData,YData,ZData_list(:,:,n),200);
        set(cont,'LineStyle','none');
        caxis(CLim_list(n,:));
                
        if i == 1
            c = colorbar('location','North','position', [ xsA, ysp + sph, sAw, subtitle_height / fig_height ], 'FontSize',subaxes_values_font_size,'TickLength', [ 0.025 0.025 ],'XAxisLocation','top');
            xlabel(c,clabel_list(n,:),'Fontsize',subtitle_font_size );
            % c.Label.String = clabel_list{n,:};
            % set(c.Label.String,'Fontsize',subtitle_font_size );
        end
        if j == ncols
            text( 1.05, 2.5, model_list(n), 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', subtitle_font_size );
        end
        if i < nrows
            set( subplot_handles(n), 'XTickLabel', [] );
        end
        if j > 1
            set( subplot_handles(n), 'YtickLabel', [] );
        end
        set( subplot_handles(n), 'TickDir', 'out', 'TickLength', [ 0.02, 0.025 ] );%, 'LineWidth', 2 );
   end
end

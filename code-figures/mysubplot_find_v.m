function subplot_handles = mysubplot_find_v( fig_width, fig_height, title, xaxis, yaxis, subtitles_list, XData, YData_list, Y2Data_list )

set(0, 'Units', 'Centimeters' );
% 
% This function creats multiple subplots in a grid but with only 1 big
% x-axis and 1 big y-axis. There are multiple versions of this function,
% adapted to every single type of figure I'm interested in producing
% 
% Input:
%   - fig_width: width of figure in cm
%   - fig_height: height of figure in cm
%   - title: big title for the entire figure (if there is one)
%   - xaxis: label for the x-axis (one large one at the bottom)
%   - yaxis: label for the y-axis (one large one at the left)
%   - subtitles_list: list of subtitles given to each panel
%   - XData: common values for the x-axes
%   - YData_list: y-values for one of the curves (v)
%   - Y2Data_list: y-values for the other curve (vAH)
% 
% In particular, this code has been created when working on the paper:
% Pellis et al (2019), Nature Communications
% 
% Update: 11/10/2019

cov = 1; % Coverage = shrinking factor to decide how much of the actual drawing of the main plot covers the area dedicated to it
subcov = 1; % subCoverage = shrinking factor to decide how much of the actual drawing of each subplot covers the subarea dedicated to it
xcorr = 0.2; % Horizontal correction of the position of the x-axis label
[ n_subplots, ~ ] = size( subtitles_list );
[ ~, ~, check ] = size( YData_list );
if n_subplots ~= check
    error('Error, number of subtitles different from number of plots');
end

%%%%%%%%%%%%%%% Figure attributes
% The figure is made with a lot of rectangles ("placeholders") where to put
% all figure elements. However, then they are adapted manually, so all is a
% bit hacky...
title_height = 0; % put 0 if you don't want a title
title_font_size = 18;
title_ypos = 1/3; % Title y position in the rectangle devoted to it

xlabel_height = 2;
ylabel_width = 2;
axes_font_size = 18;
xlabel_ypos = 1/3;
ylabel_xpos = 1/2;

vllegend_width = 0; % space for vertical left legend
vrlegend_width = 1; % space for vertical right legend
htlegend_height = 0; % space for horizontal top legend
hblegend_height = 0.5; % space for horizontal bottom legend

%%%%%%%%%%%%%%% Subplots attributes
subtitle_height = 1;
subtitle_font_size = 15;
subtitle_ypos = 1/3; % Title y position in the rectangle devoted to it

xsublabel_height = 0;
ysublabel_width = 0.5;
subaxes_font_size = 12;
xsublabel_ypos = 0;
ysublabel_xpos = 0;

vlsublegend_width = 0;
vrsublegend_width = 0;
htsublegend_height = 0;
hbsublegend_height = 0;


screen = get(0,'ScreenSize');
screen_width = screen(3);
screen_height = screen(4);

xLeft = ( screen_width - fig_width ); yTop = ( screen_height - fig_height ) / 3;

switch n_subplots
    case 1
        nrows = 1;
        ncols = 1;
    case 2
        nrows = 1;
        ncols = 2;
    case 3
        nrows = 2;
        ncols = 2;
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
        nrows = 2;
        ncols = 4;
    case 8
        nrows = 2;
        ncols = 4;
    case 9
        nrows = 3;
        ncols = 3;
    case 10
        nrows = 3;
        ncols = 4;
    case 11
        nrows = 3;
        ncols = 4;
    case 12
        nrows = 3;
        ncols = 4;
    otherwise
        warning('Too many subplots: not plotting anything...');
end

figure
clf;
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
text( ( vllegend_width + ylabel_width * ( 1 - ylabel_xpos ) ) / fig_width, ( yArea + Area_height / 2.25 ) / fig_height, yaxis, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', axes_font_size );
text( ( xArea + vlsublegend_width + ysublabel_width - xcorr + Area_width / 2 ) / fig_width, ( hblegend_height + xlabel_height * xlabel_ypos ) / fig_height, xaxis, 'HorizontalAlignment', 'center', 'FontSize', axes_font_size);
text( ( xArea + vlsublegend_width + ysublabel_width - xcorr + Area_width / 2 ) / fig_width, ( yArea + Area_height + htlegend_height + title_height * title_ypos ) / fig_height, title, 'HorizontalAlignment', 'center', 'FontSize', title_font_size);
% colormap( jet );

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
        xsubtitle = xsubplot + subplot_width / 2;
        ysubtitle = ysubplot + hbsublegend_height + xsublabel_height + subArea_height + htsublegend_height + subtitle_height * subtitle_ypos;
        xsubxlabel = xsubplot + vlsublegend_width + ysublabel_width + subArea_width / 2;
        ysubxlabel = ysubplot + hbsublegend_height + xsublabel_height * xsublabel_ypos;
        xsubylabel = xsubplot + vlsublegend_width + ysublabel_width * ysublabel_xpos;
        ysubylabel = ysubplot + hbsublegend_height + xsublabel_height + subArea_height / 2;
        
        % text( xsubylabel / fig_width, ysubylabel / fig_height, 'y subaxis', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', subaxes_font_size );
        % text( xsubxlabel / fig_width, ysubxlabel / fig_height, 'x subaxis', 'HorizontalAlignment', 'center', 'FontSize', subaxes_font_size);
        text( xsubtitle / fig_width, ysubtitle / fig_height, subtitles_list(n,:), 'HorizontalAlignment', 'center', 'FontSize', subtitle_font_size);
              
        xsA = xsubArea / fig_width;
        ysA = ysubArea / fig_height;
        
        sahandle(n) = axes( 'position', [ xsp, ysp, spw, sph ], 'visible', 'off' );
        subplot_handles(n) = axes( 'position', [ xsA, ysA, sAw, sAh ], 'Fontsize', subaxes_font_size );
        set(gca,'ColorOrder',[ 0 0 1; 1 0 0 ]); % Set the first line in blue (adults) and the second in red (children)
        hold on;
        plot(XData,YData_list(:,:,n),'Linewidth',2);
        plot([0,1],[Y2Data_list(1,n),Y2Data_list(1,n)],'b-.','Linewidth',2);
        plot([0,1],[Y2Data_list(2,n),Y2Data_list(2,n)],'r-.','Linewidth',2);
        box on;
        if i < nrows
            set( subplot_handles(n), 'XTickLabel', [] );
        else
            set( subplot_handles(n), 'XTick', 0:0.2:1 );
        end
        if j > 1
            set( subplot_handles(n), 'YtickLabel', [] );
        end
        if ( i == 2 ) && ( j == 1 )
            legend('v^A_a','v^A_c','v^{AH}_a','v^{AH}_c','Location','SouthWest');
        end
   end
end

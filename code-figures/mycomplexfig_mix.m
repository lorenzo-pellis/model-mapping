function subplot_handles = mycomplexfig_mix( D, L, T )
% 
% This function creats multiple subplots in a grid but with only 1 big
% x-axis and 1 big y-axis. This is a variant of the more general function
% "mycomplexfig.m", because creating heatmaps for 3 subplots and one square
% grid creates some problems for which the solution needs to be tailored
% specifically.
% 
% Input:
%   - D: a structure containing figure Data
%   - L: a structure for Layout information for the figure structure
%   - T: a structure for Text (title(s), axe labels, legend text, etc.)
% 
% In particular, this code has been created when working on the paper:
% Pellis, L. et al (2020), Nature Communications
% 
% Update: 26-01-2020

show_layout = false; show_plots = true;
% show_layout = true; show_plots = false; % This line is here if you want to see the skeleton of the figure (with no plots)

screen = get(0,'ScreenSize');
screen_width = screen(3);
screen_height = screen(4);

xLeft = ( screen_width - L.fig_width ) / 2; yTop = ( screen_height - L.fig_height ) / 2;

figure
clf;
set(gcf,'Units','centimeters')
set(gcf,'Position',[xLeft yTop L.fig_width L.fig_height])

% outax = axes( 'position', [ 0 0 1 1 ], 'visible', 'off', 'Color', 'none', 'XLimMode', 'manual', 'YLimMode', 'manual', 'XLim', [ 0 1 ], 'YLim', [ 0 1 ] );
outax = axes( 'position', [ 0 0 1 1 ], 'visible', 'off' );
xlim( [ 0 1 ] )
ylim( [ 0 1 ] )

% Area = part of the figure where the suplots are drawn (excluding big
% title, big legends and big axes labels)
xArea = ( L.ylabel_width + L.vllegend_width );
Area_width = ( 1 - xArea - L.vrlegend_width );
yArea = ( L.xlabel_height + L.hblegend_height );
Area_height = ( 1 - yArea - L.title_height - L.htlegend_height );

if show_layout
    axes( 'position', [ xArea yArea Area_width Area_height ], 'Xtick', [], 'XTickLabel', [], 'Ytick', [], 'YTickLabel', [], 'box', 'on' );
    axes( outax )
end

% Title
xtitle = xArea;
ytitle = yArea + Area_height;
xctitle = xArea + Area_width / 2 + T.title_xcorr;
yctitle = yArea + Area_height + L.title_height * L.title_ypos + T.title_ycorr;
if show_layout
    rectangle( 'position', [ xtitle ytitle Area_width L.title_height ]);
%     viscircles( [ xctitle, yctitle ], 0.005, 'EdgeColor', 'k' );
end
if ~isempty(T.title)
    text( xctitle, yctitle, T.title, 'HorizontalAlignment', 'center', 'FontSize', T.title_font_size_sc);
end

% xaxis
xxaxis = xArea;
yxaxis = L.hblegend_height;
if show_layout
    rectangle( 'position', [ xxaxis yxaxis Area_width L.xlabel_height ]);
end
xcxaxis = xArea + Area_width / 2 + T.labelx_xcorr ;
ycxaxis = L.hblegend_height + L.xlabel_height * L.xlabel_ypos + T.labelx_ycorr ;
text( xcxaxis, ycxaxis, T.labelx, 'HorizontalAlignment', 'center', 'FontSize', T.axes_font_size_sc);

% yaxis
xyaxis = L.vllegend_width;
yyaxis = yArea;
if show_layout
    rectangle( 'position', [ xyaxis yyaxis L.ylabel_width Area_height ]);
end
xcyaxis = L.vllegend_width + L.ylabel_width * ( 1 - L.ylabel_xpos ) + T.labely_xcorr ;
ycyaxis = yArea + Area_height / 2 + T.labely_ycorr ;
text( xcyaxis, ycyaxis, T.labely, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', T.axes_font_size_sc );

% Figure letter
xlfigletter = xcyaxis - T.axes_font_size * 0.0352777778 / L.fig_width;
ycfigletter = yctitle;
if ~isempty(T.figletter)
    text( xlfigletter, ycfigletter, [T.figletter,')'], 'HorizontalAlignment', 'left', 'FontSize', T.title_font_size_sc );
end

% subpanel = area dedicated to each subplot (including title, axes labels
% and legends)
subpanel_width = ( Area_width - ( L.ncols - 1 ) * L.vseparator ) / L.ncols;
subpanel_height = ( Area_height - ( L.nrows - 1 ) * L.hseparator ) / L.nrows;

% subArea = the are dedicated to each drawn part of the subplot (i.e. the 
% subpanel, once you remove the subplot title, axes labels and legends)
subArea_width = subpanel_width - L.ysublabel_width - L.vlsublegend_width - L.vrsublegend_width;
subArea_height = subpanel_height - L.subtitle_height - L.xsublabel_height - L.hbsublegend_height - L.htsublegend_height;

lx = length(D.X);
ly = length(D.Y);

cell_width = subArea_width / lx;
cell_height = subArea_height / ly;

n = 0;
subplot_handles = zeros(L.nrows,L.ncols);
for i = 1:L.nrows
    for j = 1:L.ncols
        n = n+1;
        axes( outax );
        xsubpanel = xArea + ( j - 1 ) * ( subpanel_width + L.vseparator ); % / fig_width;
        ysubpanel = yArea + Area_height - i * subpanel_height - ( i - 1 ) * L.hseparator; % / fig_width;

        if show_layout
            rectangle( 'position', [ xsubpanel ysubpanel subpanel_width subpanel_height ]);
        end
        
        xsubArea = xsubpanel + L.ysublabel_width + L.vlsublegend_width;
        ysubArea = ysubpanel + L.xsublabel_height + L.hbsublegend_height;

        if show_layout
            rectangle( 'position', [ xsubArea ysubArea subArea_width subArea_height ]);
        end

        % Subtitle
        xsubtitle = xsubArea;
        ysubtitle = ysubArea + subArea_height + L.htsublegend_height;
        if show_layout
            rectangle( 'position', [ xsubtitle ysubtitle subArea_width L.subtitle_height ]);
        end        
        xcsubtitle = xsubtitle + subArea_width / 2 + T.subtitles_xcorr;
        ycsubtitle = ysubtitle + L.subtitle_height * L.subtitle_ypos + T.subtitles_ycorr;
        if ~isempty( T.subtitles )
            text( xcsubtitle, ycsubtitle, T.subtitles{i,j}, 'HorizontalAlignment', 'center', 'FontSize', T.subtitle_font_size_sc);
        end
        nsubletter = double(T.first_subletter) + n - 1;
        xrsubletter = xsubtitle + T.subletter_xcorr;
        ycsubletter = ycsubtitle + T.subletter_ycorr;
        if ~isempty( T.first_subletter )
            text( xrsubletter, ycsubletter, [char(nsubletter),')'], 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', T.subtitle_font_size_sc );
        end
            
        if n < L.nrows * L.ncols
            xsubplot = xsubArea + L.lsubmargin + cell_width / 2;
            ysubplot = ysubArea + L.bsubmargin + cell_height / 2;
            subplot_width = subArea_width - L.lsubmargin - L.rsubmargin - cell_width;
            subplot_height = subArea_height - L.bsubmargin - L.tsubmargin - cell_height;

            xsubcbarbox = xsubplot;
            ysubcbarbox = ysubArea + subArea_height;
            if show_layout
                rectangle( 'position', [ xsubcbarbox ysubcbarbox subplot_width L.htsublegend_height ]);
            end    
            
            subplot_handles(i,j) = axes( 'position', [ xsubplot ysubplot subplot_width subplot_height ], 'FontSize', T.subaxes_font_size_sc );


            %%%%%%%%%%%%%%%%% This is where the real plot occurs %%%%%%%%%%%%%%
            if show_plots
                colormap( 'jet' );
                [ XX, YY ] = meshgrid(D.X,D.Y);
                [~,cont] = contourf(XX,YY,D.Z{i,j},200);
                set(cont,'LineStyle','none');
                if ~isempty( D.Clim )
                    caxis( D.Clim );
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offsetcoeff = L.screen_scale * 0.025;

            set( subplot_handles(i,j), 'YTick', 1:0.5:4, 'FontSize', T.subaxes_font_size_sc );
            if i < L.nrows
                set( subplot_handles(i,j), 'XTickLabel', [] );
            else
                my_format_ticks( subplot_handles(i,j),[],[],[],[],[],[],offsetcoeff,[] );
                % my_format_ticks( subplot_handles(i,j),[],[],[],[],[],[],offsetcoeff,offsetcoeff );
            end
            if j > 1
                set( subplot_handles(i,j), 'YtickLabel', [] );
            % else
            %     my_format_ticks( subplot_handles(i,j),[],[],[],[],[],[],offsetcoeff,offsetcoeff );
            end
            set( subplot_handles(i,j), 'TickDir', 'out', 'TickLength', [ 0.02, 0.025 ] );


            %%%%%% Subcolorbar
            if L.subcolorbar
                xsubcbar = xsubcbarbox;
                ysubcbar = ysubcbarbox + L.subcolorbar_bottom_limit_in_box * L.htsublegend_height;
                subcbar_width = subplot_width;
                subcbar_height = L.subcolorbar_fraction_of_box * L.htsublegend_height;
                colorbar('Location','North',...
                    'Position', [ xsubcbar ysubcbar subcbar_width subcbar_height ],...
                    'Limits', D.Clim, 'FontSize', T.subaxes_font_size_sc, 'XAxisLocation', 'top' );
            end
            
        else % Last figures (bottom-right) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            xsubplot = xsubArea + L.lsubmargin;
            ysubplot = ysubArea + L.bsubmargin;
            subplot_width = subArea_width - L.lsubmargin - L.rsubmargin;
            subplot_height = subArea_height - L.bsubmargin - L.tsubmargin;

            xsubcbarbox = xsubplot;
            ysubcbarbox = ysubArea + subArea_height;
            if show_layout
                rectangle( 'position', [ xsubcbarbox ysubcbarbox subplot_width L.htsublegend_height ]);
            end
            
            subplot_handles(i,j) = axes( 'position', [ xsubplot ysubplot subplot_width subplot_height ], 'FontSize', T.subaxes_font_size_sc );


            %%%%%%%%%%%%%%%%% This is where the real plot occurs %%%%%%%%%%%%%%
            if show_plots               
                imagesc( D.X, D.Y, D.Z{i,j} )
                colormap(gca,[ 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0 ]);
                set(gca,'YDir','normal')
                hold on;
                % Make the grid
                limx = get( gca, 'XLim' );
                limy = get( gca, 'YLim' );
                xstep = D.X(2) - D.X(1);
                ystep = D.Y(2) - D.Y(1);
                for ix = 1:(lx-1)
                    xx = D.X(ix) + xstep / 2;
                    plot( [ xx xx ], limy, 'k', 'LineWidth', L.screen_scale * 0.1 );
                end
                for iy = 1:(ly-1)
                    yy = D.Y(iy) + ystep / 2;
                    plot( limx, [ yy yy ], 'k', 'LineWidth', L.screen_scale * 0.1 );
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if i < L.nrows
                set( subplot_handles(i,j), 'XTickLabel', [] );
            end
            if j > 1
                set( subplot_handles(i,j), 'YtickLabel', [] );
            end
            set( subplot_handles(i,j), 'TickDir', 'out', 'TickLength', [ 0.02, 0.025 ], 'FontSize', T.subaxes_font_size_sc );


            %%%%%% Subcolorbar
            if L.subcolorbar
                smallcorr = 0.002;
                xsubcbar = xsubcbarbox;
                ysubcbar = ysubcbarbox + L.subcolorbar_bottom_limit_in_box * L.htsublegend_height + smallcorr;
                subcbar_width = subplot_width;
                subcbar_height = L.subcolorbar_fraction_of_box * L.htsublegend_height - smallcorr;
                caxis( [ -0.5 4.5 ] );
                colorbar( 'Location', 'North',...
                    'Position', [ xsubcbar ysubcbar subcbar_width subcbar_height ],...
                    'XLim', [ -0.5 4.5 ], 'FontSize', T.legend_font_size_sc,...
                    'XAxisLocation', 'top', 'TickLength', 0,...
                    'XTick', 0:4, 'XTickLabel', T.clabels );
            end
        end
%         cbfreeze;
%         freezeColors;        
    end
end
% axes( outax )
% % labels = { 'Random mixing', 'Age', 'Either', 'Households', 'Both' };
% xlegendbox = xArea;
% ylegendbox = 0;
% legendbox_width = Area_width;
% legendbox_height = L.hblegend_height;
% if show_layout
%     rectangle( 'position', [ xlegendbox ylegendbox legendbox_width legendbox_height ] );
% end
% xcbar = xArea;
% ycbar = L.colorbar_bottom_limit_in_box * legendbox_height;
% cbar_width = Area_width;
% cbar_height = L.colorbar_fraction_of_box * legendbox_height;
% caxis( D.Clim );
% if isempty( T.clabels )
%     ch = colorbar( 'peer', outax, 'Location', 'South',...
%         'Position', [ xcbar ycbar cbar_width cbar_height ],...
%         'XLim', D.Clim, 'FontSize', T.legend_font_size_sc, 'XAxisLocation', 'bottom' );
% else
%     ch = colorbar( 'peer', outax, 'Location', 'South',...
%         'Position', [ xcbar ycbar cbar_width cbar_height ],...
%         'XLim', D.Clim, 'FontSize', T.legend_font_size_sc, 'XAxisLocation', 'bottom',...
%         'TickLength', [0 0], 'XTick', 0:4, 'XTickLabel', T.clabels );
% end
e = 0;

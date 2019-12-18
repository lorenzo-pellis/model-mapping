function subplot_handles = mycomplexfig( plottype, D, L, T )
% 
% This function creats multiple subplots in a grid but with only 1 big
% x-axis and 1 big y-axis. There are multiple versions of similar functions
% (for different specific plots), but this is one of the most general
% functions of this sort that I used
% 
% Input:
%   - plottype: which actual plot I want to draw, choosing from:
%       - heatmap: normal heatmap
%       - OAR_over: Overall Acceptance Region plot, with lines overlaid
%       - OAR_grid: Overall Acceptance Region plot, with cell grid
%       - dyn: epidemic dynamics
%       - ROT_table: table for the Rule Of Thumb
%       - ROT_line: line for the Rule of Thumb
%       - line_fit: line fit
%       - metaROT: meta-Rule Of Thumb
%   - D: a structure containing figure Data
%   - L: a structure for Layout information for the figure structure
%   - T: a structure for Text (title(s), axe labels, legend text, etc.)
% 
% In particular, this code has been created when working on the paper:
% Pellis et al (2019), Nature Communications
% 
% Update: 13/10/2019

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

% outax = outer axes handle, for the axes along the border of the actual figure
outax = axes( 'position', [ 0 0 1 1 ], 'visible', 'off' );
xlim( [ 0 1 ] )
ylim( [ 0 1 ] )

% Area = part of the figure where the suplots are drawn (excluding big
% title, big legends and big axes labels)
xArea = ( L.ylabel_width + L.vllegend_width );
Area_width = ( 1 - xArea - L.vrlegend_width - L.y2label_width - L.rowtitle_width );
yArea = ( L.xlabel_height + L.hblegend_height );
Area_height = ( 1 - yArea - L.title_height - L.htlegend_height - L.coltitle_height );

if show_layout
    axes( 'position', [ xArea yArea Area_width Area_height ], 'Xtick', [], 'XTickLabel', [], 'Ytick', [], 'YTickLabel', [], 'box', 'on' );
    axes( outax )
end

% Title
xtitle = xArea;
ytitle = yArea + Area_height;
xctitle = xArea + Area_width / 2 + T.title_xcorr;
yctitle = yArea + Area_height + L.coltitle_height + L.title_height * L.title_ypos + T.title_ycorr;
if show_layout && ( L.title_height > 0 )
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
if ~isempty(T.labelx)
    text( xcxaxis, ycxaxis, T.labelx, 'HorizontalAlignment', 'center', 'FontSize', T.axes_font_size_sc);
end

% yaxis
xyaxis = L.vllegend_width;
yyaxis = yArea;
if show_layout
    rectangle( 'position', [ xyaxis yyaxis L.ylabel_width Area_height ]);
end
xcyaxis = L.vllegend_width + L.ylabel_width * ( 1 - L.ylabel_xpos ) + T.labely_xcorr ;
ycyaxis = yArea + Area_height / 2 + T.labely_ycorr ;
if ~isempty(T.labely)
    text( xcyaxis, ycyaxis, T.labely, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', T.axes_font_size_sc );
end
% y2axis
if L.y2axis
    xy2axis = L.vllegend_width + L.y2label_width + Area_width;
    yy2axis = yArea;
    if show_layout
        rectangle( 'position', [ xy2axis yy2axis L.y2label_width Area_height ]);
    end
    xcy2axis = xy2axis + L.y2label_width * ( 1 - L.y2label_xpos ) + T.labely2_xcorr ;
    ycy2axis = yy2axis + Area_height / 2 + T.labely2_ycorr ;
    text( xcy2axis, ycy2axis, T.labely2, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', T.axes_font_size_sc );
end

% Figure letter
xlfigletter = xcyaxis - T.axes_font_size * 0.0352777778 / L.fig_width + T.figletter_xcorr;
if L.title_height > 0
    ycfigletter = yctitle + T.figletter_ycorr;
else
    if L.coltitle_height > 0
        ycfigletter = yctitle - L.coltitle_height / 2 + T.figletter_ycorr;
    else
        ycfigletter = yctitle - L.subtitle_height / 2 + T.figletter_ycorr;
    end
end
if ~isempty(T.figletter)
    text( xlfigletter, ycfigletter, [T.figletter,')'], 'HorizontalAlignment', 'left', 'FontSize', T.title_font_size_sc );
end

if ~isempty( D.colormap )
    colormap( D.colormap );
end

% subpanel = area dedicated to each subplot (including title, axes labels
% and legends)
subpanel_width = ( Area_width - ( L.ncols - 1 ) * L.vseparator ) / L.ncols;
subpanel_height = ( Area_height - ( L.nrows - 1 ) * L.hseparator ) / L.nrows;

% subArea = the are dedicated to each drawn part of the subplot (i.e. the 
% subpanel, once you remove the subplot title, axes labels and legends)
subArea_width = subpanel_width - L.ysublabel_width - L.vlsublegend_width - L.vrsublegend_width;
subArea_height = subpanel_height - L.subtitle_height - L.xsublabel_height - L.hbsublegend_height - L.htsublegend_height;

n = 0;
subplot_handles = zeros(L.nrows,L.ncols);
for i = 1:L.nrows
    for j = 1:L.ncols
        n = n+1;
        axes( outax );
        xsubpanel = xArea + ( j - 1 ) * ( subpanel_width + L.vseparator );
        ysubpanel = yArea + Area_height - i * subpanel_height - ( i - 1 ) * L.hseparator;

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
        if ~isempty( T.subtitles{i,j} )
            text( xcsubtitle, ycsubtitle, T.subtitles{i,j}, 'HorizontalAlignment', 'center', 'FontSize', T.subtitle_font_size_sc);
        end
        
        % Column title
        if i == 1
            xcoltitle = xsubArea;
            ycoltitle = ysubArea + subArea_height + L.htsublegend_height;
            if show_layout
                rectangle( 'position', [ xcoltitle ycoltitle subArea_width L.coltitle_height ]);
            end        
            xccoltitle = xcoltitle + subArea_width / 2 + T.coltitles_xcorr;
            yccoltitle = ycoltitle + L.coltitle_height * L.coltitle_ypos + T.coltitles_ycorr;
            if ~isempty( T.coltitles{j} )
                text( xccoltitle, yccoltitle, T.coltitles{j}, 'HorizontalAlignment', 'center', 'FontSize', T.subtitle_font_size_sc);
            end
        end

        % Row title
        if j == L.ncols
            xrowtitle = xsubArea + subArea_width + L.vlsublegend_width;
            yrowtitle = ysubArea;
            if show_layout
                rectangle( 'position', [ xrowtitle yrowtitle L.rowtitle_width subArea_height ]);
            end        
            xcrowtitle = xrowtitle + L.rowtitle_width * L.rowtitle_xpos + T.rowtitles_xcorr;
            ycrowtitle = yrowtitle + subArea_height / 2 + T.rowtitles_ycorr;
            if ~isempty( T.rowtitles{i} )
                text( xcrowtitle, ycrowtitle, T.rowtitles{i}, 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontSize', T.subtitle_font_size_sc);
            end
        end
        
        
        if strcmp(plottype,'dyn') % Add colored square
            sqw = L.screen_scale * 0.007;
            sqh = sqw;
            xsq = xsubArea + subArea_width - sqw - 0.05;
            ysq = ycsubtitle - sqh / 2 + 0.005;
            listorder = [ 2 1 4 3 ];
            rectangle( 'position', [ xsq ysq sqw sqh ], 'FaceColor', D.colorlist(listorder(n),:), 'LineWidth', 1 );
        end
        nsubletter = double(T.first_subletter) + n - 1;
        xrsubletter = xsubtitle + T.subletter_xcorr;
        ycsubletter = ycsubtitle + T.subletter_ycorr;
        if ~isempty( T.first_subletter )
            text( xrsubletter, ycsubletter, [char(nsubletter),')'], 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', T.subtitle_font_size_sc );
        end
        
        xsubplot = xsubArea + L.lsubmargin;
        ysubplot = ysubArea + L.bsubmargin;
        subplot_width = subArea_width - L.lsubmargin - L.rsubmargin;
        subplot_height = subArea_height - L.bsubmargin - L.tsubmargin;
        subplot_handles(i,j) = axes( 'position', [ xsubplot ysubplot subplot_width subplot_height ], 'FontSize', T.subaxes_font_size_sc );
        
        
        %%%%%%%%%%%%%%%%% This is where the real plot occurs %%%%%%%%%%%%%%
        if show_plots
            if strcmp(plottype,'heatmap')
                [~,cont] = contourf(D.X,D.Y,D.Z{i,j},200);
                set(cont,'LineStyle','none');
                if ~isempty( D.Clim )
                    caxis( D.Clim );
                end
                if ~isempty( D.Zover )
                    if numel( D.Zcont ) == 1
                        D.Zcont = [ D.Zcont D.Zcont ];
                    end
                    hold on;
                    contour(D.X(:,2:end),D.Y(:,2:end),D.Zover{i,j}(:,2:end),D.Zcont,'color','k','Linestyle','--','Linewidth',L.screen_scale*0.5,'ShowText','off');
                    plot(D.X(:,1)+0.00001,D.Y(:,1),'color','k','Linestyle','--','Linewidth',L.screen_scale*0.5);
                end
                if ~isempty(D.Xticks)
                    set(gca,'XTick',D.Xticks);
                end
                if ~isempty(D.Yticks)
                    set(gca,'YTick',D.Yticks);
                end
            elseif strcmp(plottype,'OAR_over')
                imagesc( D.X, D.Y, D.Z{i,j} )
                set(gca,'YDir','normal')  
                if ~isempty( D.Zover )
                    for l = 1:D.nZover
                        currZover = D.Zover{l};
                        if numel( D.Zcont{l} ) == 1
                            D.Zcont{l} = [ D.Zcont{l} D.Zcont{l} ];
                        end
                        hold on;
                        [CC,hh] = contour(D.X,D.Y,currZover{i,j},D.Zcont{l},'color','k','Linewidth',L.screen_scale * D.Zlinewidth{l},'ShowText',D.Zshowtext{l},'LabelSpacing', 200);
%                         cl = clabel( CC, hh );
                    end
                end
                caxis( D.Clim );
                if ~isempty(D.Xticks)
                    set(gca,'XTick',D.Xticks);
                end
                if ~isempty(D.Yticks)
                    set(gca,'YTick',D.Yticks);
                end
                set( gca, 'FontSize', T.subaxes_font_size_sc )
            elseif strcmp(plottype,'OAR_grid')
                imagesc( D.X, D.Y, D.Z{i,j} )
                set(gca,'YDir','normal')
                hold on;
                % Make the grid
                limx = get( gca, 'XLim' );
                limy = get( gca, 'YLim' );
                lx = length(D.X);
                ly = length(D.Y);
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
                caxis( D.Clim );
                if ~isempty(D.Xticks)
                    set(gca,'XTick',D.Xticks);
                end
            elseif strcmp(plottype,'dyn')
                Y1 = D.Y{i,j}(:,1:4);
                Y2 = [ D.Y{i,j}(:,5:8), 100*ones(41,1) ];
                [axyy,h1,h2] = plotyy(D.X,Y1,D.X,Y2);
                for l1 = 1:4
                    set( h1(l1), 'Color', D.colorlist(l1,:), 'Linestyle', '-', 'Linewidth', L.screen_scale * 0.7 );
                end
                for l2 = 1:4
                    set( h2(l2), 'Color', D.colorlist(l2,:), 'Linestyle', '--', 'Linewidth', L.screen_scale * 0.7 );
                end
                set( h2(5), 'LineStyle','-.','Color','k', 'LineWidth', L.screen_scale * 0.5 );
                set( axyy(1), 'XLim',[ min(D.X),max(D.X) ], 'YLim', [ 0 7 ], 'YTick', 0:7, 'TickLength', [0 0], 'YGrid', 'on', 'Color', 0.9 * ones(1,3) );
                set( axyy(2), 'XLim',[ min(D.X),max(D.X) ], 'YLim', [ 0 120 ], 'YTick', 0:20:120, 'TickDir','out' );
                if n == L.nrows * L.ncols % Add legend
                    if ~isempty(T.legentries)
                        legend( T.legentries );
                    end
                end
                if ~isempty(D.Xticks)
                    set(gca,'XTick',D.Xticks);
                end
            elseif strcmp(plottype,'ROT_table')
                my_colored_tablestr( D.X, D.Y, D.Z{1}{i,j}, T.table_incells_font_size_sc, T.table_outcells_font_size_sc, [], [], true, D.Z{2}{i,j}, D.Clim, ['R_0'], D.colormap, [] ); %['   R_0 to';'  from   ']
            elseif strcmp(plottype,'ROT_line')
                hold on;
                whichk = D.Z{i,j};
                plot(D.X{i,j}(~whichk),D.Y{i,j}(~whichk),'o','color',0.6*ones(1,3),'markersize',8,'linewidth',0.1);
                plot(D.X{i,j}(whichk),D.Y{i,j}(whichk),'ko','markersize',8,'linewidth',0.1);
                xlim(D.limx)
                ylim(D.limy)
                plot([0 4.1],[D.Zover{1}{i,j}(1),D.Zover{1}{i,j}(2)*4+D.Zover{1}{i,j}(1)],'--','color',0.6*ones(1,3),'Linewidth',1.5);
                plot([0 4.1],[D.Zover{2}{i,j}(1),D.Zover{2}{i,j}(2)*4+D.Zover{2}{i,j}(1)],'k--','Linewidth',1.5);
                plot([0 4.1],[D.Zover{3}{i,j}(1),D.Zover{3}{i,j}(2)*4+D.Zover{3}{i,j}(1)],'k','Linewidth',2.5);
                if L.legend
                    if i*j == 1
                        legend(T.legentries,'Location','NorthWest','fontsize',T.legend_font_size_sc)
                    end
                end
            elseif strcmp(plottype,'line_fit')
                hold on;
                plot(D.X(D.Z{i,j}),D.Y(D.Z{i,j}),'ko','markersize',6);
                xlim(D.limx)
                ylim(D.limy)
                set(gca,'box','on');
                if ~isempty(D.Xticks)
                    set(gca,'XTick',D.Xticks);
                end
                plot([0 4.1],[D.Zover{i,j}(1),D.Zover{i,j}(2)*4+D.Zover{i,j}(1)],'k--');
                str1 = sprintf('% 5.1f',D.Zover{i,j}(2));
                str2 = sprintf('% 5.1f',D.Zover{i,j}(1));
                text(0.2,85,['\beta_1 = ',str1;'\beta_0 = ',str2]);
            elseif strcmp(plottype,'metaROT')
                hold on;
                xlim(D.limx)
                markerlist = { 'o', '^', 's' };
                for c = 1:3
                    plot(D.X(1,c),D.Y(1,c),markerlist{c},'markeredgecolor','k','markerfacecolor','w','markersize',6,'visible','on');
                end
                for e = 1:3
                    plot(D.X(e,1:3),D.Y(e,1:3),D.colorlist{e},'linewidth',3);
                end
                for e = 1:3
                    for c = 1:3
                        plot(D.X(e,c),D.Y(e,c),markerlist{c},'markeredgecolor',D.colorlist{e},'markerfacecolor','w','markersize',10,'linewidth',2);
                    end
                end
                legend(T.legentries,'Location','NorthEast')
                set(gca,'box','on');
            end
            set( gca, 'FontSize', T.subaxes_font_size_sc ) 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(plottype,'dyn') && show_plots
            if i < L.nrows
                set( axyy(1), 'XTickLabel', [] );
                set( axyy(2), 'XTickLabel', [] );
            else
                set( axyy(2), 'XTickLabel', [] );     
            end
            if j > 1
                set( axyy(1), 'YTickLabel', [] );
                set( axyy(2), 'FontSize', T.subaxes_font_size_sc );            
            end
            if j < L.ncols
                set( axyy(2), 'YTickLabel', [] );
            end
        else
            if i < L.nrows
                set( subplot_handles(i,j), 'XTickLabel', [] );
            end
            if j > 1
                set( subplot_handles(i,j), 'YtickLabel', [] );
            end
            set( subplot_handles(i,j), 'TickDir', 'out', 'TickLength', [ 0.02, 0.025 ] );
        end
        
        %%%%%% Subcolorbar
        if L.subcolorbar
            xsubcbar = xsubplot;
            ysubcbar = ysubplot + subplot_height + L.subcolorbar_bottom_limit_in_box * L.htsublegend_height;
            subcbar_width = subplot_width;
            subcbar_height = L.subcolorbar_fraction_of_box * L.htsublegend_height;
            ch = colorbar( 'peer', outax, 'Location', 'North',...
                'Position', [ xsubcbar ysubcbar subcbar_width subcbar_height ],...
                'XLim', D.Climleg, 'FontSize', T.subaxes_font_size_sc, 'XAxisLocation', 'top' );
            caxis( outax, D.Clim )
        end
   end
end
if L.colorbar
    axes( outax )
    xlegendbox = xArea;
    ylegendbox = 0;
    legendbox_width = Area_width;
    legendbox_height = L.hblegend_height;
    if show_layout
        rectangle( 'position', [ xlegendbox ylegendbox legendbox_width legendbox_height ] );
    end
    xcbar = xArea;
    ycbar = L.colorbar_bottom_limit_in_box * legendbox_height;
    cbar_width = Area_width;
    cbar_height = L.colorbar_fraction_of_box * legendbox_height;
    caxis( D.Clim );
    if isempty( T.clabels )
        ch = colorbar( 'peer', outax, 'Location', 'South',...
            'Position', [ xcbar ycbar cbar_width cbar_height ],...
            'XLim', D.Climleg, 'FontSize', T.legend_font_size_sc, 'XAxisLocation', 'bottom' );
    else
        ch = colorbar( 'peer', outax, 'Location', 'South',...
            'Position', [ xcbar ycbar cbar_width cbar_height ],...
            'XLim', D.Climleg, 'FontSize', T.legend_font_size_sc, 'XAxisLocation', 'bottom',...
            'TickLength', 0, 'XTick', 0:4, 'XTickLabel', T.clabels );
    end
    if ~isempty( T.cbar_string )
        ch.Label.String = T.cbar_string;
        if ~isempty( T.cbar_string_font_size )
            ch.Label.FontSize = T.cbar_string_font_size_sc;
        else
            ch.Label.FontSize = T.legend_font_size_sc;
        end
    end
            
end

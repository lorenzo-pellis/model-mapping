function my_colored_tablestr( b, a, T, infont, outfont, interp, grid, box, C, limc, toplefttext, cmap, nancolor )
% 
% This function creates a table and colours cells based ont he value in
% them. This function is used in the paper:
% Pellis, L. et al (2020), Nature Communications
%
% Inputs:
% Required:
%   - b: the values for each column (i.e. the x-axis of the table)
%   - a: the values for each row (i.e. the y-axis of the table)
%   - T: the table itself, as a cell of strings
% Optional:
%   - infont: font size of text inside cells
%   - outfont: font size of text outside cells
%   - interp: interpreter (tex,latex) used by Matlab to produce text
%   - grid: draw the grid only on the inside of the table or also in the
%   title rows and columns ('in','out')
%   - box: draw a box around the table? (true/false)
%   - C: matrix of colours for each table cell
%   - limc: limits for the colorbar
%   - toplefttext: text to put in the top-left corner of the table
%   - cmap: colormap
%   - nancolor: colour to be used for cells that contain NaN
% 
% Update: 28-06-2019

if ( nargin < 13 ) || isempty(nancolor)
    nancolor = [ 0 0 0 ];
end
if ( nargin < 12 ) || isempty(cmap)
    cmap = 'hot';
end
if nargin < 11
    toplefttext = [];
end
if nargin < 10 || isempty(limc)
    limc = [ min(min(T)),max(max(T)) ];
end
if nargin < 9 || isempty(C)
    C = ( T - limc(1) ) / ( limc(2) - limc(1) );
    C(C<=0) = 0.001;
    C(C>1) = 1;
end
if nargin < 8 || isempty(box)
    box = true;
end
if nargin < 7 || isempty(grid)
    grid = 'out';
end
if nargin < 6 || isempty(interp)
    interp = 'tex';
end
if nargin < 5 || isempty(outfont)
    font = 12;
end
if nargin < 4 || isempty(infont)
    font = 10;
end
if isempty(a) || isempty(b)
    if ~isempty(toplefttext)
        disp('Ignoring text in top left of table!');
        toplefttext = [];
    end
end

[ la, lb ] = size(T);
% la = number of lines
% lb = number of columns
% a = left column titles
% b = top row titles

TC = NaN((1+la),(1+lb),3); % Table for the colours (array)
TM = cell((1+la),(1+lb)); % Table for the matrix

TM{1,1} = toplefttext;
TC(1,1,:) = 1;
if ~isempty(b)
    for ib = 1:lb
        TM{1,ib+1} = num2str(b(ib));
        TC(1,ib+1,:) = 1;
    end
end
if ~isempty(a)
    for ia = 1:la
        TM{ia+1,1} = num2str(a(ia));
        TC(ia+1,1,:) = 1;
    end
end

map = colormap(cmap);
Cscaled = ( C - limc(1) ) / ( limc(2) - limc(1) );
for ia = 1:la
    for ib = 1:lb
        if isnan(C(ia,ib))
            TM{ia+1,ib+1} = '';
%             assert(strcmp(T{ia,ib},''));
            TC(ia+1,ib+1,:) = nancolor;
        else
            TM{ia+1,ib+1} = T{ia,ib};
            TC(ia+1,ib+1,:) = reshape(map(ceil(64*Cscaled(ia,ib)),:),[1 1 3]);
        end
    end
end


% figure(1)
% clf;
hold off;
im = imagesc(TC);
hold on;
if ~isempty(toplefttext)
    text(1,1,toplefttext,'HorizontalAlignment','center','interpreter',interp,'FontSize',outfont);
end
if isempty(a)
    limx = [ 1.5, lb + 1.5 ];
    for ia = 1:la
        plot(limx,[ia,ia]+0.5,'k');
    end
else
    limx = [ 0.5, lb + 1.5 ];
    for ia = 1:la
        text(1,ia+1,TM{ia+1,1},'HorizontalAlignment','center','interpreter',interp,'FontSize',outfont);
        if strcmp(grid,'out')
            plot(limx,[ia,ia]+0.5,'k');
        elseif strcmp(grid,'in')
            plot([limx(1)+0.5,limx(2)],[ia,ia]+0.5,'k');
        end
    end
end
if isempty(b)
    limy = [ 1.5, la + 1.5 ];
    for ib = 1:lb
        plot([ib,ib]+0.5,limy,'k');
    end
else
    limy = [ 0.5, la + 1.5 ];
    for ib = 1:lb
        text(ib+1,1,TM{1,ib+1},'HorizontalAlignment','center','interpreter',interp,'FontSize',outfont);
        if strcmp(grid,'out')
            plot([ib,ib]+0.5,limy,'k');
        elseif strcmp(grid,'in')
            plot([ib,ib]+0.5,[limy(1)+0.5,limy(2)],'k');
        end
    end
end
for ia = 1:la
    for ib = 1:lb
        text(ib+1,ia+1,TM{ia+1,ib+1},'HorizontalAlignment','center','interpreter',interp,'FontSize',infont);
    end
end
set(gca,'Xtick',[],'Xticklabel',[],'YTick',[],'YTickLabel',[])
if box
    set(gca,'box','on');
else
    set(gca,'box','off');
end
xlim(limx)
ylim(limy)




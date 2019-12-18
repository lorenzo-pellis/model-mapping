function D = set_mycomplexfig_data( D )
% 
% This is a hard-made functions written in support of "mycomplexfig", a
% function I use to make plots that have multiple panels showing variants
% of the same output. I mostly use it to create the figures used in the
% model mapping paper:
% 
% Reference: Pellis et al (2019), Nature Communications
% 
% Update: 22/06/2019

%%%%%% Data
if ~isfield(D,'X')
    D.X = [ 0 1 ];
end
if ~isfield(D,'Y')
    D.Y = [ 0 1 ];
end
if ~isfield(D,'Z')
    D.Z = [];
end
if ~isfield(D,'nZover')
    D.nZover = 0;
end
if ~isfield(D,'Zover')
    if D.nZover == 0
        D.Zover = [];
    else
        D.Zover = cell(1,D.nZover);
        for i = 1:D.nZover
            D.Zover{i} = [];
        end
    end
end
if ~isfield(D,'Zcont')
    if D.nZover == 0
        D.Zcont = [];
    else
        D.Zcont = cell(1,D.nZover);
        for i = 1:D.nZover
            D.Zcont{i} = [];
        end
    end
end
if ~isfield(D,'Zshowtext')
    if D.nZover == 0
        D.Zshowtext = [];
    else
        D.Zshowtext = cell(1,D.nZover);
        for i = 1:D.nZover
            D.Zshowtext{i} = 'off';
        end
    end
end
if ~isfield(D,'limx')
    D.limx = [];
end
if ~isfield(D,'limy')
    D.limy = [];
end
if ~isfield(D,'limz')
    D.limz = [];
end
if ~isfield(D,'Clim')
    D.Clim = [];
end
if ~isfield(D,'Climleg')
    D.Climleg = D.Clim;
end

% Tick labels
if ~isfield(D,'Xticks')
    D.Xticks = [];
end
if ~isfield(D,'Yticks')
    D.Yticks = [];
end
if ~isfield(D,'Zticks')
    D.Zticks = [];
end


% Colormap
if ~isfield(D,'colormap')
    D.colormap = 'jet';
end
if ~isfield(D,'colorlist')
    D.colorlist = 'k';
end


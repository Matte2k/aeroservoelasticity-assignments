function plt = plotGrid(ax, Coord, dotsize, style)
%PLOT GRID - plot reference system vectors
%
%
    if nargin == 2
        plt = plot3(ax, Coord(1), Coord(2), Coord(3), '.', 'MarkerSize', 30);
    elseif nargin == 3 && isnumeric(dotsize)
        plt = plot3(ax, Coord(1), Coord(2), Coord(3), '.', 'MarkerSize', dotsize);
    elseif nargin == 3 && isnumeric(dotsize) == 0
        plt = plot3(ax, Coord(1), Coord(2), Coord(3), '.', 'MarkerSize', 60);
    else
        plt = plot3(ax, Coord(1), Coord(2), Coord(3), style, 'MarkerSize', dotsize);
    end
    
    view(3)
    ax.Visible = 'off';
    ax.Clipping = 'off';
    axis equal

return

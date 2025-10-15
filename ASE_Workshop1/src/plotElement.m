function handle = plotElement(ax, nodePos, connectivity, handle)
%PLOT ELEMENT
% Function to draw a 3D mesh given nodes and a connectivity matrix
% Inputs:
%   nodes - N-by-3 matrix, where N is the number of nodes, each row is [x,y,z]
%   connectivity - M-by-2 matrix, where each row defines an edge between two nodes

if isempty(handle)
    handleNew = [];
end

% Extract the number of faces from the connectivity matrix
numFaces = size(connectivity, 1);

% Loop over each face and plot
for i = 1:numFaces
    % Get the indices of the four nodes that form the current quad face
    node1Idx = connectivity(i, 1);
    node2Idx = connectivity(i, 2);
    node3Idx = connectivity(i, 3);
    node4Idx = connectivity(i, 4);
    
    % Get the coordinates of the four nodes
    node1 = nodePos(node1Idx, :);
    node2 = nodePos(node2Idx, :);
    node3 = nodePos(node3Idx, :);
    node4 = nodePos(node4Idx, :);

    if isempty(handle)
        % Plot the quad face as a patch (a filled polygon)
        handleNew = [handleNew, fill3(ax,[node1(1), node2(1), node3(1), node4(1)], ...  % X coordinates
                                         [node1(2), node2(2), node3(2), node4(2)], ...  % Y coordinates
                                         [node1(3), node2(3), node3(3), node4(3)], ...  % Z coordinates
                                         'b', 'FaceAlpha', 0.3, 'EdgeColor', 'k')];    % Set face color and edge color
    else
        set(handle(i), 'Vertices', [[node1(1); node2(1); node3(1); node4(1)], ...
                                    [node1(2); node2(2); node3(2); node4(2)], ...
                                    [node1(3); node2(3); node3(3); node4(3)]]);
    end
end

if isempty(handle)
    handle = handleNew;
end

return

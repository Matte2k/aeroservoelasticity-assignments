function handle = plotElement(ax, nodePos, connectivity, handle)

% Function to draw a 3D mesh given nodes and a connectivity matrix
% Inputs:
%   nodes - N-by-3 matrix, where N is the number of nodes, each row is [x,y,z]
%   connectivity - M-by-2 matrix, where each row defines an edge between two nodes

if isempty(handle)
    handleNew = [];
end



if isempty(handle)
    handleNew = trisurf(connectivity, nodePos(:,1), nodePos(:,2), nodePos(:,3));

    set(handleNew, 'parent', ax);
    
    set(handleNew, 'facecolor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'k')
else
    set(handle, 'Vertices', nodePos);
end

if isempty(handle)
    handle = handleNew;
end

return

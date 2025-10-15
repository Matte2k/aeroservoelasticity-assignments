function plotFrame(ax, frame, scale)
%PLOT FRAME - frame plot for figure
%
%
    orig = frame.orig;
    xaxis = scale*frame.xaxis;
    yaxis = scale*frame.yaxis;
    zaxis = scale*frame.zaxis;
    
    hold on;
    axis square;    axis equal;
    %grid minor
    
    qx = quiver3(ax, orig(1), orig(2), orig(3), xaxis(1), xaxis(2), xaxis(3), 'r');
    qx.MaxHeadSize = 0.5;
    qy = quiver3(ax, orig(1), orig(2), orig(3), yaxis(1), yaxis(2), yaxis(3), 'g');
    qy.MaxHeadSize = 0.5;
    qz = quiver3(ax, orig(1), orig(2), orig(3), zaxis(1), zaxis(2), zaxis(3), 'b');
    qz.MaxHeadSize = 0.5;

return

function [fig] = modePlot(plotModel,time,zoomFactor,startMode,modeShapes,eigenvalues)
%MODE PLOT
%   Function creating a figure of the deformed model based on
%   modeShapes and eigenvalues of the considered model
%
%   SYNTAX:
%       [fig] = modePlot(plotModel, time, zoomFactor, startMode, modeShapes, eigenvalues)
%
%   INPUT:
%       plotModel,   struct: model struct containing grid data
%       time;        double: time instant of the deformated picture
%       zoomFactor,  double: amplification factor of the displacement
%       startMode,   double: starting mode to plot in the animation
%       modeShapes,  double: modeShape in modal coordinate of the model
%       eigenvalues, double: eigenvalues of the considered system
%
%   OUTPUT:
%       fig, figure: plot of the deformated model
%
%

    % Optional input
    if nargin == 5
        complex = false;
    else
        complex = true;
    end
    
    timeRad = time.*(pi);
    
    % Create new figure to overlap the two deformated model
    figName = ['displacements mode ', num2str(startMode)];
    fig = figure(Name=figName);
    tl = tiledlayout(1,length(time));
        tl.TileSpacing = 'none';
        tl.Padding = 'none';
    
    for i = 1:length(time) 
        % compute dispacements from modeShapes and eventually eigs
        if complex
            displacement = 2 * real(modeShapes(:,startMode) * exp(1i*imag(eigenvalues(startMode))*timeRad(i))) * zoomFactor;
        else
            displacement = modeShapes(:,startMode)*sin(timeRad(i))*zoomFactor;
        end

        [~, ~, figUndeformed] = meshVisualization(plotModel, []);
        [~, ~, figDeformed]   = meshVisualization(plotModel, displacement);
    
        ax = nexttile(tl);
        ax.Visible = "on";
        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.ZColor = 'none';
        ax.Clipping = "off";
        view(30, 25);
        hold(ax, "on")
    
        % Recover deformed and undeformed model plot
        fill3_objs1 = findobj(figUndeformed, 'Type', 'patch');
        fill3_objs2 = findobj(figDeformed, 'Type', 'patch');
                
        % Plot the undeformed model
        for j = 1:length(fill3_objs1)
            x1 = get(fill3_objs1(j), 'XData');
            y1 = get(fill3_objs1(j), 'YData');
            z1 = get(fill3_objs1(j), 'ZData');
            c1 = get(fill3_objs1(j), 'FaceColor');
            fill3(ax, x1, y1, z1, c1, 'EdgeColor', [.7 .7 .7], 'LineWidth',0.5);
        end

        % Plot the deformed model
        for j = 1:length(fill3_objs2)
            x2 = get(fill3_objs2(j), 'XData');
            y2 = get(fill3_objs2(j), 'YData');
            z2 = get(fill3_objs2(j), 'ZData');
            c2 = get(fill3_objs2(j), 'FaceColor');
            fill3(ax, x2, y2, z2, c2, 'LineWidth',0.8);
        end
    
        % Plot basic reference frame
        frame.orig  = [6; 0; 0];
        frame.xaxis = [1; 0; 0];
        frame.yaxis = [0; 1; 0];
        frame.zaxis = [0; 0; 1];
        plotFrame(ax, frame, 1)
        
        % set FOV of the plot
        xlim([4.75 9.25]);    ylim([0 5]);    zlim([-1 2]);

        % % set title of the current visualization
        % plotTitle = ['Time ', strrep(rats(time(i)),' ',''),'$\pi$'];
        % title(ax,plotTitle);
        
        % close building figures
        close(figUndeformed);
        close(figDeformed);
    end

    hold off;

end


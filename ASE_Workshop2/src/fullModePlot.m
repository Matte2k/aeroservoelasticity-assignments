function [fig] = fullModePlot(deformedModel,time,zoomFactor,startMode,modeShapes,eigenvalues)
%FULL MODE PLOT
%   Function creating a figure of the deformed model based on
%   modeShapes and eigenvalues of the considered model
%
%   SYNTAX:
%       [fig] = fullModePlot(deformedModel,staticModel,time,zoomFactor,startMode,modeShapes,eigenvalues)
%
%   INPUT:
%       deformedModel, struct: model struct containing deforming grid data
%       staticModel,   struct: model struct containing static grid data
%       time;          double: time instant of the deformated picture
%       zoomFactor,    double: amplification factor of the displacement
%       startMode,     double: starting mode to plot in the animation
%       modeShapes,    double: modeShape in modal coordinate of the model
%       eigenvalues,   double: eigenvalues of the considered system
%
%   OUTPUT:
%       fig, figure: plot of the deformated model
%
%

    % Optional input
    if nargin < 6 || isempty(eigenvalues)
        complex = false;
    else
        complex = true;
    end

    if ~iscell(deformedModel)
        deformedModel = {deformedModel};
    end

    if ~iscell(modeShapes)
        modeShapes = {modeShapes};
    end

    if complex
        [~, iSort] = sort(abs(imag(eigenvalues)));
        eigenvalues = eigenvalues(iSort);
        for iModel = 1:length(modeShapes)
            modeShapes{iModel} = modeShapes{iModel}(:, iSort);
        end
    end

    % Create new figure to overlap the two deformated model
    figName = ['displacements mode ', num2str(startMode)];
    fig = figure(Name=figName);
    tl = tiledlayout(1,length(time));
        tl.TileSpacing = 'none';
        tl.Padding = 'none';

    for j = 1:length(time) 

        ax = nexttile(tl);
        ax.Visible = "on";
        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.ZColor = 'none';
        ax.Clipping = "off";
        view(-50, 30);
        hold(ax, "on")

        % compute dispacements from modeShapes and eventually eigs
        nModels = length(deformedModel);
        for iModel = 1:nModels
            if complex
                displacement = 2 * real(modeShapes{iModel}(:,startMode) * exp(1i*imag(eigenvalues(startMode))*time(j))) * zoomFactor;
            else
                displacement = modeShapes{iModel}(:,startMode)*sin(timeRad(j))*zoomFactor;
            end

            [axUndeformed{iModel}, handleElement{iModel}, handlePanel{iModel}, figUndeformed{iModel}] = fullMeshVisualization(deformedModel{iModel}, []);
            
            % Recover undeformed model plot
            fill3_objs1 = findobj(axUndeformed{iModel}, 'Type', 'patch');

            % Plot the undeformed model
            for i = 1:length(fill3_objs1)
                x1 = get(fill3_objs1(i), 'XData');
                y1 = get(fill3_objs1(i), 'YData');
                z1 = get(fill3_objs1(i), 'ZData');
                c1 = get(fill3_objs1(i), 'FaceColor');
                fill3(ax, x1, y1, z1, c1, 'FaceColor','none','EdgeColor', [.7 .7 .7], 'LineWidth', 0.5);
            end

            [axDeformed{iModel}, ~, ~, figDeformed{iModel}] = fullMeshVisualization(deformedModel{iModel}, displacement);

            % Recover deformed model plot
            fill3_objs2 = findobj(axDeformed{iModel}, 'Type', 'patch');

            % Plot the deformed model
            for i = 1:length(fill3_objs2)
                x2 = get(fill3_objs2(i), 'XData');
                y2 = get(fill3_objs2(i), 'YData');
                z2 = get(fill3_objs2(i), 'ZData');
                c2 = get(fill3_objs2(i), 'FaceColor');
                fill3(ax, x2, y2, z2, c2, 'FaceColor','none','EdgeColor', [0  0  0],'LineWidth', 0.8);
            end
            
            close(figUndeformed{iModel});
            close(figDeformed{iModel});
        end   
    
        % Plot basic reference frame
        frame.orig  = [6; 0; 0];
        frame.xaxis = [1; 0; 0];
        frame.yaxis = [0; 1; 0];
        frame.zaxis = [0; 0; 1];
        plotFrame(ax, frame, 1)
        
        % set FOV of the plot
        xlim([3 9]);    ylim([0 7]);    zlim([-1 3]); 

        % % set title of the current visualization
        % plotTitle = ['Time ', strrep(rats(time(i)),' ',''),'$\pi$'];
        % title(ax,plotTitle);
        
        % close building figures

        % close(figStatic);
    end

    hold off;

end


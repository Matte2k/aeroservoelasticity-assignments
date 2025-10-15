function [handleElement,handlePanel,handleFig] = meshVisualization(model,displacements,handleElement,handlePanel)
%MESH VISUALIZATION
%   Function plotting the aircraft model and giving as output
%   handleElement, handlePanel and handleFig
%
%   SYNTAX:
%       [handleElement, handlePanel, handleFig] = meshVisualization(model, displacements, handleElement, handlePanel)
%
%   INPUT:
%       model,         struct: data organized in following fields:
%                           > model.gridPos: matrix, same order of the gridID,
%                                            three columns for x,y,z
%                           > model.element: struct with label and connectivity
%                           > model.aeroPanel: matrix, rows for the different
%                                              elements, four columns for the two IDS
%                           > model.hinge: vector with the nodes index
%                           > model.masses: vector with the nodes index
%       displacements, double: vector containing displacement info gained
%                              from mode shapes and eventually eigenvalues
%       handleElement and handlePanel: output from previous call of the
%                                      same function (ex. for animation)
%
%   OUTPUT:
%       handleElement,handlePanel: current plot state
%       handleFig,         figure: current figure variable
%
%   OPTIONAL INPUT:
%       handleElement,handlePanel: can be avoided if no animation wanted
%
%

    % Optional input
    if nargin == 2
        handleElement = [];
        handlePanel = [];
    
        % Setup window
        handleFig = figure();
        ax = axes(handleFig);
        ax.Visible = "off";
        ax.Clipping = "off";
        view(30, 30);
        hold(ax, "on")
    
        % Plot basic reference frame
        frame.orig  = [6; 0; 0];
        frame.xaxis = [1; 0; 0];
        frame.yaxis = [0; 1; 0];
        frame.zaxis = [0; 0; 1];
        plotFrame(ax, frame, 1)
        xlim([4.75 9.25]);    ylim([0 5]);    zlim([-1 2]);
    
    else
        ax = [];
    end
    
    
    %%% Apply displacements    
    if ~isempty(displacements)
        model.gridPos = model.gridPos + [displacements(1:6:end), displacements(2:6:end), displacements(3:6:end)];
    end
    
    %%% Plot Model
    % Plot cbar, we use the same function, with a degenerate panel element
    if isfield(model,"elements")
        connectivity = [];
        for i = 1:length(model.elements)
            connectivity = [connectivity;
                model.elements(i).Nodes(1), model.elements(i).Nodes(1),...
                model.elements(i).Nodes(2), model.elements(i).Nodes(2)];
        end
        handleElement = plotElement(ax, model.gridPos, connectivity, handleElement);
    end
    
    % Plot aeroPanel
    if isfield(model,"aeroPanel")
        handlePanel = plotElement(ax, model.gridPos, model.aeroPanel, handlePanel);
    end
    
    % Plot hinge
    if isfield(model,"hinge")
        for iElement = 1:size(model.hinge, 1)
            plotGrid(ax, model.gridPos(model.hinge(iElement)==model.gridID,:), 30, '.g')
        end
    end
    
    % Plot masses
    if isfield(model,"masses")
        for iElement = 1:size(model.masses, 1)
            plotGrid(ax, model.gridPos(model.masses(iElement)==model.gridID,:), 30, 'dr');
        end
    end

end
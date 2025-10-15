function [handleAxes,handleElement,handlePanel,handleFig] = fullMeshVisualization(model,displacements,handleAxes,handleElement,handlePanel)
%FULL MESH VISUALIZATION
%   Function plotting the aircraft model and giving as output
%   handleAxes, handleElement, handlePanel and handleFig
%
%   SYNTAX:
%       [handleAxes,handleElement,handlePanel,handleFig] = fullMeshVisualization(model,displacements,handleAxes,handleElement,handlePanel)
%
%   INPUT:
%       model,         struct: data organized in following fields:
%                               > model.gridPos: matrix, same order of the gridID,
%                                                three columns for x,y,z
%                               > model.element: struct with label and connectivity
%                               > model.aeroPanel: matrix, rows for the different
%                                                  elements, four columns for the two IDS
%                               > model.hinge: vector with the nodes index
%                               > model.masses: vector with the nodes index
%       handleAxes,      axes: plot to be used as base for the current mesh visualization
%       displacements, double: vector containing displacement info gained
%                              from mode shapes and eventually eigenvalues
%       handleElement and handlePanel: output from previous call of the
%                                      same function (ex. for animation)
%
%   OUTPUT:
%       handleAxes           axes: current axes created
%       handleElement,handlePanel: current plot state
%       handleFig,         figure: current figure variable
%
%   OPTIONAL INPUT:
%       handleAxes: can be avoided if no overlap with other mesh visualization is wanted
%       handleElement,handlePanel: can be avoided if no animation wanted
%
%

    % Optional input
    if nargin < 3 || isempty(handleAxes)    
        % Setup window
        handleFig = figure();
        handleAxes = axes(handleFig);
        handleAxes.Visible = "off";
        handleAxes.Clipping = "off";
        view(-40, 30);
        hold(handleAxes, "on")
            
        % Plot basic reference frame
        frame.orig  = [6; 0; 0];
        frame.xaxis = [1; 0; 0];
        frame.yaxis = [0; 1; 0];
        frame.zaxis = [0; 0; 1];
        plotFrame(handleAxes, frame, 1)
        xlim([4 9]);    ylim([0 7]);    zlim([-1 3]);    
    end
    
    if nargin < 4
        handleElement = [];
    end
    
    if nargin < 5
        handlePanel = [];
    end
    
    
    %%% Apply displacements
    if ~isempty(displacements)
        model.gridPos = model.gridPos + [displacements(1:6:end), displacements(2:6:end), displacements(3:6:end)];
    end
    
    %%% Plot Model
    % Plot cbar, we use the same function, with a degenerate panel element
    if isfield(model,"elements")
        connectivity = zeros(length(model.elements), 4);
        for i = 1:length(model.elements)
            if length(model.elements(i).Nodes) == 2
                connectivity(i,:) = [ ...
                    model.elements(i).Nodes(1), model.elements(i).Nodes(1),...
                    model.elements(i).Nodes(2), model.elements(i).Nodes(2)];
            else
                connectivity(i, :) = [ ...
                    model.elements(i).Nodes(1), model.elements(i).Nodes(2),...
                    model.elements(i).Nodes(3), model.elements(i).Nodes(4)];
            end
        end
        handleElement = plotElement(handleAxes, model.gridPos, connectivity, handleElement);
    end
    
    % Plot aeroPanel
    if isfield(model,"aeroPanel")
        handlePanel = plotElement(handleAxes, model.gridPos, model.aeroPanel, handlePanel);
    end
    
    % Plot hinge
    if isfield(model,"hinge")
        for iElement = 1:size(model.hinge, 1)
            plotGrid(handleAxes, model.gridPos(model.hinge(iElement)==model.gridID,:), 30, '.g')
        end
    end
    
    % Plot masses
    if isfield(model,"masses")
        for iElement = 1:size(model.masses, 1)
            plotGrid(handleAxes, model.gridPos(model.masses(iElement)==model.gridID,:), 30, 'dr');
        end
    end

end

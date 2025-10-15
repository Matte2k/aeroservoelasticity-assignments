function [cmap] = graphicSettings()
%GRAPHIC SETTINGS - set grapghics options for the visual output
%   

    % Graphic settings
    set(0,'defaulttextinterpreter','latex');  
    set(0,'defaultTextFontSize', 14);
    set(0,'defaultLegendFontSize', 14);
    set(0,'defaultAxesTickLabelInterpreter','latex');  
    set(0,'defaultLegendInterpreter','latex');
    set(0,'defaultLineLineWidth',1.2);
    warning('off','MATLAB:fplot:NotVectorized');

    % Color map definition
    cmap = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
    
    % Create folder for figure export
    if ~exist('fig','dir')
        mkdir('fig');
    end
    
end
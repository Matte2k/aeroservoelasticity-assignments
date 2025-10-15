function modeAnimation(plotModel,zoomFactor,startMode,modeShapes,eigenvalues)
%MODE ANIMATION
%   Function creating and animation of the deformed model based on
%   modeShapes and eigenvalues of the considered model
%
%   SYNTAX:
%       modeAnimation(plotModel,zoomFactor,iMode,modeShapes,eigenvalues)
%
%   INPUT:
%       plotModel,   struct: model struct containing grid data
%       zoomFactor,  double: amplification factor of the displacement
%       startMode,   double: starting mode to plot in the animation
%       modeShapes,  double: modeShape in modal coordinate of the model
%       eigenvalues, double: eigenvalues of the considered system
%
%   OUTPUT:
%       Only visual output of the animation
%
%

    % Optional input
    if nargin == 4
        complex = false;
    else
        complex = true;
    end
    
    [~, b] = size(modeShapes);
    if complex
        maxMode = b/2;
    else
        maxMode = b;
    end
    
    global stopFlag nextMode prevMode;
    stopFlag = false;
    nextMode = false;
    prevMode = false;
    
    
    % Prepare the figure and store the handles
    [handleElement, handlePanel] = meshVisualization(plotModel, []);
    
    % Create a stop button. Automatically added to the current figure
    stopButton = uicontrol('Style', 'pushbutton', 'String', 'Stop', ...
                           'Units', 'normalized', ...
                           'Position', [0.05 0.05 0.1 0.05], ...
                           'Callback', @(src, event) setStopFlag());
    
    % Create a next mode button
    nextButton = uicontrol('Style', 'pushbutton', 'String', 'Next', ...
                           'Units', 'normalized', ...
                           'Position', [0.20 0.05 0.1 0.05], ...
                           'Callback', @(src, event) setNextMode());
    
    % Create a previous mode button
    prevButton = uicontrol('Style', 'pushbutton', 'String', 'Previous', ...
                           'Units', 'normalized', ...
                           'Position', [0.35 0.05 0.1 0.05], ...
                           'Callback', @(src, event) setPrevMode());
    
    % Time initializtion
    i = 1;
    time = 2*pi/100;
    timeVec(i) = time;
    
    % The increment of displacement
    while ~stopFlag && ~nextMode && ~prevMode
        % compute displacements
        if complex
            displacement = 2 * real(modeShapes(:,startMode) * exp(1i*imag(eigenvalues(startMode))*time/1000)) * zoomFactor;
        else
            displacement = modeShapes(:,startMode)*cos(time)*zoomFactor;
        end
        [handleElement, handlePanel] = meshVisualization(plotModel, displacement, handleElement, handlePanel);
        pause(0.05)
        
        % Increment time for the next iteration
        time = time + 2*pi/100;
        i = i + 1;
        timeVec(i) = time;
    
        % Increment the iMode if required, but not more than the maximum mode
        if nextMode
            if complex
                startMode = startMode+2;
            else
                startMode = startMode+1;
            end
            nextMode = false;
            if startMode > maxMode
                startMode = maxMode;
            end
        % Decrement the iMode, but with minimum 1
        elseif prevMode
            if complex
                startMode = startMode-2;
            else
                startMode = startMode-1;
            end
            prevMode = false;
            if startMode < 1
                startMode = 1;
            end
        end
    end
    
    close(gcf)

end

% We need this because the workspaces are different for the figure and the
% running script
function setStopFlag()
    global stopFlag;
    stopFlag = true;
end

function setNextMode()
    global nextMode;
    nextMode = true;
end

function setPrevMode()
    global prevMode;
    prevMode = true;
end


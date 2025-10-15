function fullAnimateModel(iMode, modeShapes, plotModel, eigenvalues, initialScale, handleAxes)

if isstruct(plotModel)
    plotModel = {plotModel};
end

if ~iscell(modeShapes)
    modeShapes = {modeShapes};
end

if nargin < 4 || isempty(eigenvalues)
    complex = false;
    eigenvalues = [];
else
    complex = true;
end

if nargin < 5
    initialScale = 100;
end

if nargin < 6
    handleAxes = [];
end

if complex
    [~, iSort] = sort(abs(imag(eigenvalues)));
    
    eigenvalues = eigenvalues(iSort);
    
    for iModel = 1:length(modeShapes)
        modeShapes{iModel} = modeShapes{iModel}(:, iSort);
    end
    
end

[~, maxMode] = size(modeShapes{1});


global stopFlag nextMode prevMode zoomFactor pauseTime stepSize;
stopFlag = false;
nextMode = false;
prevMode = false;

zoomFactor = initialScale;
pauseTime = 0.05;
stepSize = 100;

nModels = length(plotModel);

% Prepare the figure and store the handles
handleElement = cell(nModels, 1);
handlePanel = cell(nModels, 1);
%handleAxes = [];
for iModel = 1:nModels
    [handleAxes, handleElement{iModel}, handlePanel{iModel}] = fullMeshVisualization(plotModel{iModel}, [], handleAxes);
end

Xlim = get(handleAxes, 'XLim');
YLim = get(handleAxes, 'YLim');
ZLim = get(handleAxes, 'ZLim');
set(handleAxes, 'XLimMode', 'manual');
set(handleAxes, 'YLimMode', 'manual');
set(handleAxes, 'ZLimMode', 'manual');
set(handleAxes, 'XLim', Xlim);
set(handleAxes, 'YLim', YLim);
set(handleAxes, 'ZLim', ZLim);

% Create a stop button. Automatically added to the current figure
stopButton = uicontrol('Style', 'pushbutton', 'String', 'Stop', ...
                       'Units', 'normalized', ...
                       'Position', [0.05 0.05 0.1 0.05], ...
                       'Callback', @(src, event) setStopFlag());

% Create a next mode button
stopButton = uicontrol('Style', 'pushbutton', 'String', 'Next', ...
                       'Units', 'normalized', ...
                       'Position', [0.20 0.05 0.1 0.05], ...
                       'Callback', @(src, event) setNextMode());

% Create a previous mode button
stopButton = uicontrol('Style', 'pushbutton', 'String', 'Previous', ...
                       'Units', 'normalized', ...
                       'Position', [0.35 0.05 0.1 0.05], ...
                       'Callback', @(src, event) setPrevMode());


hAnn = annotation('textbox', 'position', [0.45 0.05 0.1 0.05], ...);
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'mid', ...
        'LineStyle', 'none', 'String', 'Scale:');
stopButton = uicontrol('Style', 'edit', 'String', num2str(zoomFactor), ...
                       'Units', 'normalized', ...
                       'Position', [0.55 0.05 0.1 0.05], ...
                       'Callback', @(src, event) setScaleFactor(src));
                   

hAnn = annotation('textbox', 'position', [0.70 0.05 0.1 0.05], ...);
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'mid', ...
        'LineStyle', 'none', 'String', 'Steps/rev:');
stopButton = uicontrol('Style', 'edit', 'String', num2str(stepSize), ...
                       'Units', 'normalized', ...
                       'Position', [0.80 0.05 0.1 0.05], ...
                       'Callback', @(src, event) setStepSize(src));

                   
% The increment of displacement
step = 1/stepSize;

hAnn = [];

while ~stopFlag && ~nextMode && ~prevMode

    if ~isempty(eigenvalues)
        hAnn = putAnnotation(hAnn, iMode, eigenvalues);
    end
    
    for iModel = 1:nModels
        if complex
            displacement = 2*real( modeShapes{iModel}(:,iMode) * exp(1i*imag(eigenvalues(iMode))*step) )*zoomFactor;
        else
            displacement = modeShapes{iModel}(:,iMode)*cos(2*pi * step)*zoomFactor;
        end
        [handleAxes, handleElement{iModel}, handlePanel{iModel}] = fullMeshVisualization(plotModel{iModel}, displacement, handleAxes, handleElement{iModel}, handlePanel{iModel});
    end
    
    pause(0.05)
    step = step + 1/stepSize;

    % Increment the iMode if required, but not more than the maximum mode
    if nextMode
        if complex
            iMode = iMode+2;
        else
            iMode = iMode+1;
        end
        nextMode = false;
        if iMode > maxMode
            if complex
                iMode = maxMode - 1;
            else
                iMode = maxMode;
            end
        end
    % Decrement the iMode, but with minimum 1
    elseif prevMode
        if complex
            iMode = iMode-2;
        else
            iMode = iMode-1;
        end
        prevMode = false;
        if iMode < 1
            iMode = 1;
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

function setScaleFactor(src)
    global zoomFactor;
    zoomFactor = str2double(src.String);
end

function setStepSize(src)
    global stepSize;
    stepSize = str2double(src.String);
end

function hAnn = putAnnotation(hAnn, iMode, eigenvalues)

if isempty(hAnn)
    
    hAnn = annotation('textbox', 'position', [0, 0.8, 0.4, 0.2]);
    set(hAnn, 'HorizontalAlignment', 'left');
    set(hAnn, 'VerticalAlignment', 'top');
    set(hAnn, 'LineStyle', 'none');
%     set(hAnn, 'Color', 'none');
end

lambda = eigenvalues(iMode);
freq = imag(lambda)/2/pi;
damp = -real(lambda) / (abs(lambda) + eps());

textstring = sprintf('Mode %d: %1.3e + %1.3ei\nFreq = %1.3fHz  xi = %1.3f%%', ...
    iMode, real(lambda), imag(lambda), freq, damp);

set(hAnn, 'string', textstring);

end

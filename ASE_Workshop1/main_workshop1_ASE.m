clearvars;  close all; clc
addpath('src\');
graphicSettings();

% load tiltrotor Leonardo's data
data = load("LabData_Session1.mat");


%% INPUT

% plot flag settings
plotFlag.aircraft_model  = false;
plotFlag.original_modes  = false;
plotFlag.actStiff_modes  = false;
plotFlag.actDamp_modes   = false;
plotFlag.hingeSens       = false;
plotFlag.flutter_std     = false;
plotFlag.flutter_cont    = true;
plotFlag.flutter_sens    = true;

% plot export settings
saveFigure = false;          % export figure in .eps
time2plot = [1/8 1/4 1/2];  % fraction of pi
mode2plot = [1 2 3 4];      % mode number to plot

% animation flag settings
animFlag.original_modes = false;
animFlag.actStiff_modes = false;
animFlag.actDamp_modes  = false;
idxAnimMode = 1;            % first mode to visualize in the animation

% displacement amplification
ampFactor.original_modes = 10;
ampFactor.actStiff_modes = 10;
ampFactor.actDamp_modes  = 100;

% wing model settings
wing.sweep = deg2rad(-6.5);         % wing sweep angle    
wing.dihed = deg2rad(2);            % wing dihedral angle
wing.chord = 1.6;                   % wing chord [m]

% hinge model settings
hinge.stiff = 1e4;                  % actuator stiffness
hinge.dampingCase = 1;              % actuator damping selector
hinge.damp  = 200;                  % fixed physical damping    [dampingCase = 1]
hinge.alpha = 1e-2;                 % mass prop. damping        [dampingCase = 2]
hinge.beta  = 1e-4;                 % stiffness prop. damping   [dampingCase = 2]
hinge.modal = 0.03;                 % modal damping value xi

% actuator stiffness sensitivity settings
hingeSens.stiffVec  = logspace(0,6,1e3);    % stiffness vector to analyze
hingeSens.mode2plot = 1:5;                  % mode number to analyze

% flutter analysis flow model settings
standardFlutter   = true;           % compute flutter with 'eig'
flutter.v_limit   = [51.4, 128.6];  % limit values of flow velocity
flutter.rho       = 1.225;          % flow density for flutter analysis
flutter.dist.type = 'lin';          % type of speed values discretization
flutter.dist.elem = 200;            % number of element in the speed vector

% complete system flutter sensitivity settings
%flutterSens.stiffVec  = (0:500:25000);     % rivedere questi range
%flutterSens.dampVec   = (0:50:2000);       % rivedere questi range 
flutterSens.stiffVec  = (0:1000:25000);     % JUST FOR DEBUG
flutterSens.dampVec   = (0:250:2000);       % JUST FOR DEBUG
flutterSens.rho       = 1.225;              % flow density for flutter analysis
flutterSens.v_limit   = [40, 200];          % limit values of flow velocity
flutterSens.dist.type = 'lin';              % flow density for flutter analysis
flutterSens.dist.elem = 1e3;                % number of element in the speed vector
flutterSens.nWorker   = 6;                  % number of core to use in parallel pool

% continuation approach settings
continuationFlutter = true;         % compute flutter with 'eig'
contOpt.modes2solve = 1:2:10;       % mode to compute
contOpt.maxIter = 200;              % maximum iteration of Newton method
contOpt.toll    = 1e-9;             % tollerance of Newton method


%% TASK 1 - Aircraft plot

% Initialize MODEL data struct
model = struct;
    model.gridPos  = data.nodePos;
    model.elements = data.elements;

% Plot aircraft undeformed configuration
if plotFlag.aircraft_model == true
    [~,~,fig.aircraft] = meshVisualization(model, []);
end


%% TASK 2 - Structure deformation

% Initialize ORIGINAL SYSTEM data struct
original_system = struct;
    original_system.type = 'dynamic';
    original_system.M = data.Mhh;
    original_system.K = data.Khh;

% Compute EIGENANALYSIS of the original system
[original_system] = sysEigenanalysis(original_system);

% Reshape the EIGENMODE of the original system in STANDARD FORM
original_system.modeShape = zeros(size(data.Phi,1)*size(data.Phi,2),size(data.Phi,3));
for k = 1:size(data.Phi,3)
     original_system.modeShape(:,k) = reshape(data.Phi(:,:,k)',[size(data.Phi,1)*size(data.Phi,2),1]);
end

% Animate aircraft deformate configuration
if animFlag.original_modes == true
    modeAnimation(model, ampFactor.original_modes, idxAnimMode, original_system.modeShape);
end

% Plot aircraft deformate configuration
if plotFlag.original_modes == true
    fig.original_modes = cell(1,length(mode2plot));
    for k = 1:length(mode2plot)
        fig.original_modes{k} = modePlot(model, time2plot, ampFactor.original_modes, mode2plot(k), original_system.modeShape);
    end
end


%% TASK 3 - Flap hinge stiffness addition

% Initialize ACTUATOR data struct
actuator = struct;
    actuator.nodes = [1205, 2205];

% Compute the ADDED STIFFNESS introduced by the actuator
[actuator.K] = addActuatorTerm(hinge.stiff, wing.sweep, wing.dihed, data.Phi, data.nodeLabel, actuator.nodes);

% Build the NEW SYSTEM struct considering the added stiffness
actStiff_system = struct;
    actStiff_system.type = 'dynamic';
    actStiff_system.M = data.Mhh;
    actStiff_system.K = data.Khh + actuator.K.modal;
 
% Compute EIGENANALYSIS of the system with actuator stiffness
[actStiff_system] = sysEigenanalysis(actStiff_system);

% Compute NEW MODAL SHAPE of the system with actuator stiffness 
actStiff_system.modeShape = original_system.modeShape * actStiff_system.eigsVect_sorted;

% Animate aircraft deformate configuration
if animFlag.actStiff_modes
    modeAnimation(model, ampFactor.actStiff_modes, idxAnimMode, actStiff_system.modeShape);
end

% Plot aircraft deformate configuration
if plotFlag.actStiff_modes == true
    fig.actStiff_modes = cell(1,length(mode2plot));
    for k = 1:length(mode2plot)
        fig.actStiff_modes{k} = modePlot(model, time2plot, ampFactor.actStiff_modes, mode2plot(k), actStiff_system.modeShape);
    end
end


%% TASK 4 - Sensitivity analysis on flap hinge stiffenss

% Initialize ACTUATOR STIFFNESS SENSITIVITY study struct
actStiff_sensitivity = struct;
    actStiff_sensitivity.stiffVec  = hingeSens.stiffVec;
    actStiff_sensitivity.Kmodal    = cell(1,length(actStiff_sensitivity.stiffVec));
    actStiff_sensitivity.eigsVal   = zeros(length(data.Mhh),length(actStiff_sensitivity.stiffVec));
    actStiff_sensitivity.eigsFreq  = zeros(length(data.Mhh),length(actStiff_sensitivity.stiffVec));

% Perform SENSITIVITY STUDY on the system eigenfrequencies
for k = 1:length(actStiff_sensitivity.stiffVec)
    % compute the Khh matrix for the current iteration
    [deltaK] = addActuatorTerm(actStiff_sensitivity.stiffVec(k), wing.sweep, wing.dihed, data.Phi, data.nodeLabel, actuator.nodes);
    actStiff_sensitivity.Kmodal{k} = data.Khh + deltaK.modal;
    
    % compute and sort eigenfrequencies
    actStiff_sensitivity.eigsVal(:,k)  = eig(actStiff_sensitivity.Kmodal{k}, data.Mhh);
    actStiff_sensitivity.eigsFreq(:,k) = sort(sqrt(actStiff_sensitivity.eigsVal(:,k))./(2*pi),'ascend');
end
clear("deltaK");    % clear the buffer variable used in the cylcle

% Plot RESULTS of the sensitivity analysis on the stiffness
if plotFlag.hingeSens == true
    legendCell = cell(1,length(hingeSens.mode2plot));

    fig.hingeSens = figure(Name='Stiffness sensitivity');
    for k = 1:length(hingeSens.mode2plot)
        semilogx(actStiff_sensitivity.stiffVec, actStiff_sensitivity.eigsFreq(k,:));
        hold on;
        legendCell{k} = ['Mode ', num2str(hingeSens.mode2plot(k))];
    end
    grid minor;  axis padded;  box on;
    xlabel('Stiffness [Nm/rad]');   ylabel('Frequency [Hz]')
    legend(legendCell,Location='northoutside',Orientation='horizontal');
end


%% TASK 5 - Flap actuator damping addition

% Compute the ADDED STIFFNESS introduced by the actuator
    % CASE 1: physical damping addition
if hinge.dampingCase == 1   
    [actuator.C] = addActuatorTerm(hinge.damp, wing.sweep, wing.dihed, data.Phi, data.nodeLabel, actuator.nodes);

    % CASE 2: proportional damping addition
elseif hinge.dampingCase == 2
    actuator.C.modal = hinge.alpha*data.Mhh + hinge.beta*actStiff_system.K;

    % CASE 3: no damping addition
elseif hinge.dampingCase == 3
    actuator.C.modal = zeros(size(data.Mhh));

    % Error in case of invalid damping case value
else
    error('select a valid damping case');
end

% Build the NEW SYSTEM struct considering the added damping
actDamp_system = struct;
    actDamp_system.type = 'state_space';
    actDamp_system.M = data.Mhh;
    actDamp_system.K = actStiff_system.K;
    actDamp_system.C = actuator.C.modal;

% Build the STATE SPACE form of the new system
Minv = inv(actDamp_system.M);
actDamp_system.A = [  zeros(length(actDamp_system.M)),    eye(length(actDamp_system.M)); 
                      -Minv*actDamp_system.K,             -Minv*actDamp_system.C];   
clear('Minv');  % clear the buffer variable of the inverse mass matrix

% Compute EIGENANALYSIS of the system with actuator damping
[actDamp_system] = sysEigenanalysis(actDamp_system);

% Compute MODAL DECOMPOSITION of the dynamic system with actuator damping
actDamp_system.eigsVect_dynSystem = real(actDamp_system.eigsVect(length(actDamp_system.M)+1:end, 1:2:end));
actDamp_system.M_modal = actDamp_system.eigsVect_dynSystem' * actDamp_system.M  * actDamp_system.eigsVect_dynSystem;
actDamp_system.K_modal = actDamp_system.eigsVect_dynSystem' * actDamp_system.K  * actDamp_system.eigsVect_dynSystem;
actDamp_system.C_modal = actDamp_system.eigsVect_dynSystem' * actDamp_system.C  * actDamp_system.eigsVect_dynSystem; 

% Compute NEW MODAL SHAPE of the system with actuator stiffness 
actDamp_system.modeShape = original_system.modeShape * actDamp_system.eigsVect_sorted(1:length(actDamp_system.M),:);
  
% Animate aircraft deformate configuration
if animFlag.actDamp_modes == true
    modeAnimation(model, ampFactor.actDamp_modes, idxAnimMode, actDamp_system.modeShape, actDamp_system.eigsVal);
end

% Plot aircraft deformate configuration
if plotFlag.actDamp_modes == true
    mode2plot_ss  = mode2plot.*2 - 1;
    fig.actDamp_modes = cell(1,length(mode2plot));
    for k = 1:length(mode2plot)
        fig.actDamp_modes{k} = modePlot(model, time2plot, ampFactor.actDamp_modes, mode2plot_ss(k), actDamp_system.modeShape, actDamp_system.eigsVal);
    end
end


%% TASK 6 - Flutter analysis of the aircraft 

% Build the NEW SYSTEM struct considering the actuator free
flutter_system = struct;
    flutter_system.type = 'state_space';
    flutter_system.M = data.Mhh;
    flutter_system.K = data.Khh;
    flutter_system.C = zeros(length(data.Mhh));

% Build model of the FLOW starting from user input
flutter_system.flow = buildFlow(flutter.v_limit, flutter.dist, flutter.rho);

% Flutter analysis with STANDARD approach
if standardFlutter == true
    [flutter_system.standard, fig.flutterStd] = flutterStandard(flutter_system, data, wing.chord, plotFlag.flutter_std);
end

% Flutter analysis with CONTINUATION approach
if continuationFlutter == true
    [flutter_system.continuation, fig.flutterCont] = flutterContinuation(flutter_system, data, wing.chord, contOpt, plotFlag.flutter_cont);
end


%% TASK 7 - Parametric flutter analysis on actuator stiffness and damping

% Parallel flutter analysis
[flutter_sensitivity_speed, fig.flutterSens] = parFlutterSensitivity(flutterSens,wing,data,plotFlag.flutter_sens);


%% Construction of the modal damping 

% Build the NEW SYSTEM struct considering modal damping from user input
modalDamp_system = struct;
         modalDamp_system.type = 'state_space';
         modalDamp_system.M = data.Mhh;
         modalDamp_system.K = data.Khh;
         modalDamp_system.C = diag(2*hinge.modal.*diag(data.Mhh).*sqrt(diag(data.Khh)));

% Build model of the FLOW starting from user input
modalDamp_system.flow = buildFlow(flutter.v_limit, flutter.dist, flutter.rho);

% Flutter analysis with STANDARD approach
if standardFlutter == true
    [modalDamp_system.standard, fig.flutterStd] = flutterStandard(modalDamp_system, data, wing.chord, plotFlag.flutter_std);
end

% Flutter analysis with CONTINUATION approach
if continuationFlutter == true
    [modalDamp_system.continuation, fig.flutterModalCont] = flutterContinuation(modalDamp_system, data, wing.chord, contOpt, plotFlag.flutter_cont);
end


%% EXPORT FIGURES

% Export figures in eps format
if saveFigure == true
    exportFigure(fig,plotFlag);
end
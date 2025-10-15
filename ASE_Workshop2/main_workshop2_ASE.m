clearvars;  close all;  clc;
addpath("src/");
graphicSettings;

% Load data
airfrData = load('LabData_Session1.mat');
rotorData = load('LabData_Session2_rotor.mat');
speedData = load('LabData_Session2_rotor_allSpeed.mat');

% Animation flag
animFlag.freeRotor            = false;
animFlag.groundedRotor        = false;
animFlag.tiltRotor.airframe   = false;
animFlag.tiltRotor.rotor      = false;
animFlag.tiltRotor.total      = false;
animFlag.gimbalXV15.airframe  = false;
animFlag.gimbalXV15.rotor     = false;
animFlag.gimbalXV15.total     = false;
animFlag.reducedXV15.airframe = false;
animFlag.reducedXV15.rotor    = false;
animFlag.reducedXV15.total    = false;

% Rootlocus plot flag
plotFlag.rootLocus.airframe     = false;
plotFlag.rootLocus.grdRotor     = true;     %
plotFlag.rootLocus.tiltrotor    = true;     %
plotFlag.rootLocus.completeXV15 = true;     %    
plotFlag.rootLocus.gimbalXV15   = true;     %    
plotFlag.rootLocus.reducedXV15  = false;

% Damping ratio plot flag
plotFlag.dampRatio.completeXV15 = true;     %    
plotFlag.dampRatio.gimbalXV15   = true;     %    
plotFlag.dampRatio.reducedXV15  = false;

% Models plot flag
plotFlag.mesh.fullModel    = true;     %        
plotFlag.modes.grdRotor    = true;     %        
plotFlag.modes.tiltRotor   = true;     %         
plotFlag.modes.gimbalXV15  = true;     %         
plotFlag.modes.reducedXV15 = false;

% Plot settings
grdRotor.ampFactor   = [2.5 5.0 10 15];
grdRotor.time2plot   = [1.0 1.5 2.0];
tiltRotor.ampFactor  = [10 100 125 150];
tiltRotor.time2plot  = [1.0 1.5 2.0];
gimbalXV15.ampFactor = [50 80 80 80];
gimbalXV15.time2plot = [1.0 1.5 2.0];
mode2plot = [1 3 5 7];


%% TASK 0 : load model from session 1 and apply actuator stiffness

% Parametric info
sweep    = -6.5*pi/180;
dihedral = 2*pi/180;
damp  = 20;
stiff = 9800;

% Rotation matrix
RS = [cos(sweep) -sin(sweep) 0; sin(sweep) cos(sweep) 0; 0 0 1]';
RD = [1 0 0; 0 cos(dihedral) -sin(dihedral); 0 sin(dihedral) cos(dihedral)]';
Y  = [0 1 0];

% Load airframe system size
airframe.nModes = size(airfrData.Khh, 1);
airframe.nNodes = size(airfrData.nodePos, 1);
     
% Extract the modal shapes, building the modal matrices
airframe.modeShape = zeros(6*airframe.nNodes, airframe.nModes);
for iMode = 1:airframe.nModes
    airframe.modeShape(:, iMode) = reshape(airfrData.Phi(:,:,iMode)', [airframe.nNodes*6,1]);
end

% Matrices to extract the position, in the modal shapes, of the rotations
airframe.psi1 = zeros(3,airframe.nNodes*6);
    airframe.psi1(1,(find(airfrData.nodeLabel==1205)-1)*6+4) = 1;
    airframe.psi1(2,(find(airfrData.nodeLabel==1205)-1)*6+5) = 1;
    airframe.psi1(3,(find(airfrData.nodeLabel==1205)-1)*6+6) = 1;
airframe.psi2 = zeros(3,airframe.nNodes*6);
    airframe.psi2(1,(find(airfrData.nodeLabel==2205)-1)*6+4) = 1;
    airframe.psi2(2,(find(airfrData.nodeLabel==2205)-1)*6+5) = 1;
    airframe.psi2(3,(find(airfrData.nodeLabel==2205)-1)*6+6) = 1;

% Modal shapes at the points of interest, in the local direction of rotation
S = [Y*RD*RS*airframe.psi1*airframe.modeShape; Y*RD*RS*airframe.psi2*airframe.modeShape];

% Dynamic system AIRFRAME definition
airframe.M = airfrData.Mhh;
airframe.C = damp*S'*[1 -1]'*[1 -1]*S;
airframe.K = airfrData.Khh + stiff*S'*[1 -1]'*[1 -1]*S;

% State space system of AIRFRAME system
A = [zeros(airframe.nModes, airframe.nModes),   eye(airframe.nModes, airframe.nModes);
    -airframe.M\airframe.K,                     -airframe.M\airframe.C];

% Eigenanalysis of AIRFRAME system
[~, airframe.eigsVal] = eig(A);
airframe.eigsVal = diag(airframe.eigsVal);
airframe.eigsVal = airframe.eigsVal(imag(airframe.eigsVal)>=0);

% Plot AIRFRAME system root locus
if plotFlag.rootLocus.airframe == true
    fig.rootLocus.airframe = figure(Name='Airframe root-locus'); 
    hold on;  grid minor;  axis padded;  box on;
    axCompl = gca();
    hPlot = plot(axCompl, real(airframe.eigsVal)/2/pi, imag(airframe.eigsVal)/2/pi, 'o', 'displayname', 'Airframe modes');
    set(hPlot, 'markerfacecolor', get(hPlot, 'color'));
    set(hPlot, 'markersize', 6);
    xlabel('Re($\lambda/{2\pi}$)');     xlim([-3, 0.5]);    yline(0,'--')
    ylabel('Im($\lambda/{2\pi}$)');     ylim([-1, 7]);      xline(0,'--')
end


%% TASK 1 - Load the complete model

% Plot rotor and airframe complete system
rotorModel = struct;
    rotorModel.gridPos  = rotorData.nodePos;
    rotorModel.elements = rotorData.elements;

airfrModel = struct;
    airfrModel.gridPos  = airfrData.nodePos;
    airfrModel.elements = airfrData.elements;
    
if plotFlag.mesh.fullModel == true
    [airfrAxes,~,~,fig.mesh.fullModel] = fullMeshVisualization(airfrModel, []);    
    fullMeshVisualization(rotorModel, [], airfrAxes);
end

% Reshape the EIGENMODE of the original system in STANDARD FORM
rotorData.modeShape = zeros(size(rotorData.Phi,1)*size(rotorData.Phi,2),size(rotorData.Phi,3));
for i = 1:size(rotorData.Phi,3)
     rotorData.modeShape(:,i) = reshape(rotorData.Phi(:,:,i)',[size(rotorData.Phi,1)*size(rotorData.Phi,2),1]);
end

% Rotor model mode animation
if animFlag.freeRotor == true
    fullAnimateModel(1, rotorData.modeShape, rotorModel);
end


%% TASK 2 - Characterize grounded rotor system

% Find hub and swashplate nodes in GROUNDED ROTOR system
grdRotor.hubNodeRow   = find(rotorData.nodeLabel == 1   );
grdRotor.swashNodeRow = find(rotorData.nodeLabel == 9500);

% Build constraint matrix to compute GROUNDED ROTOR system
grdRotor.phiR(1:6 ,:) = squeeze(rotorData.Phi(grdRotor.hubNodeRow,:,:));
grdRotor.phiR(7:12,:) = squeeze(rotorData.Phi(grdRotor.swashNodeRow,:,:));

grdRotor.phiR1 = grdRotor.phiR(:, 1:12);
grdRotor.phiR2 = grdRotor.phiR(:, 13:end);

grdRotor.constrPhi = -inv(grdRotor.phiR1)*grdRotor.phiR2;
grdRotor.constrMat = [grdRotor.constrPhi; eye(17)];     % T matrix

% Dynamic system GROUNDED ROTOR definition
grdRotor.M = grdRotor.constrMat' * rotorData.M * grdRotor.constrMat;
grdRotor.C = grdRotor.constrMat' * rotorData.C * grdRotor.constrMat;
grdRotor.K = grdRotor.constrMat' * rotorData.K * grdRotor.constrMat;

% State space system of GROUNDED ROTOR system
grdRotor.Minv = inv(grdRotor.M);
grdRotor.A = [  zeros(length(grdRotor.M)),    eye(length(grdRotor.M)); 
                -grdRotor.Minv*grdRotor.K,   -grdRotor.Minv*grdRotor.C ]; 

% Eigenanalysis of GROUNDED ROTOR system
[grdRotor.eigsVec,grdRotor.eigsVal] = eig(grdRotor.A);
grdRotor.eigsVal = diag(grdRotor.eigsVal);

% Plot GROUNDED ROTOR system root locus
if plotFlag.rootLocus.grdRotor == true
    fig.rootLocus.grdRotor = figure(Name='Rotor root-locus');
    hold on;  grid minor;  axis padded;  box on;
    axCompl = gca();
    hPlot = plot(axCompl, real(grdRotor.eigsVal(1:2:end))/2/pi, imag(grdRotor.eigsVal(1:2:end))/2/pi, 'o', 'displayname', 'Rotor modes');
    set(hPlot, 'markerfacecolor', get(hPlot, 'color'));
    set(hPlot, 'markersize', 6);
    xlabel('Re($\lambda/{2\pi}$)');     xlim([-3, 0.5]);    yline(0,'--')
    ylabel('Im($\lambda/{2\pi}$)');     ylim([-1, 7]);      xline(0,'--')
end

% Compute GROUNDED ROTOR modal shape
grdRotor.modeShape = (rotorData.modeShape * grdRotor.constrMat) * grdRotor.eigsVec(1:length(grdRotor.M),:);

% Animation of the GROUNDED ROTOR system
if animFlag.groundedRotor == true
    fullAnimateModel(1, grdRotor.modeShape, rotorModel, grdRotor.eigsVal);
end

% Plot GROUNDED ROTOR deformate configuration
if plotFlag.modes.grdRotor == true
    fig.grdRotor = cell(1,length(mode2plot));
    for k = 1:length(mode2plot)
        fig.grdRotor{k} = fullModePlot(rotorModel, grdRotor.time2plot, grdRotor.ampFactor(k), mode2plot(k), grdRotor.modeShape, grdRotor.eigsVal);
    end
end


%% TASK 3 - Characterize complete model system

% Find hub and swashplate nodes in TILT ROTOR system
tiltRotor.hubNodeRow   = find(airfrData.nodeLabel == 1000007);
tiltRotor.swashNodeRow = find(airfrData.nodeLabel == 100000);

% Build constraint matrix to compute TILT ROTOR system
tiltRotor.phiA(1:6 ,:) = squeeze(airfrData.Phi(tiltRotor.hubNodeRow,:,:));
tiltRotor.phiA(7:12,:) = squeeze(airfrData.Phi(tiltRotor.swashNodeRow,:,:));

tiltRotor.constrMat = [ eye(length(airframe.M)),               zeros(length(airframe.M),17);
                        inv(grdRotor.phiR1)*tiltRotor.phiA,    grdRotor.constrPhi;
                        zeros(17,length(airframe.M)),          eye(17,17)];

% Define MIXED dynamic system matrix
Mmix = zeros(43);
Mmix(1:14,1:14) = airframe.M;
Mmix(15:end,15:end) = rotorData.M;

Cmix = zeros(43);
Cmix(1:14,1:14) = airframe.C;
Cmix(15:end,15:end) = rotorData.C;

Kmix = zeros(43);
Kmix(1:14,1:14) = airframe.K;
Kmix(15:end,15:end) = rotorData.K;

% Dynamic system TILT ROTOR definition
tiltRotor.M = tiltRotor.constrMat' * Mmix * tiltRotor.constrMat;
tiltRotor.C = tiltRotor.constrMat' * Cmix * tiltRotor.constrMat;
tiltRotor.K = tiltRotor.constrMat' * Kmix * tiltRotor.constrMat;

% State space system of TILT ROTOR system
tiltRotor.Minv = inv(tiltRotor.M);
tiltRotor.A = [ zeros(length(tiltRotor.M)),     eye(length(tiltRotor.M)); 
                -tiltRotor.Minv*tiltRotor.K,    -tiltRotor.Minv*tiltRotor.C];   % A

% Eigenanalysis of TILT ROTOR system
[tiltRotor.eigsVec,tiltRotor.eigsVal] = eig(tiltRotor.A);
tiltRotor.eigsVal = diag(tiltRotor.eigsVal);

% Plot TILT ROTOR system root locus
if plotFlag.rootLocus.tiltrotor == true
    fig.rootLocus.tiltrotor = figure(Name='Tiltrotor root-locus');
    hold on;  grid minor;  axis padded;  box on;
    axCompl = gca();
    hPlot = plot(axCompl, real(tiltRotor.eigsVal(1:2:end))/2/pi, imag(tiltRotor.eigsVal(1:2:end))/2/pi, 'o', 'displayname', 'Rotor modes');
    set(hPlot, 'markerfacecolor', get(hPlot, 'color'));
    set(hPlot, 'markersize', 6);
    xlabel('Re($\lambda/{2\pi}$)');     xlim([-3, 0.5]);    yline(0,'--')
    ylabel('Im($\lambda/{2\pi}$)');     ylim([-1, 7]);      xline(0,'--')
end

% Compute modal shape for TILT ROTOR system
tiltRotor.modeShape_airframe = (airframe.modeShape * tiltRotor.constrMat(1:14,:)) ...
                                * tiltRotor.eigsVec(1:length(tiltRotor.eigsVec)/2,:);

tiltRotor.modeShape_rotor    = (rotorData.modeShape * tiltRotor.constrMat(15:end,:)) ...
                                * tiltRotor.eigsVec(1:length(tiltRotor.eigsVec)/2,:);

% Animation of the AIRFRAME of the TILT ROTOR system
if animFlag.tiltRotor.airframe == true
    rotorAxes = fullMeshVisualization(rotorModel, []);
    fullAnimateModel(1, tiltRotor.modeShape_airframe, airfrModel, tiltRotor.eigsVal, 1, rotorAxes);
end

% Animation of the ROTOR of the TILT ROTOR system
if animFlag.tiltRotor.rotor == true
    airfrAxes = fullMeshVisualization(airfrModel, []);
    fullAnimateModel(1, tiltRotor.modeShape_rotor, rotorModel, tiltRotor.eigsVal, 1, airfrAxes);
end

% Animation of the TOTAL of the TILT ROTOR system
if animFlag.tiltRotor.total == true
    fullAnimateModel(1, {tiltRotor.modeShape_airframe,tiltRotor.modeShape_rotor}, {airfrModel,rotorModel}, tiltRotor.eigsVal, 1);
end

% Plot TOTAL of the TILT ROTOR deformate configuration
if plotFlag.modes.tiltRotor == true
    fig.tiltRotor = cell(1,length(mode2plot));
    for k = 1:length(mode2plot)
        fig.tiltRotor{k} = fullModePlot({airfrModel,rotorModel}, tiltRotor.time2plot, tiltRotor.ampFactor(k), mode2plot(k), {tiltRotor.modeShape_airframe,tiltRotor.modeShape_rotor}, tiltRotor.eigsVal);
    end
end


%% TASK 4 - Study the whirl flutter of the system

% Reshape modal shapes for the rotor at different speed
speedData.modeShape = zeros( size(speedData.Phi,1)*size(speedData.Phi,2), ...
                             size(speedData.Phi,3), size(speedData.Phi,4) );
for j = 1:size(speedData.Phi,4)
    for i = 1:size(speedData.Phi,3)
         speedData.modeShape(:,i,j) = reshape( speedData.Phi(:,:,i,j)', ...
                                               [size(speedData.Phi,1)*size(speedData.Phi,2),1] );
    end
end

% Define colormap for the following root-locus plot
cmap = winter(size(speedData.Phi,4));

%%% Complete XV15 real model whirl flutter analysis
completeXV15.constrMat = tiltRotor.constrMat;
[completeXV15] = computeWhirlflutter(airframe, speedData, completeXV15);

if any(completeXV15.unstable ~= 0)
    flutterIdx = find(completeXV15.unstable ~= 0, 1);

    completeXV15.modeShapeFlutter_airframe = (airframe.modeShape * completeXV15.constrMat(1:14,:)) ...
                                             * completeXV15.eigsVec(1:length(completeXV15.eigsVec)/2,:,flutterIdx);

    completeXV15.modeShapeFlutter_rotor    = (speedData.modeShape(:,:,flutterIdx) * completeXV15.constrMat(15:end,:)) ... 
                                             * completeXV15.eigsVec(1:length(completeXV15.eigsVec)/2,:,flutterIdx);
end

% Plot COMPLETE XV15 system root locus
if plotFlag.rootLocus.completeXV15 == true
    fig.rootLocus.completeXV15 = figure(Name='complete XV15 root-locus');
    hold on;  grid minor;  axis padded;  box on;
    axCompl = gca();
    for i = 1:size(completeXV15.eigsVal,2)
        hPlot = plot(axCompl, completeXV15.eigsVal(1:2:end,i)./(2*pi), 'o', 'displayname', 'complete XV15');
        set(hPlot, 'MarkerEdgeColor', cmap(i,:));
        set(hPlot, 'MarkerFaceColor', cmap(i,:));
        set(hPlot, 'MarkerSize', 6);
    end
    colormap(cmap);
    cbar = colorbar;
    cbar.Label.String = 'Airspeed [m/s]';
    cbar.Label.Interpreter = 'latex';
    cbar.Ticks = speedData.V;
    cbar.TickLabels = string(round(speedData.V,2));
    cbar.TickLabelInterpreter = 'latex';
    cbar.TicksMode  = 'auto';
    xlabel('Re($\lambda/{2\pi}$)');     xlim([-3, 0.5]);    yline(0,'--')
    ylabel('Im($\lambda/{2\pi}$)');     ylim([-1, 7]);      xline(0,'--')
end


%%% NO GIMBAL XV15 model whirl flutter analysis
gimbalXV15.constrMat = tiltRotor.constrMat(:,[1:11,17:end]);
[gimbalXV15] = computeWhirlflutter(airframe, speedData, gimbalXV15);

if any(gimbalXV15.unstable ~= 0)
    flutterIdx = find(gimbalXV15.unstable ~= 0, 1);
    
    gimbalXV15.modeShapeFlutter_airframe = (airframe.modeShape * gimbalXV15.constrMat(1:14,:)) ...
                                           * gimbalXV15.eigsVec(1:length(gimbalXV15.eigsVec)/2,:,flutterIdx);
    
    gimbalXV15.modeShapeFlutter_rotor    = (speedData.modeShape(:,:,flutterIdx) * gimbalXV15.constrMat(15:end,:)) ... 
                                           * gimbalXV15.eigsVec(1:length(gimbalXV15.eigsVec)/2,:,flutterIdx);
    
    % Animation of the AIRFRAME of the NO GIMBAL XV15 system at whirl flutter speed
    if animFlag.gimbalXV15.airframe == true
        rotorAxes = fullMeshVisualization(rotorModel, []);
        fullAnimateModel(1, gimbalXV15.modeShapeFlutter_airframe, airfrModel, gimbalXV15.eigsVal(:,flutterIdx));
    end

    % Animation of the ROTOR of the NO GIMBAL XV15 system at whirl flutter speed
    if animFlag.gimbalXV15.rotor == true
        airfrAxes = fullMeshVisualization(airfrModel, []);
        fullAnimateModel(1, gimbalXV15.modeShapeFlutter_rotor, rotorModel, gimbalXV15.eigsVal(:,flutterIdx), 1, airfrAxes);
    end
   
    % Animation of the TOTAL of the NO GIMBAL XV15 system
    if animFlag.gimbalXV15.total == true
        fullAnimateModel(1, {gimbalXV15.modeShape_airframe,gimbalXV15.modeShape_rotor}, {airfrModel,rotorModel}, gimbalXV15.eigsVal, 1);
    end

    % Plot TOTAL of the NO GIMBAL XV15 deformate configuration
    if plotFlag.modes.gimbalXV15 == true
        fig.gimbalXV15 = cell(1,length(mode2plot));
        for k = 1:length(mode2plot)
            fig.gimbalXV15{k} = fullModePlot({airfrModel,rotorModel}, gimbalXV15.time2plot, gimbalXV15.ampFactor(k), mode2plot(k), {gimbalXV15.modeShapeFlutter_airframe,gimbalXV15.modeShapeFlutter_rotor}, gimbalXV15.eigsVal);
        end
    end

end

% Plot NO GIMBAL XV15 system root locus
if plotFlag.rootLocus.gimbalXV15 == true
    fig.rootLocus.gimbalXV15 = figure(Name='gimbal XV15 root-locus');
    hold on;  grid minor;  axis padded;  box on;
    axCompl = gca();
    for i = 1:size(gimbalXV15.eigsVal,2)
        hPlot = plot(axCompl, gimbalXV15.eigsVal(1:2:end,i)./(2*pi), 'o', 'displayname', 'gimbal XV15');
        set(hPlot, 'MarkerEdgeColor', cmap(i,:));
        set(hPlot, 'MarkerFaceColor', cmap(i,:));
        set(hPlot, 'MarkerSize', 6);
    end
    colormap(cmap);
    cbar = colorbar;
    cbar.Label.String = 'Airspeed [m/s]';
    cbar.Label.Interpreter = 'latex';
    cbar.Ticks = speedData.V;
    cbar.TickLabels = string(round(speedData.V,2));
    cbar.TickLabelInterpreter = 'latex';
    cbar.TicksMode  = 'auto';
    xlabel('Re($\lambda/{2\pi}$)');     xlim([-3, 0.5]);    yline(0,'--')
    ylabel('Im($\lambda/{2\pi}$)');     ylim([-1, 7]);      xline(0,'--')
end


%%% NO GIMBAL REDUCED XV15 model whirl flutter analysis
reducedXV15.constrMat = tiltRotor.constrMat(:,[1:14,17:end]);
[reducedXV15] = computeWhirlflutter(airframe, speedData, reducedXV15);

if any(reducedXV15.unstable ~= 0)
    flutterIdx = find(reducedXV15.unstable ~= 0, 1);
    
    reducedXV15.modeShapeFlutter_airframe = (airframe.modeShape * reducedXV15.constrMat(1:14,:)) ...
                                            * reducedXV15.eigsVec(1:length(reducedXV15.eigsVec)/2,:,flutterIdx);

    reducedXV15.modeShapeFlutter_rotor    = (speedData.modeShape(:,:,flutterIdx) * reducedXV15.constrMat(15:end,:)) ... 
                                            * reducedXV15.eigsVec(1:length(reducedXV15.eigsVec)/2,:,flutterIdx);

    % Animation of the AIRFRAME of the NO GIMBAL XV15 system at whirl flutter speed
    if animFlag.reducedXV15.airframe == true
        rotorAxes = fullMeshVisualization(rotorModel, []);
        fullAnimateModel(1, reducedXV15.modeShapeFlutter_airframe, airfrModel, reducedXV15.eigsVal(:,flutterIdx));
    end

    % Animation of the ROTOR of the NO GIMBAL XV15 system at whirl flutter speed
    if animFlag.reducedXV15.rotor == true
        airfrAxes = fullMeshVisualization(airfrModel, []);
        fullAnimateModel(1, reducedXV15.modeShapeFlutter_rotor, rotorModel, reducedXV15.eigsVal(:,flutterIdx), 1, airfrAxes);
    end

    % Animation of the TOTAL of the NO GIMBAL XV15 system
    if animFlag.reducedXV15.total == true
        fullAnimateModel(1, {reducedXV15.modeShape_airframe,reducedXV15.modeShape_rotor}, {airfrModel,rotorModel}, reducedXV15.eigsVal, 1);
    end

    % Plot TOTAL of the TILTROTOR deformate configuration
    if plotFlag.modes.reducedXV15 == true
        fig.reducedXV15 = cell(1,length(mode2plot));
        for k = 1:length(mode2plot)
            fig.reducedXV15{k} = fullModePlot({airfrModel,rotorModel}, reducedXV15.time2plot, reducedXV15.ampFactor(k), mode2plot(k), {reducedXV15.modeShapeFlutter_airframe,reducedXV15.modeShapeFlutter_rotor}, reducedXV15.eigsVal);
        end
    end

end

% Plot NO GIMBAL REDUCED XV15 system root locus
if plotFlag.rootLocus.reducedXV15 == true
    fig.rootLocus.reducedXV15 = figure(Name='reduced XV15 root-locus');
    hold on;  grid minor;  axis padded;  box on;
    axCompl = gca();
    for i = 1:size(reducedXV15.eigsVal,2)
        hPlot = plot(axCompl, reducedXV15.eigsVal(1:2:end,i)./(2*pi), 'o', 'displayname', 'reduced XV15');
        set(hPlot, 'MarkerEdgeColor', cmap(i,:));
        set(hPlot, 'MarkerFaceColor', cmap(i,:));
        set(hPlot, 'MarkerSize', 6);
    end
    colormap(cmap);
    cbar = colorbar;
    cbar.Label.String = 'Airspeed [m/s]';
    cbar.Label.Interpreter = 'latex';
    cbar.Ticks = speedData.V;
    cbar.TickLabels = string(round(speedData.V,2));
    cbar.TickLabelInterpreter = 'latex';
    cbar.TicksMode  = 'auto';
    xlabel('Re($\lambda/{2\pi}$)');     xlim([-3, 0.5]);    yline(0,'--')
    ylabel('Im($\lambda/{2\pi}$)');     ylim([-1, 7]);      xline(0,'--')
end


%% Export figures
exportFigure(fig,plotFlag);
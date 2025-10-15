clearvars;  close all;  clc;
addpath("src/");
cmapStd = graphicSettings;

airframeNew = load("LabData_Session3.mat");
load('OutputData_Session2.mat');                % see main_workshop2.m


%% INPUT

% actuator stiffness
Kact_aileron = 9800;
Kact_swash = 1e8;

% actuator dynamic model
omega = 2 * pi * 12;
xi = 1; 

% spoint option flag
spointBuild    = false;
spointRotation = true;

% frequency sweep analysis
freqSweep.f0    = 0.1;    
freqSweep.f1    = 10;     
freqSweep.amp   = 0.2;     
freqSweep.tStep = 0.0001; 
freqSweep.t0    = 0;
freqSweep.t1    = 50;

% frequency dweel analysis
freqDwell.n  = 30;
freqDwell.t0 = 0;
freqDwell.t1 = 250/pi;
freqDwell.tStep = 0.0001;

% plot flag
plotFlag.bodeAirframe  = true;
plotFlag.bodeTiltrotor = true;
plotFlag.flutterStd    = false;
plotFlag.flutterCont   = true;
plotFlag.freqSweep     = true;
plotFlag.freqDwell     = true;
plotFlag.logDecrement  = true;

% save flag
saveFigures = true;
dispResults = true;


%% TASK 1

% Extract the modal shapes, building the modal matrices
airframeNew.modeShape = zeros(6*airframeNew.nNodes, airframeNew.nModes);
for iMode = 1:airframeNew.nModes
    airframeNew.modeShape(:, iMode) = reshape(airframeNew.Phi(:,:,iMode)', [airframeNew.nNodes*6,1]);
end

% Definition of the Aileron SPOINT
if spointBuild == true
    % Aileron SPOINT build as model 2
    if spointRotation == true
        airframeNew.sMatrixAileron = zeros(3,size(airframeNew.modeShape,1));
            airframeNew.sMatrixAileron(1,(find(airframeNew.nodeLabel==1205)-1)*6+4) = -1;
            airframeNew.sMatrixAileron(2,(find(airframeNew.nodeLabel==1205)-1)*6+5) = -1;
            airframeNew.sMatrixAileron(3,(find(airframeNew.nodeLabel==1205)-1)*6+6) = -1;
            airframeNew.sMatrixAileron(1,(find(airframeNew.nodeLabel==2205)-1)*6+4) = 1;
            airframeNew.sMatrixAileron(2,(find(airframeNew.nodeLabel==2205)-1)*6+5) = 1;
            airframeNew.sMatrixAileron(3,(find(airframeNew.nodeLabel==2205)-1)*6+6) = 1;
            airframeNew.sPointAileron = Y * RD * RS * airframeNew.sMatrixAileron;   % Applicazione matrice di rotazione
    else
        airframeNew.sMatrixAileron = zeros(1,size(airframeNew.modeShape,1));
            airframeNew.sMatrixAileron(1,(find(airframeNew.nodeLabel==1205)-1)*6+4) = -1;
            airframeNew.sMatrixAileron(1,(find(airframeNew.nodeLabel==2205)-1)*6+4) = 1;
            airframeNew.sPointAileron = airframeNew.sMatrixAileron;
    end
else
    % Aileron SPOINT model 3
    if spointRotation == true
        airframeNew.sMatrixAileron = zeros(3,size(airframeNew.modeShape,1));
            airframeNew.sMatrixAileron(1,(find(airframeNew.nodeLabel==91919)-1)*6+1) = 1;
            airframeNew.sMatrixAileron(2,(find(airframeNew.nodeLabel==91919)-1)*6+2) = 1;
            airframeNew.sMatrixAileron(3,(find(airframeNew.nodeLabel==91919)-1)*6+3) = 1;
            airframeNew.sPointAileron = [1 0 0] * RD * RS * airframeNew.sMatrixAileron;   % Applicazione matrice di rotazione 
    else
        airframeNew.sMatrixAileron = zeros(1,size(airframeNew.modeShape,1));
            airframeNew.sMatrixAileron(1,(find(airframeNew.nodeLabel==91919)-1)*6+1) = 1;
            airframeNew.sPointAileron = airframeNew.sMatrixAileron;   % Applicazione matrice di rotazione
    end
end

%%%
% DEBUG
aileronSpointEffect =  airframeNew.sPointAileron * airframeNew.modeShape;      % check if SPOINT is null
%%%

% Compute ACTUATOR added modal stiffness and damping 
airframeNew.Khha = airframeNew.modeShape' * airframeNew.sPointAileron' * Kact_aileron * airframeNew.sPointAileron * airframeNew.modeShape;
airframeNew.Ktot = airframeNew.Khh + airframeNew.Khha;
airframeNew.Chh  = zeros(size(airframeNew.Mhh));

% Compute coupling term between AIRFRAME and AILERON ACTUATOR
airframeNew.forzAct = -airframeNew.Mhh \ (airframeNew.modeShape' * airframeNew.sPointAileron' * Kact_aileron);

% Build the AIRFRAME + ACTUATOR state space system
airframeNew.A = [zeros(size(airframeNew.Mhh)),       eye(size(airframeNew.Mhh)),         zeros(size(airframeNew.Mhh,1), 1),  zeros(size(airframeNew.Mhh,1), 1);  ...
                 -airframeNew.Mhh\airframeNew.Ktot,  -airframeNew.Mhh\airframeNew.Chh,   airframeNew.forzAct,                zeros(size(airframeNew.Mhh,1), 1);  ...
                 zeros(1, size(airframeNew.Mhh,2)),  zeros(1, size(airframeNew.Mhh,2)),  zeros(1,1),                         eye(1,1);                           ...
                 zeros(1, size(airframeNew.Mhh,2)),  zeros(1, size(airframeNew.Mhh,2)),  -omega^2,                           -2*xi*omega];

airframeNew.B = [zeros(2*size(airframeNew.Mhh,1),1); 0; omega^2];

airframeNew.C = zeros(1,size(airframeNew.A,2));
airframeNew.C(1,2*size(airframeNew.Mhh)+1) = 1;

airframeNew.D = 0;

airframeNew.stateSpace  = ss(airframeNew.A, airframeNew.B, airframeNew.C, airframeNew.D);

% Compute the transfer function and bode plot of the AILERON ACTUATOR in the AIRFRAME system
airframeNew.transferFun = tf(airframeNew.stateSpace);
[fig.bodeAirframe] = ssBodePlot(airframeNew.transferFun, plotFlag.bodeAirframe);


%% TASK 2

% Find hub and swashplate nodes in TILT ROTOR system
tiltRotorNew.hubNodeRow   = find(airframeNew.nodeLabel == 1000007);
tiltRotorNew.swashNodeRow = find(airframeNew.nodeLabel == 100000);

% Build total TILT ROTOR system constraint matrix
tiltRotorNew.phiA(1:6 ,:) = squeeze(airframeNew.Phi(tiltRotorNew.hubNodeRow,:,:));
tiltRotorNew.phiA(7:12,:) = squeeze(airframeNew.Phi(tiltRotorNew.swashNodeRow,:,:));
tiltRotorNew.constrMat = [  eye(length(airframeNew.Mhh)),       zeros(length(airframeNew.Mhh),17);
                            grdRotor.phiR1\tiltRotorNew.phiA,   grdRotor.phiR2;
                            zeros(17,length(airframeNew.Mhh)),  eye(17,17) ];

% Build TILT ROTOR dynamic system from airframe and rotor modal system
Mmix = zeros(size(airframeNew.Mhh)+size(rotorData.M));
Mmix(1:13,1:13) = airframeNew.Mhh;
Mmix(14:end,14:end) = rotorData.M;

Cmix = zeros(size(airframeNew.Mhh)+size(rotorData.M));
Cmix(14:end,14:end) = rotorData.C;

Kmix = zeros(size(airframeNew.Mhh)+size(rotorData.M));
Kmix(1:13,1:13) = airframeNew.Khh;
Kmix(14:end,14:end) = rotorData.K;

% Definition of total TILT ROTOR mode shape and node label
tiltRotorNew.nodeLabel = [airframeNew.nodeLabel; rotorData.nodeLabel];
tiltRotorNew.modeShape = zeros(size(airframeNew.modeShape)+size(rotorData.modeShape));
tiltRotorNew.modeShape(1:size(airframeNew.modeShape,1),1:size(airframeNew.modeShape,2)) = airframeNew.modeShape;
tiltRotorNew.modeShape(size(airframeNew.modeShape,1)+1:end,size(airframeNew.modeShape,2)+1:end) = rotorData.modeShape;

% Dynamic system TILT ROTOR definition using CONSTRAINT MATRIX
tiltRotorNew.M = tiltRotorNew.constrMat' * Mmix * tiltRotorNew.constrMat;
tiltRotorNew.C = tiltRotorNew.constrMat' * Cmix * tiltRotorNew.constrMat;
tiltRotorNew.K = tiltRotorNew.constrMat' * Kmix * tiltRotorNew.constrMat;

% Definition of the Aileron SPOINT
if spointBuild == true
    % Aileron SPOINT build as model 2
    if spointRotation == true
        tiltRotorNew.sMatrixAileron = zeros(3,size(tiltRotorNew.modeShape,1));
            tiltRotorNew.sMatrixAileron(1,(find(tiltRotorNew.nodeLabel==1205)-1)*6+4) = -1;
            tiltRotorNew.sMatrixAileron(2,(find(tiltRotorNew.nodeLabel==1205)-1)*6+5) = -1;
            tiltRotorNew.sMatrixAileron(3,(find(tiltRotorNew.nodeLabel==1205)-1)*6+6) = -1;
            tiltRotorNew.sMatrixAileron(1,(find(tiltRotorNew.nodeLabel==2205)-1)*6+4) = 1;
            tiltRotorNew.sMatrixAileron(2,(find(tiltRotorNew.nodeLabel==2205)-1)*6+5) = 1;
            tiltRotorNew.sMatrixAileron(3,(find(tiltRotorNew.nodeLabel==2205)-1)*6+6) = 1;
            tiltRotorNew.sPointAileron = Y * RD * RS * tiltRotorNew.sMatrixAileron;   % Applicazione matrice di rotazione
    else
        tiltRotorNew.sMatrixAileron = zeros(1,size(tiltRotorNew.modeShape,1));
            tiltRotorNew.sMatrixAileron(1,(find(tiltRotorNew.nodeLabel==1205)-1)*6+4) = -1;
            tiltRotorNew.sMatrixAileron(1,(find(tiltRotorNew.nodeLabel==2205)-1)*6+4) = 1;
            tiltRotorNew.sPointAileron = tiltRotorNew.sMatrixAileron;
    end
else
    % Aileron SPOINT model 3
    if spointRotation == true
        tiltRotorNew.sMatrixAileron = zeros(3,size(tiltRotorNew.modeShape,1));
            tiltRotorNew.sMatrixAileron(1,(find(tiltRotorNew.nodeLabel==91919)-1)*6+1) = 1;
            tiltRotorNew.sMatrixAileron(2,(find(tiltRotorNew.nodeLabel==91919)-1)*6+2) = 1;
            tiltRotorNew.sMatrixAileron(3,(find(tiltRotorNew.nodeLabel==91919)-1)*6+3) = 1;
            tiltRotorNew.sPointAileron = [1 0 0] * RD * RS * tiltRotorNew.sMatrixAileron;   % Applicazione matrice di rotazione 
    else
        tiltRotorNew.sMatrixAileron = zeros(1,size(tiltRotorNew.modeShape,1));
            tiltRotorNew.sMatrixAileron(1,(find(tiltRotorNew.nodeLabel==91919)-1)*6+1) = 1;
            tiltRotorNew.sPointAileron = tiltRotorNew.sMatrixAileron;   % Applicazione matrice di rotazione
    end
end

%%%
% DEBUG
aileronSpointEffect =  tiltRotorNew.sPointAileron * tiltRotorNew.modeShape;      % check if SPOINT is null
%%%

% Definition of the SPOINT for swashplate
if spointRotation == true
    tiltRotorNew.sMatrixSwash1 = zeros(3,size(tiltRotorNew.modeShape,1));
        tiltRotorNew.sMatrixSwash1(1, (find(tiltRotorNew.nodeLabel==100100)-1)*6+1) = 1;
        tiltRotorNew.sMatrixSwash1(2, (find(tiltRotorNew.nodeLabel==100100)-1)*6+2) = 1;
        tiltRotorNew.sMatrixSwash1(3, (find(tiltRotorNew.nodeLabel==100100)-1)*6+3) = 1;
        tiltRotorNew.sPointSwash1 = [1 0 0] * RD * RS * tiltRotorNew.sMatrixSwash1;
else
    tiltRotorNew.sMatrixSwash1 = zeros(1,size(tiltRotorNew.modeShape,1));
        tiltRotorNew.sMatrixSwash1(1, (find(tiltRotorNew.nodeLabel==100100)-1)*6+1) = 1;
        tiltRotorNew.sPointSwash1 = tiltRotorNew.sMatrixSwash1;
end

if spointRotation == true
    tiltRotorNew.sMatrixSwash2 = zeros(3,size(tiltRotorNew.modeShape,1));
        tiltRotorNew.sMatrixSwash2(1, (find(tiltRotorNew.nodeLabel==100200)-1)*6+1) = 1;
        tiltRotorNew.sMatrixSwash2(2, (find(tiltRotorNew.nodeLabel==100200)-1)*6+2) = 1;
        tiltRotorNew.sMatrixSwash2(3, (find(tiltRotorNew.nodeLabel==100200)-1)*6+3) = 1;
        tiltRotorNew.sPointSwash2 = [1 0 0] * RD * RS * tiltRotorNew.sMatrixSwash2;
else
    tiltRotorNew.sMatrixSwash2 = zeros(1,size(tiltRotorNew.modeShape,1));
        tiltRotorNew.sMatrixSwash2(1, (find(tiltRotorNew.nodeLabel==100200)-1)*6+1) = 1;
        tiltRotorNew.sPointSwash2 = tiltRotorNew.sMatrixSwash2;
end

if spointRotation == true
    tiltRotorNew.sMatrixSwash3 = zeros(3,size(tiltRotorNew.modeShape,1));
        tiltRotorNew.sMatrixSwash3(1, (find(tiltRotorNew.nodeLabel==100300)-1)*6+1) = 1;
        tiltRotorNew.sMatrixSwash3(2, (find(tiltRotorNew.nodeLabel==100300)-1)*6+2) = 1;
        tiltRotorNew.sMatrixSwash3(3, (find(tiltRotorNew.nodeLabel==100300)-1)*6+3) = 1;
        tiltRotorNew.sPointSwash3 = [1 0 0] * RD * RS * tiltRotorNew.sMatrixSwash3;
else
    tiltRotorNew.sMatrixSwash3 = zeros(1,size(tiltRotorNew.modeShape,1));
        tiltRotorNew.sMatrixSwash3(1, (find(tiltRotorNew.nodeLabel==100300)-1)*6+1) = 1;
        tiltRotorNew.sPointSwash3 = tiltRotorNew.sMatrixSwash3;
end


% Compute ACTUATORS added modal stiffness
K_act_aileron = tiltRotorNew.modeShape' * tiltRotorNew.sPointAileron' * Kact_aileron * tiltRotorNew.sPointAileron * tiltRotorNew.modeShape;
K_act_swash1  = tiltRotorNew.modeShape' * tiltRotorNew.sPointSwash1'  * Kact_swash   * tiltRotorNew.sPointSwash1  * tiltRotorNew.modeShape;
K_act_swash2  = tiltRotorNew.modeShape' * tiltRotorNew.sPointSwash2'  * Kact_swash   * tiltRotorNew.sPointSwash2  * tiltRotorNew.modeShape;
K_act_swash3  = tiltRotorNew.modeShape' * tiltRotorNew.sPointSwash3'  * Kact_swash   * tiltRotorNew.sPointSwash3  * tiltRotorNew.modeShape;

K_act_tot = K_act_aileron + K_act_swash1 + K_act_swash2 + K_act_swash3;
tiltRotorNew.Ktot = tiltRotorNew.constrMat' * K_act_tot * tiltRotorNew.constrMat + tiltRotorNew.K;

% Compute coupling term between tiltrotor and actuators using CONSTRAINT MATRIX
forz_act_aileron = tiltRotorNew.constrMat' * tiltRotorNew.modeShape' * tiltRotorNew.sPointAileron' * Kact_aileron ;
forz_act_swash1  = tiltRotorNew.constrMat' * tiltRotorNew.modeShape' * tiltRotorNew.sPointSwash1'  * Kact_swash ;
forz_act_swash2  = tiltRotorNew.constrMat' * tiltRotorNew.modeShape' * tiltRotorNew.sPointSwash2'  * Kact_swash ;
forz_act_swash3  = tiltRotorNew.constrMat' * tiltRotorNew.modeShape' * tiltRotorNew.sPointSwash3'  * Kact_swash ;
tiltRotorNew.forzAct = [-tiltRotorNew.M\forz_act_aileron, -tiltRotorNew.M\forz_act_swash1, -tiltRotorNew.M\forz_act_swash2, -tiltRotorNew.M\forz_act_swash3];

% Build the TILT ROTOR + ACTUATORS state space system
tiltRotorNew.A = [zeros(30,30),                       eye(30),                         zeros(30,4),           zeros(30,4);  ...
                  -tiltRotorNew.M\tiltRotorNew.Ktot,  -tiltRotorNew.M\tiltRotorNew.C,  tiltRotorNew.forzAct,  zeros(30,4);  ...
                  zeros(4,30),                        zeros(4,30),                     zeros(4,4),            eye(4,4);     ...
                  zeros(4,30),                        zeros(4,30),                     -omega^2 * eye(4),     -2*xi*omega * eye(4)];

tiltRotorNew.B = zeros(size(tiltRotorNew.A,1),1);
tiltRotorNew.B(size(tiltRotorNew.A,1)-3) = omega^2;

tiltRotorNew.C = zeros(1,68);
tiltRotorNew.C(1,61) = 1;

tiltRotorNew.D = 0;

tiltRotorNew.stateSpace = ss(tiltRotorNew.A, tiltRotorNew.B, tiltRotorNew.C, tiltRotorNew.D);

% Compute the transfer function and bode plot of the AILERON ACTUATOR in the TILT ROTOR system
tiltRotorNew.transferFun = tf(tiltRotorNew.stateSpace);
[fig.bodeTiltrotor] = ssBodePlot(tiltRotorNew.transferFun, plotFlag.bodeTiltrotor);


%% TASK 3 - Build aeroelastic model

% Build ACTUATOR XV15 model @ every airspeed
for i = 1:length(speedData.V)
    % Build ACTUATOR XV15 dynamic system from airframe and rotor modal system
    Mmix = zeros(size(airframeNew.Mhh)+size(speedData.M(:,:,i)));
    Mmix(1:13,1:13) = airframeNew.Mhh;
    Mmix(14:end,14:end) = speedData.M(:,:,i);
    
    Cmix = zeros(size(airframeNew.Mhh)+size(speedData.M(:,:,i)));
    Cmix(14:end,14:end) = speedData.C(:,:,i);
    
    Kmix = zeros(size(airframeNew.Mhh)+size(speedData.M(:,:,i)));
    Kmix(1:13,1:13) = airframeNew.Khh;
    Kmix(14:end,14:end) = speedData.K(:,:,i);
    
    % Definition of total ACTUATOR XV15 mode shape and node label
    actXV15.nodeLabel = [airframeNew.nodeLabel; speedData.nodeLabel];
    actXV15.modeShape(:,:,i) = zeros(size(airframeNew.modeShape)+size(speedData.modeShape(:,:,i)));
    actXV15.modeShape(1:size(airframeNew.modeShape,1),1:size(airframeNew.modeShape,2),i) = airframeNew.modeShape;
    actXV15.modeShape(size(airframeNew.modeShape,1)+1:end,size(airframeNew.modeShape,2)+1:end,i) = speedData.modeShape(:,:,i);
    
    % Dynamic system ACTUATOR XV15 definition using CONSTRAINT MATRIX
    actXV15.Mhh(:,:,i) = tiltRotorNew.constrMat' * Mmix * tiltRotorNew.constrMat;
    actXV15.Chh(:,:,i) = tiltRotorNew.constrMat' * Cmix * tiltRotorNew.constrMat;
    actXV15.Khh(:,:,i) = tiltRotorNew.constrMat' * Kmix * tiltRotorNew.constrMat;
    
    % Compute ACTUATORS added modal stiffness
    K_act_aileron = actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointAileron' * Kact_aileron * tiltRotorNew.sPointAileron * actXV15.modeShape(:,:,i);
    K_act_swash1  = actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointSwash1'  * Kact_swash   * tiltRotorNew.sPointSwash1  * actXV15.modeShape(:,:,i);
    K_act_swash2  = actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointSwash2'  * Kact_swash   * tiltRotorNew.sPointSwash2  * actXV15.modeShape(:,:,i);
    K_act_swash3  = actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointSwash3'  * Kact_swash   * tiltRotorNew.sPointSwash3  * actXV15.modeShape(:,:,i);
    
    K_act_tot = K_act_aileron + K_act_swash1 + K_act_swash2 + K_act_swash3;
    actXV15.Ktot(:,:,i) = tiltRotorNew.constrMat' * K_act_tot * tiltRotorNew.constrMat + actXV15.Khh(:,:,i);
    
    % Compute coupling term between tiltrotor and actuators using CONSTRAINT MATRIX
    forz_act_aileron = tiltRotorNew.constrMat' * actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointAileron' * Kact_aileron ;
    forz_act_swash1  = tiltRotorNew.constrMat' * actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointSwash1'  * Kact_swash ;
    forz_act_swash2  = tiltRotorNew.constrMat' * actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointSwash2'  * Kact_swash ;
    forz_act_swash3  = tiltRotorNew.constrMat' * actXV15.modeShape(:,:,i)' * tiltRotorNew.sPointSwash3'  * Kact_swash ;
    actXV15.forzAct(:,:,i) = [-actXV15.Mhh(:,:,i)\forz_act_aileron, ...
                              -actXV15.Mhh(:,:,i)\forz_act_swash1,  ...
                              -actXV15.Mhh(:,:,i)\forz_act_swash2,  ...
                              -actXV15.Mhh(:,:,i)\forz_act_swash3];
    
    % Build the ACTUATOR XV15 + ACTUATORS state space system & transfer function
    actXV15.A(:,:,i) = [zeros(30,30),                             eye(30),                                 zeros(30,4),             zeros(30,4);  ...
                        -actXV15.Mhh(:,:,i)\actXV15.Ktot(:,:,i),  -actXV15.Mhh(:,:,i)\actXV15.Chh(:,:,i),  actXV15.forzAct(:,:,i),  zeros(30,4);  ...
                        zeros(4,30),                              zeros(4,30),                             zeros(4,4),              eye(4,4);     ...
                        zeros(4,30),                              zeros(4,30),                             -omega^2 * eye(4),       -2*xi*omega * eye(4)];
    
    actXV15.B(:,:,i) = zeros(size(actXV15.A(:,:,i),1),1);
    actXV15.B(size(actXV15.A(:,:,i),1)-3,i) = omega^2;
    
    actXV15.C(:,:,i) = zeros(1,68);
    actXV15.C(1,61,i) = 1;
    
    actXV15.D(i) = 0;
    
    actXV15.stateSpace{i}  = ss(actXV15.A(:,:,i), actXV15.B(:,:,i), actXV15.C(:,:,i), actXV15.D(i));
    actXV15.transferFun{i} = tf(actXV15.stateSpace{i});

    actXV15.eigsVal(:,i) = eig(actXV15.A(:,:,i));

end

% V-G and V-F plots for the complete system
[actXV15.flutterStd,  fig.flutterStd]  = whirlFlutterStandard(actXV15.A, speedData.V, plotFlag.flutterStd);
[actXV15.flutterCont, fig.flutterCont] = whirlFlutterContinuation(actXV15.A, speedData.V, plotFlag.flutterCont);


%% TASK 3 - Frequency sweep input 

% Built frequency sweep time and signal
freqSweep.time   = freqSweep.t0:freqSweep.tStep:freqSweep.t1; 
freqSweep.signal = chirp(freqSweep.time, freqSweep.f0, freqSweep.time(end), freqSweep.f1, 'linear') * freqSweep.amp; 

% Simulate the time response
[freqSweep.y, freqSweep.t_out, freqSweep.x] = lsim(actXV15.stateSpace{9}, freqSweep.signal, freqSweep.time);

% Plot the results
if plotFlag.freqSweep == true  
    fig.freqSweep = figure(Name='Frequency Sweep');
    % Input signal (chirp)
    subplot(3, 1, 1);
        hold on; grid minor; axis padded; box on;
        plot(freqSweep.time, freqSweep.signal, Color=cmapStd(1,:));
        %title('Input Signal: Frequency Sweep (Chirp)', 'Interpreter', 'Latex');
        ylabel('$u(t)$');       %xlabel('Time [sec]');       
        xlim([-1 51]);
    
    % Output of the system
    subplot(3, 1, 2);
        hold on; grid minor; axis padded; box on;
        plot(freqSweep.t_out, freqSweep.y, Color=cmapStd(2,:));
        %title('Output Signal ($y(t)$)');
        ylabel('$y(t)$');       %xlabel('Time [sec]');    
        xlim([-1 51]);
    
    % State variables
    subplot(3, 1, 3);
        hold on; grid minor; axis padded; box on;
        plot(freqSweep.t_out, freqSweep.x(:, 1), Color=cmapStd(3,:));
        plot(freqSweep.t_out, freqSweep.x(:, 2), Color=cmapStd(5,:));
        %title('State Variables');
        ylabel('States $x(t)$');        xlabel('Time [sec]');           
        legend('$x_1(t)$', '$x_2(t)$', 'Interpreter', 'Latex');
        xlim([-1 51]);
end


%% TASK 3 - Frequency dwell input

% Time at which to compute frequency
[freqSweep.maxState5, freqSweep.maxIdx] = max(freqSweep.x(:, 8));
freqSweep.maxTime = freqSweep.time(freqSweep.maxIdx);  

% Compute frequency at time t_query for linear chirp
freqSweep.maxFreq = freqSweep.f0 + (freqSweep.f1 - freqSweep.f0) * (freqSweep.maxTime / freqSweep.t1);

% Definition of the input signal
freqDwell.time   = freqDwell.t0:freqDwell.tStep:freqDwell.t1;  
freqDwell.signal = sin(2 * pi * freqDwell.time(freqDwell.time<=freqDwell.n/pi) * freqSweep.maxFreq);
freqDwell.signal = [freqDwell.signal, zeros(1,length(find(freqDwell.time>freqDwell.n/pi)))];

% Simulate the time response
[freqDwell.y, freqDwell.t_out, freqDwell.x] = lsim(actXV15.stateSpace{9}, freqDwell.signal, freqDwell.time);

% Plot the results
if plotFlag.freqDwell == true
    fig.freqDwell = figure(Name='Frequency Dwell');
        % Input signal
        subplot(3, 1, 1);
        hold on; grid minor; axis padded; box on;
        plot(freqDwell.time, freqDwell.signal, Color=cmapStd(1,:));
        %title('Input Signal');
        ylabel('$u(t)$');       %xlabel('Time (s)');  
        xlim([-1 81]);

        % Output of the system
        subplot(3, 1, 2);
        hold on; grid minor; axis padded; box on;
        plot(freqDwell.t_out, freqDwell.y, Color=cmapStd(2,:));
        %title('Output Signal ($y(t)$)');
        ylabel('$y(t)$');       %xlabel('Time (s)');
        xlim([-1 81]);
        
        % State variables
        subplot(3, 1, 3);
        hold on; grid minor; axis padded; box on;
        plot(freqDwell.t_out, freqDwell.x(:, 1), Color=cmapStd(3,:));
        plot(freqDwell.t_out, freqDwell.x(:, 2), Color=cmapStd(5,:));
        %title('State Variables');
        ylabel('States $x(t)$');        xlabel('Time (s)');
        xlim([-1 81]);
        legend('$x_{1}(t)$', '$x_{2}(t)$');
end


%% TASK 3 - Logaritmic decrement computation

% Identify successive peaks in the response
[logDec.peaks, logDec.peak_indices] = findpeaks(freqDwell.x(:,1), freqDwell.time); % Find peaks and their times
 
% Isolate the decrement line
logDec.peaks_dec = logDec.peaks(45:end);
logDec.peak_index_dec = logDec.peak_indices(45:end);

% Compute logarithmic decrement (using consecutive peaks)
logDec.x1 = logDec.peaks_dec(1:end-1);    % First peak amplitude
logDec.x2 = logDec.peaks_dec(2:end);      % Second peak amplitude
logDec.delta = log(logDec.x1 ./ logDec.x2);
logDec.deltaMean = mean(logDec.delta);

% Compute damping ratio from logarithmic decrement
logDec.dampRatio = logDec.deltaMean / sqrt(4 * pi^2 + logDec.deltaMean^2);

% Plot the logaritmic decrement results
if plotFlag.logDecrement == true
    fig.logDecrement = figure(Name='Log decremet');
        hold on; grid minor; axis padded; box on;
        plot(freqDwell.time, freqDwell.x(:,1), Color=cmapStd(1,:));
        plot(logDec.peak_index_dec, logDec.peaks_dec, Color=cmapStd(2,:), LineWidth=1.5)
        %title('Underdamped Response of a System');
        legend('System response', 'Logaritmic decrement')
        xlabel('Time (s)');         ylabel('Amplitude');
        xlim([-1 81]);
end


%% EXPORT FIGURES

if saveFigures == true
    exportFigure(fig,plotFlag);
end

if dispResults == true
    fprintf('Frequency and damping from the eigenanalysis: \n');
    fprintf('\tomega = %fHz \t xi = %f \n\n', actXV15.flutterCont.freq_Hz(8), actXV15.flutterCont.dampRatio(8,9));
    fprintf('Frequency and damping from the dweel:\n');
    fprintf('\tomega = %fHz \t xi = %f \n\n', freqSweep.maxFreq, logDec.dampRatio);
    fprintf('Analysis computed at %fm/s \n', speedData.V(9));
end
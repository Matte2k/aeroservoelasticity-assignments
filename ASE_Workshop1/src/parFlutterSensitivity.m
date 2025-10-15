function [flutterSpeed,fig] = parFlutterSensitivity(flutterSens,wing,data,plotFlag)
%PARALLEL FLUTTER SENSITIVITY - Compute flutter sensitivity analysis with parfor
%   Function computing the flutter sensitivity analysis for different
%   values of actuator stiffness and damping using a parallel approach to
%   optimize the cpu time cost. For best performance the Parallel Computing
%   Toolbox must be installed
%
%   SYNTAX:
%       [flutterSpeed,fig] = parFlutterSensitivity(flutterSens,wing,data,nWorker,plotFlag)
%
%   INPUT:
%       flutterSens,  struct: defines sensitivity analysis variables
%                              > flow characterization (see buildFlow inputs)
%                              > damping and stiffness vectors
%                              > number of cpu for parallel pool run
%       wing,         struct: wing sweep, dihedral and chord
%       data,         struct: data provided in 'LabData_Session1.mat'
%       plotFlag,       bool: set if the figures must be generated
%
%   OUTPUT: 
%       flutterSpeed,  double: matrix containing the flutter speed for the
%                              different damping (cols) and stiffness (row)
%                              '-1' value means that no flutter has been found
%       fig,           figure: rappresentation of the sensitivity analysis
%
%   OPTIONAL INPUT:
%       plotFlag: set by default at false
%
%

    % Optional input
    if nargin < 4
        plotFlag = false;
    end

    % Initialize flow for flutter sensitivity analysis
    flow = buildFlow(flutterSens.v_limit, flutterSens.dist, flutterSens.rho);

    % Define local variable to be used in parfor
    chord = wing.chord;
    sweep = wing.sweep;
    dihed = wing.dihed;
    rho   = flow.rho;
    v_vec = flow.v_vec;
    dampVec  = flutterSens.dampVec;
    stiffVec = flutterSens.stiffVec;
    Phi       = data.Phi;
    nodeLabel = data.nodeLabel;

    % Initialize output variables
    parallelOutput = zeros(1,length(dampVec));
    flutterSpeed   = zeros(length(stiffVec),length(dampVec));
    
    % Perform sensitivity analysis
    parpool('local', flutterSens.nWorker);
    for i = 1:length(stiffVec)
        parfor j = 1:length(dampVec)
            % Compute current iteration hinge stiffness and damping
            [iterC] = addActuatorTerm(dampVec(j),  sweep, dihed, Phi, nodeLabel);
            [iterK] = addActuatorTerm(stiffVec(i), sweep, dihed, Phi, nodeLabel);
    
            % Build current iteration system
            iterSystem = struct;
                iterSystem.type = 'state_space';
                iterSystem.M = data.Mhh;
                iterSystem.K = iterK.modal + data.Khh;
                iterSystem.C = iterC.modal;
    
            % Compute standard flutter analysis for each speed
            for k = 1:length(v_vec)
                [sys] = buildAeroSystem(iterSystem, data, chord, rho, v_vec(k));
                eigsVal = eig(sys.A, sys.V);
    
                % Stop the cycle if a flutter speed is found
                if any(real(eigsVal)>0)
                    parallelOutput(1,j) = v_vec(k);
                    break;
                else
                    parallelOutput(1,j) = nan;    % no flutter speed
                end
            end            

        end   
        % Save the parfor output in the final vector
        flutterSpeed(i,:) = parallelOutput;
    end



    % Plot surface
    if plotFlag == true
        [X, Y] = meshgrid(flutterSens.stiffVec, flutterSens.dampVec);
        fig = figure(Name='Flutter sensitivity');
        hold on;  grid minor;  axis padded;  box on;
        surf(X,Y,flutterSpeed')
        surf(X,Y,51.4*ones(size(X)),FaceAlpha=0.3,FaceColor='b',EdgeColor='none')
        surf(X,Y,128.6*ones(size(X)),FaceAlpha=0.3,FaceColor='r',EdgeColor='none')
        xlabel('Actuator Stiffness');
        ylabel('Actuator Damping');
        zlabel('Flutter Speed [m/s]');
        view(-40,30);

    else
        fig = [];
    end

    delete(gcp('nocreate'));
    
end


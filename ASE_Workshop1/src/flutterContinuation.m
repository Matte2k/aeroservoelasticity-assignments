function [flutter_solution,fig] = flutterContinuation(flutter_system,data,chord,options,plotFlag)
%FLUTTER CONTINUATION - compute flutter analysis with continuation approach
%   Function computing flutter analysis solving the eigenvalue problem 
%   at different flow field velocity using continuation approach
%
%   SYNTAX:
%       [flutter_solution,fig] = flutterContinuation(flutter_system,data,chord,options,plotFlag)
%
%   INPUT:
%       flutter_system, struct: state-space system struct of the considered
%                               model containin even the 'flow' struct
%       data,           struct: data provided in 'LabData_Session1.mat'
%       chord,          double: chord value [m] of the considered wing
%       options,        struct: settings for continuation approach
%                                > modes2solve: mode number to be computed
%                                > maxIter: maximum iteration for Newton method
%                                > toll: tollerance for Newton method
%       plotFlag(*)       bool: set if the figures must be generated
%
%   OUTPUT:
%       flutter_solution, struct: flutter analysis results
%       fig,              figure: plot of flutter results
%
%   OPTIONAL INPUT:
%       plotFlag: set by default as false
%
%   TIP:  just the first 28 eigenvalues (and their respective eigenvectors)
%         are intersting because, for flutter speed, we are intrested just
%         into structural modes. They are organized in matrix as follows:
%           > flutter_solution.eigsVal: rows=eigenvalue and columns=speed)
%           > flutter_solution.eigsVec: rows+cols=eigenvector and 3rd
%                                       dimenison is speed
%
%

    % Optional input
    if nargin < 5
        plotFlag = false;
    end

    if nargin < 4 || isempty(options)
            options.modes2solve = 1:2:28;
            options.maxIter = 100;
            options.toll    = 1e-9; 
    end

    % Initialize eigenvalues and eigenvector matrix
    flutter_solution.eigsVal = zeros(length(options.modes2solve), length(flutter_system.flow.v_vec));
    flutter_solution.eigsVec = zeros(length(options.modes2solve), ...
                                      length(flutter_system.M)*2+length(data.Aaero), ...
                                      length(flutter_system.flow.v_vec));
    
    dU = diff(flutter_system.flow.v_vec(1:2));     % NOTE: valid only if the speed distribution is linear
    warning("off");

    for i = 1:length(flutter_system.flow.v_vec)
    
        v_inf = flutter_system.flow.v_vec(i);
        [sys] = buildAeroSystem(flutter_system, data, chord, flutter_system.flow.rho, v_inf);
    
        % Build initiale guess for the first speed of the vector flutter_system.flow.v_vec and sort
        if i==1
            [q0_mat, lambda0_mat] = eig(sys.A, sys.V);
            [lambda0_mat, q0_mat] = sorting_the_eigs(lambda0_mat, q0_mat);
        end
    
        % Solve the two systems of the continuation approach     
        for j = 1:length(options.modes2solve)
            % Define the current mode to be tracked as the speed increases
            mode2track = options.modes2solve(j);
            
            % Define initial guess for the current mode at the i-th speed
            q_0      = q0_mat(:,mode2track);
            lambda_0 = lambda0_mat(mode2track, mode2track);
    
            % Solution of the FIRST SYSTEM with Newton-like approach
            delta_sol = ones(size(flutter_solution.eigsVec,2)+1, 1);
            iter = 1;
            while ((norm(delta_sol(2:end))>options.toll || norm(delta_sol(1))>options.toll) && iter<options.maxIter) 
                Matrix_A = [ sys.V * q_0,  (sys.V * lambda_0-sys.A);
                             zeros(1,1),   2*q_0' ];
    
                Vector_b_eig= [ -(sys.V*lambda_0-sys.A)*q_0;
                                1 - q_0'*q_0];
    
                delta_sol = Matrix_A\Vector_b_eig;
    
                q_0      = q_0 + delta_sol(2:end);
                lambda_0 = lambda_0 + delta_sol(1);
    
                iter = iter + 1;
            end
    
            % Store current mode computed at i^th speed in solution struct
            flutter_solution.eigsVal(j,i)   = lambda_0;
            flutter_solution.eigsVec(j,:,i) = q_0;
    
            % Solution of the SECOND SYSTEM to get the initial guess at the next speed
            Vector_b_cont   = [-(sys.dV_dU * lambda_0 - sys.dA_dU)*q_0; 0];
            derivatives_sol = Matrix_A\Vector_b_cont;
            q_0      = q_0 + derivatives_sol(2:end)*dU;
            lambda_0 = lambda_0 + derivatives_sol(1)*dU;
    
            % Store the initial guess for the current mode at the next speed
            lambda0_mat(mode2track,mode2track) = lambda_0;
            q0_mat(:,mode2track) = q_0;
        end    
    end
    warning("on");
    
    % Compute damping ratios and eigenfrequencies
    flutter_solution.dampRatio = -(real(flutter_solution.eigsVal)./abs(flutter_solution.eigsVal));
    flutter_solution.freq_Hz   = imag(flutter_solution.eigsVal)./(2*pi);

    % Plot of the results just computed
    if plotFlag == true
        legendCell = cell(1,length(options.modes2solve));
        fig.flutter = figure(Name='Flutter analysis with continuation');
        %fig.damping = figure(Name='Damping ratio');
        t1 = tiledlayout(1,2);
        nexttile
            hold on;  grid minor;  axis padded;  box on;
            for i = 1:length(options.modes2solve)
                plot(flutter_system.flow.v_vec, flutter_solution.dampRatio(i,:), '-')
                legendCell{i} = ['Mode ', num2str(options.modes2solve(i))];
            end
            %legend(legendCell);
            xlabel('Velocity [m/s]');      ylabel('Damping');
    
        %fig.freq = figure(Name='Eigenfrequencies');
        nexttile
            hold on;  grid minor;  axis padded;  box on;
            for i = 1:length(options.modes2solve)
                plot(flutter_system.flow.v_vec, -imag(flutter_solution.eigsVal(i,:))./(2*pi), '-')
                legendCell{i} = ['Mode ', num2str(i)];
            end
            %legend(legendCell);
            xlabel('Velocity [m/s]');      ylabel('Frequency Hz');
        leg = legend(legendCell, Orientation='horizontal');
        leg.Layout.Tile = 'south';
    else
        fig.damping   = [];
        fig.freq      = [];
        fig.flutter   = [];
    end

end


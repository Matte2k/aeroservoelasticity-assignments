function [flutter_solution,fig] = flutterStandard(flutter_system,data,chord,plotFlag)
%FLUTTER STANDARD - Compute standard flutter analysis
%   Function computing flutter analysis solving using MATLAB function 'eig'
%   the eigenvalue problem at different flow field velocity
%
%   SYNTAX:
%       [flutter_solution,fig] = flutterStandard(flutter_system,data,chord,plotFlag)
%
%   INPUT:
%       flutter_system, struct: state-space system struct of the considered
%                               model containin even the 'flow' struct
%       data,           struct: data provided in 'LabData_Session1.mat'
%       chord,          double: chord value [m] of the considered wing
%       plotFlag(*)       bool: set if the figures must be generated
%
%   OUTPUT:
%       flutter_solution, struct: flutter analysis results
%       fig,              figure: plot of flutter results
%
%   OPTIONAL INPUT:
%       plotFlag: set by default as false
%
%

    % Optional input
    if nargin < 4
        plotFlag = false;
    end

    % Initialize eigenvalues and eigenvector matrix
    flutter_solution.eigsVal = zeros(length(data.Mhh)*2+length(data.Aaero), length(flutter_system.flow.v_vec));
    
    for i = 1:length(flutter_system.flow.v_vec)
        v_inf = flutter_system.flow.v_vec(i);    
        [sys] = buildAeroSystem(flutter_system, data, chord, flutter_system.flow.rho, v_inf);
        flutter_solution.eigsVal(:,i) = eig(sys.A, sys.V);   
    end

    % Compute damping ratios and eigenfrequencies
    flutter_solution.dampRatio = -(real(flutter_solution.eigsVal)./abs(flutter_solution.eigsVal));
    flutter_solution.freq_Hz   = imag(flutter_solution.eigsVal)./(2*pi);
    
    % Plot of the results just computed
    if plotFlag == true
        fig.rootlocus = figure(Name='Root-locus');
        hold on; grid on; axis padded; box on;
        for i = 1:length(flutter_system.flow.v_vec)
            plot(real(flutter_solution.eigsVal(:,i)/2/pi), imag(flutter_solution.eigsVal(:,i))/2/pi,'o');
        end
        xlabel('Re($\lambda$)');       ylabel('Im($\lambda$)');
        xlim([-1 0.5]);     ylim([0 30]);     

        fig.damping = figure(Name='Damping ratio');
        hold on;  grid minor;  axis padded;  box on;       
        for i = 1:length(sys.A)
            plot(flutter_system.flow.v_vec, flutter_solution.dampRatio(i,:),'-*');
        end
        xlabel('Velocity [m/s]');  ylabel('Damping');
    
        fig.freq = figure(Name='Eigenfrequencies');
        hold on;  grid minor;  axis padded;  box on;
        for i = 1:length(sys.A)
            plot(flutter_system.flow.v_vec, flutter_solution.freq_Hz(i,:),'-*');
        end
        xlabel('Velocity [m/s]');  ylabel('Frequency Hz');
    else
        fig.rootlocus = [];
        fig.damping   = [];
        fig.freq      = [];
    end

end


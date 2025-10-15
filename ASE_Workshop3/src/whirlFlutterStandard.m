function [flutter,fig] = whirlFlutterStandard(A,uVec,plotFlag)
%WHIRL FLUTTER STANDARD - Compute standard whirl flutter analysis
%   Function computing flutter analysis solving using MATLAB function 'eig'
%   the eigenvalue problem at different flow field velocity
%
%   SYNTAX:
%       [flutter_solution,fig] = whirlFlutterStandard(A,uVec,plotFlag)
%
%   INPUT:
%       A,          double: state-space matrix of the system where 3rd
%                           dimension contains the data of the different flow velocity
%       uVec,       double: vector containing the flow velocity
%       plotFlag(*)   bool: set if the figures must be generated
%
%   OUTPUT:
%       flutter,   struct: flutter analysis results
%       fig,       figure: plot of flutter results
%
%   OPTIONAL INPUT:
%       plotFlag: set by default as false
%
%

    % Optional input
    if nargin < 3
        plotFlag = false;
    end

    % Initialize eigenvalues and eigenvector matrix
    flutter.eigsVal = zeros(size(A,1),size(A,3));
    for i = 1:size(A,3)
        flutter.eigsVal(:,i) = eig(A(:,:,i));
        imagPart = abs(imag(flutter.eigsVal(:,i)));
        [~, sortIdx] = sort(imagPart);
        flutter.eigsVal(:,i) = flutter.eigsVal(sortIdx,i);
    end
    flutter.eigsVal = flutter.eigsVal(1:2:end,:);

    % Compute damping ratios and eigenfrequencies
    flutter.dampRatio = -(real(flutter.eigsVal)./abs(flutter.eigsVal));
    flutter.freq_Hz   = imag(flutter.eigsVal)./(2*pi);
    
    % Plot of the results just computed
    if plotFlag == true
        fig.rootlocus = figure(Name='Root-locus');
        hold on; grid on; axis padded; box on;
        for i = 1:size(A,3)
            plot(real(flutter.eigsVal(:,i)/2/pi), imag(flutter.eigsVal(:,i))/2/pi,'o');
        end
        xlabel('Re($\lambda$)');       ylabel('Im($\lambda$)');     

        fig.damping = figure(Name='Damping ratio');
        hold on;  grid minor;  axis padded;  box on;       
        for i = 1:size(A,1)/2
            plot(uVec, flutter.dampRatio(i,:),'-*');
        end
        xlabel('Velocity [m/s]');  ylabel('Damping');
    
        fig.freq = figure(Name='Eigenfrequencies');
        hold on;  grid minor;  axis padded;  box on;
        for i = 1:size(A,1)/2
            plot(uVec, flutter.freq_Hz(i,:),'-*');
        end
        xlabel('Velocity [m/s]');  ylabel('Frequency Hz');
    else
        fig.rootlocus = [];
        fig.damping   = [];
        fig.freq      = [];
    end

end


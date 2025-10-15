function [fig,dataTF] = ssBodePlot(transferFun,plotFlag)
%SSBODEPLOT Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 2
        plotFlag = false;
    end
    
    [mag, phase, w] = bode(transferFun);
    dataTF = struct;
        dataTF.mag = squeeze(mag);         % Magnitude
        dataTF.phase = squeeze(phase);     % Phase
        dataTF.w = squeeze(w);             % Frequency (rad/s)
        dataTF.f = w / (2 * pi);           % Frequency (Hz)
    
    % Plotta il diagramma di Bode
    if plotFlag == true
        fig = figure(Name='Bode plot \delta_0');
        title('Magnitude');
        
        % Magnitude
        subplot(2, 1, 1);
        semilogx(dataTF.f, 20*log10(dataTF.mag));
        grid minor;  box on;
        xlabel('Frequenza (Hz)');   xlim([10^-1,50]);
        ylabel('Modulo (dB)');
        
        % Phase
        subplot(2, 1, 2);
        semilogx(dataTF.f, dataTF.phase);
        grid minor;  box on;
        xlabel('Frequenza (Hz)');   xlim([10^-1,50]);
        ylabel('Fase');
    else
        fig = [];
    end

end


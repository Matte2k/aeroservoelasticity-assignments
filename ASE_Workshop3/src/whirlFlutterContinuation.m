function [flutter,fig] = whirlFlutterContinuation(A,uVec,plotFlag)
%WHIRL FLUTTER CONTINUATION - compute whirl flutter analysis with continuation approach
%   Function computing flutter analysis solving the eigenvalue problem 
%   at different flow field velocity using continuation approach
%
%   SYNTAX:
%       [flutter,fig] = whirlFlutterContinuation(A,uVec,plotFlag)
%
%   INPUT:
%       A,          double: state-space matrix of the system where 3rd
%                           dimension contains the data of the different flow velocity
%       uVec,       double: vector containing the flow velocity
%       plotFlag(*)   bool: set if the figures must be generated
%
%   OUTPUT:
%       flutter,  struct: flutter analysis results
%       fig,      figure: plot of flutter results
%
%   OPTIONAL INPUT:
%       plotFlag: set by default as false
%
%
%

    % Optional input
    if nargin < 3
        plotFlag = false;
    end

    % Initialize eigenvalues and eigenvector matrix
    %flutter.eigsVal = zeros(size(A,1),size(A,3));
    %flutter.eigsVec = zeros(size(A,1),size(A,2),size(A,3));
    
        % Compute initial eigenvalues and eigenvector
        A_balanced = balance(A(:,:,1));
        [V,E] = eig(A_balanced);
        E = diag(E);

        % Pick only one of the complex conjugate couples
        E = E(1:2:end);
        V = V(:,1:2:end);

        % % Sorting eigenvalues from low freq to high freq
        % imagPart = abs(imag(E));
        % [~, sortIdx] = sort(imagPart);
        % E = E(sortIdx);
        % V = V(:,sortIdx);
        % % NOTE: they are already sorted later on

    % Pre-alloca per accoppiare i vettori propri
    eig_vals_history = zeros(length(E), length(uVec));
    eig_vals_history(:,1) = E;                              % Primo set di autovalori
    
    eig_vecs_history = V;                                   % Vettori propri del passo precedente
    
    for k = 2:length(uVec)
        % Calcola matrice a stati per la velocità corrente
        A_k = balance(A(:,:,k));
        [eig_vecs_k, eig_vals_k] = eig(A_k);
        eig_vals_k = diag(eig_vals_k);
        
        eig_vals_k = eig_vals_k(1:2:end);
        eig_vecs_k = eig_vecs_k(:,1:2:end);
        
        % Accoppia i modi basandoti sul prodotto scalare
        match_indices = zeros(length(E(:,1)), 1);
        for i = 1:length(E(:,1))
            % Calcola prodotto scalare con i vettori propri precedenti
            dot_products = abs(eig_vecs_history(:,i)' * eig_vecs_k);
            [~, match_idx] = max(dot_products);             % Trova l'accoppiamento migliore
            match_indices(i) = match_idx;
        end
        
        % Riordina autovalori e vettori propri in base all'accoppiamento
        eig_vals_history(:,k) = eig_vals_k(match_indices);
        eig_vecs_history = eig_vecs_k(:, match_indices);
    end

    flutter.eigsVal = zeros(size(A,1)/2,size(A,3));
    flutter.eigsVal = eig_vals_history;    
    %flutter.eigsVal = flutter.eigsVal_cont;

    % Compute damping ratios and eigenfrequencies
    flutter.dampRatio = -(real(flutter.eigsVal)./abs(flutter.eigsVal));
    flutter.freq_Hz   = imag(flutter.eigsVal)./(2*pi);

    % Resort the vectors to have modes in right order
    flutter.dampRatio = flip(flutter.dampRatio);
    flutter.freq_Hz   = flip(flutter.freq_Hz);

    % Manually tuning from 'flutter.freq_Hz'
    mode2plot = 5:10;   % exclude the rigid modes
    lgdCell = cell(1,length(mode2plot));
    
    % Plot of the results just computed
    if plotFlag == true
        fig.rootlocus = figure(Name='Root-locus');
        hold on; grid on; axis padded; box on;
        for i = 1:size(A,3)
            plot(real(flutter.eigsVal(:,i)/2/pi), imag(flutter.eigsVal(:,i))/2/pi,'o');
        end
        xlabel('Re($\lambda$)');       ylabel('Im($\lambda$)');     
        
        % fig.damping = figure(Name='Damping ratio');
        % hold on;  grid minor;  axis padded;  box on;       
        % for i = 1:length(mode2plot)
        %     plot(uVec, flutter.dampRatio(mode2plot(i),:),'-*');
        %     lgdCell{i} = sprintf('Mode %.0f',mode2plot(i));
        % end
        % legend(lgdCell);
        % xlabel('Velocity [m/s]');  ylabel('Damping');
        % 
        % fig.freq = figure(Name='Eigenfrequencies');
        % hold on;  grid minor;  axis padded;  box on;
        % for i = 1:length(mode2plot)
        %     plot(uVec, flutter.freq_Hz(mode2plot(i),:),'-*');         
        % end
        % legend(lgdCell);
        % xlabel('Velocity [m/s]');  ylabel('Frequency Hz');

        fig.vgvf = figure(Name='Flutter analysis');
        tiledlayout(1,2);
        nexttile
            hold on;  grid minor;  axis padded;  box on;
            for i = 1:length(mode2plot)
                plot(uVec, flutter.freq_Hz(mode2plot(i),:),'-*');         
            end
            xlabel('Velocity [m/s]');  ylabel('Frequency [Hz]');
        nexttile
            hold on;  grid minor;  axis padded;  box on;       
            for i = 1:length(mode2plot)
                plot(uVec, flutter.dampRatio(mode2plot(i),:),'-*');
                lgdCell{i} = sprintf('Mode %.0f',mode2plot(i));
            end
            xlabel('Velocity [m/s]');  ylabel('Damping');
        lgd = legend(lgdCell);
        lgd.Layout.Tile = 'south';
        lgd.Orientation = 'horizontal';
        
    else
        fig.rootlocus = [];
        % fig.damping   = [];
        % fig.freq      = [];
        fig.vgvf      = [];
    end

end


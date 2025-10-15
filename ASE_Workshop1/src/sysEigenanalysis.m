function [sys] = sysEigenanalysis(sys)
%SYSTEM EIGENANALYSIS - Compute eigenanalysis of system struct
%   Function computing the eigenanalysis of a given system struct. Based on
%   the info contained in 'sys.type' the function recognize if the
%   considered system must be analyzed in state space form or classical
%   dynimaic form.
%
%   SYNTAX:
%       [sys] = sysEigenanalysis(sys)
%
%   INPUT:
%       sys, struct: system struct that may contains the following fields:
%                    > type: can be 'dynamic' or 'state_space'
%                    > K,M: dynamic system matrix and stiffness matricies
%                    > A: state space system matrix of the system
%   
%   OUTPUT:
%       sys, struct: update strcut with computed eigenvalues and eigenvector
%                    organized in fields as follows:
%                    > eigsVect: matrix with one eigenvector on each column
%                    > eigsVal: vector with eigenvalues associated to eigsVect
%                    > freq_rad: eigenfrequencies in rad/s
%                    > freq_Hz: eigenfrequencies in Hz
%                    > xx_sorted: sorted quantity in 'ascend' order
%
%

    if isequal(sys.type,'dynamic')
        [sys.eigsVect, sys.eigsVal] = eig(sys.K, sys.M);
        sys.eigsVal = diag(sys.eigsVal);
        [sys.eigsVal_sorted, sys.sorting_index] = sort(sys.eigsVal);
        sys.eigsVect_sorted = sys.eigsVect(:,sys.sorting_index);
        sys.freq_rad_sorted = sqrt(sys.eigsVal_sorted);
        sys.freq_Hz_sorted  = sys.freq_rad_sorted/2/pi;
        
    elseif isequal(sys.type,'state_space')
        [sys.eigsVect, sys.eigsVal] = eig(sys.A);
        sys.eigsVal  = diag(sys.eigsVal);
        sys.freq_rad = imag(sys.eigsVal);
        sys.freq_Hz  = sys.freq_rad/2/pi;
        [sys.freq_rad_sorted,sys.sorting_index]=sort(abs(sys.freq_rad));
        sys.eigsVect_sorted = sys.eigsVect(:,sys.sorting_index);
        sys.freq_Hz_sorted  = sys.freq_Hz(sys.sorting_index);


    end
 
end


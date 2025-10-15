function [wFlutter] = computeWhirlflutter(airframe,speedData,wFlutter)
%COMPUTEWHIRLFLUTTER Summary of this function goes here
%   Detailed explanation goes here

    constrMat = wFlutter.constrMat;

    flutterVel_size = size(speedData.M,3);
    flutterSys_size = size(constrMat,2)*2;

    flutter_M = zeros(size(constrMat,1), size(constrMat,1));
    flutter_C = zeros(size(constrMat,1), size(constrMat,1));
    flutter_K = zeros(size(constrMat,1), size(constrMat,1));

    wFlutter.unstable = zeros(1,flutterVel_size);
    wFlutter.eigsVal  = zeros(flutterSys_size, flutterVel_size);
    wFlutter.eigsVec  = zeros(flutterSys_size, flutterSys_size, flutterVel_size);
    
    for i = 1:flutterVel_size
        flutter_M(1:14,1:14) = airframe.M;
        flutter_M(15:end,15:end) = speedData.M(:,:,i);
    
        flutter_C(1:14,1:14) = airframe.C;
        flutter_C(15:end,15:end) = speedData.C(:,:,i);
    
        flutter_K(1:14,1:14) = airframe.K;
        flutter_K(15:end,15:end) = speedData.K(:,:,i);
    
        Mnew_T = constrMat' * flutter_M * constrMat;
        Cnew_T = constrMat' * flutter_C * constrMat;
        Knew_T = constrMat' * flutter_K * constrMat;
    
        flutter_Minv = inv(Mnew_T);
        flutter_A    = [zeros(length(Mnew_T)),    eye(length(Mnew_T));
                        -flutter_Minv*Knew_T,   -flutter_Minv*Cnew_T];   % A

        [wFlutter.eigsVec(:,:,i),eigsValMatrix] = eig(flutter_A);
        wFlutter.eigsVal(:,i) = diag(eigsValMatrix);
        
        if any(real(wFlutter.eigsVal(:,i)) >= 0)
            wFlutter.unstable(i) = 1;
        end
    end


end


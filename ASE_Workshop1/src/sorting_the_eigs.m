function [lambda0_mat,q0_mat] = sorting_the_eigs(lambda0_mat,q0_mat)   

    [l,index] = sort(abs(imag(diag(lambda0_mat))));
    
    index = [index(43:end);index(1:42)];  % put aero mode at the end
    lambda0 = diag(lambda0_mat);

    lambda0_mat = diag(lambda0(index));
    q0_mat = q0_mat(index,:);

end
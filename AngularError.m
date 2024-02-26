function err = AngularError(R_est, R_gt)
    
    % Find the angle (deg) between two rotation matrices.
    
    err = abs(acosd((trace(R_est'*R_gt)-1)/2));

end
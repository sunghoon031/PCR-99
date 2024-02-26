function [R, t] = Rt_from_N_points(A, B, scale)

    % Compute rotation (R) and translation (t) from n pairs of
    % corresponding points.
    % A is [3 x n] ground truth points
    % B is [3 x n] estimated points
    
    centroid_A = mean(A,2);
    centroid_B = mean(B,2);
    
    A = A - centroid_A;
    B = B - centroid_B;

    [UU,~,VV] = svd(A*B');

    R = VV*diag([1,1,sign(det(VV*UU'))])*UU';
  
    t = centroid_B - scale*R*centroid_A;

end
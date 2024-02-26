function [s, R, t] = sRt_from_N_points(A, B)

    % Compute scale (s), rotation (R) and translation (t) from n pairs of
    % corresponding points.
    % A is [3 x n] ground truth points
    % B is [3 x n] estimated points
    
    centroid_A = mean(A,2);
    centroid_B = mean(B,2);
    
    A = A - centroid_A;
    B = B - centroid_B;

    [UU,~,VV] = svd(A*B');

    R = VV*diag([1,1,sign(det(VV*UU'))])*UU';

    
    den = 0;
    num = 0;
    for i=1:size(A,2)
        % Each increment is [3x1]'*[3x1]=[1x1]
        num = num + B(:,i)'*R*A(:,i);
        den = den + A(:,i)'*A(:,i);
    end

    s = num/den;
    
    t = centroid_B - s*R*centroid_A;

end
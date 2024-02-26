function [s,R,t] = sRt_from_3points(A,B)

    % Compute scale (s), rotation (R) and translation (t) from 3 pairs of
    % corresponding points.
    % A is [3 x n] ground truth points (n = 3)
    % B is [3 x n] estimated points (n = 3)
    
    centroid_A = mean(A,2);
    centroid_B = mean(B,2);
    
    A = A - centroid_A;
    B = B - centroid_B;
    

    v12 = B(:,2)-B(:,1);
    x_axis = v12/norm(v12);
    v13 = B(:,3)-B(:,1);
    v23 = cross(v12, v13);
    y_axis = v23/norm(v23);
    z_axis = cross(x_axis, y_axis);
    z_axis = z_axis/norm(z_axis);


    v12 = A(:,2)-A(:,1);
    x_axis_ = v12/norm(v12);
    v13 = A(:,3)-A(:,1);
    v23 = cross(v12, v13);
    y_axis_ = v23/norm(v23);
    z_axis_ = cross(x_axis_, y_axis_);
    z_axis_ = z_axis_/norm(z_axis_);
    
    R=[x_axis,y_axis,z_axis]*[x_axis_,y_axis_,z_axis_]';
    

    
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
function [R_gt, t_gt, xyz_gt, xyz_est] = SimulateknownScale(n, sigma, outlier_ratio, s_gt)

    % Generate two point clouds (ground truth and estimation).
    
    xyz_gt = rand(3, n);

    x_range = max(xyz_gt(1,:)) - min(xyz_gt(1,:));
    y_range = max(xyz_gt(2,:)) - min(xyz_gt(2,:));
    z_range = max(xyz_gt(3,:)) - min(xyz_gt(3,:));
    max_range = max([x_range, y_range, z_range]);
    xyz_gt = xyz_gt/max_range;

    x_mid = 0.5*(max(xyz_gt(1,:)) + min(xyz_gt(1,:)));
    y_mid = 0.5*(max(xyz_gt(2,:)) + min(xyz_gt(2,:)));
    z_mid = 0.5*(max(xyz_gt(3,:)) + min(xyz_gt(3,:)));
    xyz_gt = xyz_gt - [x_mid;y_mid;z_mid];      


    R_gt = RandomRotationMatrix;
    t_gt = rand(3,1)-0.5;
    t_gt = t_gt/norm(t_gt)*rand(1)*3;
    xyz_est = s_gt*R_gt*xyz_gt + t_gt;

    xyz_est = xyz_est + normrnd(0, sigma, [3, n]);

    xyz_outliers = zeros(3, n);
    c = 0;
    while (c < n)
        xyz_outlier = s_gt*sqrt(3)*(rand(3, 1)-0.5);
        if (norm(xyz_outlier) < s_gt*sqrt(3)/2)
            c = c + 1;
            xyz_outliers(:,c) = xyz_outlier;
        end
    end
    xyz_outliers = xyz_outliers + t_gt;

    n_outliers = round(n*outlier_ratio);
    idx_outliers = randperm(n, n_outliers);
    xyz_est(:,idx_outliers) = xyz_outliers(:,idx_outliers);

end
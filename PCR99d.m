function [R_final, t_final] = PCR99d(xyz_gt,xyz_est, sigma, thr1, thr2, n_hypo, scale)
 
    
    n = size(xyz_gt, 2);
        
    log_ratio_diff_mat = nan(n,n);

    log_ratio_gt = log(scale);
    
    idx_inliers = [];

    
    
    Sx_mat = zeros(5, n_hypo);
    Sy_mat = zeros(5, n_hypo);
    Sz_mat = zeros(5, n_hypo);
        
    max_max_nInliers = 0;
    c = 0;
    
    for it = 1:1e+9

        ijk = randperm(n, 3);
        i = ijk(1);
        j = ijk(2);
        k = ijk(3);

        if (isnan(log_ratio_diff_mat(i,j)))
            d_gt_ij = norm(xyz_gt(:,i) - xyz_gt(:,j));
            d_est_ij = norm(xyz_est(:,i) - xyz_est(:,j));

            log_ratio_diff_mat(i,j) = abs(log(d_est_ij/d_gt_ij)-log_ratio_gt);
            log_ratio_diff_mat(j,i) = log_ratio_diff_mat(i,j);
        end
        if (isnan(log_ratio_diff_mat(j,k)))
            d_gt_jk = norm(xyz_gt(:,j) - xyz_gt(:,k));
            d_est_jk = norm(xyz_est(:,j) - xyz_est(:,k));

            log_ratio_diff_mat(j,k) = abs(log(d_est_jk/d_gt_jk)-log_ratio_gt);
            log_ratio_diff_mat(k,j) = log_ratio_diff_mat(j,k);
        end
        if (isnan(log_ratio_diff_mat(k,i)))
            d_gt_ki = norm(xyz_gt(:,k) - xyz_gt(:,i));
            d_est_ki = norm(xyz_est(:,k) - xyz_est(:,i));

            log_ratio_diff_mat(k,i) = abs(log(d_est_ki/d_gt_ki)-log_ratio_gt);
            log_ratio_diff_mat(i,k) = log_ratio_diff_mat(k,i);
        end

          
        if (log_ratio_diff_mat(i,j)>thr1 || log_ratio_diff_mat(j,k)>thr1 || log_ratio_diff_mat(k,i)>thr1)
            continue;
        end

        c = c + 1;

        A = xyz_gt(:,[i, j, k]);
        B = xyz_est(:,[i, j, k]);


        [R, t] = Rt_from_3points(A,B, scale);

        Sx_mat(:,c) = [scale*R(1,1); scale*R(1,2); scale*R(1,3); t(1); 1];
        Sy_mat(:,c) = [scale*R(2,1); scale*R(2,2); scale*R(2,3); t(2); 1];
        Sz_mat(:,c) = [scale*R(3,1); scale*R(3,2); scale*R(3,3); t(3); 1];

        if (c==n_hypo)
            %disp('evaluate!')
            c = 0;
            Ex = [-xyz_gt' -ones(n,1) xyz_est(1,:)']*Sx_mat;
            Ey = [-xyz_gt' -ones(n,1) xyz_est(2,:)']*Sy_mat;
            Ez = [-xyz_gt' -ones(n,1) xyz_est(3,:)']*Sz_mat;

            E = Ex.^2 + Ey.^2 + Ez.^2;

            e_thr = (sigma*thr2)^2;

            nInliers = E<=e_thr;
            nInliers = sum(nInliers,1);

            [max_nInliers, idx] = max(nInliers);

            if (max_max_nInliers < max_nInliers)

                max_max_nInliers = max_nInliers;

                idx_inliers = find(E(:,idx) <= e_thr);

                if (max_max_nInliers >= max(9, round(n*0.009)))
                    break;
                end
            end
        end

    end


    if (isempty(idx_inliers))
        R_final = nan;
        t_final = nan;
        return;
    end
    
    A = xyz_gt(:,idx_inliers);
    B = xyz_est(:,idx_inliers);


    [R_final, t_final] = Rt_from_N_points(A,B, scale);
end

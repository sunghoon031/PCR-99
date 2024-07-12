function [R_final, t_final] = PCR99c(xyz_gt,xyz_est, sigma, thr1, thr2, n_hypo, scale)
 
  
    
    idx_inliers = [];
    
    
    % 1. Construct the log ratio matrix:

    n = size(xyz_gt, 2);
    
    large_error_mat = zeros(n,n);

    log_ratio_gt = log(scale);
    costs = nan(1,n);
        
    for i = 1:n
        p_gt_i = xyz_gt(:,i);
        p_est_i = xyz_est(:,i);
        
        d_gt = xyz_gt - p_gt_i;
        d_est = xyz_est - p_est_i;
        
        d_gt = sum(d_gt.^2, 1);
        d_est = sum(d_est.^2, 1);
        
        lr_mat = 0.5*log(d_est./d_gt);
        
        % 2. Score the correspondences:
        r = abs(lr_mat - log_ratio_gt);
        large_error_bool = r > thr1;
        r(large_error_bool) = thr1;
        costs(i) = nansum(r);
        
        large_error_mat(i, large_error_bool) = 1;
        
    end

    % 3. Sort the correspondences:
    [~, sort_idx] = sort(costs);
    
    xyz_est = xyz_est(:, sort_idx);
    xyz_gt = xyz_gt(:, sort_idx);

    % 4. Evaluate the samples sequentially, prioritizing those that have
    % smaller total ranking numbers:
    
    Sx_mat = zeros(5, n_hypo);
    Sy_mat = zeros(5, n_hypo);
    Sz_mat = zeros(5, n_hypo);
    
    break_loop = false;
    
    max_max_nInliers = 0;
    c = 0;
    for s = 6:n+(n-1)+(n-2)

        if (break_loop)
            break;
        end

        i_min = max(1, s-n-(n-1));
        i_max =  floor((s-3)/3);

        for i = i_min:i_max

            if (break_loop)
                break;
            end

            j_min = max(i+1, s-i-n);
            j_max = floor((s-i-1)/2);

            for j = j_min:j_max
                
                k = s - i - j;

                i_old = sort_idx(i);
                j_old = sort_idx(j);
                k_old = sort_idx(k);
                
                % Prescreening test:
                if (large_error_mat(i_old,j_old) ...
                        || large_error_mat(j_old,k_old) ...
                        || large_error_mat(k_old,i_old))
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
                            break_loop = true;
                            break;
                        end
                    end
                end
                
            end
        end
    end


    if (isempty(idx_inliers))
        R_final = nan;
        t_final = nan;
        return;
    end
    
    % 5. Recompute the transformation using the inlier set:

    A = xyz_gt(:,idx_inliers);
    B = xyz_est(:,idx_inliers);

    [R_final, t_final] = Rt_from_N_points(A,B, scale);
end

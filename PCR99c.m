function [R_final, t_final] = PCR99c(xyz_gt,xyz_est, sigma, thr1, thr2, n_hypo, scale)
 
  
    % 1. Construct the log ratio matrix:

    n = size(xyz_gt, 2);
    
    log_ratio_mat = nan(n,n);
    log_ratio_fail_mat = zeros(n,n);
    log_ratio_gt = log(scale);
    
    idx_inliers = [];
    
    for i = 1:n-1
        p_gt_i = xyz_gt(:,i);
        p_est_i = xyz_est(:,i);

        for j = i+1:n
            p_gt_j = xyz_gt(:,j);
            p_est_j = xyz_est(:,j);

            v_gt_ij = p_gt_i - p_gt_j;
            v_est_ij = p_est_i - p_est_j;

            d_gt_ij = norm(v_gt_ij);
            d_est_ij = norm(v_est_ij);

            log_ratio_mat(i,j) = log(d_est_ij/d_gt_ij);
            log_ratio_mat(j,i) = log_ratio_mat(i,j);
            
            log_ratio_diff = abs(log_ratio_mat(i,j) - log_ratio_gt);
            
            if (log_ratio_diff > thr1)
                log_ratio_fail_mat(i,j) = 1;
                log_ratio_fail_mat(j,i) = 1;
            end

        end
    end

    
    % 2. Score the correspondences:

    costs = nan(1,n);
    for i = 1:n
        lr_mat = log_ratio_mat(i,:);
        r = abs(lr_mat - log_ratio_gt);
        r(r>thr1) = thr1;
        costs(i) = nansum(r);
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
                if (log_ratio_fail_mat(i_old,j_old) ...
                        || log_ratio_fail_mat(j_old,k_old) ...
                        || log_ratio_fail_mat(k_old,i_old))
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
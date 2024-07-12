function [s_final, R_final, t_final] = PCR99a(xyz_gt,xyz_est, sigma, thr1, thr2, n_hypo)
 


    % 1. Construct the log ratio matrix:
    
    n = size(xyz_est, 2);
    
    log_ratio_mat = nan(n,n);
    
    idx_inliers = [];
        
    for i = 1:n
        p_gt_i = xyz_gt(:,i);
        p_est_i = xyz_est(:,i);

        d_gt = xyz_gt - p_gt_i;
        d_est = xyz_est - p_est_i;
        
        d_gt = sum(d_gt.^2, 1);
        d_est = sum(d_est.^2, 1);
        
        log_ratio_mat(i,:) = 0.5*log(d_est./d_gt);
    end


    % 2. Score the correspondences:
    
    min_costs = nan(1,n);
    for i = 1:n
        lr_mat = log_ratio_mat(i,:);

        max_lr = max(lr_mat);
        min_lr = min(lr_mat);
        lr_range = max_lr - min_lr;
        lr_step = lr_range/max(1, round(lr_range/0.1));
        lr_candidates = min_lr:lr_step:max_lr;

        min_cost = inf;
        for lr_c = lr_candidates
            r = abs(lr_mat - lr_c);
            r(r>thr1) = thr1;
            cost = nansum(r);
            if (cost < min_cost)
                min_cost = cost;
            end
        end
        min_costs(i) = min_cost;
    end
 

    % 3. Sort the correspondences:
    [~, sort_idx] = sort(min_costs);
    
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
                log_ratio_ij = log_ratio_mat(i_old,j_old);
                log_ratio_jk = log_ratio_mat(j_old,k_old);
                log_ratio_ki = log_ratio_mat(k_old,i_old);

                e1 = abs(log_ratio_ij - log_ratio_jk);
                e2 = abs(log_ratio_jk - log_ratio_ki);
                e3 = abs(log_ratio_ki - log_ratio_ij);

                if (e1 > thr1 || e2 > thr1 || e3 > thr1)
                    continue;
                end

                c = c + 1;

                A = xyz_gt(:,[i, j, k]);
                B = xyz_est(:,[i, j, k]);

      
                [scale, R, t] = sRt_from_3points(A,B);

                Sx_mat(:,c) = [scale*R(1,1); scale*R(1,2); scale*R(1,3); t(1); 1];
                Sy_mat(:,c) = [scale*R(2,1); scale*R(2,2); scale*R(2,3); t(2); 1];
                Sz_mat(:,c) = [scale*R(3,1); scale*R(3,2); scale*R(3,3); t(3); 1];

                if (c==n_hypo)
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
        s_final = nan;
        R_final = nan;
        t_final = nan;
        return;
    end

    
    % 5. Recompute the transformation using the inlier set:
    
    A = xyz_gt(:,idx_inliers);
    B = xyz_est(:,idx_inliers);

    [s_final, R_final, t_final] = sRt_from_N_points(A,B);
end

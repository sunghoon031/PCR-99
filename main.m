clear all; close all; clc;

% 1. PCR99a: Unknown scale, With sample ordering
% 2. PCR99b: Unknown scale, With random sampling
% 3. PCR99c: Known scale, With sample ordering
% 4. PCR99d: Known scale, With random sampling (NOT RECOMMENDED!)


%% Create an unknown-scale problem

n = 1000;
sigma = 0.01;
outlier_ratio = 0.99;

[s_gt, R_gt, t_gt, xyz_gt, xyz_est] = SimulateUnknownScale(n, sigma, outlier_ratio);

disp('Unknown-scale problem has been created!')

%% Solve the unknown-scale problem

thr1 = 0.1;
thr2 = 5;
n_hypo = 1000;

rng default
tic
[s_final, R_final, t_final] = PCR99a(xyz_gt, xyz_est, sigma, thr1, thr2, n_hypo);
time_elapsed = toc;
R_error=AngularError(R_gt,R_final);
t_error=norm(t_final - t_gt);
s_error=abs(s_final - s_gt);

disp(['  PCR99a (with sample ordering): R err = ', num2str(R_error), ' t err = ', num2str(t_error), ' s err = ', num2str(s_error), ' time = ' num2str(time_elapsed)])
 

rng default
tic
[s_final, R_final, t_final] = PCR99b(xyz_gt, xyz_est, sigma, thr1, thr2, n_hypo);
time_elapsed = toc;
R_error=AngularError(R_gt,R_final);
t_error=norm(t_final - t_gt);
s_error=abs(s_final - s_gt);
disp(['  PCR99b (with random sampling): R err = ', num2str(R_error), ' t err = ', num2str(t_error), ' s err = ', num2str(s_error), ' time = ' num2str(time_elapsed)])


%% Create a known-scale problem


n = 1000;
sigma = 0.01;
outlier_ratio = 0.99;
s_gt = 1+4*rand(1);

[R_gt, t_gt, xyz_gt, xyz_est] = SimulateknownScale(n, sigma, outlier_ratio, s_gt);

disp('Known-scale problem has been created!')

            
%% Solve the known-scale problem
     
rng default
tic
[R_final, t_final] = PCR99c(xyz_gt, xyz_est, sigma, thr1, thr2, n_hypo, s_gt);
time_elapsed = toc;
R_error=AngularError(R_gt,R_final);
t_error=norm(t_final - t_gt);
disp(['  PCR99c (with sample ordering): R err = ', num2str(R_error), ' t err = ', num2str(t_error), ' time = ' num2str(time_elapsed)])


rng default
tic
[R_final, t_final] = PCR99d(xyz_gt, xyz_est, sigma, thr1, thr2, n_hypo, s_gt);
time_elapsed = toc;
R_error=AngularError(R_gt,R_final);
t_error=norm(t_final - t_gt);
disp(['  PCR99d (with random sampling): R err = ', num2str(R_error), ' t err = ', num2str(t_error), ' time = ' num2str(time_elapsed)])







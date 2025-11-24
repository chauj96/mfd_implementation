% mean_logR = [-0.000000000004, -0.0000000005, -0.0000000030, ...
%               -0.0000000332, -0.0000005192, -0.0000014240, -0.0000030246];
% flux_error_transition = [8.122e-08, 4.062e-07, 7.313e-07, ...
%                          1.608e-05, 1.324e-05, 3.049e-06, 1.696e-06];
% distortion = [0.01, 0.05, 0.09, 0.2, 0.5, 0.7, 0.9]; % optional labels
% 
% figure('Position',[200 150 750 500]);
% loglog(abs(mean_logR), flux_error_transition, '-o', ...
%        'LineWidth',1.8, 'MarkerSize',7, 'MarkerFaceColor','auto');
% xlabel('|Mean log_{10}(R)|','FontSize',13);
% ylabel('Flux L2 error (transition point)','FontSize',13);
% title('Case 1 — Flux Error vs mean log_{10}(R)','FontSize',14);
% grid on;
% 
% text(abs(mean_logR)*1.1, flux_error_transition, ...
%     compose('%.2f', distortion), 'FontSize',11, 'Color',[0.3 0.3 0.3]);
% 
% % CASE 2
% mean_logR = [0.00002136, 0.00002201, 0.00002358, ...
%               0.00002434, 0.00002652, 0.00002971];
% flux_error_transition = [4.559e-07, 1.398e-04, 3.704e-03, ...
%                          3.280e-03, 6.976e-05, 1.234e-03];
% distortion = [0.02, 0.10, 0.30, 0.40, 0.70, 1.20];
% 
% r = 200 / 10000;
% 
% mean_logR_norm = mean_logR * r;
% 
% figure('Position',[250 150 780 520]);
% loglog(mean_logR_norm, flux_error_transition, '-o', ...
%        'LineWidth',1.8, 'MarkerSize',7, 'MarkerFaceColor','auto');
% xlabel('Normalized mean log_{10}(R)','FontSize',13);
% ylabel('Flux L2 error (transition point)','FontSize',13);
% title('Case 2 — Flux Error vs mean log_{10}(R)','FontSize',14);
% grid on;
% 
% text(mean_logR_norm * 1.01, flux_error_transition, ...
%      compose('%.2f', distortion), 'FontSize',11, 'Color',[0.3 0.3 0.3]);
% 
% %% CASE 3
% distortion = [0.02, 0.10, 0.30, 0.40, 0.70, 1.20];
% mean_logR = [0.060979213847, 0.300032998181, 0.837082048491, ...
%               1.065411683421, 1.618930282537, 2.269711488529];
% flux_error_transition = [1.401e-07, 2.679e-05, 2.451e-04, ...
%                          2.373e-04, 2.675e-05, 3.169e-04];
% 
% r = 200 / 10000;
% mean_logR_norm = mean_logR * r;
% 
% figure('Position',[250 150 780 520]);
% loglog(mean_logR_norm, flux_error_transition, '-o', ...
%        'LineWidth',1.8, 'MarkerSize',7, 'MarkerFaceColor','auto');
% xlabel('Normalized mean log_{10}(R)','FontSize',13);
% ylabel('Flux L2 error (transition point)','FontSize',13);
% title('Case 3 — Flux Error vs mean log_{10}(R)','FontSize',14);
% grid on;
% 
% text(mean_logR_norm*1.05, flux_error_transition, ...
%      compose('%.2f', distortion), 'FontSize',11, 'Color',[0.3 0.3 0.3]);


%%
% Case 1
distortion = [0.01, 0.05, 0.09, 0.20, 0.50, 0.70, 0.90];

err_1e4 = [2.996e-5, 1.498e-4, 2.699e-4, 6.022e-4, 4.787e-4, 4.213e-4, 2.770e-4];
err_1e5 = [2.996e-5, 5.095e-5, 2.954e-5, 1.608e-5, 1.324e-5, 3.049e-6, 1.696e-6];
err_1e6 = [2.727e-6, 4.062e-7, 7.313e-7, 1e-13, 1e-13, 1e-13, 1e-13];
err_1e7 = [8.122e-8, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13];

figure; hold on; grid on;
loglog(distortion, err_1e4, '-o', 'LineWidth', 1.8, 'DisplayName','tol = 1e-4');
loglog(distortion, err_1e5, '-s', 'LineWidth', 1.8, 'DisplayName','tol = 1e-5');
loglog(distortion, err_1e6, '-^', 'LineWidth', 1.8, 'DisplayName','tol = 1e-6');
loglog(distortion, err_1e7, '-d', 'LineWidth', 1.8, 'DisplayName','tol = 1e-7');

xlabel('Distortion (or mean log_{10}(R))','FontSize',12);
ylabel('Flux L_2 error','FontSize',12);
title('Case 1: Flux L_2 Error vs. Distortion for Fixed Tolerances','FontSize',13);
legend('Location','northwest');
axis tight;


% Case 2
distortion = [0.02, 0.10, 0.30, 0.40, 0.70, 1.20];

% Flux L2 errors for each tolerance
err_1e4 = [1.566e-3, 7.830e-3, 3.704e-3, 3.280e-3, 1.772e-3, 1.234e-3];
err_1e5 = [5.224e-4, 1.398e-4, 1e-13, 1e-13, 6.976e-5, 1e-13];
err_1e6 = [1.574e-5, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13];
err_1e7 = [4.559e-7, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13];

figure; hold on; grid on;
loglog(distortion, err_1e4, '-o', 'LineWidth', 1.8, 'DisplayName','tol = 1e-4');
loglog(distortion, err_1e5, '-s', 'LineWidth', 1.8, 'DisplayName','tol = 1e-5');
loglog(distortion, err_1e6, '-^', 'LineWidth', 1.8, 'DisplayName','tol = 1e-6');
loglog(distortion, err_1e7, '-d', 'LineWidth', 1.8, 'DisplayName','tol = 1e-7');

xlabel('Distortion (or mean log_{10}(R))','FontSize',12);
ylabel('Flux L_2 error','FontSize',12);
title('Case 2: Flux L_2 Error vs. Distortion for Fixed Tolerances','FontSize',13);
legend('Location','northeast');
axis tight;


% Case 3
distortion = [0.02, 0.10, 0.30, 0.40, 0.70, 1.20];

err_1e4 = [1.430e-3, 3.051e-4, 2.451e-4, 2.373e-4, 2.321e-4, 3.169e-4];
err_1e5 = [2.724e-5, 2.679e-5, 1e-13, 1e-13, 2.675e-5, 1e-13];
err_1e6 = [1.401e-7, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13];

figure; hold on; grid on;
loglog(distortion, err_1e4, '-o', 'LineWidth', 1.8, 'DisplayName','tol = 1e-4');
loglog(distortion, err_1e5, '-s', 'LineWidth', 1.8, 'DisplayName','tol = 1e-5');
loglog(distortion, err_1e6, '-^', 'LineWidth', 1.8, 'DisplayName','tol = 1e-6');

xlabel('Distortion (or mean log_{10}(R))','FontSize',12);
ylabel('Flux L_2 error','FontSize',12);
title('Case 3: Flux L_2 Error vs. Distortion (Anisotropy 1:50)','FontSize',13);
legend('Location','northeast');
axis tight;

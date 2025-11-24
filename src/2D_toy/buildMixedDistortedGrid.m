% clear; close all; clc;
% 
% % ---- parameters ----
% nx = 40; nz = 20;
% Lx = 1.0; Lz = 1.0;
% 
% % ---- base structured grid ----
% xv = linspace(0, Lx, nx+1);
% zv = linspace(0, Lz, nz+1);
% [X, Z] = meshgrid(xv, zv);
% X = X'; Z = Z';
% 
% % ---- region-based distortion (6 zones) ----
% for i = 2:nx
%     for j = 2:nz
%         x = X(i,j);
%         z = Z(i,j);
% 
%         % Determine horizontal zone
%         if x < 0.33
%             x_zone = 1;
%         elseif x < 0.66
%             x_zone = 2;
%         else
%             x_zone = 3;
%         end
% 
%         % Determine vertical zone
%         if z < 0.5
%             z_zone = 2;  % bottom
%         else
%             z_zone = 1;  % top
%         end
% 
%         % Define amplitude scaling depending on zone
%         if z_zone == 1
%             base_amp = 0.2 * (Lx/nx);   % top: mild
%         else
%             base_amp = 0.4 * (Lx/nx);   % bottom: stronger
%         end
% 
%         % Horizontal weighting (left→right more distorted)
%         amp_scale = [0.0, 0.5, 1.0]; 
%         amp = base_amp * amp_scale(x_zone);
%         lambda = 2 + (x_zone-1);   % frequency grows rightward
% 
%         % Apply distortion
%         dx = amp * sin(lambda*pi*z/Lz);
%         dz = 0.5 * amp * sin(lambda*pi*x/Lx); % optional z-distortion
% 
%         X(i,j) = X(i,j) + dx;
%         Z(i,j) = Z(i,j) + dz;
%     end
% end
% 
% % ---- visualization ----
% figure('Position',[100 100 850 500]);
% hold on; axis equal tight;
% for i = 1:nx
%     for j = 1:nz
%         v1 = [X(i,j)   Z(i,j)];
%         v2 = [X(i+1,j) Z(i+1,j)];
%         v3 = [X(i+1,j+1) Z(i+1,j+1)];
%         v4 = [X(i,j+1) Z(i,j+1)];
%         patch([v1(1) v2(1) v3(1) v4(1)], [v1(2) v2(2) v3(2) v4(2)], ...
%               'w', 'EdgeColor','k', 'FaceColor','none', 'LineWidth',0.5);
%     end
% end
% 
% % region separators
% xline(0.33,'--r','LineWidth',1.2);
% xline(0.66,'--r','LineWidth',1.2);
% yline(0.5,'--r','LineWidth',1.2);
% 
% text(0.16,0.95,'R1','FontSize',12,'Color','r','FontWeight','bold','HorizontalAlignment','center');
% text(0.5,0.95,'R2','FontSize',12,'Color','r','FontWeight','bold','HorizontalAlignment','center');
% text(0.83,0.95,'R3','FontSize',12,'Color','r','FontWeight','bold','HorizontalAlignment','center');
% text(0.16,0.45,'R4','FontSize',12,'Color','r','FontWeight','bold','HorizontalAlignment','center');
% text(0.5,0.45,'R5','FontSize',12,'Color','r','FontWeight','bold','HorizontalAlignment','center');
% text(0.83,0.45,'R6','FontSize',12,'Color','r','FontWeight','bold','HorizontalAlignment','center');
% 
% xlabel('x'); ylabel('z');
% title('2D Mixed-Distortion Structured Grid (6 regions)');


% clear; close all; clc;
% 
% % ---- parameters ----
% nx = 40; nz = 20;
% Lx = 1.0; Lz = 1.0;
% 
% % ---- base structured grid ----
% xv = linspace(0, Lx, nx+1);
% zv = linspace(0, Lz, nz+1);
% [X, Z] = meshgrid(xv, zv);
% X = X';  Z = Z';   % for consistent orientation
% 
% % ---- region-based distortion ----
% for i = 2:nx
%     for j = 2:nz
%         x = X(i,j);  z = Z(i,j);
% 
%         if x < 0.33
%             % Region 1: perfectly structured
%             dx = 0; dz = 0;
% 
%         elseif x < 0.66
%             % Region 2: mild sinusoidal distortion
%             amp = 0.25 * (Lx/nx);
%             lambda = 2;
%             dx = amp * sin(lambda*pi*z/Lz);
%             dz = 0;
% 
%         else
%             % Region 3: stronger distortion
%             amp = 0.5 * (Lx/nx);
%             lambda = 3;
%             dx = amp * sin(lambda*pi*z/Lz) * (x - 0.66) / (1 - 0.66);
%             dz = 0;
%         end
% 
%         X(i,j) = X(i,j) + dx;
%         Z(i,j) = Z(i,j) + dz;
%     end
% end
% 
% % ---- visualization ----
% figure('Position',[100 100 800 400]);
% hold on; axis equal tight;
% for i = 1:nx
%     for j = 1:nz
%         % each cell corner indices
%         v1 = [X(i,j)   Z(i,j)];
%         v2 = [X(i+1,j) Z(i+1,j)];
%         v3 = [X(i+1,j+1) Z(i+1,j+1)];
%         v4 = [X(i,j+1) Z(i,j+1)];
%         patch([v1(1) v2(1) v3(1) v4(1)], [v1(2) v2(2) v3(2) v4(2)], 'w', ...
%               'EdgeColor','k','FaceColor','none','LineWidth',0.6);
%     end
% end
% 
% % add region separators
% xline(0.33,'--r','LineWidth',1.2);
% xline(0.66,'--r','LineWidth',1.2);
% text(0.16,1.02,'Region 1 (Structured)','FontSize',11,'Color','r','HorizontalAlignment','center');
% text(0.5,1.02,'Region 2 (Mild distortion)','FontSize',11,'Color','r','HorizontalAlignment','center');
% text(0.83,1.02,'Region 3 (Strong distortion)','FontSize',11,'Color','r','HorizontalAlignment','center');
% 
% xlabel('x'); ylabel('z');
% title('Piecewise Mixed-Quality Structured Grid');

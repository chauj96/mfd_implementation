function G = scratch()
%% set boundary
tx = 1;
ty = 1;
tz = 0.4;
bdr   = [ 0, 0, 0;  ...
          tx, 0, 0;  ...
          tx, ty, 0;  ...
          0, ty, 0;  ...
          0, 0, tz;  ...
          tx, 0, tz;  ...
          tx, ty, tz;  ...
          0, ty, tz];

%% Set gridding parameters
fGs = tx/20;
rho = @(p) fGs*(1+0*p(:,1));

%% Triangulate the two curved surfaces exactly as in your original code
swap1 = [0,1;1,0]==1;
swap2 = [1,0;0,1]==1;
rectangle1 = [min(bdr(:, [1,3])); max(bdr(:,[1,3]))];
fixedPts = [rectangle1; rectangle1(swap1)'; rectangle1(swap2)'];
hd = @(p) rho(p)/fGs;
fd = @(p) drectangle(p, rectangle1(1),rectangle1(2), rectangle1(3), rectangle1(4));

[Pts,t] = distmesh2d(fd, hd, fGs, rectangle1, fixedPts, false);

t1.ConnectivityList = t;
t2.ConnectivityList = t;

faultHeight1z = @(p)   ty/6*ones(size(p,1),1) + 0.4*p(:,2);
faultHeight2z = @(p) 5*ty/6*ones(size(p,1),1) - 0.4*p(:,2);

t1.Points = [Pts(:,1), faultHeight1z(Pts), Pts(:,2)];
t2.Points = [Pts(:,1), faultHeight2z(Pts), Pts(:,2)];

faultHeight1x = @(p) p(:,2) + 0.2*(p(:,1)-0.75).^2 + 0.1*sin(2*pi/tx*p(:,1));
faultHeight2x = @(p) p(:,2) + 0.2*(p(:,1)-0.75).^2;

t1.Points(:,2) = faultHeight1x(t1.Points);
t2.Points(:,2) = faultHeight2x(t2.Points);

%% Visualize the two surfaces
% figure;
% color = get(gca,'colororder');
% hold on
% trisurf(t1.ConnectivityList,t1.Points(:,1),t1.Points(:,2),t1.Points(:,3),...
%     'FaceColor',color(1,:),'EdgeAlpha',.1,'FaceAlpha',0.7);
% trisurf(t2.ConnectivityList,t2.Points(:,1),t2.Points(:,2),t2.Points(:,3),...
%     'FaceColor',color(2,:),'EdgeAlpha',.1,'FaceAlpha',0.7);
% plotGrid(cartGrid([1 1 1],[tx ty tz]),'FaceColor','none');
% view(3); axis equal
% camlight left; camlight headlight
% title('Two curved internal surfaces')

%% Expand boundary slightly
% bdr = 1.01 * bdr;

%% Create reservoir sites
dt = 1 / 20; % mesh size
xmax = max(bdr(:,1)) - dt / 2; xmin = min(bdr(:,1)) + dt / 2;
ymax = max(bdr(:,2)) - dt / 2; ymin = min(bdr(:,2)) + dt / 2;
zmax = max(bdr(:,3)) - dt / 2; zmin = min(bdr(:,3)) + dt / 2;

xr = xmin:dt:xmax; yr = ymin:dt:ymax; zr = zmin:dt:zmax;
[X,Y,Z] = ndgrid(xr, yr, zr);
rSitesAll = [X(:), Y(:), Z(:)];

y1 = faultHeight1x([rSitesAll(:,1), rSitesAll(:,3)]);
y2 = faultHeight2x([rSitesAll(:,1), rSitesAll(:,3)]);

d1 = abs(rSitesAll(:,2) - y1);
d2 = abs(rSitesAll(:,2) - y2);

band1 = 1.5*dt;
band2 = 1.5*dt;
keep = (d1 > band1) & (d2 > band2);
rSites = rSitesAll(keep,:);

s1 = t1.Points;
s2 = t2.Points;

s1 = unique(round(s1,10),'rows');
s2 = unique(round(s2,10),'rows');

off = 0.1*dt;

s1u = s1; s1u(:,2) = s1u(:,2) + off;
s1d = s1; s1d(:,2) = s1d(:,2) - off;

s2u = s2; s2u(:,2) = s2u(:,2) + off;
s2d = s2; s2d(:,2) = s2d(:,2) - off;

tol = 1e-10;
inside = @(P) ...
    P(:,1) > 0+tol & P(:,1) < tx-tol & ...
    P(:,2) > 0+tol & P(:,2) < ty-tol & ...
    P(:,3) > 0+tol & P(:,3) < tz-tol;

s1u = s1u(inside(s1u),:);
s1d = s1d(inside(s1d),:);
s2u = s2u(inside(s2u),:);
s2d = s2d(inside(s2d),:);

midMask = (rSites(:,2) > faultHeight1x([rSites(:,1), rSites(:,3)]) + 0.2*dt) & ...
          (rSites(:,2) < faultHeight2x([rSites(:,1), rSites(:,3)]) - 0.2*dt);

midSites = rSites(midMask,:);
midSites = midSites(1:2:end,:); 

%% Merge all sites
sites = [s1; s2; s1u; s1d; s2u; s2d; rSites; midSites];
sites = unique(round(sites,10),'rows');
eps_b = 1e-6;

sites(:,1) = min(max(sites(:,1), eps_b), tx - eps_b);
sites(:,2) = min(max(sites(:,2), eps_b), ty - eps_b);
sites(:,3) = min(max(sites(:,3), eps_b), tz - eps_b);
sites = unique(round(sites,8),'rows');

fprintf('Total number of sites = %d\n', size(sites,1));

%% Build 3D PEBI mesh
G = mirroredPebi3D(sites, bdr);
G = computeGeometry(G);

fc = G.faces.centroids;

fprintf('Mesh built:\n');
fprintf('  #cells = %d\n', G.cells.num);
fprintf('  #faces = %d\n', G.faces.num);
fprintf('  #nodes = %d\n', G.nodes.num);

% figure; hold on
% plot3(rSites(:,1), rSites(:,2), rSites(:,3), '.', 'MarkerSize', 5, 'Color', color(4,:));
% plot3(s1(:,1), s1(:,2), s1(:,3), '.', 'MarkerSize', 12, 'Color', color(1,:));
% plot3(s2(:,1), s2(:,2), s2(:,3), '.', 'MarkerSize', 12, 'Color', color(2,:));
% plot3(s1u(:,1), s1u(:,2), s1u(:,3), '.', 'MarkerSize', 7, 'Color', color(1,:));
% plot3(s1d(:,1), s1d(:,2), s1d(:,3), '.', 'MarkerSize', 7, 'Color', color(1,:));
% plot3(s2u(:,1), s2u(:,2), s2u(:,3), '.', 'MarkerSize', 7, 'Color', color(2,:));
% plot3(s2d(:,1), s2d(:,2), s2d(:,3), '.', 'MarkerSize', 7, 'Color', color(2,:));
% plotGrid(cartGrid([1 1 1],[tx ty tz]),'FaceColor','none');
% axis equal
% view(-100,50); camlight left; camlight headlight
% title('Sites used to generate the mesh')

c = G.cells.centroids;

cl = c(:,2) > faultHeight2x([c(:,1), c(:,3)]);
cr = c(:,2) < faultHeight1x([c(:,1), c(:,3)]);
cm = ~cl & ~cr;

figure(); hold on
color = lines(7);
plotGrid(G, cl, ...
    'FaceColor', color(1,:), ...
    'EdgeColor', [0 0 0], ...
    'EdgeAlpha', 0.15);

plotGrid(G, cm, ...
    'FaceColor', color(2,:), ...
    'EdgeColor', [0 0 0], ...
    'EdgeAlpha', 0.15);

plotGrid(G, cr, ...
    'FaceColor', color(3,:), ...
    'EdgeColor', [0 0 0], ...
    'EdgeAlpha', 0.15);

plotGrid(G, ...
    'FaceColor', 'none', ...
    'EdgeColor', [0 0 0], ...
    'EdgeAlpha', 0.2);

trisurf(t1.ConnectivityList, ...
    t1.Points(:,1), t1.Points(:,2), t1.Points(:,3), ...
    'FaceColor', color(1,:), ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');

trisurf(t2.ConnectivityList, ...
    t2.Points(:,1), t2.Points(:,2), t2.Points(:,3), ...
    'FaceColor', color(2,:), ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');

view(120,30);
light('Position',[5 10 -1],'Style','infinite')
light('Position',[-1 -1 2],'Style','local')

axis equal off
set(gca,'zdir','normal')

title('Final mesh (clean boundary + sharp edges)')
end




% plot_topography_shapes.m
%
% Plots topography s(x,y) and gradients ds/dx, ds/dy for four shapes.
%
%   A — Circle              r = sqrt(x^2+y^2)          (current 2D model)
%   B — Square (L-inf)      r = max(|x|,|y|)           (sharp corners)
%   C — Smooth square (L-p) r = (|x|^p+|y|^p)^(1/p)   (C-inf everywhere)
%   D — Square smooth edges s = A*(1 - Hx(x)*Hy(y))    (each edge independent)
%
% Options A-C use flc2hs on a single "radius"; Option D applies one
% flc2hs per axis direction so each of the four sides is smoothed
% independently — the corners are the natural product of two 1D transitions.

%% Parameters
L       = 1.0;    % domain half-size [normalised]
offset  = 0.6;    % ring/square half-width [normalised]
epsilon = 0.05;   % flc2hs transition half-width [normalised]
A       = 1.0;    % topography amplitude
p       = 8;      % L-p exponent for smooth square (Option C)

%% Grid
N   = 500;
x1d = linspace(-L, L, N);
[X, Y] = meshgrid(x1d, x1d);

%% ── Option A: Circle ──────────────────────────────────────────────────────
r_c      = sqrt(X.^2 + Y.^2);
r_c_safe = max(r_c, 1e-12);
xp_c     = r_c - offset;
sA       = A * flc2hs_fn(xp_c, epsilon);
dAdx     = A * flc2hs_d(xp_c, epsilon) .* X ./ r_c_safe;
dAdy     = A * flc2hs_d(xp_c, epsilon) .* Y ./ r_c_safe;

%% ── Option B: Square, L-inf norm (sharp corners) ─────────────────────────
r_sq     = max(abs(X), abs(Y));
xp_sq    = r_sq - offset;
sB       = A * flc2hs_fn(xp_sq, epsilon);
gBx      = double(abs(X) >= abs(Y)) .* sign(X);  % diagonal -> x-face
gBy      = double(abs(Y) >  abs(X)) .* sign(Y);
dBdx     = A * flc2hs_d(xp_sq, epsilon) .* gBx;
dBdy     = A * flc2hs_d(xp_sq, epsilon) .* gBy;

%% ── Option C: Smooth square, L-p norm (C-inf) ────────────────────────────
r_p      = (abs(X).^p + abs(Y).^p).^(1/p);
r_p_safe = max(r_p, 1e-12);
xp_p     = r_p - offset;
sC       = A * flc2hs_fn(xp_p, epsilon);
gCx      = X .* abs(X).^(p-2) ./ r_p_safe.^(p-1);
gCy      = Y .* abs(Y).^(p-2) ./ r_p_safe.^(p-1);
dCdx     = A * flc2hs_d(xp_p, epsilon) .* gCx;
dCdy     = A * flc2hs_d(xp_p, epsilon) .* gCy;

%% ── Option D: Square with smoothed edges (product of 1D steps) ───────────
%
%  s(x,y) = A * [1 - Hx(x)*Hy(y)]
%
%  where  Hx(x) = flc2hs(offset - |x|, epsilon)  ~ 1 inside |x| < offset
%         Hy(y) = flc2hs(offset - |y|, epsilon)  ~ 1 inside |y| < offset
%
%  So Hx*Hy ~ 1 inside the square, ~ 0 outside -> s jumps from 0 to A
%  at each edge independently.  Corners are the product of two transitions.
%
%  Gradient:
%    ds/dx = -A * dHx/dx * Hy  =  A * flc2hs'(offset-|x|,eps)*sign(x) * Hy(y)
%    ds/dy = -A * Hx * dHy/dy  =  A * flc2hs'(offset-|y|,eps)*sign(y) * Hx(x)

Hx      = flc2hs_fn(offset - abs(X), epsilon);   % ~1 inside |x|<offset
Hy      = flc2hs_fn(offset - abs(Y), epsilon);   % ~1 inside |y|<offset
dHxdx   = -flc2hs_d(offset - abs(X), epsilon) .* sign(X);  % d(Hx)/dx
dHydy   = -flc2hs_d(offset - abs(Y), epsilon) .* sign(Y);  % d(Hy)/dy

sD      = A * (1 - Hx .* Hy);
dDdx    = A * (-dHxdx .* Hy);     % = +A * flc2hs'(...)*sign(x) * Hy
dDdy    = A * (-Hx .* dHydy);     % = +A * Hx * flc2hs'(...)*sign(y)

%% ── Figure 1: 2D maps for all four shapes ────────────────────────────────
shapes = {sA,dAdx,dAdy; sB,dBdx,dBdy; sC,dCdx,dCdy; sD,dDdx,dDdy};
rlabel = {
  'A — Circle  (L^2)',
  'B — Square  (L^\infty, sharp corners)',
  'C — Smooth square  (L^p, p=8)',
  'D — Square + smoothed edges  (product of 1D steps)'
};
clabel = {'s(x,y)', '\partials/\partialx', '\partials/\partialy'};

clims_s = [0, A];
dpk     = A / epsilon * 15/16;          % peak value of flc2hs derivative
clims_d = [-dpk, dpk] * 0.55;

figure('Name','Topography — 2D maps','Position',[30 30 1250 920]);
for row = 1:4
    for col = 1:3
        ax = subplot(4, 3, (row-1)*3 + col);
        data = shapes{row,col};
        imagesc(x1d, x1d, data);
        axis equal tight;  set(gca,'YDir','normal');
        colorbar;
        if col == 1
            clim(clims_s);  colormap(ax,'parula');
        else
            clim(clims_d);  colormap(ax,'coolwarm');
        end
        if row == 1,  title(clabel{col},'FontSize',11,'FontWeight','bold');  end
        if col == 1,  ylabel(rlabel{row},'FontSize',9,'FontWeight','bold');  end
        xlabel('x');
    end
end
sgtitle(sprintf('Topography shapes  |  offset=%.2g,  \\epsilon=%.2g', ...
        offset, epsilon), 'FontSize', 13);

%% ── Figure 2: 1D cross-sections along y=0 and diagonal y=x ──────────────
[~, iy0] = min(abs(x1d));            % index for y = 0
xd       = x1d(x1d >= 0);           % diagonal: first quadrant, y = x

% Evaluate Option D on diagonal manually
Hx_d    = flc2hs_fn(offset - abs(xd), epsilon);
Hy_d    = flc2hs_fn(offset - abs(xd), epsilon);   % y = x so same
dHx_d   = -flc2hs_d(offset - abs(xd), epsilon) .* sign(xd);
sD_diag = A * (1 - Hx_d .* Hy_d);
dDdx_d  = A * (-dHx_d .* Hy_d);

sVec    = {sA(iy0,:), sB(iy0,:), sC(iy0,:), sD(iy0,:)};
dVec    = {dAdx(iy0,:), dBdx(iy0,:), dCdx(iy0,:), dDdx(iy0,:)};
cols_   = {'b-','r--','g:','m-.'};
lw_     = [2.0, 2.0, 2.5, 2.0];
names   = {'Circle','Square L^\infty','Smooth sq (L^p)','Smooth edges (D)'};

figure('Name','1D cross-sections','Position',[30 980 1100 420]);

subplot(1,3,1);
hold on;
for k = 1:4
    plot(x1d, sVec{k}, cols_{k},'LineWidth',lw_(k),'DisplayName',names{k});
end
xlabel('x');  ylabel('s(x, 0)');
title('s(x,y=0)  —  same as 1D axisym profile');
legend('Location','northwest');  grid on;  hold off;

subplot(1,3,2);
hold on;
for k = 1:4
    plot(x1d, dVec{k}, cols_{k},'LineWidth',lw_(k),'DisplayName',names{k});
end
xlabel('x');  ylabel('\partials/\partialx  at  y=0');
title('Derivative along y = 0');
legend;  grid on;  hold off;

subplot(1,3,3);
% Diagonal cross-section (y = x)
r_c_d   = sqrt(2)*xd;
sA_d    = A * flc2hs_fn(r_c_d - offset, epsilon);
sB_d    = A * flc2hs_fn(xd - offset, epsilon);   % max(|x|,|x|) = |x|
sC_d    = A * flc2hs_fn((2*abs(xd).^p).^(1/p) - offset, epsilon);
hold on;
plot(xd, sA_d,   'b-', 'LineWidth',2.0, 'DisplayName','Circle');
plot(xd, sB_d,   'r--','LineWidth',2.0, 'DisplayName','Square L^\infty');
plot(xd, sC_d,   'g:', 'LineWidth',2.5, 'DisplayName','Smooth sq (L^p)');
plot(xd, sD_diag,'m-.','LineWidth',2.0, 'DisplayName','Smooth edges (D)');
xlabel('x  (along y = x)');  ylabel('s(x,x)');
title('Shape along diagonal y = x');
legend('Location','northwest');  grid on;  hold off;

sgtitle('Cross-section comparisons','FontSize',13);

%% ── Figure 3: 3D surfaces ────────────────────────────────────────────────
figure('Name','3D surfaces','Position',[1300 30 1100 600]);
surf3 = {sA, sB, sC, sD};
for k = 1:4
    subplot(1,4,k);
    surf(X, Y, surf3{k}, 'EdgeColor','none');
    axis tight;  view([-40,35]);
    clim(clims_s);  colormap('parula');
    xlabel('x');  ylabel('y');  zlabel('s');
    title(rlabel{k},'FontSize',9);
end
sgtitle('3D topography surfaces','FontSize',13);

%% ── Local functions ───────────────────────────────────────────────────────

function y = flc2hs_fn(x, scale)
% C2-smooth Heaviside matching COMSOL flc2hs(x, scale).
% 0 for x < -scale,  1 for x > +scale,
% 5th-order polynomial (zero 1st and 2nd derivatives at endpoints).
    u = x / scale;
    y = zeros(size(x));
    y(u > 1) = 1;
    m = (u >= -1) & (u <= 1);
    y(m) = 0.5 + (15/16)*u(m) - (5/8)*u(m).^3 + (3/16)*u(m).^5;
end

function y = flc2hs_d(x, scale)
% Derivative of flc2hs_fn w.r.t. x:
% (15/16/scale)*(1-(x/scale)^2)^2  inside the band, zero outside.
    u = x / scale;
    y = zeros(size(x));
    m = (u >= -1) & (u <= 1);
    y(m) = (15/16/scale) * (1 - u(m).^2).^2;
end

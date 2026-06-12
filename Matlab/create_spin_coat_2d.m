% create_spin_coat_2d.m
% Creates a 2D Cartesian spin coating model based on
% spin_coat_stillwagon_topography_solvent_evaporation.mph (1D axisymmetric).
%
% Variables: [hh, P, c]  (film thickness, curvature pressure, concentration)
% Eq 0 (hh): d(hh)/dt + div((hh³/(3*eta))*(∇P - ρω²r)) = -F(c)
%            [1D original: eta*d(hh)/dt = ... divided by eta(c)]
% Eq 1 (P):  div(-gamma*∇(hh+s)) + P = 0   [quasi-static, da=0]
% Eq 2 (c):  dc/dt + D*div(hh*∇c)/hh_safe = f3/hh_safe
%            [1D original: hh*dc/dt = ... divided by hh_safe]
%
% da is constant (1,0,1) after normalization — avoids COMSOL API restriction
% that rejects variable expressions in the da matrix.

addpath('C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
try mphstart; catch; end
import com.comsol.model.*
import com.comsol.model.util.*

%% Create model
model = ModelUtil.create('Model');
model.modelNode.create('comp1', true);

%% Parameters
p = model.param;
p.set('rot',     '2000[rpm]',       'Rotation speed');
p.set('w',       'rot*2*pi[rad]',   'Angular velocity');
p.set('rho',     '1000[kg/m^3]',    'Fluid density');
p.set('h_0',     '50e-6',           'Initial film thickness (m)');
p.set('A',       '1e-6[m]',         'Topography amplitude');
p.set('offset',  '90e-3',           'Topography radial position (m)');
p.set('epsilon', '5.4e-6',          'Topography transition half-width (m)');
p.set('gamma',   '0.04[N/m]',       'Surface tension');
p.set('Fo',      '10e-8[m/s]',      'Evaporation flux coefficient');
p.set('c_max',   '0.99',            'Maximum solvent volume fraction');
p.set('D',       '10e-10[m^2/s]',   'Diffusion coefficient');
p.set('tau',     '0',               'Surface shear stress (Pa)');
p.set('L',       '100e-3',          'Domain half-size (m)');
p.set('h_min',   '1e-9[m]',         'Minimum film thickness (regularization for eq2 normalization)');

%% Geometry: square domain [-L, L] x [-L, L]
geom = model.geom.create('geom1', 2);
r1   = geom.feature.create('r1', 'Rectangle');
r1.set('size', {'2*L', '2*L'});
r1.set('pos',  {'-L',  '-L'});
geom.run;

%% Functions (exact forms from 1D model)
fa = model.func.create('fa1', 'Analytic');
fa.set('funcname', 'F');
fa.set('args',    {'c'});
fa.set('expr',    'Fo*(1-c/c_max)');
fa.set('argunit', {'1'});
fa.set('fununit', 'm/s');

fb = model.func.create('fa2', 'Analytic');
fb.set('funcname', 'eta');
fb.set('args',    {'c'});
fb.set('expr',    '2[mPa*s] + 1[Pa*s]*c');
fb.set('argunit', {'1'});
fb.set('fununit', 'Pa*s');

%% Variables
v = model.variable.create('var1');
v.model('comp1');
v.set('r_pos',    'sqrt(x^2+y^2)',                             'Radial coordinate');
v.set('r_safe',   'max(r_pos, 1e-12[m])',                      'Regularised r');
v.set('x_prime',  'r_pos - offset',                            'Radial shift');
v.set('s',        'A*flc2hs(x_prime, epsilon)',                 'Substrate topography');
v.set('dsdx',     'A*flc2hs(x_prime, epsilon, 1)*x/r_safe',    'ds/dx');
v.set('dsdy',     'A*flc2hs(x_prime, epsilon, 1)*y/r_safe',    'ds/dy');
v.set('hh_safe',  'max(hh, h_min)',                             'Regularised film thickness');
v.set('eta_safe', 'max(eta(c), 1e-6[Pa*s])',                   'Regularised viscosity');

%% Physics: General Form PDE with 3 components [hh, P, c]
phys = model.physics.create('g', 'GeneralFormPDE', 'geom1');
phys.field('dimensionless').field('hh');
phys.field('dimensionless').component({'hh','P','c'});

gfeq = phys.feature('gfeq1');

% Set all 9 da entries at once — COMSOL validates the full matrix on every
% setIndex call, so any unset default can trigger "not a scalar". Using set()
% with a complete constant cell array avoids that.
%   da = diag(1, 0, 1): hh evol (normalized) | P algebraic | c evol (normalized)
gfeq.set('da', {'1','0','0'; '0','0','0'; '0','0','1'});

% --- Eq 0: d(hh)/dt + div((hh³/(3*eta))*(∇P - ρω²r)) = -F(c) ---
% Derived from: eta*(d(hh)/dt) + div(hh³/3*(∇P-ρω²r)) = -F*eta  [÷ eta(c)]
gfeq.setIndex('Ga', '(hh^3/(3*eta_safe))*(-Px + rho*w^2*x)', 0, 0);
gfeq.setIndex('Ga', '(hh^3/(3*eta_safe))*(-Py + rho*w^2*y)', 0, 1);
gfeq.setIndex('f',  '-F(c)', 0);

% --- Eq 1: div(-gamma*∇(hh+s)) + P = 0   (P = gamma*Laplacian(hh+s)) ---
gfeq.setIndex('Ga', '-gamma*(hhx + dsdx)', 1, 0);
gfeq.setIndex('Ga', '-gamma*(hhy + dsdy)', 1, 1);
gfeq.setIndex('f',  'P', 1);

% --- Eq 2: dc/dt - div(D*(hh/hh_safe)*∇c) = f3/hh_safe + D*(hhx*cx+hhy*cy)/hh_safe ---
% Derived from: hh*dc/dt - div(D*hh*∇c) = f3  [÷ hh_safe]
% Ga captures the Laplacian part (hh/hh_safe * D*Δc);
% the grad(hh)*grad(c)/hh_safe correction is added to f.
gfeq.setIndex('Ga', '-D*(hh/hh_safe)*cx', 2, 0);
gfeq.setIndex('Ga', '-D*(hh/hh_safe)*cy', 2, 1);
f3 = ['(c*F(c)' ...
      ' - (hh^2*tau/(2*eta_safe))*c' ...
      ' - (hh^3/(3*eta_safe))*((Px - rho*w^2*x)*cx + (Py - rho*w^2*y)*cy))/hh_safe' ...
      ' + D*(hhx*cx + hhy*cy)/hh_safe'];
gfeq.setIndex('f', f3, 2);

%% Initial conditions
init = phys.feature('init1');
init.setIndex('u', 'h_0',  0);   % hh
init.setIndex('u', '0',    1);   % P  (relaxes instantly to curvature)
init.setIndex('u', '0.01', 2);   % c  (matches 1D model)

%% Boundary conditions
% Outer edges: allow centrifugal outflow for hh; zero flux (sealed) for P and c.
% n.Ga_hh = -(hh³/(3*eta)) * ρω² * (x*nx + y*ny)  [centrifugal, no pressure term]
flx = phys.feature.create('flux1', 'FluxBoundary', 1);
flx.selection.all;
flx.setIndex('g', '0',                                                    0);  % P
flx.setIndex('g', '-(hh^3/(3*eta_safe))*rho*w^2*(x*nx + y*ny)',          1);  % hh
flx.setIndex('g', '0',                                                    2);  % c

%% Mesh
mesh = model.mesh.create('mesh1', 'geom1');
mesh.feature('size').set('hauto', 4);  % fine mesh
mesh.feature.create('ftri1', 'FreeTri');
mesh.run;

%% Study
std   = model.study.create('std1');
tnode = std.feature.create('time', 'Transient');
tnode.set('tlist', 'range(0,1,100)');
tnode.set('rtol',  '1e-3');

%% Save
save_path = 'C:\Users\20251195\OneDrive - TU Eindhoven\Documents\COMSOL\Project\spin_coat_2d.mph';
model.save(save_path);
fprintf('Saved: %s\n', save_path);

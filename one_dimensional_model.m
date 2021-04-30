function [Y_mod, Z_mod, dY_mod_dy, d2Y_mod_dy2, dZ_mod_dy, d2Z_mod_dy2] ...
    = one_dimensional_model(f, rho, l, w, t, s11E, eps33T, d31)
% Calculate the admittance (Y_mod) and impedance (Z_mod)
% in the one-dimensional vibration model and their
% derivatives and Hessian matrices (with respect to y).

N = size(f, 1); % Number of frequency data points

% Define a non-dimensional frequency to simplify the following expressions.
f_star = pi * f * l * sqrt(rho * s11E); 

% Admittance (model) (S)
Y_mod = 2i * pi * f * l * w / t .* (eps33T - d31^2 / s11E ...
    + d31^2 / s11E * tan(f_star) ./ f_star);

% Impedance (model) (Ohm)
Z_mod = Y_mod .^ (-1); 

if nargout == 2
    return; % Return when only Y_mod and Z_mod are needed.
end

% Partial derivative of Y_mod with respect to s11E
dY_mod_ds11E = 1i * pi * f * l * w * d31^2 / (t * s11E^2) ...
    .* (2 - 3 * tan(f_star) ./ f_star + cos(f_star).^(-2));

% Partial derivative of Y_mod with respect to eps33T
dY_mod_deps33T = 2i * pi * f * l * w / t;

% Partial derivative of Y_mod with respect to d31
dY_mod_dd31 = 4i * pi * f * l * w * d31 / (t * s11E) ...
    .* (-1 + tan(f_star) ./ f_star);

dY_mod_dy = zeros(N,6); % Gradient of the Y_mod with respect to y
dY_mod_dy(:,1) = dY_mod_ds11E;
dY_mod_dy(:,2) = 1i * dY_mod_ds11E;
dY_mod_dy(:,3) = dY_mod_deps33T;
dY_mod_dy(:,4) = 1i * dY_mod_deps33T;
dY_mod_dy(:,5) = dY_mod_dd31;
dY_mod_dy(:,6) = 1i * dY_mod_dd31;

% Second-order partial derivative with respect to s11E
d2Y_mod_ds11E2 = 1i * pi * f * l * w * d31^2 / (2 * t * s11E^3) ...
    .* (-8 + 15 * tan(f_star) ./ f_star ...
    - 7 * cos(f_star).^(-2) + 2 * f_star .* sin(f_star) ./ cos(f_star).^3);

% Second-order partial derivative with respect to s11E and eps33T
d2Y_mod_ds11E_deps33T = zeros(N,1);

% Second-order partial derivative with respect to s11E and d31
d2Y_mod_ds11E_dd31 = 2i * pi * f * l * w * d31 / (t * s11E^2) ...
    .* (2 - 3 * tan(f_star) ./ f_star + cos(f_star).^(-2));

% Second-order partial derivative with respect to eps33T
d2Y_mod_deps33T2 = zeros(N,1);

% Second-order partial derivative with respect to eps33T and d31
d2Y_mod_deps33T_dd31 = zeros(N,1);

% Second-order partial derivative with respect to d31
d2Y_mod_dd312 =  4i * pi * f * l * w / (t * s11E) ...
    .* (-1 + tan(f_star) ./ f_star);

% Hessian matrix of Y_mod with respect to y
% Note that those d2Y_mod_dy2(:,j,k) with j>k are not used in the
% calculation of d2E_dy2 due to symmetry and are therefore omitted here.
d2Y_mod_dy2 = zeros(N,6,6);
d2Y_mod_dy2(:,1,1) = d2Y_mod_ds11E2;
d2Y_mod_dy2(:,1,2) = 1i * d2Y_mod_ds11E2;
d2Y_mod_dy2(:,1,3) = d2Y_mod_ds11E_deps33T;
d2Y_mod_dy2(:,1,4) = 1i * d2Y_mod_ds11E_deps33T;
d2Y_mod_dy2(:,1,5) = d2Y_mod_ds11E_dd31;
d2Y_mod_dy2(:,1,6) = 1i * d2Y_mod_ds11E_dd31;

d2Y_mod_dy2(:,2,2) = -d2Y_mod_ds11E2;
d2Y_mod_dy2(:,2,3) = 1i * d2Y_mod_ds11E_deps33T;
d2Y_mod_dy2(:,2,4) = -d2Y_mod_ds11E_deps33T;
d2Y_mod_dy2(:,2,5) = 1i * d2Y_mod_ds11E_dd31;
d2Y_mod_dy2(:,2,6) = -d2Y_mod_ds11E_dd31;

d2Y_mod_dy2(:,3,3) = d2Y_mod_deps33T2;
d2Y_mod_dy2(:,3,4) = 1i * d2Y_mod_deps33T2;
d2Y_mod_dy2(:,3,5) = d2Y_mod_deps33T_dd31;
d2Y_mod_dy2(:,3,6) = 1i * d2Y_mod_deps33T_dd31;

d2Y_mod_dy2(:,4,4) = -d2Y_mod_deps33T2;
d2Y_mod_dy2(:,4,5) = 1i * d2Y_mod_deps33T_dd31;
d2Y_mod_dy2(:,4,6) = -d2Y_mod_deps33T_dd31;

d2Y_mod_dy2(:,5,5) = d2Y_mod_dd312;
d2Y_mod_dy2(:,5,6) = 1i * d2Y_mod_dd312;

d2Y_mod_dy2(:,6,6) = -d2Y_mod_dd312;

dZ_mod_dy = zeros(N,6); % Gradient of Z_mod with respect to y
for j = 1 : 6
    dZ_mod_dy(:,j) = -dY_mod_dy(:,j) ./ Y_mod.^2;
end

% Hessian matrix of Z_mod with respect to y
% Note that those d2Z_mod_dy2(:,j,k) with j>k are not used in the
% calculation of d2E_dy2 due to symmetry and are therefore omitted here.
d2Z_mod_dy2 = zeros(N,6,6);
for k = 1 : 6
    for j = 1 : k
        d2Z_mod_dy2(:,j,k) = -d2Y_mod_dy2(:,j,k) ./ Y_mod.^2 ...
            + 2 * dY_mod_dy(:,j) .* dY_mod_dy(:,k) ./ Y_mod.^3;
    end
end
end


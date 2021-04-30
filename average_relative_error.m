function [E, dE_dy, d2E_dy2] = average_relative_error(Y_exp, Z_exp, ...
    Y_mod, Z_mod, dY_mod_dy, dZ_mod_dy, d2Y_mod_dy2, d2Z_mod_dy2)
% Calculate the average relative error (E) of the admittance and
% impedance in the one-dimensional vibration model and its
% derivative and Hessian matrix (with respect to y).

N = size(Y_exp, 1); % Number of frequency data points

G_exp = real(Y_exp); % Conductance (experiment) (S)
B_exp = imag(Y_exp); % Susceptance (experiment) (S)
R_exp = real(Z_exp); % Resistance (experiment) (Ohm)
X_exp = imag(Z_exp); % Reactance (experiment) (Ohm)

G_mod = real(Y_mod); % Conductance (model) (S)
B_mod = imag(Y_mod); % Susceptance (model) (S)
R_mod = real(Z_mod); % Resistance (model) (Ohm)
X_mod = imag(Z_mod); % Reactance (model) (Ohm)

% Relative errors of G_mod, B_mod, R_mod, and X_mod
E_G = sqrt(sum(((G_mod - G_exp) ./ G_exp).^2) / N);
E_B = sqrt(sum(((B_mod - B_exp) ./ B_exp).^2) / N);
E_R = sqrt(sum(((R_mod - R_exp) ./ R_exp).^2) / N);
E_X = sqrt(sum(((X_mod - X_exp) ./ X_exp).^2) / N);

% Average relative error
E = 0.25 * (sum(E_G) + sum(E_B) + sum(E_R) + sum(E_X));

if nargout == 1
    return; % Return when only the average relative error is needed.
end

coef_G = 1 / E_G * G_exp.^(-2);
coef_B = 1 / E_B * B_exp.^(-2);
coef_R = 1 / E_R * R_exp.^(-2);
coef_X = 1 / E_X * X_exp.^(-2);

% Gradients of E_G, E_B, E_R, and E_X with respect to y
dEG_dy = 1 / N * (coef_G' * ((G_mod - G_exp) .* real(dY_mod_dy)))';
dEB_dy = 1 / N * (coef_B' * ((B_mod - B_exp) .* imag(dY_mod_dy)))';
dER_dy = 1 / N * (coef_R' * ((R_mod - R_exp) .* real(dZ_mod_dy)))';
dEX_dy = 1 / N * (coef_X' * ((X_mod - X_exp) .* imag(dZ_mod_dy)))';

% Gradient of the average relative error with respect to y
dE_dy = 0.25 * (dEG_dy + dEB_dy + dER_dy + dEX_dy);

% Hessian matrix of the average relative error with respect to y
d2E_dy2 = zeros(6,6);

% Only those d2E_dy2(j,k) with j<=k need to be calculated due to symmetry.
for k = 1 : 6
    for j = 1 : k
        d2E_dy2(j,k) = -0.25 * (dEG_dy(j) * dEG_dy(k) / E_G ...
            + dEB_dy(j) * dEB_dy(k) / E_B ...
            + dER_dy(j) * dER_dy(k) / E_R ...
            + dEX_dy(j) * dEX_dy(k) / E_X) ...
            + 0.25 / N * sum( ...
            coef_G .* (real(dY_mod_dy(:,j)) .* real(dY_mod_dy(:,k)) ...
            + (G_mod - G_exp) .* real(d2Y_mod_dy2(:,j,k))) ...
            + coef_B .* (imag(dY_mod_dy(:,j)) .* imag(dY_mod_dy(:,k)) ...
            + (B_mod - B_exp) .* imag(d2Y_mod_dy2(:,j,k))) ...
            + coef_R .* (real(dZ_mod_dy(:,j)) .* real(dZ_mod_dy(:,k)) ...
            + (R_mod - R_exp) .* real(d2Z_mod_dy2(:,j,k))) ...
            + coef_X .* (imag(dZ_mod_dy(:,j)) .* imag(dZ_mod_dy(:,k)) ...
            + (X_mod - X_exp) .* imag(d2Z_mod_dy2(:,j,k))));
    end
end

% Calculate the elements below the diagnonal of the Hessian matrix
% according to symmetry.
for k = 1 : 6
    for j = k+1 : 6
        d2E_dy2(j,k) = d2E_dy2(k,j);
    end
end
end


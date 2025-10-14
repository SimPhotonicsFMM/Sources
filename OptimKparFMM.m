function [Kx_opt, TabRelError] = OptimKparFMM(index, geom, lambda, theta, inc, varargin)

% OptimKparFMM
%   Finds the optimal value of Kx (parallel wave vector component) for a given Ky (default Ky = 0).
%   This is achieved iteratively by minimizing the relative error.
%
% Syntax
%   [Kx_opt, TabRelError] = OptimKparFMM(index, geom, lambda, theta, inc, varargin)
%   [Kx_opt, TabRelError] = OptimKparFMM(index, geom, lambda, theta, inc, 'Kx', InitialKx)
%
% Parameters:
%   index   : Refractive indices of the layers.
%   geom    : Geometry of the system (e.g., layer thicknesses).
%   lambda  : Wavelength of the incident light.
%   theta   : Incident angle (radians).
%   inc     : Incident polarization state (+1: from top, -1: from bottom).
%   varargin: Additional parameters, including:
%       'Kx'    : Initial value of Kx (required).
%       'SymY'  : Optional parameter for symmetry considerations in Y.
%
% Outputs:
%   Kx_opt      : Optimal value of Kx.
%   TabRelError : Relative error at each iteration (for convergence monitoring).
%
% Note:
%   For further details, see the description in "Spectrum.m".
%
% Author: Mondher Besbes (LCF / CNRS / IOGS)
% Date  : 01 October 2024
%

%% Parameters and Initialization
% Maximum number of iterations and convergence threshold
IterMax = 30;  % Maximum iterations
RelTol = 1e-5; % Relative error tolerance

% Initialize relative error tracking
TabRelError = [];
TabRelError(1:2) = 1; % Initial placeholder values

% Parse optional input parameters
varin = varargin;
paramNames = varin(1:2:end);

% Test geom or Mesh
if isfield(geom(1),'Cn'), Mesh = geom; end

% Determine polarization symmetry
PolSymIndex = find(ismember(paramNames, 'SymY'));
if isempty(PolSymIndex) 
    Pol = 2; % Default: consider all polarization components
else
    Pol = varin{2 * PolSymIndex};
    varin{2 * PolSymIndex} = 2; % Adjust to ensure valid input
end

% Retrieve initial Kx value
KxIndex = find(ismember(paramNames, 'Kx'));
if isempty(KxIndex)
    error('Initial value for Kx is required. Use ''Kx'' parameter.');
end

Kxi = varin{2 * KxIndex};
k0 = Kxi * [0.99, 1, 1.01]; % Initial range around the provided Kx value
fk = zeros(size(k0)); % Initialize function values for the range

%% Initial Function Evaluation
for kw = 1:2
    Kx = k0(kw);
    varin{2 * KxIndex} = Kx;
    if isfield(geom(1),'Cn')
        s = Spectrum(index, Mesh, lambda, theta, inc, varin); % Compute field
    else
        s = Spectrum(index, geom, lambda, theta, inc, varin); % Compute field
    end
    [E,H] = CalculFieldFMM(s,0,0,0);
    % Evaluate the function based on polarization
    switch Pol
        case 0
            fk(kw) = 1 / sum(E(1:3)); % TM polarization
            %fk(kw) = 1 / sum(s.CoefD(:,1));
        case 1
            %fk(kw) = 1 / sum(E(4:6)); % TE polarization
            fk(kw) = 1 / sum(s.CoefD(:,2));
        otherwise
            %fk(kw) = 1 / sum(E(:));   % Combined polarizations
            fk(kw) = 1 / sum(s.CoefD(:));
    end
end

%% Iterative Optimization
for iter = 1:IterMax
    % Update Kx and compute fields
    Kx = k0(3);
    varin{2 * KxIndex} = Kx;
    if isfield(geom(1),'Cn')
        s = Spectrum(index, Mesh, lambda, theta, inc, varin); % Compute field
    else
        s = Spectrum(index, geom, lambda, theta, inc, varin); % Compute field
    end
    [E,H] = CalculFieldFMM(s,0,0,0);
    % Update function values based on polarization
    switch Pol
        case 0
            fk(3) = 1 / sum(E(1:3));
            %fk(3) = 1 / sum(s.CoefD(:,1));
        case 1
            %fk(3) = 1 / sum(E(4:6));
            fk(3) = 1 / sum(s.CoefD(:,2));
        otherwise
            %fk(3) = 1 / sum(E(:));
            fk(3) = 1 / sum(s.CoefD(:));
    end

    % Update k0 and calculate relative error
    [k0, fk, RelError] = OptimK(k0, fk);

    % Display progress
    if iter == 1
        disp('Relative Error at each iteration:');
    end
    fprintf('Iteration %d: Relative Error = %.3e\n', iter, RelError);

    % Store relative error for monitoring
    TabRelError(iter + 2) = RelError;

    % Check convergence
    if RelError < RelTol
        fprintf('Convergence achieved in %d iterations.\n', iter);
        break;
    end
end

% Return optimal value of Kx
Kx_opt = k0(3);

end
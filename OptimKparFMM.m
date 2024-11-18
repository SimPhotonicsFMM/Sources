function [Kx_opt,TabRelError,s] = OptimKparFMM(index,geom,lambda,theta,inc,varargin)

% OptimKparFMM
%   Search mode of K_paralal (Kx) for a giving of Ky (per default Ky = 0)
%
% Syntaxe
%   [Kx_opt,TabRelError] = OptimKparFMM(index,geom,lambda,theta,inc,varargin)
%   [Kx_opt,TabRelError] = OptimKparFMM(index,geom,lambda,theta,inc,'Kx',ValKx)
%
% Description
%   See the descirption of "Spectrum.m"
%
%   Kx_opt   : Optimal valur of Kx
%   TabRelError : Relative error at each iteration

% Date of the latest version : 01 October 2024
% Author : Mondher Besbes (LCF / CNRS / IOGS)



%
%%
% Initial calculation
% -------------------

IterMax = 30;
RelTol = 1e-9;
TabRelError = [];
TabRelError(1:2) = 1;
%
varin = varargin;
f = varin(1:2:end);
P = find(ismember(f,'Kx'));
if isempty(P), error('Give an intial value for Kx'); end
%
Kxi = varin{2*P};
k0 = Kxi*[.99 1 1.01];
fk = zeros(size(k0));
%
for kw = 1:2
    Kx = k0(kw);
    %
    varin{2*P} = Kx;
    %
    s = Spectrum(index,geom,lambda,theta,inc,varin);    %
    %
    [~,P0] = max(abs(s.CoefD(:)));
    fk(kw) = 1/s.CoefD(P0);
end

% Iterative computation

for iter = 1:IterMax 
    %
    Kx = k0(3);
    %
    varin{2*P} = Kx;
    %
    s = Spectrum(index,geom,lambda,theta,inc,varin);    %
    %
    [~,P0] = max(abs(s.CoefD(:)));
    fk(3) = 1/s.CoefD(P0);
    %
    [k0,fk,RelError] = OptimK(k0,fk);
    %
    if iter == 1, disp('RelError'); end 
    disp(RelError)
    TabRelError(iter+2) = RelError;
    if RelError < RelTol, disp('Number of iterations'); disp(iter); break; end
    %return
end

Kx_opt = k0(3);

end
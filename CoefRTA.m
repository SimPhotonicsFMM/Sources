function [CoefR0,CoefT0,CoefR,CoefT,CoefA] = CoefRTA(r,t,k)

% CoefRTA
%   Extraction of diffractive order 
%       2D grating : (TM, TM-TE, TE-TM, TE) 
%
% Syntax
%   [CoefR0,CoefT0,CoefR,CoefT,CoefA] = CoefRTA(r,t); % k = 0
%   [CoefR0,CoefT0] = CoefRTA(r,t,k); % order k
%   [CoefR,CoefT] = CoefRTA(r,t);     % Total R & T 
%   [CoefR,CoefT,CoefA] = CoefRTA(r,t);    
%
% Description

%   r : Reflectivity of different orders
%   t : Transmittivity of different orders
%   k : diffractive order
%
%   CoefR0 : Reflectivity in k-order
%   CoefT0 : Transmittivity in k-order
%   CoefR  : Total Reflectivity 
%   CoefT  : Total Transmittivity
%   CoefA  : Absorptivity 

% Date of the latest version : 14 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

if nargin == 2, k = 0; end

switch size(r,2)
    case 1  % réseau 1D
        CoefR = sum(abs(r(1:size(r,1),1))); % TME
        CoefT = sum(abs(t(1:size(r,1),1)));
        CoefR0 = r((size(r,1)-1)/2+1+k,1); 
        CoefT0 = t((size(r,1)-1)/2+1+k,1);
        %
        CoefA = 1-CoefR-CoefT;

    case 2 % réseau 2D
        CoefR(1) = sum(abs(r(size(r,1)/2+1:end,1))); % TMM
        CoefT(1) = sum(abs(t(size(r,1)/2+1:end,1)));
        CoefR0(1) = r(end-(size(r,1)/2-1)/2+k,1); 
        CoefT0(1) = t(end-(size(r,1)/2-1)/2+k,1);
        %
        CoefR(2) = sum(abs(r(1:size(r,1)/2,1))); % TME
        CoefT(2) = sum(abs(t(1:size(r,1)/2,1)));
        CoefR0(2) = r((size(r,1)/2-1)/2+1+k,1); 
        CoefT0(2) = t((size(r,1)/2-1)/2+1+k,1);
        %
        CoefA(1) = 1-CoefR(1)-CoefR(2)-CoefT(1)-CoefT(2);
        %
        CoefR(3) = sum(abs(r(size(r,1)/2+1:end,2))); % TEM
        CoefT(3) = sum(abs(t(size(r,1)/2+1:end,2)));
        CoefR0(3) = r(end-(size(r,1)/2-1)/2+k,2); 
        CoefT0(3) = t(end-(size(r,1)/2-1)/2+k,2); 
        %
        CoefR(4) = sum(abs(r(1:size(r,1)/2,2))); % TEE
        CoefT(4) = sum(abs(t(1:size(r,1)/2,2)));
        CoefR0(4) = r((size(r,1)/2-1)/2+1+k,2); 
        CoefT0(4) = t((size(r,1)/2-1)/2+1+k,2); 
        %
        CoefA(2) = 1-CoefR(3)-CoefR(4)-CoefT(3)-CoefT(4);
end

if nargin == 2 && nargout <= 3
    CoefR0 = CoefR; CoefT0 = CoefT; CoefR = CoefA; 
end

end


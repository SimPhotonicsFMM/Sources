function [r,t,CoefD,R,T] = CalculFMM(Data,Mesh,Phys,N)

% CalculFMM 
%   Calculation of diffraction efficiency 
%
% Syntax
%   [r,t,CoefD,R,T] = CalculFMM(Data,Mesh);
%   [r,t,CoefD,R,T] = CalculFMM(Data,Mesh,Phys);
%   [r,t,CoefD,R,T] = CalculFMM(Data,Mesh,Phys,N);
%
% Description
%   Data : Data of the priblem
%   Mesh : Mesh of the structure (CoorN, Cn, Nsd,...CoorA, Ca, ExtAr, )
%   Phys : Material properties (Epsr, Mur, ...,Lambda0, K0, Omega, TypePol)
%
%   r : Reflectivity
%   t : Transmittivity
%   CoefD : Complex diffraction coeficients
%   R,T : Fresnel coefficients
%
% Example : Gold grating+substrat
%   ld = .85;
%   dx = .2; lix = .1; dy = .2; liy = .1; h = .03;
%   Mesh = MeshLayer(dx,lix,dy,liy,h,2,2,2); 
%   Data = SetData('Lambda0',ld,'Theta0',0,'Phi0',0,'ChampInc',-1,...
%       'mx',5,'my',5,'nh',1.333,'nb',1.7,'Indice',[IndexVal('Au',ld) 1.333]);
%   Phys = CaractMat(Mesh,Data);
%   [r,t] = CalculFMM(Data,Mesh,Phys);

% Date of the latest version : 13 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

if nargin < 3, Phys = CaractMat(Mesh,Data); end
if nargin < 4, N = 1; end 
%
Sb = CalculMatS(Data,Mesh,Phys,-1); % milieu bas
Sh = CalculMatS(Data,Mesh,Phys,+1); % milieu haut
%
if length(Mesh) == 1 && max(diff(Mesh.CoorN(:,end)))<2*eps
    [r,t,CoefD,R,T] = CalculCoefRT(Sb,Sh); % Cas d'un dioptre
else
    MatS = CalculMatS(Data,Mesh,Phys); % Matrices S des différentes couches
    if N > 1, MatS = ProdMatS(MatS,N); end
    % Calcul de R et T
    [r,t,CoefD,R,T] = CalculCoefRT(Sb,MatS,Sh);
end

%
end
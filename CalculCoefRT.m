function [r,t,CoefD,R,T] = CalculCoefRT(Sb,MatS,Sh)

% CalculCoefRT
%   Calculation of diffraction efficiency by a product of S-matrices
%
% Syntax
%   [r,t,CoefD,R,T] = CalculCoefRT(Sb,MatS,Sh);
%   [r,t,CoefD,R,T] = CalculCoefRT(Sb,Sh);
%
% Description
%   Sb   : Substrat S-matrix
%   MatS : S-Matrices of different layers
%   Sh   : Upstrat S-matrix
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
%   Data = SetData('Lambda0',ld,'Theta0',0,'Phi0',0,'ChampInc',-1,'TypePol',2,...
%       'mx',5,'my',5,'nh',1.33,'nb',1.7,'Indice',[IndexVal('Au',ld) 1.33]);
%   Phys = CaractMat(Mesh,Data);
%   Sb = CalculMatS(Data,Mesh,Phys,-1); % milieu bas
%   Sh = CalculMatS(Data,Mesh,Phys,+1); % milieu haut
%   MatS = CalculMatS(Data,Mesh,Phys); 
%   [r,t,CoefD,R,T] = CalculCoefRT(Sb,MatS,Sh);
%

% Date of the latest version : 13 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

m = size(Sb{3},1);

if nargin == 2 
    S = Sb;
    Sh = MatS;
    MatD = ProdMatS(S,Sh);
else    
    %S = ProdMatS(Sb,ProdMatS(MatS));
    S = Sb; for k=1:size(MatS,1), S = ProdMatS(S,MatS(k,:)); end
    MatD = ProdMatS(S,Sh);
end



[Pdi,Pd,Pui,Pu] = deal(Sb{7},Sb{8},Sh{7},Sh{8});
%
%MatD = ProdMatS(S,Sh);
%
CoefD = MatD{1};
if isempty(CoefD), CoefD = MatD{2}; end

Dim = size(CoefD,2);
%
if size(MatS{1},1) == 2*m
    [r,t] = deal(zeros(m,Dim));
    [R,T] = deal(zeros(m,Dim));
    %
    rt = abs(CoefD).^2;
    
    t(Pu,1) = rt(1:length(Pu),1)';
    r(Pd,1) = rt(length(Pu)+1:end,1)';
    %
    T(Pu,1) = CoefD(1:length(Pu),1).';
    R(Pd,1) = CoefD(length(Pu)+1:end,1).';
else
    %
    rt = abs(CoefD).^2;
    
    t = sum(rt(1:length(Pu),1));
    r = sum(rt(length(Pu)+1:end,1));
    %
    T = CoefD(1:length(Pu),1);
    R = CoefD(length(Pu)+1:end,1);

end


%
if size(rt,2) == 2 
    t(Pu,2) = rt(1:length(Pu),2)'; 
    r(Pd,2) = rt(length(Pu)+1:end,2)';
    %
    T(Pu,2) = CoefD(1:length(Pu),2).'; 
    R(Pd,2) = CoefD(length(Pu)+1:end,2).';

end

if isempty(Pdi), [r,t] = deal(t,r); [R,T] = deal(T,R); end  % cas incidence en bas 

return

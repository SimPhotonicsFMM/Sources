function MatS = CalculMatS(Data,Mesh,Phys,Coef)

% CalculMatS 
%   Calculation of S-Matrix MMF or HYB model
%
% Syntax
%   MatS = CalculMatS(Data,Mesh,Phys);
%   MatS = CalculMatS(Data,Mesh,Phys,Coef);
%
% Description
%   Data : Data of the priblem
%   Mesh : Mesh of the structure (CoorN, Cn, Nsd,...CoorA, Ca, ExtAr, )
%   Phys : Material properties (Epsr, Mur, ...,Lambda0, K0, Omega, TypePol)
%   Coef : +1 Upstrat S-Matrix ; -1 Substrat S-Matrix
%
%   MatS : S-matrix   [Eh;Hb] = MatS*[Eb;Hh]
%
% Example : Gold grating+substrat
%   ld = .85;
%   dx = .2; lix = .1; dy = .2; liy = .1; h = .03;
%   Mesh = MeshLayer(dx,lix,dy,liy,h,2,2,2); 
%   Data = SetData('Lambda0',ld,'Theta0',0,'Phi0',0,'ChampInc',-1,'TypePol',2,...
%       'mx',5,'my',5,'nh',1.333,'nb',1.7,'Indice',[IndexVal('Au',ld) 1.333]);
%   Phys = CaractMat(Mesh,Data);
%   Sb = CalculMatS(Data,Mesh,Phys,-1); % Substrat
%   Sh = CalculMatS(Data,Mesh,Phys,+1); % Upstrat
%   MatS = CalculMatS(Data,Mesh,Phys); 

% Date of the latest version : 13 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)


%
if nargin == 4 
    % Coef : -1 bas, +1 haut
    if size(Mesh(1).CoorN,2) == 2
        MatS = CalculMatS2D(Data(1),Mesh(1),Phys(1),Coef);
    else
        MatS = CalculMatS3D(Data(1),Mesh(1),Phys(1),Coef);
    end
else
   for kc = 1:length(Mesh)
        if size(Mesh(kc).CoorN,2) == 2
            MatS(kc,:) = CalculMatS2D(Data(kc),Mesh(kc),Phys(kc));
        else
            MatS(kc,:) = CalculMatS3D(Data(kc),Mesh(kc),Phys(kc));
        end
    end
end

end
function [Mesh,index] = MeshGrating1D(Mesh,TabIndex,Flag)

% From 2D Grating to 1D Grating

if nargin == 2, Flag = 0; end
if Flag ==0
[Mesh.CoorN(:,2), Mesh.CoorN(:,3)] = deal(Mesh.CoorN(:,3),Mesh.CoorN(:,2));
[Mesh.CoorA(:,2), Mesh.CoorA(:,3)] = deal(Mesh.CoorA(:,3),Mesh.CoorA(:,2));
[Mesh.CoorF(:,2), Mesh.CoorF(:,3)] = deal(Mesh.CoorF(:,3),Mesh.CoorF(:,2));
[Mesh.CoorV(:,2), Mesh.CoorV(:,3)] = deal(Mesh.CoorV(:,3),Mesh.CoorV(:,2));
Mesh = MoveMesh(utilMesh2(Mesh),[0 -max(Mesh.CoorN(:,2))/2 0]);
Mesh.zv = Mesh.yv;
Mesh.yv = [];
end
%
% Transform 3D mesh into slices
%
Mesh0 = Mesh;
clear Mesh;
%
z = [min(Mesh0.CoorN(:,3)) ;unique(Mesh0.CoorV(:,3))];
%
Nz = length(z);
%Mesh(Nz) = struct([]);
for k = 1:Nz-1
    Pe = find(Mesh0.CoorV(:,3)>z(k)+2*eps & Mesh0.CoorV(:,3)<=z(k+1)+2*eps);
    Mesh1(k) = ExtractMesh(utilMesh(Mesh0),Pe);
end

for k = 1:length(Mesh1), Mesh1(k).xv = [];  Mesh1(k).zv = []; end
%
%%
index = cell(1,length(Mesh1)+2);
if iscell(TabIndex)
    index{1,1} = TabIndex{1}; % supstrat
    index{1,end} = TabIndex{end}; % substrat
else
    index{1,1} = TabIndex(1); % supstrat
    index{1,end} = TabIndex(end); % substrat
end

%

for k = 1:length(Mesh1)
    NumSD = unique(Mesh1(k).Nsd);
    if iscell(TabIndex)
        for kd = 1:length(NumSD), index{1,end-k}{kd} = TabIndex{NumSD(kd)+1}; end
    else
        index{1,end-k} = TabIndex(NumSD+1); 
    end
end
%%
%% 
% Update number of facets according to subdomain number
Mesh = Mesh1;

for k = 1:length(Mesh1)
    NumSD = unique(Mesh1(k).Nsd);
    %
    for km = 1:length(NumSD)
        Pe1 = Mesh1(k).Nsd == NumSD(km);
        Mesh(k).Nsd(Pe1) = km;
        %Pf = Mesh1(k).TabNsdF{1,km}~=0;
        %Mesh(k).TabNsdF{1,km}(Pf) = km; 
    end 
    Mesh(k) = NumFacet(utilMesh(Mesh(k)));
end

% Concatenate adjacent films
%[Mesh,index] = OptimFreeForm(Mesh,index);
% Simplify mesh
%Mesh = SimplifyMesh(Mesh);
end
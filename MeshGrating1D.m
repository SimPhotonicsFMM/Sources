function Mesh = MeshGrating1D(Mesh)

[Mesh.CoorN(:,2), Mesh.CoorN(:,3)] = deal(Mesh.CoorN(:,3),Mesh.CoorN(:,2));
[Mesh.CoorA(:,2), Mesh.CoorA(:,3)] = deal(Mesh.CoorA(:,3),Mesh.CoorA(:,2));
[Mesh.CoorF(:,2), Mesh.CoorF(:,3)] = deal(Mesh.CoorF(:,3),Mesh.CoorF(:,2));
[Mesh.CoorV(:,2), Mesh.CoorV(:,3)] = deal(Mesh.CoorV(:,3),Mesh.CoorV(:,2));
Mesh = MoveMesh(utilMesh2(Mesh),[0 -max(Mesh.CoorN(:,2))/2 0]);
Mesh.zv = Mesh.yv;
Mesh.yv = [];
%
Mesh0 = Mesh;
clear Mesh;
%
z = [min(Mesh0.CoorN(:,3)) ;unique(Mesh0.CoorV(:,3))];
%
Nz = length(z);
%Mesh(Nz) = struct([]);
for k = 1:Nz-1
    Pe = find(Mesh0.CoorV(:,3)>z(k) & Mesh0.CoorV(:,3)<=z(k+1)+eps);
    Mesh(k) = ExtractMesh(utilMesh(Mesh0),Pe);
end

end
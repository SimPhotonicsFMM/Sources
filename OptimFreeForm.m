function [Mesh2,index2] = OptimFreeForm(Mesh,index)


dx = max(Mesh(1).CoorN(:,1))-min(Mesh(1).CoorN(:,1));
dy = max(Mesh(1).CoorN(:,2))-min(Mesh(1).CoorN(:,2));
%
P = []; 
for k = 1:length(Mesh)
    if all(Mesh(k).Nsd == 1), P = [P k]; end 
end

if isempty(P), Mesh2 = Mesh; index2 = index; return; end

dP = [0 find(diff(P) ~= 1) length(P)];

Mesh2 = Mesh(1:P(1)-1);
index0 = fliplr(index(2:end-1));
index2 = index0(1:P(1)-1);
%
for k = 1:length(dP)-1
    k1 = 1+dP(k);   
    h = max(Mesh(P(k1)).CoorN(:,3))-min(Mesh(P(k1)).CoorN(:,3));
    for km = k1+1:dP(k+1)
        h=h+max(Mesh(P(km)).CoorN(:,3))-min(Mesh(P(km)).CoorN(:,3));
    end
    geom0 = SetGeom('dx',dx,'dy',dy,'hc',h);
    Meshk = MeshLayer(geom0); Meshk.zv = [];
    if P(k1)>1
        Meshk = MoveMesh(utilMesh(Meshk),[0 0 max(Mesh2(end).CoorN(:,3))]);
    end
    Mesh2 = [Mesh2 Meshk];
    index2 = [index2 index0(P(k1))];
    if k<length(dP)-1
        Mesh2 = [Mesh2 Mesh(P(dP(k+1))+1:P(dP(k+1)+1)-1)];
        index2 = [index2 index0(P(dP(k+1))+1:P(dP(k+1)+1)-1)];
    else
        Mesh2 = [Mesh2 Mesh(P(end)+1:end)];
        index2 = [index2 index0(P(end)+1:end)];
    end
end

index2 = [index(1) fliplr(index2) index(end)];

end

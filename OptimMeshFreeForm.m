function Mesh1 = OptimMeshFreeForm(Mesh)

% Optimisation du maillage sur une couche contenant uniquement une seule
% inclusion

Mesh1 = Mesh;
for k = 1:length(Mesh)
    if max(Mesh(k).Nsd) == 2 && length(unique(Mesh(k).CoorV(:,2))) > 3
        Dx = max(Mesh(k).CoorN(:,1))-min(Mesh(k).CoorN(:,1));
        Dy = max(Mesh(k).CoorN(:,2))-min(Mesh(k).CoorN(:,2));
        h = max(Mesh(k).CoorN(:,3))-min(Mesh(k).CoorN(:,3));
        npx = length(unique(Mesh(k).CoorN(:,1)));
        npy = length(unique(Mesh(k).CoorN(:,2)));
        %
        Pf = find(Mesh(k).TabNsdF{1}==1 & Mesh(k).TabNsdF{2}==2);
        if length(Pf) > 8
            
        xv = Mesh(k).CoorF(Pf,1); xv = xv';
        yv = Mesh(k).CoorF(Pf,2); yv = yv';
        r = sqrt(xv.^2+yv.^2);
        th = atan2(yv,xv);
        [th,P] = sort(th);
        xv = r(P).*cos(th);
        yv = r(P).*sin(th);
        %
        ie = inpolygon(Mesh(k).CoorV(:,1),Mesh(k).CoorV(:,2),xv,yv);
        Num = unique(Mesh(k).Nsd(ie));
        %
        geom = SetGeom('dx',Dx,'dy',Dy,'xv',xv,'yv',yv,'npx',npx,'npy',npy,...
                       'npz',1,'hc',h,'NumSD',1);
        Meshk = MeshLayer(geom);
        if isfield(Mesh(k),'zv'), Meshk.zv = Mesh(k).zv; end
        Meshk = MoveMesh(utilMesh(Meshk),[0 0 min(Mesh(k).CoorN(:,3))]);
        if Num ~= 1
            ie = Meshk.Nsd == 1;
            Meshk.Nsd(ie) = 2;
            Meshk.Nsd(~ie) = 1;
            Meshk = NumFacet(utilMesh(Meshk));
        end
        Mesh1(k) = Meshk;
        end
    end
end
%
end

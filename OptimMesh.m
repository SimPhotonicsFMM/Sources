function [x0,y0] = OptimMesh(x0,y0,xv,yv)

a1 = max(x0(:))-min(x0(:));
a2 = max(y0(:))-min(y0(:));
%
[x,y] = meshgrid(x0,y0);
P1 = inpolygon(x(:),y(:),xv(1,:),yv(1,:));
x0 = unique2(util(x(P1)),1e-12);
y0 = unique2(util(y(P1)),1e-12);
%
for k = 2:size(xv,1)
    P1 = inpolygon(x(:),y(:),xv(k,:),yv(k,:));
    x0 = [x0(:) ; unique2(util(x(P1)),1e-12)];
    y0 = [y0(:) ; unique2(util(y(P1)),1e-12)];
end
%
x0 = sort([0; -a1/2; +a1/2; x0]);
x0 = unique2(util(sort([x0(:); min(xv,[],2); max(xv,[],2)])),1e-10);

y0 = sort([0; -a2/2; +a2/2; y0]);
y0 = unique2(util(sort([y0(:); min(yv,[],2); max(yv,[],2)])),1e-10);

% Mesh 2D

[x,y] = meshgrid(x0,y0);

Cn = [];
for i = 1:length(y0)-1
    i0 = (1:length(x0)-1) + (i-1)*length(x0);
    i1 = i0 + length(x0);
    Cn = [Cn; [i0' i0'+1 i1'+1 i1']];
end
    
CoorN = [];
j = 1:length(x0);
for i = 1:length(y0)
    CoorN = [CoorN; [x(i,j)' y(i,j)']];
end

% Numéro des sous-domaines

Nsd = ones(size(Cn,1),1);    % domaine extrême

[x,y] = deal(CoorN(:,1),CoorN(:,2));
[X,Y] = deal(x(Cn),y(Cn));
[xg,yg] = deal(sum(X,2)/4,sum(Y,2)/4);



% Recherche des numéros des arêtes ...
[Ca,ExtAr,CoorA,CoorF] = RechercheArete(Cn,CoorN);
TypeElmt = 'Quadrangle';

% Structure Mesh
Mesh = struct('Cn',Cn,'Ca',Ca,'CoorN',CoorN,'CoorA',CoorA,'CoorF',CoorF,...
              'Nsd',Nsd,'ExtAr',ExtAr,'TypeElmt',TypeElmt);

figure, VisuMesh(Mesh)


end

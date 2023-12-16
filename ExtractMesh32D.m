function Mesh1 = ExtractMesh32D(Mesh,Valk,Geom,ChoixAffiche)


if iscell(Valk)
    Pf = Valk{1};
    Valk = Valk{2};
else
    if Valk<0,
        Pf = abs(Mesh.CoorF(:,-Valk)-min(Mesh.CoorF(:,-Valk)))<max(abs(Mesh.CoorF(:)))/1e5;
    else
        Pf = abs(Mesh.CoorF(:,Valk)-max(Mesh.CoorF(:,Valk)))<max(abs(Mesh.CoorF(:)))/1e5;
    end
end

switch abs(Valk)
    case 1,  Dim = [2 3];
    case 2,  Dim = [1 3];
    case 3,  Dim = [1 2];
end
%

Pn = unique(Mesh.ExtFa(Pf,:));
Pa = unique(abs(Mesh.ExtFaAr(Pf,:)));

%
Cn = Mesh.ExtFa(Pf,:);
Ca = Mesh.ExtFaAr(Pf,:);
%
Mesh1.CoorN =  Mesh.CoorN(Pn,Dim);
P1 = zeros(size(Mesh.CoorN,1),1);
for k = 1:length(Pn), P1(Pn(k)) = k; end  % compresser la numérotation des noeuds
Mesh1.Cn = P1(Cn);
ExtAr = Mesh.ExtAr(Pa,:);
Mesh1.ExtAr = P1(ExtAr);

%
Mesh1.CoorA =  Mesh.CoorA(Pa,Dim);
P1 = zeros(size(Mesh.CoorA,1),1);
for k = 1:length(Pa), P1(Pa(k)) = k; end  % compresser la numérotation des arêtes
Mesh1.Ca = P1(abs(Ca)).*double(sign(Ca));

%
Mesh1.CoorF =  Mesh.CoorF(Pf,Dim);
Mesh1.Nsd = Mesh.NsdF(Pf);
if isfield(Mesh,'Surf'), Mesh1.Surf =  Mesh.Surf(Pf); end

%

P0 = Mesh1.Nsd == 0;
Mesh1.Nsd(P0) = max(Mesh1.Nsd)+1;

if ~isempty(Geom),
    Mesh1 = AddSymMesh(Geom,Mesh1,ChoixAffiche);
end



return
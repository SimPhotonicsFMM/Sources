function geom = AddGeom(geom1,geom2)

geom = geom1;

geom.hc = [geom1.hc geom2.hc];
geom.mn = [geom1.mn; geom2.mn];
geom.ab = [geom1.ab; geom2.ab];
geom.Dep = [geom1.Dep; geom2.Dep];
geom.npx = [geom1.npx; geom2.npx];
geom.npy = [geom1.npy; geom2.npy];
geom.NumSD = [geom1.NumSD; geom2.NumSD];

end
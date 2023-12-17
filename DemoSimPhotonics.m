function DemoSimPhotonics()

clc
type readme.txt

%% Example of photonic crystal
geom = SetGeom('dx',0.2,'dy',0.2,'hc',0.05,'mn',4,'ab',[.03 .03]);
%
index = {1 , {1 3.51}, 1.};
% Incident Plane Wave
lambda = 0.7;%linspace(0.4,0.8,41);  % Wavelength(µm)
theta = linspace(0,89.9,31)*pi/180;                   % Incident angle (rd)
inc = +1;                       % inc = -1: from down, +1: from top
% Spectrum Calculation
[R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc,'mx',5,'my',5);
%
% Field calculation
s = Spectrum(index,geom,lambda(1),theta(1),inc,'mx',5,'my',5);
%
[x,y,z] = deal(linspace(-geom.dx/2,geom.dx/2,100),...
               linspace(-geom.dy/2,geom.dy/2,101),...
               geom.hc/2); % (xoy)
% 
[E,H] = CalculFieldFMM(s,x,y,z);         

% Plot Field distribution
figure('Position',[400 400 1100 300]), 
subplot(131), VisuMesh(s.Mesh), colorbar off
subplot(132), hold on 
PlotCourbe(theta*180/pi,R_te*100,'','','','--r')
PlotCourbe(theta*180/pi,T_te*100,'','','','--g')
PlotCoefRTA(lambda,theta,R_tm,T_tm), legend('Rte','Tte','Rtm','Ttm')
subplot(133),hold on
VisuFieldFMM(abs(E(:,3)),x,y,z), axis tight,title('abs(Ez)') 
VisuContour(s.Mesh,3) % Nanoparticle

end
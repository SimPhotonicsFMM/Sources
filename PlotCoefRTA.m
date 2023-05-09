function PlotCoefRTA(TabLambda0,TabTheta0,R,T,A)

% PlotCoefRTA
%   Plot Reflectivity and/or Transmittivity and/or Absorptivity versus
%   wavelength and angle of incidence (curve or map)
%
% Syntax
%   PlotCoefRTA(TabLambda0,TabTheta0,R)
%   PlotCoefRTA(TabLambda0,TabTheta0,R,T)
%   PlotCoefRTA(TabLambda0,TabTheta0,R,T,A)
%
% Description
%   TabLambda0 : Wavelength (scalar or array)
%   TabTheta0 : Angle of incidence (scalar or array)
%   R : Reflectivity 
%   T : Transmittivity
%   A : Absorptivity
%
% Example
%   lambda = linspace(0.4,0.8,41); theta = linspace(50,70,61)*pi/180;  
%   geom = 0.05; index = {1.33 IndexVal('Au') 1.7}; inc = -1;
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc);
%   figure, PlotCoefRTA(lambda,theta,R_tm,T_tm,1-R_tm-T_tm)
%   figure, PlotCoefRTA(lambda,theta,R_te,T_te,1-R_te-T_te)
%   figure, PlotCoefRTA(lambda(end),theta,R_tm(end,:))

% Date of the latest version : 31 January 2022
% Author : Mondher Besbes (LCF / CNRS / IOGS)

[Xl,Xa] = ndgrid(TabLambda0,TabTheta0);
if nargin == 4, A = nan(size(R)); end
if nargin == 3, [T,A] = deal(nan(size(R))); end

set(gcf,'Color','w')
hold on

    
%
if length(TabLambda0) == 1
    StrTitle = ['\lambda = ' num2str(TabLambda0) ' µm'];
    plot(TabTheta0*180/pi,R*100,'LineWidth',2)
    %PlotCourbe(TabTheta0*180/pi,R*100,'\Theta (°)',' (%)',StrTitle,'r')
    if nargin > 3
        plot(TabTheta0*180/pi,T*100,'LineWidth',2)
        %PlotCourbe(TabTheta0*180/pi,T*100,'\Theta (°)','R , T (%)',StrTitle,'g')
        legend('R (%)', 'T (%)','Location','best');
        if nargin > 4
            plot(TabTheta0*180/pi,A*100,'LineWidth',2)
            %PlotCourbe(TabTheta0*180/pi,A*100,'\Theta (°)','R , T , A (%)',StrTitle,'b')
            legend('R (%)', 'T (%)', 'A (%)','Location','best'); 
        end
    end
    xlabel('\Theta (°)','FontWeight','bold','FontSize',12),
    ylabel('(%)','FontWeight','bold','FontSize',12),
    title(StrTitle,'FontWeight','bold','FontSize',12)

    ylim([0 100]), 

elseif length(TabTheta0) == 1
    StrTitle = ['\theta = ' num2str(TabTheta0*180/pi) '°'];
    plot(TabLambda0,R*100,'LineWidth',2)
    %PlotCourbe(TabLambda0,R*100,'\lambda (µm)',' (%)',StrTitle,'r')
    if nargin > 3
        plot(TabLambda0,T*100,'LineWidth',2)
        %PlotCourbe(TabLambda0,T*100,'\lambda (µm)','R , T (%)',StrTitle,'g')
        legend('R (%)', 'T (%)','Location','best');
        if nargin > 4
            plot(TabLambda0,A*100,'LineWidth',2)
            %PlotCourbe(TabLambda0,A*100,'\lambda (µm)','R , T , A (%)',StrTitle,'b')
            legend('R (%)', 'T (%)', 'A (%)','Location','best')
        end
    end
    xlabel('\lambda (µm)','FontWeight','bold','FontSize',12),
    ylabel('(%)','FontWeight','bold','FontSize',12),
    title(StrTitle,'FontWeight','bold','FontSize',12)
    ylim([0 100]), 
else
     if nargin > 3
        subplot(1,3,1), pcolor(Xl',Xa'*180/pi,R'*100)
        xlabel('\lambda (µm)','FontWeight','bold','FontSize',12),
        ylabel('\Theta (°)','FontWeight','bold','FontSize',12),
        title('R (%)','FontWeight','bold','FontSize',12)
        shading interp; colormap jet; colorbar, axis tight
        set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12)
        %
        subplot(1,3,2), pcolor(Xl',Xa'*180/pi,T'*100)
        xlabel('\lambda (µm)','FontWeight','bold','FontSize',12),
        ylabel('\Theta (°)','FontWeight','bold','FontSize',12),
        title('T (%)','FontWeight','bold','FontSize',12)
        shading interp; colormap jet; colorbar, axis tight
        set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12)
        %
        subplot(1,3,3), pcolor(Xl',Xa'*180/pi,A'*100)
        xlabel('\lambda (µm)','FontWeight','bold','FontSize',12),
        ylabel('\Theta (°)','FontWeight','bold','FontSize',12),
        title('A (%)','FontWeight','bold','FontSize',12)
        shading interp; colormap jet; colorbar, axis tight
    set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12)
     else
        pcolor(Xa*180/pi,Xl,R*100)
        ylabel('\lambda (µm)','FontWeight','bold','FontSize',12),
        xlabel('\Theta (°)','FontWeight','bold','FontSize',12),
        title('R (%)','FontWeight','bold','FontSize',12)
        shading interp; colormap jet; colorbar, axis tight
        set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12)
     end
         
end
box on
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12)

end
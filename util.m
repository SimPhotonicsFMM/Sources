classdef util

% util
%   class of utility methods 
%
% Properties
%   Value : scalar, vector or cell array
%
% Methods
%   max2 : as "max" Matlab function with input argument array or cell
%       y = max2(obj)
%
%   Find : as "find" Matlab function with two output arguments P and ~P 
%       [P1,P2] = Find(obj)
%
%   FlipUD : as "fipud" Matlab function with geom as input argument
%       geom = FlipUD(obj)
%
%   unique2 : as "unique" Matlab function with input argument integer or real
%       y = unique2(obj,tol)
%
%   MinMax : Minimum ans maximum elements of an array 
%       [xm,xM] = MinMax(obj)
%       [xm,xM,ym,yM] = MinMax(obj)
%       [xm,xM,ym,yM,zm,zM] = MinMax(obj)
%
%   linspace2 : as "linspace" Matlab function with extra parameters
%           Alfa0 : regular or quadratic step
%           Ext : 0: without d1 and d2 , 1: with only d1 , 2 : with only d2
%       y = linspace2(obj1, obj2, n, Alfa0, Ext)


% Date of the latest version : 14 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)    

    properties
        Value 
    end
    methods
        function obj = util(val)
            if nargin == 1
                obj.Value = val;
            end
        end
        %
        % max of array or cell
        function y = max2(obj)  
            x = [obj.Value];
            if iscell(x)
                xnum = cell2mat(x);
                y = max(xnum(:));
            else
                y = max(x(:));
            end
        end
        %
        % unique elements with tolerance
        function [y,P] = unique2(obj,tol) 
            x = [obj.Value];
            if nargin == 1, tol = 0; end
            if tol == 0
                y = unique(x);
            else
                [~,P] = unique(round(x/tol));
                y = x(P);
            end
        end 
        %
        % min and max of coordinates
        function [xm,xM,ym,yM,zm,zM] = MinMax(obj) 
            Coor = [obj.Value];
            TabMin = min(Coor);
            TabMax = max(Coor);
            switch nargout
                case 2
                    xm = TabMin;
                    xM = TabMax;
                case 4
                    [xm,ym] = deal(TabMin(1),TabMin(2));
                    [xM,yM] = deal(TabMax(1),TabMax(2));
                case 6
                    [xm,ym,zm] = deal(TabMin(1),TabMin(2),TabMin(3));
                    [xM,yM,zM] = deal(TabMax(1),TabMax(2),TabMax(3));
            end
        end
        %
        % 
        function [P1,P2] = Find(obj)
            a = [obj.Value];
            if isempty(a), P1=[]; P2=[]; return, end
            
            P1 = find(a);
            P2=find(~a); 
        
        end
        %
        function y = linspace2(obj1, obj2, n, Alfa0, Ext)

            % linspace2
            %   générer une liste de valeurs entre d1 et d2 avec options 
            %
            % Syntaxe
            %   y = linspace2(d1, d2, n, Alfa0, Ext);
            %   y = linspace2(d1, d2, n, Alfa0);
            %   y = linspace2(d1, d2, n);  % identique à linspace
            %   y = linspace2(d1, d2);  % n = 100
            %
            % Description
            %   d1 : valeur de départ
            %   d2 : valuer d'arrivée
            %   n   : nombre de points
            %   Alfa0 : pour pas variable sinon = 1
            %   Ext : prise en compte ou non des valeurs extrêmes
            %           0: sans d1 et d2 , 1: avec d1 sans d2 , 2 : sans d1 avec d2
            %
            %   y : tableau de valeurs généré
            
            % Date de la dernière version : 25 janvier 2019
            % Auteur : Mondher Besbes (LCF / CNRS / IOGS)

            d1 = [obj1.Value];
            d2 = [obj2.Value];
            
            if nargin == 2
                n = 100;
            end
            
            n = floor(n);
            
            if Alfa0 == 1,
                y = linspace(d1, d2, n);
            else
            N = n-1;
            Alfa = 0; 
            for in = 1:N, Alfa = Alfa + sqrt(Alfa0^(in-1)); end
            
            d = (d2-d1)/Alfa;
            
            y(1) = d1+d;
            for in = 2:n-1, y(in) = y(in-1) + sqrt(Alfa0^(in-1))*d; end
            y = [d1 y];
            end
            
            if nargin == 5   % Prise en compte des extrémités
                switch Ext
                    case 0, y = y(2:end-1);
                    case 1, y = y(1:end-1);
                    case 2, y = y(2:end);
                end
            end
            
        end
        %
        % Flip layers in up/down direction
        function geom1 = FlipUD(obj)
            geom0 = [obj.Value];
            FieldData = fieldnames(geom0);
            geom1 = geom0;
            for k = 1:length(FieldData)
                if iscell(geom0.(FieldData{k}))
                    geom1.(FieldData{k}) = geom0.(FieldData{k})(end:-1:1);
                end
            end
            geom1.hc = geom0.hc(end:-1:1);
            
        end
        % Toplitz
        function y = toeplitz2(obj)
            x = obj.Value;
            x = x(:);
            n = floor((length(x)+1)/2);
            y = zeros(n,n);
            for ii=1:n; y(:,ii) = x(n-ii+1:2*n-ii); end
            y = sparse(y);
        end


    end
end


function Data = InterpIndex(Data0,ld)

% InterpIndex
%   Update the value of refractice index (Data.Indice, nb, nh) 
%
% Syntax
%   Data = InterpIndex(Data,ld);
%   Data = InterpIndex(Data);
%
% Description
%   ld   : wavelength (scaler or array)
%   Data : E/S Data of the problem or Array of refractive indices
%
% Example
%   Data = SetData('Indice',{IndexVal('Au') 1.5},'nb',IndexVal('BK7'));
%   Data0 = InterpIndex(Data,0.5);
%

% Date of the latest version : 06 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

if iscell(Data0)
    if isstruct(Data0{1})
        for k = 1:length(Data0)
            Data(k) = InterpIndex(Data0{k},ld);
        end
        return
    else
        Data.TabIndice = Data0{1};
    end
end
%
if isstruct(Data0)
    Data = Data0;
    if length(Data0)>1
        for k = 1:length(Data0)
            Data(k) = InterpIndex(Data0(k));
        end
        return
    end
else
    Data.TabIndice = Data0; 
end
%    
if nargin == 1
    if isstruct(Data0), ld = Data.Lambda0; else, ld =0; end
else
    Data.Lambda0 = ld;
end
%
if iscell(Data.TabIndice)
    for k = 1:length(Data.TabIndice)
        if isa(Data.TabIndice{k},'function_handle')
            Data.Indice(k) = Data.TabIndice{k}(ld);
        elseif isa(Data.TabIndice{k},'double')
            if numel(Data.TabIndice{k}) == 1
                Data.Indice(k) = Data.TabIndice{k}(1);
            else
                Data.Indice(k) = interp1(Data.TabIndice{k}(:,1),...
                                        Data.TabIndice{k}(:,2),ld);
            end
        end
    end
else
    if isa(Data.TabIndice,'function_handle')
            Data.Indice = Data.TabIndice(ld);
    else
        if size(Data.TabIndice,2) == length(Data.Indice)+1
            for k = 1:length(Data.Indice)
                Data.Indice(k) = interp1(Data.TabIndice(:,1),...
                                         Data.TabIndice(:,k+1),ld);
            end
            if size(Data.TabIndice,2) == length(Data.Indice)+2    
                Data.Indice = Data.TabIndice(1,1:end-2);
                warning('Indice égal à la 1ère ligne de TabIndice')
            else
                Data.Indice = Data.TabIndice(1,:);
                warning('Indice égal à la 1ère ligne de TabIndice')
            end       
        end
    end
end

if isa(Data.nb,'function_handle'), Data.nb = Data.nb(ld); end
if isa(Data.nh,'function_handle'), Data.nh = Data.nh(ld); end

if ~isstruct(Data0), Data = Data.Indice; end
return
    
        
        
        
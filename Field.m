function varargout = Field(index,geom,lambda,theta,inc,varargin)


varin = varargin;
f = varin(1:2:end);
P = find(ismember(f,'Field'));
dim = length(varin);
if isempty(P) || varin{2*P} == 0
    varin{dim+1} = 'Field';
    varin{dim+2} = 1;
%     error('For Field Calculation, set parameter ''Field'' equal to 1 ' )
end
%
if length(lambda)>1 || length(theta)>1
    error('Calculation allowed for only one value of lambda and theta')
end
%
s = Spectrum(index,geom,lambda,theta,inc,varin);

[x,y,z] = deal(geom.x,geom.y,geom.z);

if isfield(s.Data,'Nsub') && sum(cell2mat({s.Data.Nsub})) ~= 0
    [E,H] = CalculFieldFD_FMM(s,x,y,z);
else
    [E,H] = CalculFieldFMM(s,x,y,z);
end

if nargout == 4
    varargout = {E(:,1:3) H(:,1:3) E(:,4:6) H(:,4:6)};
elseif nargout == 2
    varargout = {E H};
else
    varargout = {E};
end

end

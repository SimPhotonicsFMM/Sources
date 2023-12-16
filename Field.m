function varargout = Field(index,geom,lambda,theta,inc,varargin)

varin = varargin;
s = Spectrum(index,geom,lambda,theta,inc,varin);

[x,y,z] = deal(geom.x,geom.y,geom.z);

if isfield(s.Data,'Nsub')
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

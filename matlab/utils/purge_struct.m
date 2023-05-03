% Removes fields and subfields of a struct that contain
% simulink-incompatible data types. 
% 
% This is a processing step for the specific use case of passing a struct 
% as a parameter to a Simulink Matlab Function. Simulink does not accept 
% structs that have fields containing: cells, strings, empty data, or 
% timeseries. Of course this loses data, so be careful in usage!
% 
% For best results, may need to apply twice, ie: 
%       s = purge_struct(purge_struct(s))


function s = purge_struct(s)

fns = fieldnames(s);
for i = 1:length(fns)     
  fn = fns{i};
  if isstruct(s.(fn))
    if isempty(fieldnames(s.(fn)))
      s = rmfield(s, fn);
    else
      s.(fn) = purge_struct(s.(fn));  % recursive layer
    end
  elseif iscell(s.(fn)) ||  ischar(s.(fn)) || isempty(s.(fn)) || isa(s.(fn), 'timeseries')
    s = rmfield(s, fn); 
  end
end

end









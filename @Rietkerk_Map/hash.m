% Mon  6 Dec 09:43:57 CET 2021
% Karl Kästner, Berlin
%
%% has the model parameters for filename generation
%
function [key_val, key_str] = hash(obj,rt)
	hashfield_C = obj.hashfield_C;
	if (~isempty(rt.innersolver) && (~isnumeric(rt.innersolver)))
		hashfield_C{end+1} = 'innersolver';
	end

	key_str = keyfun('',rt,hashfield_C,'');
	key_val = hash_str(key_str);
 
	function key_str = keyfun(key_str,s,field_C,prefix)
		for idx=1:length(field_C)
			val = s.(field_C{idx});
			if (isstruct(val))
				key_str = keyfun(key_str,val,fieldnames(val),[prefix,'.',field_C{idx}]);
			else
				if (isstr(val))
					val_s = val;
				elseif (isa(val,'function_handle'))
					val_s = func2str(val);
				elseif (iscell(val))
					for cdx=1:length(val)
						if (isstr(val{cdx}))
							if (1==cdx)
								val_s = val{cdx};
							else
								val_s = [val_s,' ',val{cdx}];
							end
						else
							error('not yet implmented')
						end
					end % for cdx
				else
					if (length(val)<3)
						val_s = sprintf('%g',val);
					else
						% TODO vector hash?
						val_s = 'vector';
					end
				end
				if (isempty(prefix))
					%key_str = [key_str, sprintf('%s=%g;',field_C{idx},val)];
					key_str = [key_str, sprintf('%s=%s;',field_C{idx},val_s)];
					%val)];
				else
					key_str = [key_str,sprintf('%s=%s;',[prefix,'.',field_C{idx}],val_s)];
					%key_str = [key_str,sprintf('%s=%g;',[prefix,'.',field_C{idx}],val)];
				end
			end
		end
	end % keyfun
%	if (~strcmp(func2str(rt.solver),func2str(@ode23)))
%		key_str = [key_str,'-solver-',func2str(rt.solver)];
%	end
end % hash


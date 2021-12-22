% Mon  6 Dec 09:43:57 CET 2021
function [key_val, key_str] = hash(obj,rt)
	key_str = keyfun('',rt,obj.hashfield_C,'');
	key_val = hash_str(key_str);
 
	function key_str = keyfun(key_str,s,field_C,prefix)
		for idx=1:length(field_C)
			val = s.(field_C{idx});
			if (isstruct(val))
				key_str = keyfun(key_str,val,fieldnames(val),[prefix,'.',field_C{idx}]);
			else
				if (isempty(prefix))
					key_str = [key_str, sprintf('%s=%g;',field_C{idx},s.(field_C{idx}))];
				else
					key_str = [key_str,sprintf('%s=%g;',[prefix,'.',field_C{idx}],s.(field_C{idx}))];
				end
			end
		end
	end % keyfun
end % hash


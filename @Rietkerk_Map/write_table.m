% Sat 11 Dec 11:28:49 CET 2021
% write hashtable as human readable csv
function write_table(obj)

	if (~isempty(obj.map))
		filename = [obj.path_str,filesep,obj.map_str];

		% write hashtable as mat-file
		map = obj.map;
		save([obj.path_str,filesep,obj.map_str],'map');

		fid = fopen([filename(1:end-4),'.csv'],'w');
		if (fid <= 0)
			error('cannot open file for writing');
		end

		% print the header
		fprintf(fid,'hash;');
		print_names(obj.map(0),obj.hashfield_C,'');
		fprintf(fid,'\n');

		% print entries with values
		key_C = obj.map.keys;
		for idx=1:length(key_C)
			rk   = obj.map(key_C{idx});
			hash = obj.hash(rk);
			fprintf(fid,'%d;',key_C{idx}); %hash);
			print_values(rk,obj.hashfield_C);
			fprintf(fid,'\n');
		end % for idx
	
		fclose(fid);
	end % if ~isempty(obj.map)
	
	function print_names(s,field_C,prefix)
		for idx=1:length(field_C)
			val = s.(field_C{idx});
			if (isstruct(val))
				print_names(val,fieldnames(val),[prefix,'.',field_C{idx}]);
			else
				if (isempty(prefix))
					fprintf(fid,'%s;',field_C{idx});
				else
					fprintf(fid,'%s;',[prefix,'.',field_C{idx}]);
				end
			end
		end
	end % print_names

	function print_values(s,field_C)
		for idx=1:length(field_C)
			val = s.(field_C{idx});
			if (isstruct(val))
				print_values(val,fieldnames(val));
			else
				fprintf(fid,'%g;',val);
			end
		end
	end % print_values
end % write_table


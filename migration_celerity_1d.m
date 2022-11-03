% 2021-05-31 21:59:40.374630793 +0200
%
%% estimate migration celerity of a travelling wave
%
% f(tdx,idx)
% 	tdx : time index
%	idx : space index
%
% function [c,cme] = migration_celerity_1d(f,dt,dx)
function [c,cme] = migration_celerity_1d(f,dt,dx)
	nt = size(f,1);
	% assume circular boundary conditions
	df_dx = diff(mid([f,f(:,1)],1),[],2)/dx;
	df_dt = diff(mid([f,f(:,1)],2),[],1)/dt;
	c = zeros(nt-1,1);
	cme = zeros(nt-1,1);
	for idx=1:nt-1
		c(idx) = df_dx(idx,:)' \ df_dt(idx,:)';
		cme(idx) = median(df_dt(idx,:)./df_dx(idx,:));
	end
end




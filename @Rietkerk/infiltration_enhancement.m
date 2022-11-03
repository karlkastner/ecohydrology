% 2021-10-27 19:38:39.504517593 +0200
% Karl Kästner, Berlin
%
%% infiltration enhancement of the Rietkerk model
%
function ie = infiltration_enhancement(obj,b)
	ie = (b + obj.p.kb.*obj.p.w0)./(b+obj.p.kb);
end

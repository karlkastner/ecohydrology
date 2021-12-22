function ie = infiltration_enhancement(obj,b)
	ie = (b + obj.p.kb.*obj.p.w0)./(b+obj.p.kb);
end

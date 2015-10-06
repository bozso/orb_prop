function ellipse = CalcEllipse (time, tau, e, a, n)
	M = n * (time - tau);

	for iii = 1 : numel(time)
		[E(iii), sinv(iii), cosv(iii)] = newtonm_scv (e, M(iii));
	end

	r = a * (1 - e * cos(E));
	ellipse(:,1) = r .* cosv;
	ellipse(:,2) = r .* sinv;
end

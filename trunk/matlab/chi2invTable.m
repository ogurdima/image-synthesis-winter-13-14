


for p=0.01:0.01:0.99
	for v = 1:1:128
		%disp(['p=', num2str(p), ' v=', num2str(v), ': ', num2str(chi2inv(p,v))]);
		disp(['chi2inv[', num2str(p*100), '][', num2str(v), '] = ', num2str(chi2inv(p,v)), ';']);
	end
end


function [whole, frac] = whole_and_frac(num)
	whole = floor(num);
	frac = num - whole;
end

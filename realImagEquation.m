function Fout = realImagEquation(y, K, p)

x = y(1) + 1i*y(2);

F = dispersionRelation(x, K, p);
scale = max(1, K^2);
if ~isfinite(real(F)) || ~isfinite(imag(F))
    Fout = [1e30; 1e30];
else
    Fout = [real(F); imag(F)] / scale;
end
end
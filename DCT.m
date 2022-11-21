function D = DCT(n)
% DCT-dictionary function 

D = zeros(n,n);
for i = 1:n
    e      = zeros(1,n);
    e(i)   = 1;
    D(:,i) = idct2(e);
end

end
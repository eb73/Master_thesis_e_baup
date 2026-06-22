function y = rownorms(A)

% returns vector with Frobenious norms of rows of A

m = size(A,2);

y = abs(A(:,1)).^2;
for k = 2:m
    y = y + abs(A(:,k)).^2;
end
y = sqrt(y);
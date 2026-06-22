function [mep] = arma11(y)

N = length(y);
M = N-1;

y1 = y(1:N-1);
y2 = y(2:N);

I = eye(M);
T = [zeros(1,M); eye(M-1) zeros(M-1,1)] + [zeros(M-1,1) eye(M-1); zeros(1,M)];

A1 = [y2' zeros(1,M) y1' 0; zeros(1,M) y2' zeros(1,N); I zeros(M,2*M) y1; zeros(M,M) I T zeros(M,1); zeros(M,2*M) I y2];
A2 = [y1' zeros(1,2*M+1); zeros(1,M) y1' zeros(1,N); zeros(2*M,3*M+1); zeros(M,3*M) y1];
A3 = [zeros(2,3*M+1); T zeros(M,2*M+1); zeros(M,M) T 2*I zeros(M,1); zeros(M,2*M) T zeros(M,1)];
A4 = zeros(3*M+2,3*M+1);
A5 = zeros(3*M+2,3*M+1);
A6 = [zeros(2,3*M+1); I zeros(M,2*M+1); zeros(M,M) I zeros(M,M+1); zeros(M,2*M) I zeros(M,1)];

mat = {A1,A2,A3,A4,A5,A6};
mep = mepstruct(mat,2,2);
end

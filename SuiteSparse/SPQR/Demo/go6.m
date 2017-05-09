clear
Prob = UFget (620)
A = Prob.A' ;

[C,R,E,B,X, err] = spqr_gpu2 (2, A) ;

[C3,R3,E3,B3,X3, err3] = spqr_gpu3 (2, A) ;

[m n] = size (A) ;

S = A (:,E) ;

[C2, R2] = qr (S, B, 0) ;

% normalize the diagonal of R
e = spdiags (sign (full (diag (R))), 0, n, n) ;
R = e*R ;

e = spdiags (sign (full (diag (R2))), 0, n, n) ;
R2 = e*R2 ;

[i j x] = find (R) ;
ok = find (abs (x) > 1e-12) ;
R = sparse (i(ok), j (ok), x (ok), n, n) ;

[i j x] = find (R2) ;
ok = find (abs (x) > 1e-12) ;
R2 = sparse (i(ok), j (ok), x (ok), n, n) ;

e = R-R2 ;

s = abs (e) > 1e-10 ;
% nnz (R)
% nnz (R2)
% nnz(s)

figure (1) ; 
cspy (s)

% R.*s
% R2.*s

[i j x] = find (R) ;
[i j x2] = find (R2) ;

x  = sort (abs (x)) ;
x2 = sort (abs (x2)) ;
% figure (2)
% plot ( 1:length(x), log10 (x),'ro', 1:length(x2),log10(x2),  'go') ;


[i j x] = find (s) ;
fprintf ('i messed %d to %d, j mess %d to %d\n', min (i), max (i), ...
    min (j), max (j)) ;

% Useful functions for solving specific linear systems
function [G,A,B,C] = BlockThomas(A,B,C,G)
% This function solves AX=G for X using the Thomas algorithm, where A = tridiag{a,b,c}.
% The solution X is returned on exit, and {if requested}, the three diagonals of A are
% replaced by the m_ij and U.
% See also ThomasLU. Verify with ThomasTest.
%B= diag(B); C = diag(B,1); A = diag(B,-1);% Takes tri diagonal matrix and decomposes to tridiag{a,b,c}
%[G,a,b,c] = Thomas{a,b,c,G,n}
n = length(B);

for j = 1:n-1,   % FORWARD SWEEP
    [X,Amod] = Gauss(B{j}.',-A{j+1}.',length(B{n}));
    A{j+1} = X.';               % This code is just a simplified version of
    B{j+1}   = B{j+1}   + A{j+1}*C{j};        % Gauss.m {using a different notation for the
    G{j+1} = G{j+1,:} + A{j+1}*G{j};      % matrix and RHS vector}, to which the reader
end                                          % is referred for explanatory comments.

[X,Amod] = Gauss(B{n},G{n},length(B{n}));
G{n} = X;
                      
for i = n-1:-1:1, % BACK SUBSTITUTION
    [X1,Amod] = Gauss(B{i},(G{i} - C{i}*G{i+1}),length(B{n}));
    G{i} = X1;    
end
end % end of block thomas


function [G] = BlockThomasLU(A,B,C,G)
%
% This function uses the LU decomposition returned in Block Thomas
% 
% 
n = length(B);

for j = 1:n-1,
    G{j+1} = G{j+1,:} + A{j+1}*G{j}; % FORWARD SUBSTITUTION
end                                  
[X,Amod] = Gauss(B{n},G{n},length(B{n}));
G{n} = X;
for i = n-1:-1:1,
   [X1,Amod] = Gauss(B{i},(G{i} - C{i}*G{i+1}),length(B{n}));  % BACK SUBSTITUTION
    G{i} = X1;           
end
end % function Block ThomasLU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test code for Block Thomas.

disp('Now testing Thomas & ThomasLU on a random tridiagonal matrix.')
n=3; m=4;   % note that m is number of blocks and n is size of each block
BlockA = {};
BlockF = {}; ThomasA{1,1} = zeros(n);
for i = 1:m
    for j = 1:m  %this iterates through and and creates each block which is filled with nxn random vectors
        if i == j
            F=randn(n,1);
            Gnew=randn(n,1);
            Gnew1{j,1} = Gnew;
            BlockF{j,1} = F;
            a=randn(n,1); b=randn(n,1); c=randn(n,1); a(1)=0; c(n)=0;
            A = diag(a(2:n),-1)+diag(b(1:n),0)+diag(c(1:n-1),1);
            BlockA{i,j} = A;
            ThomasB{1,j} = A;
        elseif i == j + 1
            b=randn(n,1);
            B = diag(b(1:n),0);
            BlockA{i,j} = B;
            ThomasA{1,j+1} = B;
        elseif i == j - 1
            b=randn(n,1);
            C = diag(b(1:n),0);
            BlockA{i,j} = C;
            ThomasC{1,j-1} = C;
        else
           C = zeros(n,n);
           BlockA{i,j} = C;
        end
    end
end
ThomasC{1,m} = zeros(n);

[X,amod,bmod,cmod] = BlockThomas(ThomasA,ThomasB,ThomasC,BlockF);  %runs blockThomas calling our test random blocks
A_times_X_from_Thomas=AA*cell2mat(X), cell2mat(BlockF),

AA = cell2mat(BlockA); F = cell2mat(BlockF); x = AA\F
error=norm((AA*cell2mat(X))-F) % checks norm of our A*B subtracted from origional b values to create error

[Xnew]=BlockThomasLU(amod,bmod,cmod,Gnew1); % tests out Thomas block leveraging the LU decomposition
A_times_Xnew_from_ThomasLU=AA*cell2mat(Xnew), cell2mat(Gnew1), error_LU=norm(AA*cell2mat(Xnew)-cell2mat(Gnew1))

% end script Block Thomas Test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,y] = NLeggedThomas(A,B,C,G,D,n,p,e)

for i = 1:n
    a = A{i}; b = B{i}; c = C{i}; g = G{i};
    for j = 1:p-1,                               % FORWARD SWEEP
        a(j+1)   = - a(j+1) / b(j);               % This code is just a simplified version of
        b(j+1)   = b(j+1)   + a(j+1)*c(j);        % Gauss.m (using a different notation for the
        g(j+1,:) = g(j+1,:) + a(j+1)*g(j,:);      % matrix and RHS vector), to which the reader
    end
    A{i} = a; B{i}=b; C{i} = c; G{i} = g;
end

Bp = []; Cp =[]; Gp = [];
for i = 1:n
    pp = B{i}; ppp = pp(p); Bp = [Bp ppp]; %Just iterates through so we can get our Bp and Cp values
    cc = C{i}; ccc = cc(p); Cp = [Cp; ccc];
    gg = G{i}; ggg = gg(p); Gp = [Gp; ggg];
    
end
Gp = [Gp; e];

Solver = diag(Bp(1:n),0); %creates matrix with diagonal Bp
Solver = [Solver Cp];       %adds on a column of Cp
Solver = [Solver ; D];       %adds on a row of D

[Xp,A] = Gauss(Solver,Gp,n+1);   %run gauss of our combo matrix

y = Xp(n+1);
                      % BACK SUBSTITUTION
for i = 1:n
     b = B{i}; c = C{i}; g = G{i};
     g(p+1) = y;
    g(p) = Xp(i); 
    for j = p-1:-1:1,
        g(j) = (g(j) - c(j) * g(j+1) ) / b(j);
    end
    X{i} = g;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NLeggedThomasTest
clear
disp('Now testing NLeggedThomasTest.')
n = 5; % number of simultaneous matrix
p=12; m=1;     % matrix size p x p+1
AA_fat = {}; A={}; B={}; C={};
G_cell = {};
for i=1:n
   
    a=randn(p,1); b=randn(p,1); c=randn(p,1); G=randn(p,m); a(1)=0;
    A{i} = a; B{i} = b; C{i} = c;
    A_full=diag(a(2:p),-1)+diag(b(1:p),0)+diag(c(1:p-1),1);
    A_fat = [A_full zeros(p,1)]; A_fat(p,p+1) = c(p);
    AA_fat{i} = A_fat;
    G_cell{i} = G;
    
    A_Test{i} = ones(p,1); B_Test{i} = -2*ones(p,1); C_Test{i} = ones(p,1);
    G_Test{i} = [randn zeros(1,p-1)]';
end
D = randn(1,n+1);
%D = ones(1,n); D = [D, -2];
e = randn;
%e=0;
A_Test;

[X] = NLeggedThomas(A,B,C,G_cell,D,n,p,e);
X
%[X,y] = NLeggedThomas(A_Test,B_Test,C_Test,G_Test,D,n,p,e);
XX = cell2mat(X);
%e
%y

for i=1:n
    plot(XX(:,i)); hold on
end % end of test script for N legged thomas test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, error_V, resi] = NewAlgorithm(f,df,ddf,x,e)
% New fuction for Newton-Raphson method

tol=1e-10; residual=2*tol;
X = [x];
resi = [];

while (residual>tol)
    
    if abs(e) > .01
        x =double( x - (df(x)/ddf(x)) + sqrt((df(x)^2 / ddf(x)^2) - 2*(f(x)/ddf(x))));
        X = [X; x];
    else
        x =double( x - (f(x) / df(x)) - .5*(f(x)^2*ddf(x) / df(x)^3));
        X = [X; x];
    end
    
   residual=norm(double(f(x)));
   resi = [resi; residual];
end
error_V = [];
for i=1:max(size(X))-1
    
    error = norm(X(i) - X(max(size(X))));
    error_V = [error_V; error];
end
   
end % end of new function for Newton-Raphson method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NewAlgorithm_Test
%script for new  Newton-Raphson test.

syms x
%f(x) = 2*x^3+x^2+1
f(x) = 256*x^9 - 576*x^7 + 432*x^5 - 120*x^3 + 9*x;  %input function to test

df = diff(f,x);     % f'
ddf = diff(df,x);   %f''

ep(x) = (-2*f*ddf)/df;
x=.9                    % Guess value
e = ep(x)               % epilon value

[X, error_V, resi] = NewAlgorithm(f,df,ddf,x,e)     % new Newton-Raphson

[X_NR, error_V_NR, resi_NR] = NewRap(f,df,ddf,x,e)  %runs Newton-Raphson

figure 
plot(log(error_V),'r')
hold on
plot(log(error_V_NR))
title ' log of error, guess x = .9, red is new fnx'

figure 
plot(log(resi), 'r')
hold on
plot(log(resi_NR))
title ' log of residual, guess x = .9, red is new fnx'

%end of script for new  Newton-Raphson test.


function [X, error_V, resi] = NewRap(f,df,ddf,x,e)
% edited simplifid and wrote my own code for Newton-Raphson test.

tol=1e-10; residual=2*tol;
X = [x];
resi = [];

while (residual>tol)
    
        x =double( x - f(x)/df(x));
        X = [X; x];
    
   residual=norm(double(f(x)));
   resi = [resi; residual];
end
error_V = [];
for i=1:max(size(X))-1
    
    error = norm(X(i) - X(max(size(X))));
    error_V = [error_V; error];
end
   
end % end of code for my Newton-Raphson test.
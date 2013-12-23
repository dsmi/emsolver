function test_gpof
% Test the Generalized Pencil of Function Method implementation.

ns = 101; % Number of samples
x1 = 0;   % Beginning of the interval
x2 = 10;  % End of the interval
dx = (x2-x1)/(ns - 1); % Sampling period
x = linspace(x1, x2, ns);
y = cos(x)+cos(2*x)+cos(4*x)+cos(8*x);

[s,b]=gpof(y,dx);

s_test = [ 0+8j; 0-8j; 0+4j; 0-4j; 0+2j; 0-2j; 0+j; 0-j ];
b_test = repmat(0.5, 8, 1);

[ss, ix] = sort(s);
[sts, ixt] = sort(s_test);
assertEquals(s_test(ixt), s(ix), 1e-11);
assertEquals(b_test(ixt), b(ix), 1e-11);

yt = sum(repmat(b,1,ns).*exp(s*x),1);
assertEquals(y, yt, 1e-11);


ns = 101; % Number of samples
x1 = 0;   % Beginning of the interval
x2 = 10;  % End of the interval
dx = (x2-x1)/(ns - 1); % Sampling period
x = linspace(x1, x2, ns);
y = sin(x)+ sin(3*x)+ sin(7*x);

[s,b]=gpof(y,dx);

s_test = [ 0+7j; 0-7j; 0+3j; 0-3j; 0+1j; 0-1j ];
b_test = [ 0-0.5j; 0+0.5j; 0-0.5j; 0+0.5j; 0-0.5j; 0+0.5j ];

[ss, ix] = sort(s);
[sts, ixt] = sort(s_test);
assertEquals(s_test(ixt), s(ix), 1e-11);
assertEquals(b_test(ixt), b(ix), 1e-11);

yt = sum(repmat(b,1,ns).*exp(s*x),1);
assertEquals(y, yt, 1e-11);

ns = 10; % Number of samples
x1 = 0;  % Beginning of the interval
x2 = 9;  % End of the interval
dx = (x2-x1)/(ns - 1); % Sampling period
x = linspace(x1, x2, ns);
y = [ 1 2 3 43 5 6 7 8 9 10 ];

[s,b]=gpof(y,dx);

s_test = [ -0.20412043-3.7232028E-17j; -1.031164+1.5920029j; -1.031164-1.5920029j ];
b_test = [ 26.65605+0.0j; -11.704902+32.226086j; -11.704902-32.226086j ];

[ss, ix] = sort(s);
[sts, ixt] = sort(s_test);
assertEquals(s_test(ixt), s(ix), 1e-6);
assertEquals(b_test(ixt), b(ix), 1e-6);

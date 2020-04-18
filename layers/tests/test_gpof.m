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
y = sin(x/10);

[s,b]=gpof(y,dx);

yt = sum(repmat(b,1,ns).*exp(s*x),1);

yt = sum(repmat(b,1,ns).*exp(s*x),1);
assertEquals(y, yt, 1e-11);


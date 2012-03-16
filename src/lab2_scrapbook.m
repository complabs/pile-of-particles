%% Lab2 scrap-book
% Contains various fragments of code for various tests
%
% Filename: lab2_scrapbook.m
% Date:     2012-02-28
% Author:   Mikica B Kocic

%=========================================================================================
% Fourier transform of the particle height

% load MAT-file with results first

L = 8192;
NFFT = 8192; % Next power of 2 from length of signal
Fs = 1 / params.h_soft;   % Sampling frequency
f = Fs / 2 * linspace( 0, 1, NFFT/2 + 1 );

figure;
ph = loglog( NaN, NaN );

% title( ['Single-Sided Amplitude Spectrum of {\itz}_1 at {\itt} = ', ...
%    num2str(i/Fs) , ' s'] )

th = title( '', 'Interpreter', 'none' );
xlabel('Frequency {\itf} /Hz', 'FontName', 'Times', 'FontSize', 12 );
ylabel('| {\itz}_1({\itf}) | /m', 'FontName', 'Times', 'FontSize', 12 );

ylim([ 10^-5, 0.05 ]);
xlim([ 0.2, 1000 ]);

for i = 1:5:track.n-NFFT

    % Plot single-sided amplitude spectrum
    signal = track.X( i + (0:NFFT-1), 1 );
    Y = fft( signal, NFFT ) / L;
    set( ph, 'XData', f, 'YData', 2 * abs( Y( 1:NFFT/2+1 ) ) );

    % Update title & display
    set( th, 'String', [ num2str( i / 1000 ), ' s' ] );
    drawnow;
    
end

return

%=========================================================================================

% Semi-implicit Euler vs Leapfrog

x =   0.    ; % position, m
v =  10.    ; % velocity, m/s
m =   1.    ; % mass, kg
g = -10.    ; % acceleration of gravity, m/s^2
h =   0.01  ; % time-step size, s
n =   0     ; % time-step counter

E_0 = 0.5 * v^2 - m * g * x; % Initial total energy

fprintf( 'Initial: n = %3d, x_0 = %4g m, v_0 = %3g m/s, E_0 = %3g J\n', n, x, v, E_0 );

while x >= 0   % Solve eqs until x drops bellow 0
    n = n + 1;
    
    if false
        % SEMI-IMPLICIT EULER:
        v = v + h * g;
        x = x + h * v;
    else
        % LEAPFROG:
        v = v + g * h/2;
        x = x + v * h;
        v = v + g * h/2;
    end
end

E_n = 0.5 * v^2 - m * g * x; % Final total energy
fprintf( 'Final:   n = %3d, x_n = %4g m, v_n = %3g m/s, E_n = %3g J\n', n, x, v, E_n );

err_E = - m * g^2 * h^2 * n / 2; % Solved E_n - E_0 from difference equation
fprintf( 'E_n - E_0 = %g, Solved Err_E = %g\n', E_n - E_0, err_E );

return

%=========================================================================================

% slower:

    E_k_n = 0.5 * ( ( M(:) .* V(:) )' * V(:) );   % Kinetic energy
    
% faster

    E_k_n = 0.5 * sum( sum( ( M .* V ) .* V ) );   % Kinetic energy

%=========================================================================================

clear all
load envx 
[top,bot] = envelope( x );

close all
figure
hold on
plot( 1:n, x, 'b' )
plot( 1:n, top, 'r' )
plot( 1:n, bot, 'g' )

%=========================================================================================

X1= rand( 10000, 3 );
X2 = rand( 10000, 3 );
v = rand( 1, 3 );

%=========================================================================================

tic
for i = 1:10000
    dVp = repmat( v, size(X1,1), 1 );
end
toc
% Elapsed time is 4.677915 seconds.

tic
XZ = zeros( size(X1) );
for i = 1:10000
    dVp = bsxfun( @minus, v, XZ );
end
toc
% Elapsed time is 2.710573 seconds.

tic
for i = 1:10000
    dVp = bsxfun( @minus, v, zeros( size(X1) ) );
end
toc
% Elapsed time is 2.745042 seconds.

clear XZ  % clearing memory helps

tic
for i = 1:10000
    dVp = bsxfun( @minus, v, zeros( size(X1) ) );
end
toc
% Elapsed time is 1.109519 seconds. <<<<<<<< Fastest

%=========================================================================================

tic
for i = 1:10000
    Y = bsxfun( @times, X1, X2 );
end
toc
% Elapsed time is 1.249867 seconds.

tic
for i = 1:10000
    Y = X1 .* X2;
end
toc
% Elapsed time is 0.776563 seconds. <<<<<<<<< Faster

return

%=========================================================================================
% Euclidean distance matrix
% Assume X and Y are an m-by-p matrix representing m points in p-dimensional space. Then, 
% to compute the m-by-m distance matrix D where D(i,j) is the Euclidean distance X(i,:) 
% between X(j,:), use:
%
Dij = repmat( permute( X, [1 3 2] ), [1 m 1] ) ...
    - repmat( permute( X, [3 1 2] ), [m 1 1] );

D = sqrt( sum( abs(Dij).^2, 3 ) );     % distances (= norm)

Nij = Dij ./ repmat( D, [ 1 1 3 ] );   % unit vectors

%=========================================================================================
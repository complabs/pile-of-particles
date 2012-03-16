%% Lab2: Script that analyzes results from simulations 2.3.*
%
% Note: Load MAT-file with the results first!
%
% Filename: lab2_sim23.m
% Date:     2012-03-14
% Author:   Mikica B Kocic

np = 10;       % Analyze lowest 10 particles
dp = 1:4:np;   % Display results for even-numbered particles

% Calculate avg particle diameter as difference between two particle centers
% Particle diameter of particle with height B is: ( C - B ) / 2 + ( B - A ) / 2,
% where A and C are heights of particles bellow and above particle B, respectively.
% This means that particle B diameter is ( C - A ) / 2.
Diameter = ( track.X(:,3:np ) - track.X(:,1:np-2 ) ) / 2;  % for particles 2, 3...
Diameter = [ track.X(:,1) - (-0.3), Diameter(:,:) ]; % for particle #1

% Calculate oscilations of the particle diameter compared to the reference linear fit
Derr = zeros( size( Diameter ) );
for i = dp
    % Linear fit of the diameter of the particle as a function of time
    k = polyfit( track.T, Diameter(:,i), 1 );
    % Calculate 'error' from the linear fit
    Derr(:,i) = Diameter(:,i) - ( k(1) * track.T + k(2) );
end

% We have our data. Now, create a figure

figure( 'Name', 'Lab2: 2.3.*', ...
    'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
    'PaperSize', [ 17, 12 ], 'PaperPosition', [ 0, 0, 17, 12 ] ...
);

% Display particle diameter in the upper subplot

axTop = subplot( 211 );
title( 'Particle Diameter', 'FontSize', 11 );
hold on, grid on
ylabel( '{\itD} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
ylim([ 0.29, 0.61 ]); % Expected diameter is between r and 2r

plot( track.T, Diameter(:,dp) );  % --------- Plot diameter

% Use particle numbers as legend
legend( cellfun( @num2str, num2cell(dp), 'uniformoutput', 0 ), ...
    'FontSize', 8, 'EdgeColor', 'w', 'Location', 'SouthWest' );

% Display particle diameter oscillations in the lower subplot

axBot = subplot( 212 );
title( 'Diameter Oscillations', 'FontSize', 11 );
hold on, grid on
ylabel( '\Delta{\itD} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
ylim([ -0.031, 0.031 ]);

plot( track.T, Derr(:,dp) );  % --------- Plot diameter oscillations

% Linked x-axes with the common label (time)
xlabel( '{\itt} /{\rms}', 'FontSize', 12, 'FontName', 'Times' );
linkaxes( [ axTop, axBot ], 'x' );

% Export figure as EPS
print( gcf, '-r300', '-depsc2', [ 'lab2_fig23_', params.ID ] );

% Temporary variables (in order of appearance)
%clearvars np dp Diameter i k Derr axTop axBot


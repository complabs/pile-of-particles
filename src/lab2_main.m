%% Lab2: Colliding Particles; The Main Simulation Engine
% Simulation of a number of spherical objects bouncing in a 1-D, 2-D or 3-D box with
% either the impulse method (instantenous transfer of momenta) or the penalty method
% (linear springs) as the contact model.
%
% Filename: lab2_main.m
% Date:     2012-03-11
% Author:   Mikica B Kocic 
%
% Usage:
%   [params,track] = lab2_main
%           runs simulation using default parameters/constants.
%
%   [params,track] = lab2_main( params )
%           runs simulation using parameters/constants in 'params' structure.
%           See lab2_defs function for declaration of 'params' structure.
%           (lab2_defs serves as the constructor of the structure).
%
%   [params,track] = lab2_main( 'MAT-filename' )
%           loads 'params' structure from MAT-file and runs simulation.
%
%   varargout = lab2_main( '@fooname', x1, ..., xn )
%           evaluates the nested function, fooname, using arguments x1 through xn.
%
% Returns:
% * params: Original 'params' structure with the results (the final state variables)
%           added to the structure (so it could be used in further simulations).
% * track:  A structure with trajectory of the system and total energy track over time.
%
% Glossary:
% * Jolt:        An instantenous change of the state variable during a time-step.
% * Projection:  The component of one vector that is in the direction of the other.
% * Trajectory:  A path in a phase plot that shows how the state of a system
%                changes over time (a trace of the state variables of a system).

function varargout = lab2_main( varargin )

    if nargin && ischar( varargin{1} ) && varargin{1}(1) == '@'

        % Nested function callback
        callback = eval( varargin{1} );
        if nargout
            [varargout{1:nargout}] = callback( varargin{2:nargin} );
        else
            callback( varargin{2:nargin} );
        end
        return

    elseif ~nargin

        % Demo mode, when called without arguments
        params = lab2_defs ();
        params.Title = 'Lab2: Demo';
        params.N_p = 1;
        params.T_emit = 1;
        params.PAsSpheres = true;
        params.PlotContacts = true;
        params.PlotEnergy = true;
        params.PlotMomentum = true;
        params.PlotHeight = true;
        params.t_f = 10;

    elseif ischar( varargin{1} )

        % If argument is a string, assume that it is a MAT-filename and try to load
        % 'params' structure from it
        try
            params = load( varargin{1}, 'params' );
            params = params.params;
        catch err
            disp( getReport( err ) );
        end

    else
        % Assume the first argument is 'params' structure
        params = varargin{1};
    end

    assert( isstruct( params ) && isfield( params, 'Title' ), ...
        [ '\nArgument to function %s must be a valid ''params'' structure\n' ...
          'or a MAT-filename that contains a ''params'' structure;\n', ...
          'See <a href="lab2_defs">lab2_defs</a> for the definition of ''params''.' ...
          ], ...
        mfilename );

%=========================================================================================
%% Initialization
% Create the data structures for the simulation and assigning them initial data.

    % Global flag used to cancel main simulation loop
    global lab2_running
    lab2_running = true;  % The simulation will run while lab2_running is true

    %%% Dump running parameters ----------------------------------------------------------

    fprintf( '\n---------------------- lab2_main\n\n' );
    params.ID = datestr( now, 'YYYYmmdd_HHMMSS' ); % Timestamp the simulation

    disp( 'Parameters:' ), disp( params );
    disp( 'Enumerations:' ), disp( params.enum );

    %%% Nickname some parametric definitions and constants to gain speed -----------------

    enum    = params.enum;     % Enumeration constants
    N_dim   = params.N_dim;    % Number of space dimensions
    N_p     = params.N_p;      % Number of particles
    T_emit  = params.T_emit;   % Emitter period between two particles, s (may be Inf)
    g_n     = params.g_n;      % Standard acceleration of Gravity (intensity), m/s^2
    k_s     = params.k_s;      % Spring coefficient
    k_d     = params.k_d;      % Spring damping factor
    k_e     = params.k_e;      % Impulse collision restitution (0 == max dissipation)
    k_p     = params.k_p;      % Position projection fraction (1 == max projection)
    t_0     = params.t_0;      % Initial time, s
    t_f     = params.t_f;      % Final (simulation) time, s
    imodel  = params.imodel;   % Interaction model
    stepper = params.stepper;  % Integrator method
    h_soft  = params.h_soft;   % Time-step used in linear spring model, s
    h_hard  = params.h_hard;   % Time-step used in impulse model, s

    %%% Assert the number of space dimensions --------------------------------------------

    N_dim = round( N_dim );
    assert( 1 <= N_dim && N_dim <= 3, ...
        '\n%s: The number of space dimensions must be either 1, 2 or 3', mfilename );

    % Define a 'zero' spatial row-vector

    ZeroVec = zeros( 1, N_dim );

    % Setup acceleration of gravity as a vector and check number of dimensions
    % Assume gravity along the vertical axis (always denoted as Z-axis).

    switch N_dim
        case 1,  g = -g_n;
        case 2,  g = [ 0, -g_n ];
        case 3,  g = [ 0, 0, -g_n ];
    end

    % Initialize 'L' as spatial vector holding dimensions of the system (box)
    % NOTE 1) The height is always the last dimension in L, i.e. L(end).
    % NOTE 2) Box always starts at [0,0,0] and ends at L
    % NOTE 3) If L(i) is negative, it is considered as 'non-limiting' box side
    %         i.e. the box side without the surface that collides with particles

    switch N_dim
        case 1,  L = [ params.BoxHeight ];
        case 2,  L = [ params.BoxWidth, params.BoxHeight ];
        case 3,  L = [ params.BoxWidth, params.BoxDepth, params.BoxHeight ];
    end

    %%% Initialize surfaces --------------------------------------------------------------

    % A surface is defined by some point on the surface XSFC, and 
    % the unit vector normal to the surface NSFC. Note that norm(NSFC) must be 1.
    % The following surfaces enclose the system in a box given by sizes in 'L'.
    % Note that negative values in L will not be considered as 'limiting' surfaces.

    XSFC = [ zeros( N_dim ); diag( L ) ];    % Points on boundary surfaces
    NSFC = [ eye( N_dim ); -eye( N_dim ) ];  % Unit vectors on boundary surfaces

    % Remove surfaces that are negative, i.e. not limiting the box
    [infRows,~] = find( XSFC < 0 );
    infRows = unique( infRows );
    XSFC(infRows,:) = [];
    NSFC(infRows,:) = [];
    clearvars infRows
    
    L = abs( L ); % Finally, make all dimensions positive
    
    %%% Initial state variables: particle mass, radius, position, velocity ---------------

    % Don't allow to add interactively particles in 3-D
    if params.ManualNp && N_dim == 3
        if nargout >= 1, varargout{1} = params;   end
        if nargout >= 2, varargout{2} = struct(); end  % empty trajectory/energy track
        lab2_running = false;
        return
    end

    % Reset number of particles, if particles should be place interactivelly
    if params.ManualNp
        N_p = 0;
    end

    % ---------------------------- Initialize particle masses
    if isfield( params, 'M' ) && all( [ N_p, N_dim ] <= size( params.M ) )
        % Use provided particle masses (reduced to our dimensions, if larger)
        M = diminishDimensions( params.M, N_p, N_dim );
    else
        % Use default particle masses
        M = params.m_p * ones( N_p, N_dim );
    end

    % ---------------------------- Initialize particle radii
    if isfield( params, 'R' ) && N_p <= size( params.R, 1 )
        % Use provided particle radii (reduced to our dimensions, if larger)
        R = params.R(1:N_p,:);
    else
        % Use default particle radii
        R = params.r_p * ones( N_p, 1 );
    end

    % ---------------------------- Initial positions
    if isfield( params, 'X' ) && all( [ N_p, N_dim ] <= size( params.X ) )
        % Use provided particle positions (reduced to our dimensions, if larger)
        X = diminishDimensions( params.X, N_p, N_dim );
    else
        X = zeros( N_p, N_dim ); % Allocate space
        
        %%% Pile up particles in a box (in a grid) with a small distances inbetween

        if N_p
            d_max = 2 * max( max( R ) ) + params.dx_emit; % Raster of the grid
            i_bound = floor( L / d_max ) - 1; % Number of particles per grid dimension
            i_np = prod( i_bound );           % Maximal number of particles in the grid
            i_cur = zeros( 1, N_dim );        % Current position in the grid, starting from 0

            % Set grid reference position (bottom-left = position of the first particle)
            if params.dx_emit > 0
                X(1,:) = params.dx_emit + ( 1 + rand( 1, N_dim ) ) * d_max/2;
            else
                X(1,:) = R(1) .* ones( 1, N_dim ); % put on the ground, if dx_emit == 0
            end
            % however, if a single particle, set it it's hight is fully random
            if N_p == 1
                X(end,:) = ( 0.5 + 0.5 * rand(1) ) * ( L(end) - d_max );
            end
        end

        % Distribute particles snapping to the grid
        for i = 1:N_p
            X(i,:) = X(1,:) + i_cur * d_max;
            
            i_cur(end) = i_cur(end) + 1; % Increment position in the grid along z-axis
            if i >= i_np
                % Grid is full; continue filling along z-axis...
                % (this also suppresses filling along non-existing dimensions)
                continue
            elseif i_cur(end) >= i_bound(end)
                % Reached z-axis bound; now continue along x-axis...
                i_cur(end) = 0;
                i_cur(1) = i_cur(1) + 1;
                if i_cur(1) >= i_bound(1)
                    % Reached x-axis bound; now continue on y-axis...
                    i_cur(1) = 0;
                    i_cur(2) = i_cur(2) + 1;
                end
            end
        end
        
        clearvars d_max i_bound i_np i_cur
    end

    % ---------------------------- Initialize particle face colors

    if isfield( params, 'faceColor' ) && N_p <= size( params.faceColor, 1 )
        % Use provided face colors (reduced to our N_p, if larger)
        faceColor = params.faceColor(1:N_p,:);
    else
        % Generate particle face colors: pastel, with particle# as hue
        % Note that 4-th component of faceColor is alpha (transparency)
        faceColor = zeros( N_p, 4 );
        for i = 1:N_p
            faceColor(i,:) = randomFaceColor( ( i - 1 ) / N_p );
        end
    end

    % ---------------------------- Initial velocities
    if isfield( params, 'V' ) && all( [ N_p, N_dim ] <= size( params.V ) )
        % Use provided particle velocities (reduced to our dimensions, if larger)
        V = diminishDimensions( params.V, N_p, N_dim );
    else
        V = zeros( N_p, N_dim ); % Allocate space
        
        if g_n == 0
            % Generate random velocities without, if no gravity
            V = bsxfun( @times, -1 + 2 * rand( N_p, N_dim ), L );
        elseif N_dim >= 2
            % Induce very small velocity disturbances in horizontal plane
            V(:,1:end-1) = -params.dv_emit + 2 * params.dv_emit * rand( N_p, N_dim-1 );
        end
    end

    %%% Initialize calculated quantities -------------------------------------------------

    F_tot  = bsxfun( @times, M, g ); % total force,
    A      = F_tot ./ M;             % and acceleration.

    VJ_tot = zeros( N_p, N_dim );    % Velocity deflections (jolt)
    XJ_tot = zeros( N_p, N_dim );    % Position projections

    %%% Initialize derived quantities ----------------------------------------------------

    % Compute the linear momentum of the system

    P = M .* V;

    % Compute the initial kinetic, potential and total energy of the system

    E_k_n = 0.5 * sum( sum( P .* V ) );  % Kinetic energy
    E_p_n = - sum( sum( X .* bsxfun( @times, M, g ) ) ); % Potential energy
    E_tot_n = E_k_n + E_p_n;  % Total energy
    
    %%% Setup the integrator variables ---------------------------------------------------

    assert( stepper == enum.ForwardEuler      ...
         || stepper == enum.SemiImplicitEuler ...
         || stepper == enum.Leapfrog, ...
        '\nUnknown integration method; stepper = %d', stepper );

    if stepper == enum.ForwardEuler && imodel == enum.ImpulseModel
        throw( MException( 'lab2_main:feim', ...
            'Forward Euler integrator is not supported for the Impulse Model' ...
        ) );
    end

    % Initial simulation time and time-step counter
    t = t_0;

    % Time-step depends on collision model:
    switch imodel
        case enum.SpringModel,  h = h_soft;  % Time-step size for linear spring model
        case enum.ImpulseModel, h = h_hard;  % Time-step size for impulse model
    end
    assert( h > 0 );
    params.h = h;  % Remember selection

    % Number of time steps must be integer
    NT = ceil( ( t_f - t_0 ) / h );

    %%% Allocate vectors/matrices for traced variables -----------------------------------

    track = struct( 'n', 1 ); % Structure that will contain traced system variables
                              % where track.n = 1 corresponds to t_0

    % Do we need to track particle trajectory?
    traceXP = nargout >= 2 || params.SaveData || params.PlotPath ...
              || params.PlotContacts || params.PlotMomentum || params.PlotHeight;

    % Alocate matrices holding position and linear momentum trajectories
    if traceXP
        % Elapsed time: T(i) = t_0 + i * h
        track.T = []; % Just declare the field--it will be defined after the main loop.

        % Position trajectory in phase-space
        track.X = zeros( [ NT + 1, N_p, N_dim ] );
        track.X( track.n, :, : ) = X; % Initial positions

        % Linear momentum trajectory in phase-space
        track.P = zeros( [ NT + 1, N_p, N_dim ] );
        track.P( track.n, :, : ) = P; % Initial linear momenta
        
        % Number of collisions per time-step where:
        % 1st column contains particle-to-particle collision count and 
        % 2nd column contains particle-to-surface collision count
        track.colc = zeros( NT + 1, 2 );
    end

    % Do we need to track energy of the system over time?
    traceE  = params.SaveData || params.PlotEnergy;

    % Allocate vectors holding energy trace and track initial values
    if traceE
        track.E_k   = zeros( NT + 1, 1 );   track.E_k  ( track.n ) = E_k_n;
        track.E_p   = zeros( NT + 1, 1 );   track.E_p  ( track.n ) = E_p_n;
        track.E_tot = zeros( NT + 1, 1 );   track.E_tot( track.n ) = E_tot_n;
    end

%=========================================================================================
%% Initialization of the Animation, if enabled
% Initialize the main figure that will be used for the animation, together with graphics 
% objects representing particles and their trajectories.

    % Get the size of the screen
    screen.Rect = get( 0, 'ScreenSize' ); % gets [ left, bottom, width, height ]
    screen.L = screen.Rect(1);  screen.B = screen.Rect(2);
    screen.W = screen.Rect(3);  screen.H = screen.Rect(4);

    if params.Animate || params.ManualNp

        % Calculate paper dimensions for the figure (height and width) which
        % depends whether the picture is flat or 3-D
        switch N_dim
            case 1,  screen.box_W = 2 * params.LTick;
                     screen.box_H = L(1);
            case 2,  screen.box_W = L(1);
                     screen.box_H = L(2);
            case 3,  screen.box_W = ( L(1) + L(2) ) * 0.86;
                     screen.box_H = L(3) + max( L(1), L(2) ) * 0.5;
        end

        % Normalize paper height/width to max. allowed figure width (in centimeters)
        screen.box_maxWH = max( screen.box_W, screen.box_H );
        screen.pap_W = 3 + ( params.MaxFigSize - 3 ) * screen.box_W / screen.box_maxWH;
        screen.pap_H = 3 + ( params.MaxFigSize - 3 ) * screen.box_H / screen.box_maxWH;
        screen.pap_maxWH = max( screen.pap_W, screen.pap_H );

        % Printer to screen conversion factor (pixels per centimeter), so the whole 
        % printer figure fits in 1/2 of the screen
        screen.dpcm = min( screen.H, screen.W ) / screen.pap_maxWH / 2;

        % Configure screen size and position (position is centered with slight offset)
        screen.fig_W = screen.pap_W * screen.dpcm; % Figure width on screen
        screen.fig_H = screen.pap_H * screen.dpcm; % Figure height on screen
        screen.fig_L = ( screen.W - screen.fig_W ) * 0.3; % Figure left position
        screen.fig_B = ( screen.H - screen.fig_H ) * 0.6; % Figure bottom position

        % Finally, open the main figure
        mainFig = figure( ...
            'Name', [ 'Colliding Particles (', params.ID, ')' ], ... % window title
            'FileName', [ 'lab2_fig1_', params.ID, '.fig' ], ... % name of the FIG-file 
            'DoubleBuffer', 'on', ... % flicker-free rendering
            'Position', [ screen.fig_L, screen.fig_B, screen.fig_W, screen.fig_H ], ...
            'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
            'PaperSize', [ screen.pap_W, screen.pap_H ], ...
            'PaperPosition', [ 0, 0, screen.pap_W, screen.pap_H ] ...
        );

        % Add the 'stop simulation' menu item, if not requested to be interactive
        if ~params.ManualNp
            miCancel = uimenu( mainFig, 'Label', '&Stop', ...
                'Callback', 'global lab2_running; lab2_running = false;', ...
                'ForegroundColor', [ 0.7, 0, 0 ] ...
            );
        end

        % Add context menu items to export figure as EPS/EMF/PNG
        cxMenu = uicontextmenu;
        uimenu( cxMenu, 'Label', 'Export as EPS', 'Callback', [ ...
            'lab2_main( ''@exportFigure'', gcbf, ''-depsc2'', ', ...
            '''lab2_fig1_', params.ID, ''' );' ...
        ]);
        uimenu( cxMenu, 'Label', 'Export as PNG', 'Callback', [ ...
            'lab2_main( ''@exportFigure'', gcbf, ''-dpng'', ', ...
            '''lab2_fig1_', params.ID, ''' );' ...
        ]);
        uimenu( cxMenu, 'Label', 'Export as EMF', 'Callback', [ ...
            'lab2_main( ''@exportFigure'', gcbf, ''-dmeta'', ', ...
            '''lab2_fig1_', params.ID, ''' );' ...
        ]);
        % Add context menu item to export annotation as EPS
        uimenu( cxMenu, 'Label', 'Export Annotation', 'Separator', 'on', ...
            'Callback', [ ...
                'lab2_main( ''@exportAnnotation'', gcbf, ', ...
                '''figCaption'', ''-depsc2'', ''', ...
                [ 'lab2_info_', params.ID ], ''' );' ...
            ] ...
        );
        set( gcf, 'UIContextMenu', cxMenu );
        set( gca, 'UIContextMenu', cxMenu );
        clearvars cxMenu

        % Axes, common part:
        axis equal;           % Sets equal scaling on all axes
        tck = params.LTick;   % Distance between ticks
        gridColor = [ 0.5, 0.8, 0.8 ];  % dark cian

        switch N_dim   % Axes, N_dim dependent part:
          case 1
            % Axes and labels
            axis( [ -tck, tck, -0.2, L(1) ] );
            ylabel( '{\itz} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
            % Ticks
            set( gca, 'XTick', -tck:tck:tck, 'YTick', 0:tck:L(1) );
            % Grid
            grid on;
          case 2
            % Axes and labels
            axis( [ 0, L(1), -0.2, L(2) ] );
            xlabel( '{\itx} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
            ylabel( '{\itz} /{\rmm}', 'FontSize', 12, 'FontName', 'Times' );
            % Ticks
            set( gca, 'XTick', 0:tck:L(1), 'YTick', 0:tck:L(2) );
            % Grid
            grid on;
            % For grid in different color than ticks do instead:
            % [a,b] = meshgrid( [ 0:tck:L(1), NaN ], [ 0:tck:L(2), NaN ] );
            % line( ...
            %     [ a(:); reshape( a', [], 1 ) ], [ b(:); reshape( b', [], 1 ) ], ...
            %     'LineStyle', ':', 'Color', gridColor ...
            % );
            % clearvars a b
          case 3
            % Axes and labels
            axis( [ 0, L(1), 0, L(2), 0, L(3) ] );
            set( gca, 'XTickLabel', '', 'YTickLabel', '', 'ZTickLabel', '' );
            % Ticks
            set( gca, 'XTick', 0:tck:L(1), 'YTick', 0:tck:L(2), 'ZTick', 0:tck:L(3) );
            set( gca, 'TickDir', 'out' );  % Ticks outside
            % Grid
            set( gca, 'XColor', gridColor, 'YColor', gridColor, 'ZColor', gridColor );
            grid on;
            box on;
            % Projection and view
            set( gca, 'Projection', 'perspective' );
            view( 146, 15 );
        end

        hold on; % Retain subsequent plots in figure

        % Draw the ground level either as horizontal thick black line or half-transparent
        % yellowish plane crossing point [0,0,0]

        switch N_dim
          case 1
            plot( [ -tck, tck ], [ 0, 0 ], '-k', 'LineWidth', 1.5 );
          case 2
            plot( [ -tck, L(1) + tck ], [ 0, 0 ], '-k', 'LineWidth', 1.5 );
          case 3
            fill( [ 0, L(1), L(1), 0, 0 ], [ 0, 0, L(2), L(2), 0 ], [ 1, 1, 0 ], ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none' );
        end

        clearvars tck gridColor   % end of temporary variable

        % Init paticle visual representation template: either a sphere or a circle

        if N_dim == 3   % Render particles as spheres always when plotting in 3D:
            params.PAsSpheres = true;
        end

        if params.PAsSpheres
            % Generates surface points for a unit sphere
            if N_dim == 3
                [ US_x, US_y, US_z ] = sphere( params.N_usf );
            else
                % Swapped Y/Z axies when plotting in 1- or 2-D
                [ US_x, US_z, US_y ] = sphere( params.N_usf );
            end
            % Also, turn-on OpenGL if ploting surfaces
            set( gcf, 'Renderer', 'OpenGL' );
        else
            % Generates edge points for a unit circle
            ang = linspace( 0, 2 * pi, params.N_uck );
            [ US_x, US_y, US_z ] = deal( cos(ang), sin(ang), zeros(size(ang)) );
            clearvars ang
        end

        % Use the light source when displaying particles as spheres

        if params.PAsSpheres
            light( 'Position', [-1 1 0.4], 'Style', 'Infinite' );
        end

        % Capture movie; prepare the new file and create an animation.

        if params.CaptureMovie
            movieObj = VideoWriter([ 'lab2_', params.ID, '.avi' ]);
            movieObj.FrameRate = params.fps;
            movieObj.Quality = 75;
            open( movieObj );
        end

        % Animation state variables

        t_nextFrame = 0;   % Schedule next animation redraw, s
        wiggleCounter = 0; % Used for wiggle stereoscopy

        % Setup screen update interval as the inverse of update frequency
        % Note that if fps is 0 than the period is Inf, which is ok.
        framePeriod = 1 / params.fps;

        %%% Initialize graphics object handles associated to particles -------------------

        OH = zeros( N_p, 1 ); % Allocate space for particle graphics object handles

        % Allocate space for trajectory plot handles
        if params.PlotPath
            trajPlot = zeros( N_p, 1 ); 
        end
        
        % Finaly, create graphical objects representing particles
        for i = 1:N_p
            initParticle_Visual( i );
        end

        % Display velocity vector and particle state variables, if not manual placing
        % mode and if AnnotInit is set or initial and final times are equal and number
        % of particles is less or equal 8).
        if ~params.ManualNp && ( params.AnnotInit || ( t_0 == t_f && params.N_p <= 8 ) )
            for i = 1:N_p
                annotateParticle( i, struct( 'm', M(i,:), 'r', R(i,:), ...
                    'x', X(i,:), 'v', V(i,:), 'color', faceColor(i,:) ) );
            end
        end

        % Display simulation status in the top region of the figure
        textInfo = text( 0.5, 1.01,  '', 'Units', 'normalized', ...
            'FontName', 'Helvetica', 'FontSize', 10, 'Color', [ 0, 0, 0.9 ], ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none', 'EdgeColor', 'none' ...
        );
    
        %%% Interactivelly add new particles ---------------------------------------------
        
        if params.ManualNp

            addNewParticles_User ();

            params.ManualNp = false;  % we're done; turn off the 'manual mode' flag
            params.N_p = N_p;         % update number of particles

            if N_p == 0 % If no particles were added
                delete( mainFig );   % close the window
            else
                % Save the particle state variables ...
                params.M_0 = M;  % initial masses
                params.R_0 = R;  % initial radii
                params.X_0 = X;  % initial positions
                params.V_0 = V;  % initial velocities
                params.M   = M;  % final masses
                params.R   = R;  % final radii
                params.X   = X;  % final positions
                params.V   = V;  % final velocities
                % ... and particle colors:
                params.faceColor_0 = faceColor;  % initial face colors
                params.faceColor   = faceColor;  % final face colors
            end

            % Save the result also as MAT-file, if requested
            if params.SaveData
                tFilename = [ 'Init_', params.ID, '.mat' ];
                save( tFilename, 'params' );
            end

            % Return results to the caller
            if nargout >= 1, varargout{1} = params; end
            if nargout >= 2, varargout{2} = track;  end
            lab2_running = false;
            return
        end
    end

%=========================================================================================
%% Sanity-check and snapshot of the initial state variables

    assert( isequal( size( M ), [ N_p, N_dim ] ) );
    assert( isequal( size( R ), [ N_p, 1 ] ) );
    assert( isequal( size( X ), [ N_p, N_dim ] ) );
    assert( isequal( size( V ), [ N_p, N_dim ] ) );
    assert( isequal( size( faceColor ), [ N_p, 4 ] ) );
    assert( all( M(:) > 0 ) );
    assert( all( R(:) >= 0 ) );

    if t_0 == 0          % Save state vars only if the initial time is a real initial time
        params.M_0 = M;  % save initial masses and ...
        params.R_0 = R;  % initial radii
        params.X_0 = X;  % initial positions
        params.V_0 = V;  % initial velocities
        params.faceColor_0 = faceColor;  % initial face colors
    end

%=========================================================================================
%% Main Simulation Loop

    tic; % Start measuring elapsed (real) time

    % Initialize next time to emit a new particle
    t_nextEmit = t_0 + T_emit - h;  

    % If not animating particles, open the progress bar
    if ~params.Animate && t_f > t_0
        t_waitbInc = ( t_f - t_0 ) / 20;   % Waitbar update each 5% of time
        t_waitbUpdate = t_0 + t_waitbInc;  % Schedule next update
        waitbh = waitbar( 0, sprintf( 't = %.3f s', t ), ...
            'Name', [ 'Running ', params.Title ], ...
            'WindowStyle', 'modal', 'Color', [ 0.925, 0.914, 0.847 ], ...
            'CreateCancelBtn', 'global lab2_running; lab2_running = false;' );
    else
        % There is no progress bar = no need for updates
        t_waitbUpdate = Inf;
    end

    % Determine number of cycles for the main simulation loop
    if stepper == enum.Leapfrog && imodel == enum.SpringModel
        % When using leapfrog integrator in the linear spring model, we are actually 
        % counting half-steps, i.e. our main loop will have 2*NT half-steps.
        % For how even and odd half-steps are handled, see 'Step Forward and Update 
        % State Variables' section bellow.
        N_steps = NT * 2;
    else
        N_steps = NT;
    end

    % Do the initial half-step for leapfrog integrator in the impulse model
    if stepper == enum.Leapfrog && imodel == enum.ImpulseModel
        V = V + A * h/2;
    end

for n = 1 : N_steps

    %-------------------------------------------------------------------------------------
    %% » Collision detection and collision response
    % Search for colliding geometries. Collision response impulse transfer, projections
    % that resolve overlapping geometries. Compute also *penalty* forces for interacting 
    % objects. Penalty forces are accumulated for interaction and later added to a total 
    % force for each object.

    %%% Initialize collision response variables

    F_tot = zeros( N_p, N_dim );     % Total force
    E_penalty = 0;                   % Start with zero total energy
    collisionCount = zeros( 1, 2 );  % Collision count

    if imodel == enum.ImpulseModel
        VJ_tot = zeros( N_p, N_dim );  % Total velocity deflection (jolt)
        XJ_tot = zeros( N_p, N_dim );  % Total position projection
    end

    %%% Find the penalty forces (in soft-core interactions) and/or the velocity
    %%% deflections (jolts) and position projections (in hard-core interactions)

    for i = 1:N_p    % for each particle (i)

        % *NOTE:* To avoid jittery behaviour in multiple-particle hard-core interaction, 
        % we should iterate this loop until velocity deflections and position projections
        % (jolts) settle down (which is slightly above the purpose of this lab).
        % Alternative is to use position projections factor k_p (in range between 0 and 1)
        % which softens overlappings over time (see usage of k_p bellow).

        %---------------------------------------------------------------------------------
        %% »» Displacement and distance from other particles (j) to current particle (i)

        dXp = bsxfun( @minus, X(i,:), X );     % dXp(i,j) = X(i) - X(j)
        norm_dXp = sqrt( sum( dXp .^ 2, 2 ) ); % norm_dXp(i,j) = |dXp(i,j)| (sum per row)

        %%% Penetration depth from other particles (overlap)

        dX = R(i) + R - norm_dXp;        % dX(i,j) = R(i) + R(j) - |X(i) - X(j)|

        %%% Find all overlaping particles (j), excluding self-overlap (where j == i)

        overlaps = find( dX > 0 );       % It is overlap with (j) if dX(j) > 0
        overlaps( overlaps == i ) = [];  % Exclude the self-overlap (where j == i)

        if isempty( overlaps )
            switch imodel
              case enum.SpringModel
                f_part = ZeroVec;    % Total penalty force for (i) from other particles
              case enum.ImpulseModel
                vj_part = ZeroVec;   % Total velocity jolt for (i)
                xj_part = ZeroVec;   % Total position projection for (i)
            end
        else
            % Update number of collisions between particles in this time-step.
            % Since the particle interacts always with some other particle, we should 
            % count only half of the interactions per one particle.
            if traceXP
                collisionCount(1) = collisionCount(1) + length( overlaps ) / 2;
            end

            % From now on, consider only overlapping particles...
            dX = dX( overlaps, : );
            dXp = dXp( overlaps, : );
            norm_dXp = norm_dXp( overlaps, : );

            % Unit direction from j to i
            % nX(i,j) = ( X(i) - X(j) ) / |X(i) - X(j)|
            nX = bsxfun( @rdivide, dXp, norm_dXp );

            % Relative velocities from other overlapping particles
            % dVp(i,j) = V(i) - V(j)
            dVp = bsxfun( @minus, V(i,:), V(overlaps,:) );

            % Get the scalar product dVp(i,j)' * nX(i,j)
            dVp_nX = sum( dVp .* nX, 2 );

            switch imodel
              case enum.SpringModel  % ---------------------------------------------------

                % Total penalty forces from overlapping particles (sum per column)
                f_part = sum( bsxfun( @times, k_s * dX - k_d * dVp_nX, nX ), 1 );
                % Sum-up the potential energy of the penalty forces per particle.
                % Since we are summing spring energy twice (for each side of the spring
                % i.e. for each particle) we should sum only 1/2 of the 1/2 *k_s * dX^2.
                E_penalty = E_penalty + 0.25 * k_s * sum( sum( dX .^ 2 ) );

              case enum.ImpulseModel  % --------------------------------------------------

                % Suppress impulse jolt for separating particles !!!
                dVp_nX( dVp_nX > 0 ) = 0;
                % Reduced mass ratio: m(i)^-1 / ( m(i)^-1 + m(j)^-1 )
                m_ratio = M(overlaps,1) ./ bsxfun( @plus, M(i,1), M(overlaps,1) );
                % Total impulse jolt from overlapping particles (sum per column)
                vj_part = sum( bsxfun( @times, - (1 + k_e) * m_ratio .* dVp_nX, nX ), 1 );
                % Position projection from overlapping particles (multiplied by a position
                % projectsion factor k_p where k_p = 0 disables projections, and k_p = 1
                % fully enables projections while value inbetween softens overlappings)
                if k_p ~= 0
                    xj_part = k_p * sum( bsxfun( @times, m_ratio .* dX, nX ), 1 );
                else
                    xj_part = ZeroVec;
                end

            end % switch imodel ----------------------------------------------------------
        end

        %---------------------------------------------------------------------------------
        %% »» Now, proceed with the displacement and distance from surface (s)
        % dXp(i,s) = | ( X(i) - R(s) )' * n(s) |
        dXp = sum( bsxfun( @minus, X(i,:), XSFC ) .* NSFC, 2 ); % sum up per row

        %%% Penetration depth (overlap) of particle to surface (s) 
        % dX(i,s) = R(i) - dXp(i,s)
        dX = R(i) - dXp;

        %%% Find all overlaping surfaces (s)
        overlaps = find( dX > 0 );

        if isempty( overlaps )
            switch imodel
              case enum.SpringModel
                f_sfcs = ZeroVec;    % Total penalty force for (i) from surfaces
              case enum.ImpulseModel
                vj_sfcs = ZeroVec;   % Total velocity jolt for (i) from surfaces
                xj_sfcs = ZeroVec;   % Total position projection for (i) from surfaces
            end
        else
            % Update number of collisions with surfaces in this time-step
            if traceXP
                collisionCount(2) = collisionCount(2) + length( overlaps );
            end

            % From now on, consider only overlapping surfaces...
            dX = dX( overlaps, : );

            % Unit direction from surface
            nX = NSFC( overlaps, : );

            % Relative velocity of particle (i) to surface (where surfaces are stationary)
            dVp = bsxfun( @minus, V(i,:), zeros(size(dX)) );

            % Get the scalar product dVp(i,j)' * nX(i,j)
            dVp_nX = sum( dVp .* nX, 2 );

            switch imodel
              case enum.SpringModel  % ---------------------------------------------------

                % Total penalty forces from overlapping surfaces (sum per column)
                f_sfcs = sum( bsxfun( @times, k_s * dX - k_d * dVp_nX, nX ), 1 );
                % Sum-up the potential energy of the penalty forces
                E_penalty = E_penalty + 0.5 * k_s * sum( sum( dX .^ 2 ) );

              case enum.ImpulseModel  % --------------------------------------------------

                % Suppress impulse jolt when separating from the surface !!!
                dVp_nX( dVp_nX > 0 ) = 0;
                % Total impulse jolt from overlapping surfaces (sum per column)
                vj_sfcs = sum( bsxfun( @times, -( 1 + k_e ) * dVp_nX, nX ), 1 );
                % Position projection from overlappings surfaces
                if k_p ~= 0
                    xj_sfcs = k_p * sum( bsxfun( @times, dX, nX ), 1 );
                else
                    xj_sfcs = ZeroVec;
                end

            end % switch imodel ----------------------------------------------------------
        end

        %---------------------------------------------------------------------------------
        %%% Sum up the total penalty force in linear spring model, or alternativelly,
        %%% the total velocity jolt and the position projection in impulse model

        switch imodel
          case enum.SpringModel  % -------------------------------------------------------

            % Add the sum of all penalty forces from other particles and surfaces
            F_tot(i,:) = F_tot(i,:) + f_part + f_sfcs;

          case enum.ImpulseModel  % ------------------------------------------------------

            % Add the sum of all impulse jolts from other particles and surfaces
            VJ_tot(i,:) = VJ_tot(i,:) + vj_part + vj_sfcs;
            % Add the sum of all position projections from other particles and surfaces
            XJ_tot(i,:) = XJ_tot(i,:) + xj_part + xj_sfcs;

        end % switch imodel
    end

    %-------------------------------------------------------------------------------------
    %% » Compute Forces and Constraints
    % Compute external forces on the objects and internal forces (except penalty forces
    % which were computed earlier) for interacting objects. Forces are accumulated to a
    % total force for each object. Initial F_tot is non-zero only if the penalty forces
    % were previously calculated.

    %%% Add gravitational force to the total force

    F_tot = F_tot + bsxfun( @times, M, g );

    %-------------------------------------------------------------------------------------
    %% » Step Forward and Update State Variables
    % Advances the simulation data from the current point of time to the next and
    % also computes the new state variables (solving a system of equations).

    if stepper == enum.Leapfrog && imodel == enum.ImpulseModel  % ------------------------

        A = F_tot ./ M;            % Solve the acceleration (first half-step)
        V = V + A * h/2 + VJ_tot;  % Solve the velocity (first half-step)
        X = X + V * h   + XJ_tot;  % Solve the position
        % A = F_tot ./ M;          % *) Solve the acceleration (second half-step)
        V = V + A * h/2;           % Solve the velocity (second half-step)
        t = t_0 + n * h;           % The time-step (t += h produces rounding errors)

        % *) Assuming that the total force does not depend on T, X or V (e.g. where
        % the only force is gravity), the total force at this point is not changed and 
        % the acceleration does not need to be solved again for the second half-step of 
        % the leapfrog. This means that this implementation of the leapfrog integrator 
        % *does not* work with the penalty collision model where F_tot depends on X.
        % The following code deals with this remedy:

    elseif stepper == enum.Leapfrog && imodel == enum.SpringModel  % ---------------------

        % In the leapfrog integrator for the linear-spring model, one cycle of the main 
        % loop corresponds to one half-step (h/2) of the integrator. This means that
        % the even cycle number 'n' correspond to a full-step, for which X is solved 
        % and derived quantities are calculated.
        % Note also that in the linear-spring model we do not use VJ_tot/XJ_tot.

        A = F_tot ./ M;            % Solve the acceleration (each half-step)
        V = V + A * h/2;           % Solve the velocity (each half-step)

        if mod( n, 2 ) == 1        % If it's the first (odd) half-step, don't solve X 
            continue;              % i.e. just do n++
        end
                                   % If it's the second (even) half-step, i.e. full-step:
        X = X + V * h;             % Solve the position
        t = t_0 + n * h/2;         % Solve the time = full time-step = 2 half-steps

    elseif stepper == enum.SemiImplicitEuler  % ------------------------------------------

        A = F_tot ./ M;            % Solve the acceleration
        V = V + A * h + VJ_tot;    % Solve the velocity
        X = X + V * h + XJ_tot;    % Solve the position
        t = t_0 + n * h;           % The time-step (t += h produces rounding errors)

    elseif stepper == enum.ForwardEuler  % -----------------------------------------------

        A = F_tot ./ M;            % Solve the acceleration
        X = X + V * h + VJ_tot;    % Solve the position
        V = V + A * h + XJ_tot;    % Solve the velocity
        t = t_0 + n * h;           % The time-step (t += h produces rounding errors)
        
    end % stepper ------------------------------------------------------------------------

    %-------------------------------------------------------------------------------------
    %% » Update Derived Quantities
    % After the state variables are updated, we can update all derived quantities.

    %%% Calculate the linear momentum
    
    P = M .* V;    % Note: The total linear momentum of the system is: sum( P, 1 )

    %%% Compute the kinetic, potential and total energy of the system

    E_k_n   = 0.5 * sum( sum( P .* V ) );  % Kinetic energy
    E_p_n   = E_penalty - sum( sum( X .* bsxfun( @times, M, g ) ) ); % Potential energy
    E_tot_n = E_k_n + E_p_n;  % Total energy

    %-------------------------------------------------------------------------------------
    %% » Simulation I/O (handle input and output of the simulation)
    % Input can be user interaction or data streaming from another simulation or hardware 
    % in the loop. Output can be data storage for post-processing, realtime graphics 
    % rendering, signals to haptic force-feedback unit.

    %%% Trace selected variables like position and momentum trajectory, or energies

    % Update the current track-record number, where track.n = 1 corresponds to t_0
    % i.e. track.n = n + 1 corresponds to t_n
    track.n = track.n + 1;

    % Save the position and the linear-momentum trajectory
    if traceXP
        track.X   ( track.n, :, : ) = X;
        track.P   ( track.n, :, : ) = P;
        track.colc( track.n, :    ) = collisionCount;
    end

    % Save the total energy trail
    if traceE
        track.E_k  ( track.n ) = E_k_n;
        track.E_p  ( track.n ) = E_p_n;
        track.E_tot( track.n ) = E_tot_n;
    end

    %%% If scheduled, update progress bar

    if t >= t_waitbUpdate
        t_waitbUpdate = t + t_waitbInc;  % Schedule next update
        waitbar( ( t - t_0 ) / ( t_f - t_0 ), waitbh, sprintf( 't = %.3f s', t ) );
    end

    %%% Real-Time Graphics Rendering

    if params.Animate

        % Determine whether the current frame should be rendered or not.
        % It should be rendered if CPU is idle or if it's scheduled at this time.
        
        t_pause = 0; % by default, don't try to pause if animating too fast
        if params.CaptureMovie
            t_now = t; % current time = simulated time
        else
            t_now = toc;  % current time = real-time
            if params.InRealTime
                % Calculate whether to pause if rendering is too fast:
                % time to sleep = simulation elapsed time - real elapsed time
                t_pause = ( t - t_0 ) - toc; 
            end
        end

        % Plot particles and trajectories: 
        % if we are idle, we are schedule to plot, or it's the last frame

        if t_pause > 0 || t_now >= t_nextFrame || t >= t_f

            % Set main figure as current (without figure pop-up)
            set( 0, 'CurrentFigure', mainFig );

            for i = 1:N_p
                % Update particle postion, either as 3D-surface or 2D-polygon
                switch N_dim
                    case 1,  c = [ 0, X(i), 0 ];
                    case 2,  c = [ X(i,:), 0 ];
                    case 3,  c = X(i,:);
                end
                set( OH(i), 'XData', R(i) * US_x + c(1), ...
                            'YData', R(i) * US_y + c(2), ...
                            'ZData', R(i) * US_z + c(3) );
                % Update particle trajectory, either as 3D or 2D-curve.
                % Note that it's ridicilous to plot trajectory in 1-D.
                if params.PlotPath
                    if N_dim == 2
                        set( trajPlot(i), ...
                            'XData', track.X( 1:track.n, i, 1 ), ...
                            'YData', track.X( 1:track.n, i, 2 ) );
                    elseif N_dim == 3
                        set( trajPlot(i), ...
                            'XData', track.X( 1:track.n, i, 1 ), ...
                            'YData', track.X( 1:track.n, i, 2 ), ...
                            'ZData', track.X( 1:track.n, i, 3 ) );
                    end
                end
            end

            % Update status info

            set( textInfo, 'String', sprintf( 'N_p = %d, t = %.3f s', N_p, t ) );

            %%% Wiggle stereoscopy, only if plotting in 3-D and if it's enabled

            if N_dim == 3 && params.Wiggle3D
                if t_pause > 0 
                    dtheta = 0.1; % smooth orbit camera, if CPU seems free
                else
                    dtheta = 0.5; % otherwise faster, but more shaky orbit
                end
                % Wiggle camera left/right
                wiggleCounter = mod( wiggleCounter + 1, 20 );
                if wiggleCounter >= 10
                    camorbit( dtheta, 0 ); % half-time orbit camera right
                else
                    camorbit( -dtheta, 0 ); % half-time orbit camera left
                end
            end
        end

        % Recalculate whether to pause or not

        t_pause = 0; % by default, don't try to pause if animating too fast
        if params.CaptureMovie
            t_now = t; % current time = simulated time
        else
            t_now = toc;  % current time = real-time
            if params.InRealTime
                % Calculate whether to pause if rendering is too fast:
                % time to sleep = simulation elapsed time - real elapsed time
                t_pause = ( t - t_0 ) - toc; 
            end
        end

        %%% Screen update

        if t_pause > 0
            % If too fast, pause so the simulation is in real-time
            pause( t_pause );
            t_nextFrame = toc + framePeriod;
        elseif t_now >= t_nextFrame
            % Update every 'framePeriod' secs in lengthy calculations
            t_nextFrame = t_now + framePeriod; 
            if params.CaptureMovie
                writeVideo( movieObj, getframe( gca ) );
            else
                drawnow;
            end
        else
            % Skips displaying frames if time-steps are lengthy to calculate
        end
    end

    %% Stop the simulation, if the global flag is reset

    if ~lab2_running
        fprintf( '\n>>> CANCELLED at n = %d, t = %g s, t_f = %g s, NT = %d\n\n', ...
            track.n - 1, t, t_f, NT );
        break;
    end

    %% Emit a new particle, if it's scheduled

    if t >= t_nextEmit
        t_nextEmit = t_nextEmit + T_emit - h; % Setup next time to emit a new particle
        emitNewParticle_HiZ( n + 1 );
    end

end % of the Main Simulation Loop --------------------------------------------------------

    t_elapsed = toc;  % get elapsed time

%=========================================================================================
%% Post-processing
% Processing of data stored during the simulation into quantities as required
% e.g., energy, temperature, velocity fields, time-averaged force. Data is stored in
% file or made into graphs, tables or animation. The system state variables may be saved 
% and used for initialization data for another simulation.

    %%% Close progress bar, if any

    if exist( 'waitbh', 'var' )
        delete( waitbh );
    end

    %%% Create common caption/annotation for the figures (in TeX format) 

    % Prepare info about interaction model
    switch imodel
      case params.enum.ImpulseModel
        imodel_info = [ params.verb.imodel{imodel}, ...
          ': {\ite} = ',   sprintf( '%g', k_e ),    ...
          ',  {\itk}_{\fontsize{8}p} = ', sprintf( '%g', k_p ) ];
      case params.enum.SpringModel
        imodel_info = [ params.verb.imodel{imodel},      ...
          ': {\itk}_{\fontsize{8}s} = ', ...
            sprintf( '%g', k_s ), ' N/m', ...
          ',  {\itk}_{\fontsize{8}d} = ', ...
            sprintf( '%g', k_d ), ' N m^{\fontsize{7}-1} s  ' ];
    end

    % Create caption for the symulation containing: title,
    % system info, collision model and integrator algorithm
    params.figCaption = { ...
        [ '{\fontname{Arial}',                     ... % Title
          '\fontsize{10}',                         ...
          '\bf', params.Title, '}  ' ];            ...
        '{\fontsize{6} }';                         ...
        [ '\fontname{Times New Roman}',            ... % System info
          '\fontsize{11}',                         ...
          num2str( N_dim ), '-D system',           ...
          ': {\itg} = ',   num2str( g_n ), ' m/s^{\fontsize{7}2}', ...
          ',  {\itN}_{\fontsize{8}p} = ', num2str( N_p ) ];  ...
        '{\fontsize{2} }';                         ...
        imodel_info;                               ... % Interaction model
        '{\fontsize{2} }';                         ...
        [ params.verb.stepper{params.stepper},     ... % Integrator algorithm
          ': {\ith} = ',   num2str( h   ), ' s',  ...
          ',  {\itt}_{\fontsize{8}f} = ', num2str( t   ), ' s  ' ] ...
    };

    %%% Update animation figure and optionaly save the figure as PNG file and close
    %%% captured movie file

    if params.Animate
        % Be more verbose since we don't need fast updates any more
        set( textInfo, 'String', sprintf( '{\\itN}_p = %d, {\\itt} = %.3f s', N_p, t ) );
        set( textInfo, 'Interpreter', 'tex', 'Color', 'k' );
        
        % Display velocity vector and particle state variables
        if params.AnnotInit && ~params.ManualNp
            for i = 1:N_p
                annotateParticle( i, struct( 'm', M(i,:), 'r', R(i,:), ...
                    'x', X(i,:), 'v', V(i,:), 'color', faceColor(i,:) ) );
            end
        end

        % Remove 'cancel simulation' menu item
        delete( miCancel );
        drawnow;  % Flush event queue and update figure window
        set( 0, 'CurrentFigure', mainFig );  % Go back to mainFig

        % Close the movie file
        if params.CaptureMovie
            close( movieObj );            
        end

        % Annotate simulation with the caption in the top-left corner of the axes
        axlim.x = xlim;
        axlim.y = ylim;
        caph = text( axlim.x(1) + params.LTick/10, axlim.y(2) - params.LTick/20, ...
            params.figCaption, 'Interpreter', 'TeX', 'Tag', 'figCaption', ...
            'FontName', 'Helvetica', 'FontSize', 10, 'Margin', 5, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
            'BackgroundColor', 'none', 'EdgeColor', 'none', 'LineStyle', '-' ...
        );
        set( caph, 'UIContextMenu', get( gcf, 'UIContextMenu' ) );
        if N_dim == 3
            axlim.z = zlim;
            set( caph, 'Position', [ axlim.x(2), params.LTick/5, axlim.z(2) ] );
            set( caph, 'Units', 'pixels' ); % convert to pixels i.e. fix position in 2-D
        end
        clearvars axlim caph

        % Convert all 'Helvetica' fonts to 'Arial' and 'Times' to 'Times New Roman'
        set( findall( findall( gcf, '-property', 'FontName' ), ...
            'FontName', 'Helvetica' ), 'FontName', 'Arial' );
        set( findall( findall( gcf, '-property', 'FontName' ), ...
            'FontName', 'Times' ), 'FontName', 'Times New Roman' );
    end

    %% Plot energy, linear-momentum, interactions and barycenter height

    % Determine total number of subplots

    if track.n <= 1
        subplot_count = 0;
    else
        subplot_count = params.PlotContacts + params.PlotEnergy ...
                      + params.PlotMomentum + params.PlotHeight;
    end

    if subplot_count
        
        % Setup window title for the next figure
        fig_name = '';
        if params.PlotContacts
            fig_name = [ fig_name, '/Collisions' ];
        end
        if params.PlotEnergy
            fig_name = [ fig_name, '/Energy' ];
        end
        if params.PlotMomentum
            fig_name = [ fig_name, '/Momentum' ];
        end
        if params.PlotHeight
            fig_name = [ fig_name, '/Height' ];
        end
        % Finally, add simulation ID to the title and get rid of leading '/'
        fig_name = [ fig_name, ' (', params.ID, ')' ]; 
        fig_name = fig_name(2:end);

        % Calculate position (*Warning*: Reuse screen structure of the main figure!)
        screen.fig_W = screen.H * 2/3; % Width  = 2/3 of the screen height (not width)
        screen.fig_H = screen.H * 2/3; % Height = 2/3 of the screen height
        screen.fig_L = ( screen.W - screen.fig_W ) * 0.3 + 30; % Figure left position
        screen.fig_B = ( screen.H - screen.fig_H ) * 0.6 - 30; % Figure bottom position

        % Reduce height depending of the number of the subplots
        max_plots = max( 3, subplot_count );
        plot_offset = screen.fig_H * ( max_plots - subplot_count ) / max_plots;
        screen.fig_H = screen.fig_H - plot_offset;
        screen.fig_B = screen.fig_B + plot_offset;
        clearvars plot_offset

        % Setup paper size
        screen.pap_W = params.PlotWidth;
        screen.pap_H = params.SubplotHeight * subplot_count;

        % Create figure
        energyFig = figure( 'Name', fig_name, ...
            'FileName', [ 'lab2_fig2_', params.ID, '.fig' ], ... % name of the FIG-file 
            'Position', [ screen.fig_L, screen.fig_B, screen.fig_W, screen.fig_H ], ...
            'PaperPositionMode', 'manual', 'PaperUnits', 'centimeters', ... 
            'PaperSize', [ screen.pap_W, screen.pap_H ], ...
            'PaperPosition', [ 0, 0, screen.pap_W, screen.pap_H ] ...
        );

        % Allocate space for individual plots axes handles
        ax = zeros( subplot_count, 1 );
        if subplot_count == 1
            ax(1) = gca;
        end
        cur_subplot = 0;  % Track the current subplot number

        % Initialize time as abscissa used for subsequent plots
        T = t_0 + h * ( 0 : track.n-1 );
    end

    %% » Plot collision statistics per time-step

    if params.PlotContacts && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end
        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Collision Statistics', 'FontSize', 11 );
        set( gca, 'YGrid', 'on' );
        ylabel( { 'Total # of contacts'; 'per time-step' }, ...
            'FontSize', 11, 'FontName', 'Times', 'Margin', 1 );

        % Calculate total number of colissions per time-step
        totColCount = sum( track.colc( 1:track.n, : ), 2 );

        % Plot statistics for total number of collisions
        i = find( totColCount > 0 );
        if ~isempty(i)
            stem( T(i), totColCount(i), 'DisplayName', 'part-part', ...
                'Color', [ 0.4 0.8 0.8 ], 'Marker', '.', ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none' );
        end

        % Plot statistics for collisions particles to surfaces
        % This plots overlays and splits previous plot into two parts: upper part with
        % particle-particle collision count and lower part with particle-surface collision
        % count.
        i = find( track.colc( 1:track.n, 2 ) > 0 );
        if ~isempty(i)
            stem( T(i), track.colc(i,2), 'DisplayName', 'part-surface', ...
                'Color', [ 0.8 0 0 ], 'Marker', '.', ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none' );
        end

        % Finally, we have plots so we can place legend
        set( legend( 'show' ), 'Location', 'SouthEast', 'FontSize', 9 );
    end
    
    %% » Plot the kinetic, potential and total energy of the system over time

    if params.PlotEnergy && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end
        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Energy', 'FontSize', 11 );
        grid on;
        ylabel( '{\itE} /\fontsize{11}{\rmJ}', ...
            'FontSize', 12, 'FontName', 'Times' );

        % Plots:
        plot( T, track.E_k  ( 1:track.n ), 'b-', 'LineWidth', 0.5 ); % blue, thin
        plot( T, track.E_p  ( 1:track.n ), 'r-', 'LineWidth', 0.5 ); % red, dot-sash, thin
        plot( T, track.E_tot( 1:track.n ), 'k-', 'LineWidth', 1.5 ); % black, thick

        % Legend:
        set( legend( '{\itE}_{\rmk}', '{\itE}_{\rmp}', '{\itE}_{\rmtot}' ), ...
            'FontSize', 10, 'FontName', 'Times', 'Location', 'SouthEast' );
    end

    %% » Plot the linear momentum of the system over time

    if params.PlotMomentum && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end
        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Linear Momentum', 'FontSize', 11 );
        grid on;
        ylabel( '{\itp} /\fontsize{11}({\rmkg \cdot m/s})', ...
            'FontSize', 12, 'FontName', 'Times' );

        % Sum-up total linear-momentum, per time-step/dimesion (as matrix NT+1 * N_dim)
        P_tot = permute( sum( track.P, 2 ), [ 1 3 2 ] );

        % Configure plot colors and legend
        switch N_dim
            case 1, ctype = { 'z'              };
                    ltype = { 'b-'             };
            case 2, ctype = { 'x',  'z'        };
                    ltype = { 'r-', 'b-'       };
            case 3, ctype = { 'x',  'y',  'z'  };
                    ltype = { 'r-', 'g-', 'b-' };
        end

        % Plot linear-momentum components:
        for i = 1:N_dim
            plot( T, P_tot( 1:track.n, i ), ltype{i}, 'DisplayName', sprintf( ...
                '\\fontname{Times New Roman}\\fontsize{11}{\\itp}_{\\it%s\\rm,tot}', ...
                ctype{i} ...
            ));
        end
        set( legend( 'show' ), 'Location', 'SouthEast' );
    end

    %% » Plot the height of the barycenter of the system over time

    if params.PlotHeight && subplot_count
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig

        % Create as a subplot, if more than one plot
        cur_subplot = cur_subplot + 1;
        if subplot_count > 1
            ax( cur_subplot ) = subplot( subplot_count, 1, cur_subplot );
        end

        hold on; % Retain subsequent plots in figure

        % Title, axes and grid:
        title( 'Barycenter Height', 'FontSize', 11 );
        grid on;
        ylabel( '{\itz}_{\rmcm} /\fontsize{11}{\rmm}', ...
            'FontSize', 12, 'FontName', 'Times' );

        % Find out the barycenter of the system over time as
        % X_cm = { Sum X(i) * M(i) } / { Sum M(i) }

        X_cm = zeros( track.n, N_dim );
        for i = 1 : track.n
            % Take particle positions at time-step i as X_i
            X_i = permute( track.X(i,:,:), [ 2 3 1 ] );
            % Find particles with NaN positions, i.e. non-existing particles at the time
            [NaN_particles,~] = find( isnan( X_i ) );
            % Existing particles ar not-NaN particles:
            ix = setdiff( 1:size(X_i,1), NaN_particles );
            % Calculate the center of mass for existing particles at the time:
            X_cm(i,:) = sum(  X_i(ix,:) .* M(ix,:), 1 ) ./ sum( M(ix,:), 1 );
        end
        clearvars X_i NaN_particles ix

        % Get the height of the barycenter and the upper envelope of the height
        height = X_cm(:,end);
        uenv = envelope( height );

        % Plots:
        plot( T, height, 'b-', 'LineWidth', 0.5 ); % height: blue, thick
        plot( T, uenv,   'r:', 'LineWidth', 0.5 ); % envelope: red, thin
        set( legend( 'height', 'envelope' ), 'Location', 'SouthEast', 'FontSize', 9 );
        
        clearvars X_cm height uenv
    end

    %% » Link axes and optionally save figure as EPS file

    if subplot_count

        % Label horizontal axis as time and link axes over time
        xlabel( '{\itt} /{\rms}', 'FontSize', 12, 'FontName', 'Times' );
        if subplot_count > 1
            linkaxes( flipdim( ax, 1 ), 'x' );
        end

        % Add context menu to export figure
        cxMenu = uicontextmenu;
        uimenu( cxMenu, 'Label', 'Export as EPS', 'Callback', [ ...
            'lab2_main( ''@exportFigure'', gcbf, ''-depsc2'', ', ...
            '''lab2_fig2_', params.ID, ''' );' ...
        ]);
        uimenu( cxMenu, 'Label', 'Export as PNG', 'Callback', [ ...
            'lab2_main( ''@exportFigure'', gcbf, ''-dpng'', ', ...
            '''lab2_fig2_', params.ID, ''' );' ...
        ]);
        uimenu( cxMenu, 'Label', 'Export as EMF', 'Callback', [ ...
            'lab2_main( ''@exportFigure'', gcbf, ''-dmeta'', ', ...
            '''lab2_fig2_', params.ID, ''' );' ...
        ]);
        set( gcf, 'UIContextMenu', cxMenu );
        set( ax, 'UIContextMenu', cxMenu );

        % Convert all 'Helvetica' fonts to 'Arial' and 'Times' to 'Times New Roman'
        set( findall( findall( gcf, '-property', 'FontName' ), ...
            'FontName', 'Helvetica' ), 'FontName', 'Arial' );
        set( findall( findall( gcf, '-property', 'FontName' ), ...
            'FontName', 'Times' ), 'FontName', 'Times New Roman' );

        clearvars cxMenu ax cur_subplot
        
        drawnow; % Flush event queue and update figure window
        set( 0, 'CurrentFigure', energyFig ); % Go back to energyFig
    end

    %% » Update 'params' structure with the final values of the state variables
    % The structure will be returned to the caller and can be used as an input for 
    % further simulation(s).

    params.N_p = N_p;  % Number of particles (may have changed if emitting new ones)
    params.M   = M;    % Final masses
    params.R   = R;    % Final radii
    params.X   = X;    % Final positions
    params.V   = V;    % Final velocities
    params.faceColor = faceColor;  % Save final particle face colors
    
    params.NT  = n;    % Number of iterated time-steps
    params.t_f = t;    % Simulation time
    params.ElapsedTime = t_elapsed;  % Elapsed real-time

    %% » Save the final values of the state variables and the traced variables in MAT-file

    if traceXP
        % Generate the time-step vector; note that 'n' is the final time-step
        % and remove non-recorded data from position and linear momentum trace
        track.T = ( t_0 + h * ( 0 : track.n-1 ) )';

        % Get rid of time-steps that are not simulated
        track.X    ( track.n+1:end, :, : ) = [];
        track.P    ( track.n+1:end, :, : ) = [];
        track.colc ( track.n+1:end, : )    = [];
    end
    
    if traceE
        % Get rid of time-steps that are not simulated
        track.E_k   ( track.n+1:end ) = [];
        track.E_p   ( track.n+1:end ) = [];
        track.E_tot ( track.n+1:end ) = [];
    end
        
    if params.SaveData
        tFilename = [ 'Trace_', params.ID, '.mat' ];
        save( tFilename, 'params', 'track' );
    end

    if ~nargout
        disp( 'Results:' ), disp( params );
        disp( 'Track:' ), disp( track );
    end

    % Return results to the caller and indicate that the simulation is done
    if nargout >= 1, varargout{1} = params; end
    if nargout >= 2, varargout{2} = track;  end

    lab2_running = false;
    return

%=========================================================================================
%% NESTED FUNCTIONS ----------------------------------------------------------------------

    %-------------------------------------------------------------------------------------
    %% diminishDimensions : Diminish dimensions of a matrix of spatial vectors
    
    function x = diminishDimensions( x, m, n )

        % Reduce number of rows (N_p), if needed
        if m < size(x,1)
            x( m+1:end, : ) = [];
        end

        % Reduce number of columns (N_dim), if needed. 
        if n < size(x,2)
            % Note that the last column is the height, which should be kept always
            switch n
                case 1, x( :, 1:end-1 ) = []; % Remove X and opt. Y from (X,[Y,],Z)
                case 2, x( :, 2 ) = []; % Remove Y from (X,Y,Z)
            end
        end
        
    end % function diminishDimensions

    %-------------------------------------------------------------------------------------
    %% envelope( sig, threshold ) : Find upper and lower envelopes of a given signal
    % Usage:
    %   upperenv = envelope( sig )
    %   upperenv = envelope( sig, threshold )
    %   [upperenv,lowerenv] = envelope( sig )
    %   [upperenv,lowerenv] = envelope( sig, threshold )
    %
    % Arguments:
    % * sig    : vector of input signal (either as 1xN or Nx1 matrix)
    % * threshold : ignore signal variations bellow this value; default 1e-4
    %
    % Returns:
    % * upperenv : upper envelope of the input signal
    % * lowerenv : lower envelope of the input signal

    function varargout = envelope( sig, threshold )

        if nargin < 2
            threshold = 1e-4;
        end

        % Find out the first derivative of the signal
        delta = diff( sig );

        % Flatten signal variations bellow the threshold
        delta( abs(delta) < threshold ) = 0;
        delta = sign( delta );

        % Remove all regions where signal remains constant
        for xi = 2 : length( delta )
            if delta( xi ) == 0
                delta( xi ) = delta( xi - 1 );
            end
        end

        % Find out the second derivative
        delta = diff( delta );

        % Determine local maximum and minimum points
        upper_ind = find( delta < 0 ) + 1;  % maximum if f''(x) < 0
        lower_ind = find( delta > 0 ) + 1;  % minimum if f''(x) > 0

        first = 1;
        last = length(sig);
        xi = first : last;
        
        if length(xi) <= 1
            if nargout >= 1
                varargout{1} = sig;
            end
            if nargout >= 2
                varargout{2} = sig;
            end
            return
        end

        if length( upper_ind ) < 2
            upper_ind = [ first, upper_ind(:), last ];
        end
        if length( lower_ind ) < 2
            lower_ind = [ first, lower_ind(:), last ];
        end

        if nargout >= 1
            varargout{1} = interp1( upper_ind, sig(upper_ind), xi, 'linear', 'extrap' );
        end
        if nargout >= 2
            varargout{2} = interp1( lower_ind, sig(lower_ind), xi, 'linear', 'extrap' );
        end

    end % function envelope

    %-------------------------------------------------------------------------------------
    %% initParticle_StateVars : Initialize state variables for a new particle
    % Initializes state, calculated and derived variables, and optionally a phase-space
    % trajectory.
    %
    % Arguments:
    % * i:   particle number
    % * nc:  current time-step number (n = 1 corresponds to t_0)
    % * pp:  structure with the particle state variables as fields:
    %        'm' as mass, 'r' as radius, 'x' as position and 'v' as velocity

    function initParticle_StateVars( i, nc, pp )

        %%% Initialize state variables
        M(i,:) = pp.m * ones( 1, N_dim );
        R(i,:) = pp.r;
        X(i,:) = pp.x;
        V(i,:) = pp.v;

        %%% Initialize calculated and derived variables
        F_tot  (i,:) = M(i,:) .* g;
        A      (i,:) = F_tot(i,:) ./ M(i,:);
        VJ_tot (i,:) = ZeroVec;
        XJ_tot (i,:) = ZeroVec;
        P      (i,:) = M(i,:) .* V(i,:);

        %%% Initialize trajectory in phase space
        if traceXP
            track.X(1:nc-1,i,:) = NaN;  track.X(nc,i,:) = X(i,:);
            track.P(1:nc-1,i,:) = 0;    track.P(nc,i,:) = P(i,:);
        end
    end

    %-------------------------------------------------------------------------------------
    %% initParticle_Visual : Initialize graphics object handles associated to particles
    % Creates graphic handles for displaying particle and its position trajectory.
    % Arguments:
    % * i:  particle number

    function initParticle_Visual( i )

        % Make the figure mainFig to current (i.e. target for graphics output),
        % but do not change its visibility or stacking with respect to other figures
        set( 0, 'CurrentFigure', mainFig );
            
        if params.PAsSpheres
            % The position of the particle in figure
            switch N_dim
                case 1,  c = [ 0, X(i), 0 ];
                case 2,  c = [ X(i,:), 0 ];
                case 3,  c = X(i,:);
            end
            % Create the graphic object, a 3D surface, representing the particle
            % (with edge color only in 3D)
            OH(i) = surf( ...
                R(i) * US_x + c(1), R(i) * US_y + c(2), R(i) * US_z + c(3), ...
                'DisplayName', sprintf( 'Particle %d', i ), ...
                'EdgeLighting', 'gouraud', 'EdgeColor', 'none', ...
                'FaceLighting', 'gouraud', 'FaceColor', faceColor(i,1:3) ...
            );
            if N_dim == 3
                set( OH(i), 'EdgeColor', faceColor(i,1:3) );
            end
        else
            % The position of the particle in figure
            switch N_dim
                case 1,  c = [ 0, X(i) ];
                case 2,  c = X(i,:);
                case 3,  c = [ X(i,1), X(i,3) ];
            end
            % Create the graphic object, a 2D fill, representing the particle
            OH(i) = fill( R(i) * US_x + c(1), R(i) * US_y + c(2), ...
                '', 'DisplayName', sprintf( 'Particle %d', i ), ...
                'EdgeColor', faceColor(i,1:3) .* 0.8, ...
                'FaceColor',  faceColor(i,1:3), 'FaceAlpha', faceColor(i,4) ...
            );
        end

        % Create an empty trajectory plot object for particle (i)
        % with the same but slightly darker color as the particle
        if params.PlotPath
            trajPlot(i) = plot( NaN, NaN, '-', 'Color', faceColor(i,1:3) .* 0.8 );
        end

    end % function initParticle_Visual

    %-------------------------------------------------------------------------------------
    %% randomFaceColor : Generate a random color to be used as particle face color
    % Returns vector with 4 components: first 3 components are RGB color components 
    % and the 4-th is alpha (transparency)

    function color = randomFaceColor( hue )
        
        if params.PAsSpheres
            % random hue with max value, but reduced saturation
            color = [ hue, 0.4 + 0.3 * rand, 1 ];
            alpha = 1;
        else
            % random hue with reduced value and saturation (pastell colors)
            color = [ hue, 0.5 + 0.3 * rand, 0.7 + 0.3 * rand ]; 
            alpha = 0.8;
        end

        % Return RGB + alpha
        color =  [ hsv2rgb( color ), alpha ];
    end

    %-------------------------------------------------------------------------------------
    %% emitNewParticle_HiZ : Emit new particle
    % Emits a new particle placing it on the top of the particle with the highest
    % position.
    % Arguments:
    % * nc: current time-step number
    
    function emitNewParticle_HiZ( nc )
       
        % Find particle with maximum height (the top-most particle)
        [~,zi] = max( X(:,end) );

        % Calculate initial state variables
        if isempty(zi)
            % If there is no the top-mostparticle, generate random position for
            % a new particle inside the box (with some margin)
            x_0 = 1.1 * params.r_p + rand( 1, N_dim ) .* ( L - 2.2 * params.r_p );
            v_0 = ZeroVec;
        else
            % Copy the state variables from the top-most particle
            x_0 = X(zi,:);
            v_0 = V(zi,:);
            % Place a new particle above the top-most particle, having the horizontal
            % initial velocity the same the top-most particle but slightly disturbed. 
            % Note that the last index 'end' always denotes the vertical coordinate.
            x_0(end) = x_0(end) + 2 * R(zi) + params.dx_emit;
            v_0(1:end-1) = v_0(1:end-1) + params.dv_emit * ( -1 + 2 * rand(1,N_dim-1) );
        end

        % Now, add particle to the system...
        N_p = N_p + 1;
        
        initParticle_StateVars( N_p, nc, ...
            struct( 'm', params.m_p, 'r', params.r_p, 'x', x_0, 'v', v_0 ) );
        
        faceColor(N_p,:) = randomFaceColor( rand ); % init new face color

        if params.Animate
            initParticle_Visual( N_p );
        end

    end % function emitNewParticle

    %-------------------------------------------------------------------------------------
    %% addNewParticles_User : Emit new particles interactivelly
    % Allow user to specify mass, radius, position and velocity of particles
    % interactivelly.
    % Arguments:
    % * nc: current time-step number (optional; default 1)
    
    function addNewParticles_User( nc )

        if ~nargin
            nc = 1;
        end

        % Figure header
        header = annotation( 'textbox', [ 0, 0, 1, 0.995 ], ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Margin', 5, ...
            'TextColor', [ 0.7, 0, 0 ], 'BackgroundColor', [ 0.925, 0.914, 0.847 ], ...
            'LineStyle', '-', 'EdgeColor', [ 1, 1, 1 ], 'FitBoxToText', 'on' );

        hue = 0;

        while true

            ip = N_p + 1;

            set( header, 'String', [ 'Click ''{\bfOK}'' to place the particle ', ...
                '{\color{gray}... or ...} Click ''{\bfCancel}'' to finish and exit' ] ...
            );

            % Ask user for particle mass and radius
            answer = inputdlg( ...
                { [ '\fontsize{10}Enter particle mass: {{\itm}_', ... % question 1
                    num2str(ip), ' / kg =}' ],...
                  [ '\fontsize{10}Enter particle radius: {{\itr}_', ... % question 2
                    num2str(ip), ' / m =}' ] ...
                  }, ...
                [ 'Add particle ', num2str(ip) ], ... % title
                [ 1, 40; 1, 40 ], ... % size 2 x 40 characters
                { num2str( params.m_p ), num2str( params.r_p ) }, ... % def. values
                struct( 'Resize', 'off', 'WindowsStyle', 'modal', ... % properties
                        'Interpreter', 'TeX' ) ...
            );

            if isempty( answer )  % Quit loop if user selects 'Cancel'
                break;
            end

            % Parse mass and radius
            pp.m = str2double( answer{1} ); % parse mass
            pp.r = str2double( answer{2} ); % parse radius

            % Remember particle mass and radius as the default
            params.m_p = pp.m;
            params.r_p = pp.r;

            % ... check for errors (dont' report them, just ask for a new particle)
            if isnan(pp.m) || isinf(pp.m) || pp.m <= 0
                continue
            elseif isnan(pp.r) || isinf(pp.r) || pp.r <= 0 || pp.r >= min(L)/2
                continue
            end

            % Get the particle position
            set( header, 'String', ...
                'Select the particle position {\fontname{Times New Roman}{\bfx}}' );
            drawnow; figure( mainFig ); % Update figure window pop-up our figure
            pp.x = ginput( 1 );

            % Remove the x-component of the position when in 1-D
            if N_dim == 1, pp.x(1) = []; end

            % Add a new particle to the system
            N_p = N_p + 1;
            pp.v = ZeroVec; % with zero initial velocity
            pp.color = randomFaceColor( hue ); % and new face color
            faceColor(N_p,:) = pp.color;
            initParticle_StateVars( N_p, nc, pp );
            initParticle_Visual( N_p );

            % Get the the head of the velocity vector
            set( header, 'String', ...
                'Place the head of the velocity vector {\fontname{Times New Roman}{\bfv}}' );
            drawnow; figure( mainFig ); % Update figure window pop-up our figure
            pp.v = ginput( 1 );

            % Remove the x-component of the velocity when in 1-D
            if N_dim == 1, pp.v(1) = []; end

            % Set the velocity = difference between head and tail of the vector
            % and reinitialize the particle state variables since we changed velocity
            pp.v = pp.v - pp.x;
            initParticle_StateVars( N_p, nc, pp );
            
            % Display velocity vector and particle state variables
            annotateParticle( N_p, pp );

            % Get next hue from the color wheel (leaping forward among base colors)
            hue = hue + 0.31;
            if ( hue > 1 )
                hue = hue - 1;
            end
        end

        delete( header );
        
    end % function emitNewParticle_User

    %-------------------------------------------------------------------------------------
    %% annotateParticle( id, pp )
    % Displays velocity vector and state variables of the particle
    % Arguments:
    % * id: particle#
    % * pp:  structure with the particle state variables as fields:
    %        'm' as mass, 'r' as radius, 'x' as position, 'v' as velocity and
    %        'color' as face color
    
    function annotateParticle( id, pp )
        
        % Calculate tail and head points for the particle velocity
        switch N_dim
            case 1
                vecX = [ 0, 0 ];
                vecY = [ pp.x(1), pp.x(1) + pp.v(1) ];
            case 2
                vecX = [ pp.x(1), pp.x(1) + pp.v(1) ];
                vecY = [ pp.x(2), pp.x(2) + pp.v(2) ];
            case 3
                return % Annotation is not supported in 3-D
        end

        % Plot the velocity vector as line 'o----.'
        color = pp.color(1:3);
        plot( vecX, vecY, ...
            '-', 'Color', color * 0.8 );
        plot( vecX(2), vecY(2), '.', ...
            'Color', color * 0.8, 'MarkerFaceColor', color );

        % Annotate vector with the the particle state variables
        idstr = num2str( id );
        text( vecX(1) - pp.r/2, vecY(1) - pp.r, {
                [ '{\itm}_', idstr, ' =', sprintf( ' %g', pp.m(1) ), ' {\rmkg}'  ]; ...
                [ '{\itr}_', idstr, ' =', sprintf( ' %g', pp.r(1) ), ' {\rmm}'   ]; ...
                [ '{\bfx}_', idstr, ' = [', sprintf( ' %.2g', pp.x ), ' ] {\rmm}'   ]; ...
                [ '{\bfv}_', idstr, ' = [', sprintf( ' %.2g', pp.v ), ' ] {\rmm/s}' ]  ...
            }, ...
            'Interpreter', 'tex', 'FontName', 'Times', 'FontSize', 9, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left' ...
        );

    end % function annotateParticle

    %-------------------------------------------------------------------------------------
    %% exportFigure( figh, dtype, filename ) : Export figure to file
    % Saves figure to a file like 'print' function, but, in case of vector output
    % formats, makes all graphic objects opaque during printing to avoid that output 
    % vector files are rendered with Z-buffer or OpenGL, in which case, vector files 
    % will contain uncompressered rasterized graphics.

    function exportFigure( figh, dtype, filename ) %#ok<DEFNU>

        switch lower( dtype )
            case '-dmeta',  isVector = true;
            case '-deps',   isVector = true;
            case '-depsc',  isVector = true;
            case '-deps2',  isVector = true;
            case '-depsc2', isVector = true;
            case '-dill',   isVector = true;
            case '-dpdf',   isVector = true;
            case '-dsvg',   isVector = true;
            otherwise,      isVector = false;
        end

        if ~isVector
            % Export rasterized output in 300 dpi
            print( figh, '-r300', dtype, filename );
        else
            % When exporting in vector format, make all graphic objects opaque.
            % Find all transparent objects ...
            objh = findobj( figh, '-property', 'FaceAlpha' );
            objh = findobj( objh, '-not', 'FaceAlpha', 1 );
            % ... and remember their alpha channel
            objv = get( objh, 'FaceAlpha' );

            % Make objects opaque
            set( objh, 'FaceAlpha', 1 );

            % Render the output file with the given format
            print( figh, '-r300', dtype, filename );

            % Restore objects' transparencies
            if length(objh) == 1
                set( objh, 'FaceAlpha', objv );
            else
                set( objh, {'FaceAlpha'}, objv );
            end
        end

        % Fix PS fonts for EPS files
        if length( dtype ) >= 5 && strcmp( dtype(1:5), '-deps' )
            fixPsFonts( filename );
        end

    end % function exportFigure

    %-------------------------------------------------------------------------------------
    %% fixPsFonts( filename ) : Fix Postscript fonts
    % Fixes MATLAB's choice of base PS fonts, after using print command to export
    % figures to PS files; replaces 'Helvetica' with 'Arial' and 'Times-Roman' 
    % with 'Times New Roman'.
    % See http://www.mathworks.se/help/techdoc/creating_plots/f3-103191.html#f3-96850)
    % for font support in PS in MATLAB.

    function fixPsFonts( filename )

        if length( filename ) < 4 || ~strcmpi( filename(end-3:end), '.eps' )
            filename = [ filename, '.eps' ];
        end

        % Read in the EPS file
        fid = fopen( filename );
        ff = fread( fid, '*char' )';   % ff = char(fread(fid))';
        fclose( fid );   

        % Replace MATLAB fonts used in EPS to Windows
        % See allowed fontnames in %windir%\fonts\AdobeFnt.lst
        %
        conv_list = { ...
            '/Helvetica-BoldOblique',  '/Arial-BoldItalicMT'; ...
            '/Helvetica-Bold',         '/Arial-BoldMT'; ...
            '/Helvetica-Oblique',      '/Arial-ItalicMT'; ...
            '/Helvetica',              '/ArialMT'; ...
            '/Times-BoldItalic',       '/TimesNewRomanPS-BoldItalicMT'; ...
            '/Times-Bold',             '/TimesNewRomanPS-BoldMT'; ...
            '/Times-Italic',           '/TimesNewRomanPS-ItalicMT'; ...
            '/Times-Roman',            '/TimesNewRomanPSMT'; ...
        };
   
        for k = 1 : size( conv_list, 1 )
            ff = strrep( ff, conv_list{k,1}, conv_list{k,2} );
        end

        % Rerite the file with new contents
        fid = fopen( filename, 'w' );
        fprintf( fid, '%s', ff );
        fclose( fid );

    end % function fixPsFonts

    %-------------------------------------------------------------------------------------
    %% exportAnnotation : Copy textbox annotation to clipboard
    % Saves tagged annotation of the current figure to EPS file and deletes object.
    % Arguments:
    % * figh: figure holding annotation
    % * tag: tag of the annotation text
    % * dtype: output format type (for print function)
    % * filename: filename to save annotation

    function exportAnnotation( figh, tag, dtype, filename ) %#ok<DEFNU>

        % Parameter: margins added around the annotation text box (in pixels)
        sz.margin = [ 5, 5; 20, 5 ]; % (1,:) = left, bottom; (2,:) = right, top

        % Get annotaiton object by tag
        info.origh = findall( figh, 'Tag', tag );
        if isempty( info.origh )
            return;
        end

        % Create invisible figure with white background with invisible axes
        invh = figure( 'Visible', 'off', 'Color', 'w', ...
            'Position', [ 100, 100, 100, 100 ] );

        try
            axis( [ 0, 1, 0, 1 ] );
            set( gca, 'Visible', 'off', ...
                'ActivePositionProperty', 'position', ...
                'Units', 'normalized', 'Position', [ 0, 0, 1, 1 ] );

            % Create text holding annotation with the same text and interpeter as
            % the tagged object
            info.texth = text( 0, 0, get( info.origh, 'String' ), ...
                'FontName', 'Arial', 'FontSize', 10, ...
                'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left', ...
                'Interpreter', get( info.origh, 'Interpreter' ), 'Margin', 1 );

            % Get extent of the text and size of the figure
            sz.extent = get( info.texth, 'Extent' );
            sz.fig = get( invh, 'Position' );

            % Convert extent to hold width/height of the text in pixel units
            sz.extent = sz.fig(3:4) .* sz.extent(3:4);

            % Set figure size = text size + margins
            sz.fig(3:4) = sz.extent + sum( sz.margin, 1 ); % add margines 

            % Convert margin units to data space units
            sz.margin = bsxfun( @rdivide, sz.margin, sz.extent );

            % Reconfigure position of the figure both on screen and paper accordingly
            set( invh, 'Position', sz.fig, ...
                'PaperPositionMode', 'manual', 'PaperUnits', 'points', ... 
                'PaperSize', sz.fig(3:4), 'PaperPosition', [ 0, 0, sz.fig(3:4) ] ...
            );

            % Move the annotation text also keeping requested margins
            set( info.texth, 'Units', 'data', ...
                'Position', [ sz.margin(1,1), 1 - sz.margin(2,2) ] ...
            );

            % Render the figure to an output file
            print( invh, '-r300', dtype, filename );
            
            % Fix PS fonts for EPS files
            if length( dtype ) >= 5 && strcmp( dtype(1:5), '-deps' )
                fixPsFonts( filename );
            end

            % Finally, delete original object
            delete( info.origh );

        catch err
            disp( getReport( err ) );
        end

        close( invh );
    end

%=========================================================================================
end % function lab2_main

%% Lab2: Colliding Particles; Simulation Parameters Constructor
% Constructs the default 'params' structure which contains properties that define the 
% mathematical model and simulation constants, e.g., gravity, friction, mass and 
% geometry of the objects, visualization flags etc.
%
% Filename: lab2_defs.m
% Date:     2012-03-11
% Author:   Mikica B Kocic 
%
% Usage:
%    params = lab2_defs;
%    params = lab2_defs( 'configId' );
%    [~,configIds] = lab2_defs;
%
% Arguments:
%    configIds: list of recognized configuration profiles

function varargout = lab2_defs( varargin )

    params = struct();

    %%% Enumeration constants
    
    % Interaction model enumeration
    enum.SpringModel  = 1;   % Linear spring model
    enum.ImpulseModel = 2;   % Collision impulse model

    % Integrator method enumeration
    enum.ForwardEuler      = 1;   % Forward Euler integrator
    enum.SemiImplicitEuler = 2;   % Semi-implicit Euler integrator
    enum.Leapfrog          = 3;   % Leapfrog integrator
    
    % Enumeration captions
    verb.imodel  = { 'Linear-spring model'; 'Impulse collision model' };
    verb.stepper = { 'Forward Euler'; 'Semi-implicit Euler'; 'Leapfrog' };
    
    %%% Simulation description
    
    params.Title = 'Lab2: Default';
    params.ID = '';

    %%% Physical constants

    params.N_dim = 2;     % Number of space dimesnsions
    params.g_n   = 9.81;  % Standard acceleration of Gravity (intensity), m/s^2

    %%% System of particle parameters

    params.N_p       = 1;     % Number of particles
    params.T_emit    = Inf;   % Emitter period (for new particles), s; may be Inf
    params.ManualNp  = false; % Add new particles interactivelly and quit (!)
    params.dx_emit   = 0.01;  % Emitter: position gap between particles, m
    params.dv_emit   = 0.01;  % Emitter: horizontal velocity disturbance, fraction
    params.m_p       = 1.0;   % Particle mass, kg
    params.r_p       = 0.3;   % Particle radius, m
    params.BoxHeight = -5;    % Box width, m  (negative means: without limiting surface)
    params.BoxWidth  =  5;    % Box width, m
    params.BoxDepth  =  5;    % Box width, m

    %%% Integration parameters

    params.stepper = enum.Leapfrog;  % Used stepper

    params.t_0   =  0.0;  % Initial time, s
    params.t_f   = 10.0;  % Final time, s
    params.t_inc =  5.0;  % Final time increment, when continuing simulation, s

    %%% Collision model physical parameters
    
    params.imodel = enum.ImpulseModel ; % Used interaction model

    params.h_soft = 0.001; % Time-step used in linear spring model, s
    params.k_s    = 4000;  % Spring coefficient, N/m
    params.k_d    = 20;    % Spring damping factor, N/(m/s)

    params.h_hard = 0.01;  % Time-step used in impulse model, s
    params.k_e    = 0.8;   % Impulse collision restitution (0 == max dissipation)
    params.k_p    = 0.3;   % Position projection fraction (1 == max projection)

    %%% Results control flags

    params.SaveData      = false ; % Save simulation data and traces in MAT file
    params.PlotContacts  = false ; % Plot number of collisions (contacts) per time-step
    params.PlotEnergy    = false ; % Plot kinetic, potential and total energy over time
    params.PlotHeight    = false ; % Plot envelope of the first particle
    params.PlotMomentum  = false ; % Plot linear momentum of the system

    %%% Visualization control flags

    params.Animate      = true  ; % Animate particles
    params.InRealTime   = true  ; % Try animation in real-time, if possible
    params.PlotPath     = false ; % Plot particle trajectories
    params.PAsSpheres   = false ; % Display particles as spheres; othewise circles
    params.Wiggle3D     = false ; % Use 'Wiggle Stereoscopy' to convay 3-D depth
    params.AnnotInit    = false ; % Annotate the initial particle state variables
    params.AnnotFinal   = false ; % Annotate the final particle state variables

    %%% Animation parameters

    params.fps   = 25;  % Animation frame rate, Hz
    params.N_usf = 24;  % Number of faces on a sphere
    params.N_uck = 48;  % Number of segments in a circle
    params.LTick = 1;   % Distance between ticks, m

    %%% Saving figures and capturing animation flags

    params.CaptureMovie  = false ; % Capture animation as movie
    params.MaxFigSize    = 17    ; % Max width or height of saved figure, in centimeters
    params.PlotWidth     = 17    ; % Width of the plot, in centimeters
    params.SubplotHeight = 5     ; % Height of the subplot (when stacked), in centimeters

    %%% Used enumeration contants
    params.enum = enum;
    params.verb = verb;

    if nargout >= 1
        varargout{1} = params;
    end

    %% Return configuartion menu as the second argument

    if nargout >= 2
        % Note: Use strtrim( sscanf( caption, ' (%[^)])' ) ) to extract config id 
        % from the caption
        
        configIds = { ...
            '<html><em>Select test profile, then click ''Start''...</em></html>'; ...
            '  (Demo) - Demo simulation (piling balls)'; ...
            '  (Default) - Default parameters'; ...
            '<html><strong>Energy conservation:</strong></html>'; ...
            '  (1.2.a) - Energy plot, Impulse model, e = 1, Semi-implicit Euler'; ...
            '  (1.2.b) - Energy plot, Impulse model, e = 1, Leapfrog'; ...
            '  (1.2.c) - Energy plot, Spring model, k_d = 0, Semi-implicit Euler'; ...
            '  (1.2.d) - Energy plot, Spring model, k_d = 0, Leapfrog'; ...
            '<html><strong>Linear momentum conservation:</strong></html>'; ...
            '  (1.3.a) - Billiard balls, Impulse model, e = 1, Leapfrog'; ...
            '  (1.3.b) - Billiard balls, Impulse model, e < 1, Leapfrog'; ...
            '  (2.1.a) - Billiard balls, Spring model, k_d = 0, Semi-implicit Euler'; ...
            '  (2.1.b) - Billiard balls, Spring model, k_d > 0, Semi-implicit Euler'; ...
            '<html><strong>Stability of the method:</strong></html>'; ...
            '  (2.2.a) - Piling particles, Impulse model, e < 1, Leapfrog, h = 0.01'; ...
            '  (2.2.b) - Piling particles, Impulse model, e < 1, Leapfrog, h = 0.005'; ...
            '  (2.3.a) - Piling particles, Spring model, k_s = 4000, k_d > 0, Semi-implicit Euler'; ...
            '  (2.3.b) - Piling particles, Spring model, k_s = 8000, k_d > 0, Semi-implicit Euler'; ...
            '  (2.3.c) - Piling particles, Spring model, k_s = 4000, k_d = 0, Semi-implicit Euler'; ...
            '  (2.3.d) - Piling particles, Spring model, k_s = 8000, k_d = 0, Semi-implicit Euler'; ...
            '  (2.3.e) - Piling particles, Spring model, k_s = 40000, k_d > 0, Semi-implicit Euler'; ...
            '<html><strong>2D and 3D models:</strong></html>'; ...
            '  (2.4.a) - Piling 2D, Impulse model, Leapfrog'; ...
            '  (2.4.b) - Piling 2D, Spring model, Semi-implicit Euler'; ...
            '  (2.5.a) - Piling 3D, Impulse model, Leapfrog'; ...
            '  (2.5.b) - Piling 3D, Spring model, Semi-implicit Euler'; ...
        };

        % Just in case; remove html tags from captions if not using java.swing
        if ~usejava( 'swing' )
            for i = 1 : length( configIds )
                configIds{i} = regexprep( configIds{i}, '<(.|\s)*?>', '' );
            end
        end

        varargout{2} = configIds;
    end

    %% Quit if configuration profile is not specified

    if ~nargin
        return
    end

    %% Customize params structure depending on the selected configuration profile

    configId = varargin{1};
    params.Title = [ 'Lab2: ', configId ];

    switch configId
        case 'Default'
            % use default values as provided by the lab2_main

        case 'Demo' % 
            % Particles
            params.N_p        = 5;
            params.T_emit     = 1;
            % Stepper
            params.t_f        = 10;
            params.t_inc      = 10;
            % Flags
            params.PAsSpheres = true;

        case '1.2.a' % Energy plot, Immpel model, e = 1, Semi-implicit Euler
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -5;
            % Particles
            params.N_p        = 1;
            params.X          = 2.3; 
            params.V          = 0;
            % Model
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 1;
            params.k_p        = 0.01;
            % Stepper
            params.t_f        = 50;
            params.t_inc      = 10;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = false;
            params.PlotEnergy = true;
            params.PlotHeight = true;

        case '1.2.b' % Energy plot, Impulse model, e = 1, Leapfrog
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -5;
            % Particles
            params.N_p        = 1;
            params.X          = 2.3; 
            params.V          = 0;
            % Model
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 1;
            params.k_p        = 0;
            % Stepper
            params.t_f        = 50;
            params.t_inc      = 10;
            params.stepper    = enum.Leapfrog;
            % Flags
            params.Animate    = false;
            params.PlotEnergy = true;
            params.PlotHeight = true;

        case '1.2.c' % Energy plot, Spring model, k_d = 0, Semi-implicit Euler
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -5;
            % Particles
            params.N_p        = 1;
            params.X          = 2.3; 
            params.V          = 0;
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 4000;
            params.k_d        = 0;
            % Stepper
            params.t_f        = 40;
            params.t_inc      = 10;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = false;
            params.PlotEnergy = true;
            params.PlotHeight = true;

        case '1.2.d' % Energy plot, Spring model, k_d = 0, Leapfrog
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -5;
            % Particles
            params.N_p        = 1;
            params.X          = 2.3; 
            params.V          = 0;
            % Flags
            params.imodel     = enum.SpringModel;
            params.k_s        = 4000;
            params.k_d        = 0;
            % Flags
            params.t_f        = 40;
            params.t_inc      = 10;
            params.stepper    = enum.Leapfrog;
            % Flags
            params.Animate    = false;
            params.PlotEnergy = true;
            params.PlotHeight = true;

        case '1.3.a' % Billiard balls, Impulse model, e = 1, Leapfrog
            % System dimensions
            params.N_dim      = 2;
            params.BoxHeight  = 2;
            params.BoxWidth   = 7;
            params.g_n        = 0;
            % Particles
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 1;
            params.k_p        = 0.1;
            % Particles
            params.t_f        = 5;
            params.t_inc      = 5;
            params.stepper    = enum.Leapfrog;
            % Flags
            params.Animate    = true;
            params.PlotContacts = true;
            params.PlotEnergy   = true;
            params.PlotMomentum = true;
            % Particles
            % 5 billiard balls in 2-D with different masses and different velocities:
            params.N_p   = 5;
            params.M = [   5,   5;  2,    2;   3,   3;   4,   4;   5,   5 ];
            params.R = [      0.3;      0.3;      0.3;      0.3;      0.3 ];
            params.X = [ 0.5, 0.3;  2,  0.3; 2.6, 0.3; 3.2, 0.3; 6.4, 0.3 ]; 
            params.V = [   1,   0;  0,    0;   0,   0;   0,   0;  -2,   0 ];

        case '1.3.b' % Billiard balls, Impulse model, e < 1, Leapfrog
            % System dimensions
            params.N_dim      = 2;
            params.BoxHeight  = 2;
            params.BoxWidth   = 7;
            params.g_n        = 0;
            % Model
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 0.8;
            params.k_p        = 0.1;
            % Stepper
            params.t_f        = 5;
            params.t_inc      = 5;
            params.stepper    = enum.Leapfrog;
            % Flags
            params.Animate      = true;
            params.PlotContacts = true;
            params.PlotEnergy   = true;
            params.PlotMomentum = true;
            % Particles
            % 5 billiard balls in 2-D with different masses and different velocities:
            params.N_p   = 5;
            params.M = [   5,   5;  2,    2;   3,   3;   4,   4;   5,   5 ];
            params.R = [      0.3;      0.3;      0.3;      0.3;      0.3 ];
            params.X = [ 0.5, 0.3;  2,  0.3; 2.6, 0.3; 3.2, 0.3; 6.4, 0.3 ]; 
            params.V = [   1,   0;  0,    0;   0,   0;   0,   0;  -2,   0 ];

        case '2.1.a' % Billiard balls, Spring model, k_d = 0, Semi-Implicit Euler
            % System dimensions
            params.N_dim      = 2;
            params.BoxHeight  = 2;
            params.BoxWidth   = 7;
            params.g_n        = 0;
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 4000;
            params.k_d        = 0;
            % Stepper
            params.t_f        = 5;
            params.t_inc      = 5;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate      = true;
            params.PlotContacts = true;
            params.PlotEnergy   = true;
            params.PlotMomentum = true;
            % Particles
            % 5 billiard balls in 2-D with different masses and different velocities:
            params.N_p   = 5;
            params.M = [   5,   5;  2,    2;   3,   3;   4,   4;   5,   5 ];
            params.R = [      0.3;      0.3;      0.3;      0.3;      0.3 ];
            params.X = [ 0.5, 0.3;  2,  0.3; 2.6, 0.3; 3.2, 0.3; 6.4, 0.3 ]; 
            params.V = [   1,   0;  0,    0;   0,   0;   0,   0;  -2,   0 ];

        case '2.1.b' % Billiard balls, Impulse model, k_d > 1, Semi-implicit Euler
            % System dimensions
            params.N_dim      = 2;
            params.BoxHeight  = 2;
            params.BoxWidth   = 7;
            params.g_n        = 0;
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 4000;
            params.k_d        = 50;
            % Stepper
            params.t_f        = 5;
            params.t_inc      = 5;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = true;
            params.PlotContacts = true;
            params.PlotEnergy   = true;
            params.PlotMomentum = true;
            % Particles
            % 5 billiard balls in 2-D with different masses and different velocities:
            params.N_p   = 5;
            params.M = [   5,   5;  2,    2;   3,   3;   4,   4;   5,   5 ];
            params.R = [      0.3;      0.3;      0.3;      0.3;      0.3 ];
            params.X = [ 0.5, 0.3;  2,  0.3; 2.6, 0.3; 3.2, 0.3; 6.4, 0.3 ]; 
            params.V = [   1,   0;  0,    0;   0,   0;   0,   0;  -2,   0 ];

        case '2.2.a' % Piling particles, Impulse model, Leapfrog, h = 0.01
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -10;
            params.LTick      = 2;
            % Particles
            params.N_p        = 1;
            params.T_emit     = 1; % emit new particle every second
            params.dx_emit    = 0; % no spacing between particles
            params.M          = [ 1, 1 ];
            params.R          = 0.3;
            params.X          = 0.3;
            params.V          = 0.0;
            % Model
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 0.9;
            params.k_p        = 0.1;
            % Stepper
            params.t_f        = 20;
            params.t_inc      = 1;
            params.stepper    = enum.Leapfrog;
            params.h_hard     = 0.01;
            % Flags
            params.Animate    = true;
            params.PAsSpheres = true;
            params.PlotEnergy = true;

        case '2.2.b' % Piling particles, Impulse model, Leapfrog, h = 0.005
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -10;
            params.LTick      = 2;
            % Particles
            params.N_p        = 15;
            params.T_emit     = 1; % emit new particle every second
            params.dx_emit    = 0; % no spacing between particles
            % Model
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 0.9;
            params.k_p        = 0.1;
            % Stepper
            params.t_f        = 10;
            params.t_inc      = 1;
            params.stepper    = enum.Leapfrog;
            params.h_hard     = 0.005;
            % Flags
            params.Animate    = true;
            params.PAsSpheres = true;
            params.PlotEnergy = true;

        case '2.3.a' % Piling particles, Spring model, k_s = 4000, Semi-implicit Euler
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -10;
            params.LTick      = 2;
            % Particles
            params.N_p        = 10;
            params.T_emit     = 0.1; % emit new particle every 100 ms
            params.dx_emit    = 0; % no spacing between particles
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 4000;
            params.k_d        = 40;
            % Stepper
            params.t_f        = 12;
            params.t_inc      = 1;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = true;
            params.PAsSpheres = true;
            params.SaveData   = true;

        case '2.3.b' % Piling particles, Spring model, k_s = 8000, Semi-implicit Euler
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -10;
            params.LTick      = 2;
            % Particles
            params.N_p        = 10;
            params.T_emit     = 0.1; % emit new particle every 100 ms
            params.dx_emit    = 0; % no spacing between particles
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 8000;
            params.k_d        = 40;
            % Stepper
            params.t_f        = 12;
            params.t_inc      = 1;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = true;
            params.PAsSpheres = true;
            params.SaveData   = true;

        case '2.3.c' % Piling particles, Spring model, k_s = 4000, k_d = 0, 
                     % Semi-implicit Euler
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -10;
            params.LTick      = 2;
            % Particles
            params.N_p        = 10;
            params.T_emit     = 0.1; % emit new particle every 100 ms
            params.dx_emit    = 0; % no spacing between particles
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 4000;
            params.k_d        = 0;
            % Stepper
            params.t_f        = 12;
            params.t_inc      = 1;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = true;
            params.PAsSpheres = true;
            params.SaveData   = true;

        case '2.3.d' % Piling particles, Spring model, k_s = 8000, k_d = 0, 
                     % Semi-implicit Euler
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -10;
            params.LTick      = 2;
            % Particles
            params.N_p        = 10;
            params.T_emit     = 0.1; % emit new particle every 100 ms
            params.dx_emit    = 0; % no spacing between particles
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 8000;
            params.k_d        = 0;
            % Stepper
            params.t_f        = 12;
            params.t_inc      = 1;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = true;
            params.PAsSpheres = true;
            params.SaveData   = true;

        case '2.3.e' % Piling particles, Spring model, k_s = 40000, Semi-implicit Euler
            % System dimensions
            params.N_dim      = 1;
            params.BoxHeight  = -10;
            params.LTick      = 2;
            % Particles
            params.N_p        = 10;
            params.T_emit     = 0.1; % emit new particle every 100 ms
            params.dx_emit    = 0; % no spacing between particles
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 40000;
            params.k_d        = 40;
            % Stepper
            params.t_f        = 12;
            params.t_inc      = 1;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.Animate    = true;
            params.PAsSpheres = true;
            params.SaveData   = true;

        case '2.4.a' % 2D, Impulse model, Leapfrog
            % System dimensions
            params.N_dim      = 2;
            % Particles
            params.N_p        = 1;
            params.X          = [ 3, 1 ];
            params.V          = [ -0.1, 1 + rand ];
            params.T_emit     = 0.2; % emit new particle every 200 ms
            params.dv_emit    = 0.5; % variate horizontal velocity +/- 0.5 m
            % Model 
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 0.5;
            params.k_p        = 0.4;
            % Stepper
            params.t_f        = 10;
            params.t_inc      = 1;
            params.stepper    = enum.Leapfrog;
            params.h_hard     = 0.002;
            % Flags
            params.PAsSpheres = true;
            params.MaxFigSize = 12; % Fig dimensions, in centimeters

        case '2.4.b' % 2D, Spring model, Semi-implicit Euler'
            % System dimensions
            params.N_dim      = 2;
            % Particles
            params.N_p        = 1;
            params.X          = [ 3, 1 ];
            params.V          = [ -0.1, 1 + rand ];
            params.T_emit     = 0.2; % emit new particle every 200 ms
            params.dv_emit    = 0.5; % variate horizontal velocity +/- 0.5 m
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 8000;
            params.k_d        = 50;
            % Stepper
            params.t_f        = 10;
            params.t_inc      = 1;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.PAsSpheres = true;
            params.MaxFigSize = 12; % Fig dimensions, in centimeters
            
        case '2.5.a' % 3D, Impulse model, Leapfrog
            % System dimensions
            params.N_dim      = 3;
            % Particles
            params.N_p        = 1;
            params.X          = [ 2, 2, 1 ];
            params.V          = [ -0.1, 0.1, 1 + rand ];
            params.T_emit     = 0.2; % emit new particle every 200 ms
            params.dv_emit    = 0.5; % variate horizontal velocity +/- 0.5 m
            % Model
            params.imodel     = enum.ImpulseModel;
            params.k_e        = 0.5;
            params.k_p        = 0.4;
            % Stepper
            params.t_f        = 70;
            params.t_inc      = 1;
            params.stepper    = enum.Leapfrog;
            params.h_hard     = 0.002;
            % Flags
            params.MaxFigSize = 12; % Fig dimensions, in centimeters

        case '2.5.b' % 3D, Spring model, Semi-implicit Euler'
            % System dimensions
            params.N_dim      = 3;
            % Particles
            params.N_p        = 1;
            params.X          = [ 2, 2, 1 ];
            params.V          = [ -0.1, 0.1, 1 + rand ];
            params.T_emit     = 0.2; % emit new particle every 200 ms
            params.dv_emit    = 0.5; % variate horizontal velocity +/- 0.5 m
            params.dx_emit    = 0;
            % Model
            params.imodel     = enum.SpringModel;
            params.k_s        = 8000;
            params.k_d        = 100;
            % Stepper
            params.t_f        = 70;
            params.t_inc      = 1;
            params.stepper    = enum.SemiImplicitEuler;
            % Flags
            params.MaxFigSize = 12; % Fig dimensions, in centimeters
            
        otherwise
            error( 'lab1_defs:invarg', 'Unercognized configuration ID %s', configId );
    end

    if nargout >= 1
        varargout{1} = params;
    end

end % function lab2_defs

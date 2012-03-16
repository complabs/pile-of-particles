%% LAB2: Colliding Particles; Front-End Application (gui code for lab2.fig) 
% Allows the user to interactivelly configure various parameters and run according
% simulation or set of simulations. Simulations can be paused/continued and parameters 
% changed inbetween. Further, the simulation state variables and trace data (e.g.
% trajectories or energies) can be saved as snapshots in MAT-files and later reused 
% to start new simulations or to perform further data analysis.
%
% Filename: lab2.m
% Date:     2012-03-11
% Author:   Mikica B Kocic
%
%#ok<*DEFNU,*INUSL,*INUSD>

function varargout = lab2( varargin )

    % ------------------------------------------------------------------------------------
    % Begin initialization code - DO NOT EDIT
    % ------------------------------------------------------------------------------------
    gui_Singleton = 1;

    gui_State = struct( 'gui_Name',       mfilename, ...
                        'gui_Singleton',  gui_Singleton, ...
                        'gui_OpeningFcn', @lab2_OpeningFcn, ...
                        'gui_OutputFcn',  @lab2_OutputFcn, ...
                        'gui_LayoutFcn',  [] , ...
                        'gui_Callback',   [] );

    if nargin && ischar( varargin{1} ) && varargin{1}(1) ~= '-'
        gui_State.gui_Callback = str2func( varargin{1} );
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn( gui_State, varargin{:} );
    else
        gui_mainfcn( gui_State, varargin{:} );
    end
    % ------------------------------------------------------------------------------------
    % End initialization code - DO NOT EDIT
    % ------------------------------------------------------------------------------------

% ========================================================================================
%% Nested functions ----------------------------------------------------------------------

% ----------------------------------------------------------------------------------------
%% - initialize_GUI( hObject, handles )
% Initialzes values and strings of the GUI components from params structure

function initialize_GUI( hObject, handles )

    global lab2_running   % Global flag used to enable simulation loop in lab2_main

    enum = handles.params.enum;
    imodel = handles.params.imodel;

    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.combo_N_dim,     'Value',           handles.params.N_dim         );
    set( handles.edit_g_n,        'String', num2str( handles.params.g_n          ));
    set( handles.edit_N_p,        'String', num2str( handles.params.N_p          ));
    set( handles.edit_T_emit,     'String', num2str( handles.params.T_emit       ));
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.edit_m_p,        'String', num2str( handles.params.m_p          ));
    set( handles.edit_r_p,        'String', num2str( handles.params.r_p          ));
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.edit_BoxHeight,  'String', num2str( handles.params.BoxHeight    ));
    set( handles.edit_BoxWidth,   'String', num2str( handles.params.BoxWidth     ));
    set( handles.edit_BoxDepth,   'String', num2str( handles.params.BoxDepth     ));
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.combo_stepper,   'Value',           handles.params.stepper       );
    set( handles.edit_t_0,        'String', num2str( handles.params.t_0          ));
    set( handles.edit_t_f,        'String', num2str( handles.params.t_f          ));
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.rb_SpringModel,  'Value',           imodel == enum.SpringModel   );
    set( handles.rb_ImpulseModel, 'Value',           imodel == enum.ImpulseModel  );
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.edit_h_soft,     'String', num2str( handles.params.h_soft       ));
    set( handles.edit_k_s,        'String', num2str( handles.params.k_s          ));
    set( handles.edit_k_d,        'String', num2str( handles.params.k_d          ));
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.edit_h_hard,     'String', num2str( handles.params.h_hard       ));
    set( handles.edit_k_e,        'String', num2str( handles.params.k_e          ));
    set( handles.edit_k_p,        'String', num2str( handles.params.k_p          ));
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.cb_Animate,      'Value',           handles.params.Animate       );
    set( handles.cb_InRealTime,   'Value',           handles.params.InRealTime    );
    set( handles.cb_PlotPath,     'Value',           handles.params.PlotPath      );
    set( handles.cb_PAsSpheres,   'Value',           handles.params.PAsSpheres    );
    set( handles.cb_Wiggle3D,     'Value',           handles.params.Wiggle3D      );
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.edit_fps,        'String', num2str( handles.params.fps          ));
    %    ------------------------ --------  -------- -------------------------- ---
    set( handles.cb_SaveData,     'Value',           handles.params.SaveData      );
    set( handles.cb_PlotContacts, 'Value',           handles.params.PlotContacts  );
    set( handles.cb_PlotEnergy,   'Value',           handles.params.PlotEnergy    );
    set( handles.cb_PlotHeight,   'Value',           handles.params.PlotHeight    );
    set( handles.cb_PlotMomentum, 'Value',           handles.params.PlotMomentum  );
    set( handles.cb_CaptureMovie, 'Value',           handles.params.CaptureMovie  );
    %    ------------------------ --------  -------- -------------------------- ---

    % Reset button caption in java, if using java.swing
    if usejava( 'swing' )
        set( handles.button_Reset, 'String', ...
            '<html>Reset <em>t</em><small><sub>0</sub></small></html>' );
    end
            
    % Setup button Reset depending on whether lab2_main is running or not
    if lab2_running
        set( handles.button_Reset, 'Enable', 'off' );
    else
        set( handles.button_Reset, 'Enable', 'on' );
    end

    % Setup button Start depending on wheter lab2_main is running or not 
    % and if not running on value of t_0
    if lab2_running
        set( handles.button_Start, 'ForegroundColor', [ 0, 0, 0 ] );
        set( handles.button_Start, 'FontWeight', 'normal' );
        set( handles.button_Start, 'String', 'Stop' );
    else
        set( handles.button_Start, 'ForegroundColor', [ 0, 0, 0 ] );
        set( handles.button_Start, 'FontWeight', 'normal' );
        if handles.params.t_0 == 0
            set( handles.edit_t_0, 'ForegroundColor', [ 0, 0, 0 ] );
            set( handles.button_Start, 'String', 'Start', ...
                'ForegroundColor', [ 0, 0, 0 ], 'FontWeight', 'normal' );
            set( handles.button_Reset, 'Enable', 'off' );
        else
            set( handles.edit_t_0, 'ForegroundColor', [ 0, 0, 0.8 ] );
            set( handles.button_Start, 'String', 'Continue', ...
                'ForegroundColor', [ 0, 0, 0.8 ], 'FontWeight', 'bold' );
        end
    end

    guidata( hObject, handles ); % Update handles structure

% ----------------------------------------------------------------------------------------
%% - guidata_UpdateDimensions( hObject, handles, N_dim )
% Enables/disables box dimensions depending on configured number of space dimesnsions.

function guidata_UpdateDimensions( hObject, handles, N_dim )

    enable = { 'on', 'off', 'off' ; ...  % N_dim == 1: H
               'on', 'on',  'off' ; ...  % N_dim == 2: H W
               'on', 'on',  'on'  };     % N_dim == 3: H W D

    set( handles.text_BoxHeight,  'Enable', enable{ N_dim, 1 } );
    set( handles.edit_BoxHeight,  'Enable', enable{ N_dim, 1 } );
    set( handles.textu_BoxHeight, 'Enable', enable{ N_dim, 1 } );
    set( handles.text_BoxWidth,   'Enable', enable{ N_dim, 2 } );
    set( handles.edit_BoxWidth,   'Enable', enable{ N_dim, 2 } );
    set( handles.textu_BoxWidth,  'Enable', enable{ N_dim, 2 } );
    set( handles.text_BoxDepth,   'Enable', enable{ N_dim, 3 } );
    set( handles.edit_BoxDepth,   'Enable', enable{ N_dim, 3 } );
    set( handles.textu_BoxDepth,  'Enable', enable{ N_dim, 3 } );

    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
%% - parseEditBox( hObject, handles, field, str2value, isValid, errMsg )
% Parses and validates string from editbox and stores value in 'params' field.
%
% * hObject:   edit box handle
% * handles:   structure with handles and user data (see GUIDATA)
% * field:     name of the field in handles.params structure
% * str2Value: function handle that converts string to param value
% * isError:   function handle that validates param value
% * errMsg:    error message displayed if validation fails

function parseEditBox( hObject, handles, field, str2value, isValid, errMsg )

    val = str2value( get( hObject, 'String' ) );

    if isValid( val )
        set( hObject, 'String', val );
        guidata( hObject, handles ) % Update handles structure
        handles.params.(field) = val;
    else % error:
        uicontrol( hObject );
        errordlg( errMsg, 'Lab2: Validation Error' );
        set( hObject, 'String', handles.params.(field) );
    end
    
    guidata( hObject, handles ) % Update handles structure
    
% ========================================================================================
%% App Event Handlers --------------------------------------------------------------------

% ----------------------------------------------------------------------------------------
%% - lab2_OpeningFcn 
% Executes just before lab2 is made visible.
%
% * hObject:    handle to figure
% * eventdata:  reserved - to be defined in a future version of MATLAB
% * handles:    structure with handles and user data (see GUIDATA)
% * varargin:   command line arguments to lab2 (see VARARGIN)

function lab2_OpeningFcn( hObject, eventdata, handles, varargin )

    % Choose default command line output for lab2
    handles.output = hObject;

    global lab2_running   % Global flag used to enable simulation loop in lab2_main
    lab2_running = false; % Stop any running simulations
    
    if length(varargin) < 1
        % Load default parameters and update values of GUI components
        handles.params = lab2_defs ();
    else
        % Load parameters from MAT-file (saved as the result of previous simulation)
        fname = varargin{1};
        if isequal( fname, '-' )
            % filename '-' means to get filename from the user interactivelly
            [fname,pname] = uigetfile( '*.mat','Enter data file' );
            if fname == 0 % user responded with 'Stop'
                delete( handles.appFigure ); 
                return; % quit
            end
            fname = [ pname, fname ];
        elseif fname(1) == '-'
            fname = fname(2:end);
        end
        % Load 'params' structure from specified file
        disp([ 'Loading file ''', fname, '''' ]);
        try
            data = load( fname, 'params' );
        catch %#ok<CTCH>
            error( 'lab2:invf', '\nFile ''%s'' not found.', fname );
        end
        % Verify that we got valid params structure ...
        if isfield( data, 'params' ) ...
            && isstruct( data.params ) ...
            && isfield( data.params, 'Title' )
        % then
            handles.params = data.params;
        else
            % Quit as we failed to find params structure
            error( 'lab2:invp', ...
                '\nMAT-file ''%s'' does not contain valid ''params'' structure', fname );
        end
        disp([ 'Loaded: ', handles.params.Title ]);
    end

    % Add new invisible axes (in case someone draws to our figure)
    axes; set( gca, 'Visible', 'off' );

    % Setup captions for combo boxes
    set( handles.combo_Configuration, 'String', testSuite_Captions );
    set( handles.combo_stepper, 'String', handles.params.verb.stepper );

    % Delayed initialization of GUI components using single-shot timer (note that 
    % at this time in OpeningFcn, uicontrol objects are not still created, so changing
    % their properties will result in error)
    start( timer( 'ExecutionMode', 'singleShot', ...
        'TimerFcn', { ...
            @( obj, event, handles ) ...
                combo_Configuration_Callback( ...
                    handles.combo_Configuration, 0, handles ...
                ), ...
            handles ...
        }...
    ) );

% ----------------------------------------------------------------------------------------
%% - lab2_OutputFcn
% Outputs from this function are returned to the command line.
%
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function varargout = lab2_OutputFcn( hObject, eventdata, handles ) 

    % Get default command line output from handles structure
    if isfield( handles, 'output' )
        varargout{1} = handles.output;
    end

% ========================================================================================
%% Button Event Handlers -----------------------------------------------------------------

% ----------------------------------------------------------------------------------------
%% - button_ManualNp Callback

function button_ManualNp_Callback( hObject, eventdata, handles )

    global lab2_running   % Global flag used to enable simulation loop in lab2_main

    if lab2_running
        % Stop any running simulations
        lab2_running = false; 
        return;
    elseif handles.params.t_f < handles.params.t_0
        % Nothing to simulate
        return;
    end

    % Setup mode to manually add particles and 'Start' simulation
    handles.params.ManualNp = true;
    button_Start_Callback( hObject, eventdata, handles );
    
% ----------------------------------------------------------------------------------------
%% - button_Start Callback

function button_Start_Callback( hObject, eventdata, handles )

    global lab2_running   % Global flag used to enable simulation loop in lab2_main

    if lab2_running
        % Stop any running simulations
        lab2_running = false; 
        return;
    elseif handles.params.t_f < handles.params.t_0
        % Nothing to simulate
        return;
    end

    % Enable lab2_main to run
    lab2_running = true;
    
    % Remember simulation final time (to check if simulation was not completed)
    t_f = handles.params.t_f; 

    if handles.params.ManualNp
        if handles.params.N_dim == 3
            handles.params.N_dim = 2; % works only in 1- or 2-D
        end
        t_inc = handles.params.t_f;
        handles.params.N_p = 0;
        handles.params.t_0 = 0;
        handles.params.t_f = 0;
        handles.params.T_emit = Inf;
    else
        t_inc = handles.params.t_inc; % increment final time, as user prefers
    end

    % Update button status to reflect that we are running application
    set( handles.button_Start, 'String', 'Stop' );
    set( handles.button_Start, 'FontWeight', 'bold' );
    set( handles.button_Start, 'ForegroundColor', [ 0.8, 0, 0 ] );
    set( handles.button_Reset, 'Enable', 'off' );
    
    guidata( hObject, handles ) % Update handles structure 
    drawnow; % and refresh window before proceeding with simulation

    % Start simulation
    try
        handles.params = lab2_main( handles.params );
    catch err
        % Change global flag indicating that lab2_main is not running
        lab2_running = false;

        % Update current params from lab2_main result data
        drawnow; set( 0, 'CurrentFigure', handles.appFigure );
        initialize_GUI( handles.appFigure, handles );
        panel_Simulation_SelectionChangeFcn( handles.appFigure, 0, handles );

        % Report recognized exceptions separately
        if strcmp( err.identifier, 'MATLAB:hg:dt_conv:Matrix_to_HObject:BadHandle' ) ...
        || strcmp( err.identifier, 'lab2_main:feim' )
            disp( '>>> Simulation interrupted' );
            disp([ '>>> ', err.message, ':' ]);
            if ~isempty( err.stack )
                disp( err.stack(1) );
            end
            return
        end

        % If not recognized, reissue the exception
        rethrow( err );
    end

    % Adjust simulation time
    if handles.params.t_f < t_f
        % Simulation was stopped prematurelly; prepare to continue:
        handles.params.t_0 = handles.params.t_f;
        handles.params.t_f = t_f;
    else
        % Simulation ended; increment the final time for a given amount
        handles.params.t_0 = handles.params.t_f;
        handles.params.t_f = handles.params.t_0 + t_inc;
    end

    % Update current params from lab2_main result data
    drawnow; set( 0, 'CurrentFigure', handles.appFigure );
    initialize_GUI( handles.appFigure, handles );
    panel_Simulation_SelectionChangeFcn( handles.appFigure, 0, handles );

    % Popup this application
    drawnow;
    figure( handles.appFigure );

% ----------------------------------------------------------------------------------------
%% - button_Reset Callback

function button_Reset_Callback( hObject, eventdata, handles )

    global lab2_running   % Global flag used to enable simulation loop in lab2_main
    if lab2_running
        return  % don't reset while simulation is running
    end

    % Reset initial time back to 0 and remove saved state variables
    handles.params.t_0 = 0; % Reset initial time

    % Reload saved state variables
    orig_N_p = handles.params.N_p;

    if isfield( handles.params, 'M_0' )
        handles.params.M = handles.params.M_0;
        handles.params.N_p = min( handles.params.N_p, size( handles.params.M, 1 ) );
    elseif isfield( handles.params, 'M' )
        handles.params = rmfield( handles.params, 'M' );
    end
    
    if isfield( handles.params, 'R_0' )
        handles.params.R = handles.params.R_0;
        handles.params.N_p = min( handles.params.N_p, size( handles.params.R, 1 ) );
    elseif isfield( handles.params, 'R' )
        handles.params = rmfield( handles.params, 'R' );
    end

    if isfield( handles.params, 'X_0' )
        handles.params.X = handles.params.X_0;
        handles.params.N_p = min( handles.params.N_p, size( handles.params.X, 1 ) );
    elseif isfield( handles.params, 'X' )
        handles.params = rmfield( handles.params, 'X' );
    end

    if isfield( handles.params, 'V_0' )
        handles.params.V = handles.params.V_0;
        handles.params.N_p = min( handles.params.N_p, size( handles.params.V, 1 ) );
    elseif isfield( handles.params, 'V' )
        handles.params = rmfield( handles.params, 'V' );
    end

    if isfield( handles.params, 'faceColor_0' )
        handles.params.faceColor = handles.params.faceColor_0;
    elseif isfield( handles.params, 'faceColor' )
        handles.params = rmfield( handles.params, 'faceColor' );
    end

    % Update N_p
    if handles.params.N_p == 0
        handles.params.N_p = orig_N_p;
    end
    set( handles.edit_N_p,'String', num2str( handles.params.N_p ) );

    % Reconfigure buttons and initial time info
    set( handles.edit_t_0, 'String', num2str( handles.params.t_0 ), ...
        'ForegroundColor', [ 0, 0, 0 ] );
    set( handles.button_Start, 'String', 'Start', ...
        'ForegroundColor', [ 0, 0, 0 ], 'FontWeight', 'normal' );
    set( handles.button_Reset, 'Enable', 'off' );

    guidata( hObject, handles ); % Update handles structure

% ----------------------------------------------------------------------------------------
%% - button_Cascade Callback
% Cascade existing figures so that they don't directly overlap

function button_Cascade_Callback( hObject, eventdata, handles )

    % Find and sort all figures in order of creation
    figs = sort( findobj( 0, 'Type', 'figure' ) );

    % Remove our figure from the list
    figs( figs == handles.appFigure ) = []; 

    if length( figs ) <= 1
        return; % Nothing to do...
    end

    % Hide all figures
    set( figs, 'Visible', 'off' )

    % Cascade image positions for d_horz/d_vert amount
    d_horz = 30;
    d_vert = 30;
    for n = 2:length(figs)
        pPrev = get( figs(n-1), 'Position' ); % gets the rectangle of the previous figure
        pos = get( figs(n), 'Position' ); % gets the rectangle of the current figure
        % Now, move the current figure (change only position) relative to the previous
        pos(1:2) = pPrev(1:2) + [ d_horz, pPrev(4) - d_vert ] - [ 0, pos(4) ];
        set( figs(n), 'Position',  pos );
    end

    % Make all figures in order and visible
    set( 0, {'CurrentFigure'}, mat2cell( figs(n), 1 ) );
    set( figs, 'Visible', 'on' )

    % Popup our window
    set( handles.appFigure, 'Visible', 'on' );
    figure( handles.appFigure );
    
    % global ffigs; ffigs = figs; % debug

% ----------------------------------------------------------------------------------------
%% - button_ExportAnnotation Callback

function button_ExportAnnotation_Callback( hObject, eventdata, handles )

    if isfield( handles.params.figCaption )
        exportLaTeX( [ 'lab2_info_', handles.params.ID ], handles.params.figCaption );
    end

% ========================================================================================
%% EditBox Event Handlers  ---------------------------------------------------------------
% Perform parsing and validation of input data and the update of the respective fields
% in 'params' structure.

% ----------------------------------------------------------------------------------------
function edit_g_n_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'g_n', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) ), ... is valid?
        'Acceleration of gravity must be a real number.' ... error if not
    );

% ----------------------------------------------------------------------------------------
function edit_N_p_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'N_p', ...
        @(str) round( str2double( str ) ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 ), ... is valid?
        'Number of particles must be either a positive integer or zero' ... msg if not
    );

% ----------------------------------------------------------------------------------------
function edit_T_emit_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'T_emit', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && val > 0 ), ... is valid?
        'Emitter period must be either Inf or a positive real number' ... msg if not
    );

% ----------------------------------------------------------------------------------------
function edit_m_p_Callback( hObject, eventdata, handles ) 

    parseEditBox( hObject, handles, 'm_p', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val > 0 ), ... is valid?
        'Particle mass must be a real number greater than zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_r_p_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'r_p', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val > 0 ), ... is valid?
        'Particle radius must be a real number greater than zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_BoxHeight_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'BoxHeight', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val ~= 0 ), ... is valid?
        'Box height must be a real number not equal to zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_BoxWidth_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'BoxWidth', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val ~= 0 ), ... is valid?
        'Box width must be a real number not equal to zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_BoxDepth_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'BoxDepth', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val ~= 0 ), ... is valid?
        'Box depth must be a real number not equal to zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_t_0_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 't_0', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 ), ... is valid?
        'Initial time must be a real number greater or equal than zero' ... if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_t_f_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 't_f', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 ), ... is valid?
        'Final time must be a real number greater or equal than zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_h_soft_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'h_soft', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val > 0 ), ... is valid?
        'Time-step must be a real number greater than zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_k_s_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'k_s', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 ), ... is valid?
        'Spring coefficent must be either a positive real number or zero' ... if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_k_d_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'k_d', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 ), ... is valid?
        'Spring dumping coefficent must be either a positive real number or zero' ...
    );

% ----------------------------------------------------------------------------------------
function edit_h_hard_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'h_hard', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val > 0 ), ... is valid?
        'Time-step must be a real number greater than zero' ... msg if not valid
    );

% ----------------------------------------------------------------------------------------
function edit_k_e_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'k_e', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 && val <= 1 ), ... is valid?
        'Coefficient of restitution must be a real number between 0 and 1 ' ... msg if not
    );

% ----------------------------------------------------------------------------------------
function edit_k_p_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'k_p', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 && val <= 1 ), ... is valid?
        'Position projection factor must be a real number between 0 and 1 ' ... msg if not
    );

% ----------------------------------------------------------------------------------------
function edit_fps_Callback( hObject, eventdata, handles )

    parseEditBox( hObject, handles, 'fps', ...
        @(str) str2double( str ), ... convert str to value
        @(val)( ~isnan(val) && ~isinf(val) && val >= 0 ), ... is valid?
        'Frames per second must be a positive real number or 0' ... msg if not valid
    );

% ========================================================================================
%% CheckBox Event Handlers ---------------------------------------------------------------
% Update respective fields in params structure with new checkbox values.

% ----------------------------------------------------------------------------------------
function cb_ManualNp_Callback( hObject, eventdata, handles )

    handles.params.ManualNp = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

function cb_Animate_Callback( hObject, eventdata, handles )

    handles.params.Animate = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_InRealTime_Callback( hObject, eventdata, handles )

    handles.params.InRealTime = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_PlotPath_Callback( hObject, eventdata, handles )

    handles.params.PlotPath = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_PAsSpheres_Callback( hObject, eventdata, handles )

    handles.params.PAsSpheres = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_Wiggle3D_Callback( hObject, eventdata, handles )

    handles.params.Wiggle3D = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_SaveData_Callback( hObject, eventdata, handles )

    handles.params.SaveData = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_PlotEnergy_Callback( hObject, eventdata, handles )

    handles.params.PlotEnergy = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_PlotHeight_Callback( hObject, eventdata, handles )

    handles.params.PlotHeight = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_PlotMomentum_Callback( hObject, eventdata, handles )

    handles.params.PlotMomentum = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_PlotContacts_Callback( hObject, eventdata, handles )

    handles.params.PlotContacts = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
function cb_CaptureMovie_Callback( hObject, eventdata, handles )

    handles.params.CaptureMovie = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ========================================================================================
%% ComboBox Event Handlers ---------------------------------------------------------------

% ----------------------------------------------------------------------------------------
%% - testSuite_Captions
% Returns captions for combo_Configuration uicontrol

function captions = testSuite_Captions

    [~,captions] = lab2_defs;   % get list of recognized configurations

% ----------------------------------------------------------------------------------------
%% - combo_Configuration Callback
% Executes on selection change in combo_Configuration.
%
% * hObject    handle to combo_Configuration (see GCBO)
% * eventdata  reserved - to be defined in a future version of MATLAB
% * handles    structure with handles and user data (see GUIDATA)

function combo_Configuration_Callback( hObject, eventdata, handles )

    id = get( hObject, 'Value' );
    info = cellstr( get( hObject, 'String' ) );

    configId = strtrim( sscanf( info{id}, ' (%[^)])' ) );
    if isempty( configId )
        % try next in the list
        id = id + 1;
        configId = strtrim( sscanf( info{id}, ' (%[^)])' ) );
        if isempty( configId )
            return
        elseif id ~= 2
            set( hObject, 'Value', id );
        end
    end

    handles.params = lab2_defs( configId );
    handles.params.Description = strtrim( info{id} );

    initialize_GUI( handles.appFigure, handles );
    panel_Simulation_SelectionChangeFcn( handles.appFigure, 0, handles );

% ----------------------------------------------------------------------------------------
%% - combo_N_dim Callback
% Executes on selection change in combo_N_dim.
%
% * hObject    handle to combo_N_dim (see GCBO)
% * eventdata  reserved - to be defined in a future version of MATLAB
% * handles    structure with handles and user data (see GUIDATA)

function combo_N_dim_Callback( hObject, eventdata, handles )

    handles.params.N_dim = get( hObject, 'Value' );

    handles.params.PAsSpheres = handles.params.N_dim == 3;
    set( handles.cb_PAsSpheres, 'Value', handles.params.PAsSpheres );

    % Update gui to reflect new N_dim
    guidata_UpdateDimensions( hObject, handles, handles.params.N_dim );
    
% ----------------------------------------------------------------------------------------
%% - combo_stepper Callback
% Executes on selection change in combo_stepper.
%
% * hObject    handle to combo_stepper (see GCBO)
% * eventdata  reserved - to be defined in a future version of MATLAB
% * handles    structure with handles and user data (see GUIDATA)

function combo_stepper_Callback( hObject, eventdata, handles )

    handles.params.stepper = get( hObject, 'Value' );
    guidata( hObject, handles ) % Update handles structure

% ----------------------------------------------------------------------------------------
%% - panel_Simulation SelectionChangeFcn
% Executes when selected object is changed in panel_Simulation.
%
% * hObject    handle to combo_stepper (see GCBO)
% * eventdata:  structure with the following fields (see UIBUTTONGROUP)
%   EventName:  string 'SelectionChanged' (read only)
%   OldValue:   handle of the previously selected object or empty if none was selected
%   NewValue:   handle of the currently selected object
% * handles     structure with handles and user data (see GUIDATA)

function panel_Simulation_SelectionChangeFcn( hObject, eventdata, handles )

    % Colors that will be used for enabled/disabled controls
    enabled_color = get( handles.panel_Simulation, 'ForegroundColor' );
    disabled_color = [ 0.6 0.55 0.55 ];

    springModel = get( handles.rb_SpringModel, 'Value' );
    
    % Setup flags and colors depending on mode
    if springModel
        sm_enable = 'on';   sm_color = enabled_color;
        im_enable = 'off';  im_color = disabled_color;
    else
        sm_enable = 'off';  sm_color = disabled_color;
        im_enable = 'on';   im_color = enabled_color;
    end

    % Apply flags/colors to Spring Model panel and its children
    set( handles.rb_SpringModel, 'ForegroundColor', sm_color );
    for i = get( handles.panel_SpringModel, 'Children' )
        set( i, 'Enable', sm_enable );
    end
    
    % Apply flags/colors to Impulse Model panel and its children
    set( handles.rb_ImpulseModel, 'ForegroundColor', im_color );
    for i = get( handles.panel_ImpulseModel, 'Children' )
        set( i, 'Enable', im_enable );
    end
    
    % Parse interaction model and stepper
    enum = handles.params.enum;
    if springModel
        handles.params.imodel  = enum.SpringModel;
    else
        handles.params.imodel  = enum.ImpulseModel;
    end

    % Reconfigure also dimensions
    guidata_UpdateDimensions( hObject, handles, handles.params.N_dim );

% ========================================================================================
%% CreateFcn Event Handler (for all objects) ---------------------------------------------
% Executes during object creation, after setting all properties.

function edit_CreateFcn( hObject, eventdata, handles )

    if ispc && isequal( get( hObject, 'BackgroundColor' ), ...
                        get( 0, 'defaultUicontrolBackgroundColor' ) )
    % then            
        set( hObject, 'BackgroundColor', 'white' );
    end

function combo_Configuration_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function combo_N_dim_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_g_n_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_N_p_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_T_emit_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_m_p_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_r_p_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_BoxHeight_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_BoxWidth_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_BoxDepth_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function combo_stepper_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_t_0_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_t_f_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_h_soft_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_k_s_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_k_d_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_h_hard_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_k_e_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_k_p_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

function edit_fps_CreateFcn( hObject, eventdata, handles )

    edit_CreateFcn( hObject, eventdata, handles );

% ========================================================================================

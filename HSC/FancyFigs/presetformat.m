%PRESETFORMAT applies pre-defined formats for the BAFIGEXPRO suite.
%The different styles available here correspond to various aspects of a 
%figure (font, dimension, color, margin, etc). Within each aspect, a number
%of mutually exclusive styles are available. The various aspects can (and
%should) be combined. The styles are simply sets of consistent options of
%FORMATFIG and FIG2PUBLIC. PRESETFORMAT is accessed from FIG2PUBLIC via the
%'base_config' parameter.
%
%Available pre-defined formats:
%
%   'article'       Defines font size, marker size and line width. Smaller
%       OR          for the 'article' mode and larger for the 'slides'
%   'slides'        mode. Figure titles are removed. For 'slides', the axes
%                   dimensions are adapted to avoid label cropping.
%
%   'small'         Defines figure width. 'small' corresponds to the half
%       OR          of the text width of a standard A4 page, while 'large'
%   'large'         corresponds to the full text width. For 'small', the
%                   axes dimensions are adapted to avoid label clipping.
%
%   'color'         Defines the figure color. 'color' doesn't change
%       OR          anything. 'bwtl' maps the solid line colors to line
%   'bwtl'          styles, sets all the text and lines to black and 
%       OR          changes the colormap to get a b/w printer friendly
%   'bw'            set of colors. 'bw' is equivalent to 'bwtl' but a gray
%                   colormap is used and the 'eps_format' is changed to
%                   export in gray scale.
%
%   'tight'         Decreases the figure margins while keeping one or both
%       OR          figure dimension untouched (see FORMATFIG for details).
%   'tightx'        The 'eps_loose' option is also activated to obtain
%       OR          a perfect match between the formatted figure and the
%   'tighty'        exported one.
%
%   'aligned'       If given, the axes alignment is conserved during
%                   figure resizing. Internal labels of multiple subplots
%                   might be cropped in that case.
%
%   'log'           If given, the axis ticks are locked. Particularly 
%                   useful to avoid bad re-ticking of log axis.
%
%   'heavy'         If given, a bitmap renderer is used instead of the 
%                   default vectorial renderer. Useful for large files.
%
%   'fixinside'     If given, the text boxes are forced to remain inside
%                   their axes at figure resizing. In addition, legends
%                   are repositioned to the 'best' location.
%
%   'old'           If given, the original formatting is used for 'article'
%                   and 'slides'.
%
%USAGE:
%
%   When calling FIG2PUBLIC, give the 'base_config' parameter a combination
%   of the pre-defined formats described above. For example:
%
%   Publications:
%       ..., 'base_config', 'article small bwtl tightx aligned', ...
%
%   Reports:
%       ..., 'base_config', 'article large bwtl tightx aligned', ...
%
%   Slides:
%       ..., 'base_config', 'slides large color tightx aligned', ...
%
%   NB: The format strings can also be juxtaposed ('articlecolor') or 
%       separated by '-' or '_' ('article_color-tight').
%
%
%  Last modification: Jonathan Rossel, 15.12.2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Modif log: (see file)


%Inputs: 
%   formatfigopt, expfigopt     option structures defined in fig2public
%   bc                          value of 'base_config' parameter
%
%Outputs:
%   updated option structures, using the pre-defined styles.

function [ formatfigopt, expfigopt ] = presetformat( formatfigopt, expfigopt, bc );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inits

%Structure of switches
%Use 2 levels. Switches in the 2nd level are mutually exclusive.
%Styles are applied in the order of the fields of this structure.
%The user only refers to 2nd level switch names.
%All the 2nd level switches must be different. The first match
%in the user input is kept (i.e. 'tightx' is kept, and not 'tight').
swistyle = struct( ...
    'base', struct( ...                 %Overall formatting
        'article', false, ...
        'slides', false ), ...
    'width', struct( ...                %Set absolute figure width
        'small', false, ...
        'large', false ), ...
    'colormode', struct( ...            %Figure color
        'color', false, ...
        'bwtl', false, ...
        'bw', false ), ...
    'margin', struct( ...               %Margin handling and loose export
        'tightx', false, ...
        'tighty', false, ...
        'tight', false ), ...
    'alignment', struct( ...            %Alignement conservation (active position)
        'aligned', false ), ...
    'logaxis', struct( ...              %Tick locking
        'log', false ), ...
    'heavyfigs', struct( ...            %Bitmap renderer
        'heavy', false ), ...
    'innerelems', struct( ...           %Legends and texts
        'fixinside', false ), ...
    'oldsettings', struct( ...          %To keep backward compatibility
        'old', false ) ...
    );

%version
verstr = version;
majorver = str2num(verstr(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parse and check input

level1fn = fieldnames( swistyle );

for ii = 1 : length( level1fn )
    
    %Loop on first level
    
    %2nd level field names    
    level2fn = fieldnames( swistyle.( level1fn{ ii } ) );
    
    %Parse the input
    
    for jj = 1 : length( level2fn )
        
        strind = strfind( lower( bc ), level2fn{ jj } ); %index of field name in user input    
        lstr = length( level2fn{ jj } ); %field name length
            
        if length( strind ) == 1
            
            %this switch has been given (once) by the user. Activate it.
            swistyle.( level1fn{ ii } ).( level2fn{ jj } ) = true;
            
            %remove the switch from the user input, to avoid secondary detections (e.g. for 'tightx' and 'tight' cases)
            bc( strind : strind + lstr - 1 ) = [];
            
        elseif length( strind ) > 1
        
            error( [ 'The switch ''', level2fn{ jj }, ''' has been detected more than once in ''base_config''. Please give it only once.' ] );
            
        end
        
    end
    
    %Check mutual exclusion
    
    swilevel2 = cell2mat( struct2cell( swistyle.( level1fn{ ii } ) ) );
    
    if sum( swilevel2 ) > 1, 
        
        indmut = find( swilevel2 );
        error( [ level2fn{ indmut(1) }, ' and ', level2fn{ indmut(2) }, ' are mutually exclusive and cannot be simultaneously activated.' ] );
        
    end
    
end

%Remove '-' or '_' or ' ' for final check
bctest = strrep( strrep( strrep( bc, ' ', '' ), '-', '' ), '_', '' );

if length( bctest ) ~= 0

    error( [ 'The parsing of ''base_config'' is incomplete. Extra characters remain: ''', bc, '''' ] );
    
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Apply style presets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Base format

if swistyle.base.article
    
    %Article format
    
    formatfigopt.markermode = 'fixed';
    formatfigopt.markersize = 4;
    formatfigopt.fontmode = 'fixed';
    formatfigopt.fontsize = 10;
    formatfigopt.linemode = 'fixed';
    formatfigopt.linewidth = 1;
    formatfigopt.rmtitle = true;
    
    if swistyle.oldsettings.old
        
        %Old settings
        formatfigopt.markersize = 6;
        formatfigopt.rmtitle = false;
    
    end

elseif swistyle.base.slides

    %Slides format
    
    formatfigopt.markermode = 'fixed';
    formatfigopt.markersize = 10;
    formatfigopt.fontmode = 'fixed';
    formatfigopt.fontsize = 16;
    formatfigopt.linemode = 'fixed';
    formatfigopt.linewidth = 2;
    formatfigopt.adaptaxes = true;
    formatfigopt.rmtitle = true;
    
    if majorver < 7,
        formatfigopt.matchcolorbars = true;
    end
    
    if swistyle.oldsettings.old
        
        %Old settings
        formatfigopt.keeptextin = true;
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Width format

if swistyle.width.small
    
    %Small format: half text width of standard A4
    %i.e. with 2.5cm margins: ( 21 - 2*2.5 ) / 2 = 8
    formatfigopt.width = 8;
    
    %Text cropping is highly probable in that case:
    formatfigopt.adaptaxes = true;
    if majorver < 7,
        formatfigopt.matchcolorbars = true;
    end

elseif swistyle.width.large

    %Large format: text width of standard A4
    %i.e. with 2.5cm margins: 21 - 2*2.5 = 16
    formatfigopt.width = 16; %NB: default screen size of Matlab figures is 15.8
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Color format

if swistyle.colormode.color
    
    %Color format
    %Matches FORMATFIG's defaults
    
elseif swistyle.colormode.bwtl

    %bwtl format
    %Text and line objects are converted to black
    %Solid lines are mapped from colors to linestyle
    %Colormap is changed to be more bw printing friendly
    
    %%%%%%%%%
    %Change the colormap, to obtain a final result that is more
    %compatible with bw printing while remaining in colors

    %based on hot, but cut the blackest region so that
    %black lines remain distinguishable
    %Numbering is based on columns in colormap
    blacklim = 0.1;
    f1 = 0; %slope factor
    f3 = 1;
    a = ( 1/(1+f3) + 1 + ( 1 - blacklim ) / (1-f1) ) / 64; %desired increment
    a3 = (1+f3)*a;
    a2 = a;
    a1 = (1-f1)*a;
    varipart3 = linspace( 1, a3, round( 1 / a3 ) )';
    varipart2 = linspace( 1, a2, round( 1 / a2 ) )';
    varipart1 = linspace( 1, 0, round( 1 / a1 ) )';
    nv3 = length( varipart3 );
    nv2 = length( varipart2 );
    nv1 = length( varipart1 );
    cm = zeros( 64, 3 );
    cm( 1:nv3, 3 ) = varipart3;
    cm( 1 : nv3+nv2, 2 ) = [ ones( nv3, 1 ); varipart2 ];
    cm( :, 1 ) = [ ones( nv3+nv2, 1 ); varipart1( 1 : 64 - (nv3+nv2) ) ];
    
    %%%%%%%%%
    %Settings
    
    formatfigopt.color = 'bwtl';
    formatfigopt.colormap = cm;
    formatfigopt.stylemap = 'bw';

elseif swistyle.colormode.bw

    %bw format
    %Text and line objects are converted to black
    %Solid lines are mapped from colors to linestyle
    %All other objects are converted to gray scales, using the 'gray' colormap
    %Output eps format is modified from colors to gray scale.
    
    formatfigopt.color = 'bw';
    formatfigopt.stylemap = 'bw'; 
    formatfigopt.colormap = 'gray';
    expfigopt.eps_format = '-deps2';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Margin format

if swistyle.margin.tightx
    
    %Tightx format
    formatfigopt.tight = 'x';
    expfigopt.eps_loose = true;

elseif swistyle.margin.tighty

    %Tighty format
    formatfigopt.tight = 'y';
    expfigopt.eps_loose = true;
    
elseif swistyle.margin.tight

    %Tight format
    formatfigopt.tight = 'xy';
    expfigopt.eps_loose = true;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alignment format

if swistyle.alignment.aligned
    
    %Use 'position' as reference so that axes remain aligned while
    %resizing
    formatfigopt.activeposition = 'position';
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%log format

if swistyle.logaxis.log
    
    %lock the ticks
    formatfigopt.lockticks = true;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%heavyfigs format

if swistyle.heavyfigs.heavy
    
    %heavy figure format
    expfigopt.eps_renderer = 'zbuffer';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%innerelems format

if swistyle.innerelems.fixinside
    
    %innerelems format
    formatfigopt.movelegto = 'best';
    formatfigopt.keeptextin = true;
    
end
    


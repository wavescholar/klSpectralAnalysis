%FIG2PUBLIC provides an automatic method to export figures, with optional
%uniform styling. Several input/output folders can be given in a single
%call. The routine works both for windows and linux (from Matlab 6.5 to
%7.x).
%
%FIG2PUBLIC is generally used to:
%   - make uniform figure styling with batch processing. Styling includes
%       font size, marker size, line width, color and figure dimensions
%       (for a cleaner figure inclusion). BW color conversion is possible.
%   - export figures to eps and optionally pdf, jpeg, png and tiff files.
%       By using the eps file as an intermediate file for the other format,
%       the typical huge white blank frame around the figures that Matlab
%       adds when exporting in other formats is not present.
%   - provide a number of pre-defined formats (article, slides, ...)
%
%REQUIREMENTS:
%   - formatfig.m : a reduced and developed version of exportfig.m, a
%       routine developed by Ben Hinkle (@Mathworks).
%   - eps2xxx.m : a figure format converter using ghostscript.
%   - presetformat.m : library of format options.
%   - addtotex.m : a .tex generator for simpler figure inclusion.
%   - ghostscript (minimal version unknown).
%
%CALL:
%
%   FIG2PUBLIC( [ OPT_STRUCT, ] [ 'PARAM_NAME', PARAM_VALUE, ] ... )
%
%   OPT_STRUCT: structure with fields .PARAM_NAME and values .PARAM_VALUE
%
%   NB: for parameters present both in structure and list, the values of
%       the list will be used.
%
%PARAMETERS:
%
%   inp_folder      string or cell array of strings containing the absolute
%                   input folder paths. Only .fig files contained in these
%                   folders will be considered. Use the / for linux and the
%                   \ for windows. Folders starting with '.' are ignored.
%                   Default: current folder.
%
%   dsubset         regular expression used to select a subset of
%                   directories in 'inp_folder'. Useful to use a particular
%                   set of parameters for a subset of folders having a
%                   regular name. E.g. 'dsubset' set to 'rect_*' limits
%                   target folders to those starting with 'rect_'. Default:
%                   []. Compatibility with 'recursive' is limited.
%
%   fsubset         regular expression used to select a subset of files in 
%                   'inp_folder'. Useful to use a particular set of
%                   parameters for a subset of files having a regular name.
%                   E.g. 'fsubset' set to 'rect_*' limits target files to
%                   those starting with 'rect_'. Default: [].
%
%   out_folder      string or cell array of strings containing the absolute
%                   output folder paths. Warning: files in these folders
%                   having the same name as input files will be
%                   OVERWRITTEN. If a cell array is given for inp_folder, a
%                   cell array MUST also be given for out_folder. Default:
%                   current folder. If the given folder does not exist, it
%                   will be created.
%
%   recursive       true (default) / false to operate recursively on the
%                   sub-folders of inp_folder. If there are sub-folders,
%                   the sub-folder name will be used to target the output
%                   sub-folder of out_folder.
%
%   base_config     init a set of parameters with default values. These
%                   values can then be overwritten by setting parameters
%                   explicitly. 'base_config' must be set to a string
%                   combining different default format styles, e.g.:
%                       'article' or 'slides'
%                       'color', 'bwtl' or 'bw' (see PRESETFORMAT)
%                   Examples:
%                       'articlecolor' exports the figures with the
%                       default 'article' sizes, using color output.
%                   Note: optional. Leave empty to use default values of 
%                       FORMATFIG. To get a list of all available default
%                       format styles, type: <help presetformat>
%
%   displayonly     true / false (default). If true, the figure is
%                   displayed and formatting is applied. Then the routine
%                   is paused and waits for the user to press on "Return".
%                   No file is saved or exported.
%
%   exportonly      true / false (default). If true, no figure formatting
%                   is applied before exportation.
%
%   includefig      true / false (default). If true, a formatted copy of
%                   the figure is saved in the output folder as a .fig
%                   file.
%
%   latexcode       true (default) / false. If true, a tex file is added to
%                   the output folder. This file contains latex code
%                   that can be used to include the figures in the latex
%                   document. In 'recursive' mode, the sub-folder pattern
%                   is kept in the figure path to allow direct
%                   copy-pasting. 
%                   Changed to false if 'displayonly' is true.
%                   Changed to true if 'latexcodeonly' is true.
%
%   latexcodeonly   true / false (default). If true, only the .tex files
%                   will be generated. Nothing else will be done. Useful to
%                   recreate the correct .tex files if fig2public has been
%                   called a posteriori on subfolders.
%
%   eps_format      Matlab driver for the eps output. Default: '-depsc2'.
%                   Type <help print> to have a complete list of drivers.
%                   NB: if a non eps driver is given, the eps conversion
%                   will fail. However, if an empty value is given to
%                   'other_formats', this parameter can be used to call
%                   default Matlab export functions. In that case,
%                   .print_ext must be defined accordingly.
%
%   print_ext       Extension of the file returned by the call to print.
%                   Default: 'eps'. Required because print returns a file
%                   with a corrupted extension if the filename contains a
%                   dot. In addition, .eps_format allows any driver, so
%                   that the file format cannot be guessed.
%
%   eps_res         resolution of the eps output. Default: 300 dpi.
%
%   eps_renderer    Renderer used for the eps output. Default: 'painters'.
%                   If white diagonal stripes appear on pcolor plot,
%                   'zbuffer' should be used (Warning: bitmap format). To
%                   reduce file size, use an appropriate 'width' (see
%                   FORMATFIG) and a bitmap renderer ('zbuffer').
%
%   eps_loose       True / False (default). If True, the option '-loose'
%                   is used in the call to PRINT. This is useful to have
%                   an exact match between the figure displayed in 
%                   Matlab and the eps file.
%
%   gs_path         path to the gs executable. Windows users can avoid
%                   setting this parameter by defining the path in the
%                   "environment variable" setting. Linux users should not
%                   have to define it at all.
%                   For Windows users (at least in XP): go to (in Windows)
%                   Configuration pannel -> System -> Advanced ->
%                   Environment Variable -> System Variables -> "Path"
%                   Then click on "Modify" and add the path to the 
%                   gs executable in the path list.                   
%
%   orientation     see EPS2XXX.m. Default: 0.
%
%   other_formats   cell array of export format. Possible values are: 
%                   'pdf', 'jpeg', 'png' and 'tiff'. Default: {'pdf'}
%
%   other_resolution
%                   resolution of the other format output. Default: 300
%                   dpi.
%
%   dointerp        see EPS2XXX.m. Default: true.
%
%   ANY OTHER FORMATFIG parameter-value pair (use: <help formatfig> to get 
%   a list) is allowed as input of FIG2PUBLIC.
%
%
%STANDARD CALLS:
%
%   For standard publications:
%
%       fig2public('inp_folder',{'my_inp_path1','my_inp_path2',...}, ...
%           'out_folder',{'my_out_path1','my_out_path2',...}, ...
%           'base_config','article small bwtl tightx aligned')
%
%   For reports, change 'base_config' to:
%       ..., 'base_config', 'article large bwtl tightx aligned', ...
%
%   For slides:
%       ..., 'base_config', 'slides large color tightx aligned', ...
%
%SPECIAL CALL EXAMPLES
%
%   Give custom values to some FORMATFIG options (NB: last options in the
%   list superseed previous ones):
%
%       fig2public( ..., 'linemode', 'fixed', 'linewidth', 2, ...
%                        'lockticks', true, 'tight', 'x', ... );
%
%       The options above means that all lines will have a width of 2, that
%       the axes ticks will be locked (useful for log axis in general) and
%       that the figure margins will be minimized while keeping the figure
%       width untouched and the height to scale. Note that some parameters
%       are given as strings, others as boolean and others as scalars.
%
%   Advanced usage: absolute control of exported figure width:
%
%       fig2public( ..., 'width', 10, 'eps_loose', true, 'tight', 'x', ...
%                        'activeposition', 'position' );
%
%       Strictly speaking, only 'width' and 'eps_loose' are necessary for
%       absolute control of exported figure width. 'tight' minimizes the
%       margins and 'activeposition' is used to get a better margin
%       control.
%
%
%  Last modification: Jonathan Rossel, 15.12.2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Modif log: (see file)

%       22.07.08        J.Rossel        Added: create output folder if
%                                       inexistant.
%                                       Added: 'recursive' option.
%
%       29.07.08        J.Rossel        Set 'article' configuration:
%                                           -fontsize: 10 (prev. 12)
%                                           -markersize: 6 (prev. 8)
%                                           -linewidth: 0.5 (prev. 1)
%                                       Added 'dointerp' option
%
%       29.10.08        J.Rossel        Added: 'latexcode' option.
%
%       20.04.09        J.Rossel        Modif:
%                                           - No outpuf folder required if
%                                           'displayonly' is true
%                                           - Error if 'includefig' is true
%                                           and same input/output folders.
%                                           - 'latexcode' is false if
%                                           'displayonly' is true
%                                       Added: 'latexcodeonly' option
%
%       10.06.09        J.Rossel        - 'bw' option of 'base_config'
%                                       modified
%                                       - Recursive option modified to avoid
%                                       infinite loop if the output folder
%                                       is a child of the input folder
%
%       16.06.09        J.Rossel        eps_renderer option added
%
%       24.06.09        J.Rossel        - 'fsubset' option added.       
%                                       - structure input allowed.
%       
%       18.12.09        J.Rossel        'print_ext' option added to allow
%                                       correct treatment of files with
%                                       dots in their name.
%
%       11.01.10        J.Rossel        allow call without parameters.
%
%       16.06.10        J.Rossel        'dsubset' option added, 'fsubset'
%                                       limited to files.
%
%       16.08.10        J.Rossel        input formatting for formatfig
%                                       modified (intern).
%
%       12.10.10        J.Rossel        'eps_loose' added
%
%       25.11.10        J.Rossel        Code cleaning and test input param 
%                                       name.
%
%       15.12.10        J.Rossel        Default format styles externalized
%                                       in PRESETFORMAT.m

function fig2public(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Misc inits

%version
verstr = version;
majorver = str2num(verstr(1));

%determine the kind of slash that is used
kindslash = filesep;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand first argument if structure. Store varargin.

%check if first arg is a structure. Expand it if nec.
if nargin > 0,
    if isstruct( varargin{1} ),
        param_values = struct2cell( varargin{1} );
        param_names = fieldnames( varargin{1} );
        param = cell( 1, 2*length( param_names ) );
        param(1:2:end-1) = param_names;
        param(2:2:end) = param_values;
        varargin(1) = []; %remove the structure
        varargin = { param{:}, varargin{:} }; %add the structure as a list (start the param list with it to keep parameter hierarchy)
    end
else
    varargin = {};
end

%save varargin (for recursive calls)
oldvarargin = varargin;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set default parameters for local and export uses

%set default paths and local options
localopt = struct( ...
    'inp_folder', {{pwd}}, ...
    'out_folder', {{pwd}}, ...
    'dsubset', [], ...
    'fsubset', [], ...
    'exportonly', false, ...
    'includefig', false, ...
    'displayonly', false, ...
    'recursive', true, ...
    'latexcode', true, ... %false if displayonly=true, see below
    'latexcodeonly', false, ...
    'latexfigroot', []); %latexfigroot is an internal option only. It is used in recursive mode.
fieldnameslocal = fieldnames(localopt);

%use default configurations as a base for input parameters
%for formatfig and exporting
formatfigopt = formatfig; %empty call to get default structure
fieldnamesformat = fieldnames( formatfigopt );

expfigopt = struct(...
    'eps_format', '-depsc2', ...
    'eps_res', 300, ...
    'eps_renderer', 'painters', ...
    'print_ext', 'eps', ...
    'eps_loose', false, ...
    'gspath', [], ...
    'orientation', 0, ...
    'other_formats', {{'pdf'}}, ...
    'other_resolution', 300, ...
    'dointerp', true );
fieldnamesexp = fieldnames(expfigopt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check that given input parameters are allowed.

%Group all the allowed parameters in a single cell array
allowedpar = [ fieldnameslocal; fieldnamesformat; fieldnamesexp; { 'base_config' } ]; 

inppar = varargin( 1 : 2 : end ); %cell array with input parameter names

%loop on input parameter names and check their value

for ii = 1 : length( inppar )
    
    if ~any( strcmpi( allowedpar, inppar{ ii } ) )
        error( [ 'FIG2PUBLIC: ''', inppar{ ii }, ''' is not an allowed parameter.' ] );
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Apply base configurations

indbc = find( strcmpi( varargin, 'base_config' ) );

if ~isempty( indbc ),
    
    bc = varargin{ indbc( end ) + 1 }; %use 'end' in case the parameter has been given multiple times
    
    [ formatfigopt, expfigopt ] = presetformat( formatfigopt, expfigopt, bc );
    
    %remove the base configuration parameter
    varargin( indbc( end ) : indbc( end ) + 1 ) = [];
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Replace defaults by user-defined values

%fill the local and export structures with the input
%NB: for formatfig, the loop is done in formatfig itself.
ind = [];
for j=1:2:length(varargin)-1 %for figure exportation
    if any(strcmpi(fieldnamesexp,varargin{j})),
        expfigopt = subsasgn(expfigopt,struct('type','.','subs',lower(varargin{j})),varargin{j+1});
        %keep a track of the parameters set in the structure
        ind = [ind,j,j+1];
    end
end
for j=1:2:length(varargin)-1 %for local
    if any(strcmpi(fieldnameslocal,varargin{j})),
        localopt = subsasgn(localopt,struct('type','.','subs',lower(varargin{j})),varargin{j+1});
        %keep a track of the parameters set in the structure
        ind = [ind,j,j+1];
    end
end

%remove the unnecessary parameters from the varargin list
varargin(ind) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check some inputs: target folders, generation of latex code

%init variables
inp_folder = localopt.inp_folder;
out_folder = localopt.out_folder;
latexfigroot = localopt.latexfigroot;

%check out_folder
if localopt.includefig & any( strcmp( inp_folder, out_folder ) ),
    error( 'If includefig==true, the output folder must be different from the input folder');
end

%check: only's. Max 1 at a time.
tmp = [ localopt.displayonly, localopt.exportonly, localopt.latexcodeonly ];
if length( find( tmp ) ) > 1,
    error( 'Max 1 "only" option can be set to true at a time' );
end

%avoid latex code generation when displayonly
if localopt.displayonly, 
    localopt.latexcode = false; 
end

%force latexcode if latexcodeonly
if localopt.latexcodeonly,
    localopt.latexcode = true;
end

%clean input format
if ischar(inp_folder),
    inp_folder = {inp_folder};
end
if ischar(out_folder),
    out_folder = {out_folder};
end
if isempty(latexfigroot),
    latexfigroot = out_folder;
elseif ischar(latexfigroot),
    latexfigroot = {latexfigroot};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop on target folders

for jj = 1:length(inp_folder),
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check folder input format and existence. Prepare latex file if nec.
    
    %clean input format
    if ~strcmp(inp_folder{jj}(end),kindslash),
        inp_folder{jj} = [inp_folder{jj},kindslash];
    end
    if ~strcmp(out_folder{jj}(end),kindslash),
        out_folder{jj} = [out_folder{jj},kindslash];
    end
    if ~strcmp(latexfigroot{jj}(end),kindslash),
        latexfigroot{jj} = [latexfigroot{jj},kindslash];
    end
    
    %check existence of folders
    if exist(inp_folder{jj},'dir') ~= 7, error(['The input folder ',inp_folder{jj},' has not been found.']); end
    if ~localopt.displayonly,
        if exist(out_folder{jj},'dir') ~= 7,
            %create the folder
            eval(['!mkdir ', out_folder{jj}]);
        end
    end
    
    if localopt.latexcode,
        %a file containing code allowing figure inclusion in Latex will be
        %added in each output folder
        texfile = [out_folder{jj},'texinclcode.tex'];
        fidtex = fopen(texfile,'w'); %create or overwrite
    end      
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Loop on fig files
        
    %get filenames (filenames only) in current directory
	filenames = dir( [ inp_folder{jj}, localopt.fsubset ] ); %limit to fsubset
    fileisdir = [filenames(:).isdir]; %the file is actually a folder
	filenames = {filenames( ~fileisdir ).name};
    
    %process files in current directory
    
    for ii = 1:length(filenames),
        
        %check the file
        thisfile = filenames{ii};
        
        if strcmp( thisfile( 1 ), '.' ),
            %file starting with . are ignored
            continue
        end
        
        ext = thisfile(end-2:end);
        if ~strcmp(ext,'fig'),
            %the file is a file but not a fig file
            continue
        end
        
        %remove the extension
        thisfile(end-3:end) = [];
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Figure formatting
        
        if ~localopt.latexcodeonly,
            
            %open the figure
            openfig([inp_folder{jj},thisfile,'.fig']);
            
            if ~localopt.exportonly,
                %format
                formatfig(gcf,formatfigopt,varargin{:});
            end

            if localopt.displayonly,
                pause
                close
                continue
            end

            if localopt.includefig,
                %save a formatted fig file in the output folder
                saveas(gcf,[out_folder{jj},thisfile],'fig');    
            end
        end
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Latex code generation
    
        %add the code to the tex file
        if localopt.latexcode,
            tmp = [out_folder{jj},thisfile]; %complete file name, without extension
            tmp = strrep(tmp,latexfigroot{jj},''); %remove the latexfigroot from the file name
            tmp = strrep(tmp,'\','/'); %latex slash format
            addtotex(fidtex,tmp,thisfile);
            
            if localopt.latexcodeonly,
                continue
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Figure exportation
        
        %export to eps (or other format defined in fields .eps_format and .print_ext)
        printoptions = { expfigopt.eps_format, ...
                         ['-',expfigopt.eps_renderer], ...
                         ['-r',num2str(expfigopt.eps_res)], ...
                         [out_folder{jj},thisfile,'.',expfigopt.print_ext] };
        if expfigopt.eps_loose
            printoptions = { printoptions{ 1 : 3 }, '-loose', printoptions{ 4 : end } };
        end
        print( printoptions{:} ); %call print
        
        %convert to other formats
        if ~isempty(expfigopt.other_formats) & strcmpi( expfigopt.print_ext, 'eps' ),
            eps2xxx([out_folder{jj},thisfile,'.eps'],expfigopt.other_formats, ...
                expfigopt.gspath,expfigopt.orientation,expfigopt.other_resolution,expfigopt.dointerp);
        end
        
        close
        
    end %end loop on files
    
    if localopt.latexcode,
        %close the file
        fclose(fidtex);
    end      
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recursive mode
    
    if localopt.recursive,
        
        %process directories in current directory

        %get dirnames (directories only)
        dirnames = dir( [ inp_folder{jj}, localopt.dsubset ] ); %limit to dsubset
        fileisdir = [dirnames(:).isdir]; %the file is actually a folder
        dirnames = {dirnames( fileisdir ).name};

        for ii = 1:length(dirnames),

            %check the directory
            thisdir = dirnames{ii};

            if strcmp( thisdir( 1 ), '.' ),
                %. .. or folders starting with . are ignored
                continue
            end

            %call fig2public recursively
            newvarargin = oldvarargin; %don't modify the old one

            newinp = [inp_folder{jj},thisdir]; %complete path to input subfolder
            newout = [out_folder{jj},thisdir]; %complete path to output subfolder. Will be created if necessary

            %check that the new input is not the current output folder
            %in which case an infinite loop would happen. If it is,
            %skip it
            if strcmp( [ newinp, kindslash ], out_folder{ jj } ),
                continue
            end

            tmp = find(strcmpi(newvarargin,'inp_folder'));
            if isempty( tmp ),
                %default path was used, add the option at the end of
                %the list
                newvarargin( end+1 : end+2 ) = { 'inp_folder', newinp };
            else
                newvarargin{tmp+1} = newinp; %update folders
            end

            tmp = find(strcmpi(newvarargin,'out_folder'));
            if isempty( tmp ),
                %default path was used, add the option at the end of
                %the list
                newvarargin( end+1 : end+2 ) = { 'out_folder', newout };
            else
                newvarargin{tmp+1} = newout; %update folders
            end

            tmp = find(strcmpi(oldvarargin,'latexfigroot'));
            if isempty(tmp),
                %this is the external call of fig2public
                %set latexfigroot to keep the sub-folders in the fig
                %paths
                newvarargin(end+1:end+2) = {'latexfigroot',out_folder{jj}};
            end               

            %call fig2public with all the same parameters, apart from
            %input/output folders
            fig2public(newvarargin{:});
            
        end %end loop on subdir

    end %end recursive mode
   
end %end loop on target folders

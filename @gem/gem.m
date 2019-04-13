% GEM : an implementation of GMP Eigen Matrices for MATLAB
%
% gem is a MATLAB class wrapper to the GmpEigenMatrix C++ class
%
% It can be used as follows:
%   a = gem(3.5)
%   b = gem(2.5+1i)
%   c = a+b

% The precision of gem objects is set at construction. Either explicitely,
% or by using the global function gemWorkingPrecision. The default value
% is 50 digits.

% For numbers specified as strings, the precision is always set to be at least
% large enough to translate all significant digits in the provided string.

% The precision of displayed number can be adjusted through the
% gemDisplayPrecision function.

% Note that the result of some operation may contain more digits than can
% be truly garanteed... (this is the standard behavior of mpfrc++)
% Note also that even though a number may have a large precision, if it
% finishes with a string of zeros, these won't be printed (again this is
% the standard behavior of mpfrc++)

% Note: A gem object is a handle-type object. This garantees that
% matlab manages its memory cleanly (calling the destructor whenever
% needed). This also means that "y=x" performes a soft copy of the object:
% modifying x afterwards _also_ modifies y. For this reason, any procedure
% that modifies a gem object should always work on a copy! (or at least
% produce the result in a new object)

% Usage : Best usage is to always start by setting the working precision to
%   whatever precision we want the library to work. This default precision
%   enters at several stages. For instance, when computing log10(x), even
%   if x has 100 digits, the result is only precise up to the precision
%   specified throught gemWorkingPrecision.

% Constants : The following special calls return mathematical constants
%   - gem('log2')    % log2 = 0.693...
%   - gem('pi')      % Pi = 3.14...
%   - gem('e')       % e = 2.718...
%   - gem('euler')   % Euler constant = 0.577...
%   - gem('catalan') % The catalan number = 0.915965594...

classdef gem < handle
    properties (SetAccess = private, Hidden = true)
        objectIdentifier; % The identifier of the underlying C++ class instance (i.e. the integer version of the pointer to the handle_class instance containing a pointer to the object of interest)
    end

    methods
        %% Constructor
        % This function can be called to
        %  - Create a new C++ class empty instance
        %  - Create a new copy of a C++ instance (i.e. a new C++ object
        %    with same data)
        %  - Create a new C++ class instance from a matlab object
        %  - Create a gem object to encapsulate a reference to an
        %    already existing C++ class instance. For this the constructor
        %    should be called with two arguments : 'encapsulate', followed
        %    by the integer reference to the existing object to be
        %    encapsulated.
        %  - Create a C++ class instance for any of the mathematical
        %    constants above
        function this = gem(varargin)

            this.checkForBinaries;
            
            % The following variable is used in some instances to tell the
            % constructor to use the precision requested by the user rather
            % than the one it thinks is better for the considered number
            persistent forceDefaultPrecision;
            if nargin == 0
                % Without further argument we construct a new empty instance
                this.objectIdentifier = gem_mex('new');
            elseif nargin == 1
                % Then we check that the argument is of the same class
                % before copying it
                if isequal(class(varargin{1}), 'gem')
                    % If it is an object of type 'gem', we interpret
                    % this call to a constructor as a call for copying the object
                    % into a new one
                    this.objectIdentifier = gem_mex('new', varargin{1}.objectIdentifier);
                elseif isequal(class(varargin{1}), 'sgem')
                    % The we create a dense version of the provided sparse
                    % sgem object
                    this.objectIdentifier = gem_mex('full', varargin{1}.objectIdentifier);
                elseif isnumeric(varargin{1}) || islogical(varargin{1})
                    % Then we interpret this call as a call for the library to
                    % create an instance of this class from some numerical
                    % matlab data (e.g. a numerical number). We thus transfer
                    % this data to the C++ library, and let it decide how to
                    % create a class instance from it. Upon completion, we
                    % return the newly created object containing this data.
%                    if issparse(varargin{1})
%                        warning('Creating a dense gem object from a sparse matrix. Use ''sgem'' or ''gemify'' to create a sparse gem object.');
%                    end
                    if isa(varargin{1}, 'uint8') || isa(varargin{1}, 'uint16') || isa(varargin{1}, 'uint32') || isa(varargin{1}, 'uint64') ...
                        || isa(varargin{1}, 'int8') || isa(varargin{1}, 'int16') || isa(varargin{1}, 'int32') || isa(varargin{1}, 'int64')
                        % For potentially large integers, translating into
                        % strings guarantees that we don't forget any digit
                        tmp = cell(size(varargin{1}));
                        for i = 1:numel(varargin{1})
                            tmp{i} = num2str(varargin{1}(i));
                        end
                        this = gem(tmp);
                    else
                        if ~isa(varargin{1}, 'double')
                            varargin{1} = double(varargin{1});
                        end
                        this.objectIdentifier = gem_mex('newFromMatlab', full(varargin{1}), this.getWorkingPrecision);
                    end
                elseif ischar(varargin{1})
                    % We embed the string into a cell array so that the c++
                    % interface interprets it correctly.
                    
                    % First, we check if the string denotes a mathematical
                    % constant
                    switch lower(varargin{1})
                        case 'log2'
                            precision = this.getWorkingPrecision; % This is the precision at which this constant will be computed (with this line, we make sure the precision has been set in the C++ library)
                            this.objectIdentifier = gem_mex('const_log2');
                            return;
                        case 'pi'
                            precision = this.getWorkingPrecision;
                            this.objectIdentifier = gem_mex('const_pi');
                            return;
                        case 'e'
                            this = exp(gem(1));
                            return;
                        case 'euler'
                            precision = this.getWorkingPrecision;
                            this.objectIdentifier = gem_mex('const_euler');
                            return;
                        case 'catalan'
                            precision = this.getWorkingPrecision;
                            this.objectIdentifier = gem_mex('const_catalan');
                            return;
                    end

                    % If we are here, the string is not a mathematical
                    % constant, so we check that the string only contains a
                    % real part
                    if ~this.isRealString(varargin{1})
                        [rPart iPart] = this.separateRIString(varargin{1});
                        if ~this.isRealString(rPart) || ~this.isRealString(iPart) % To prevent an infinite loop
                            error('The complex string was not correctly separated.');
                        end

                        % We create a gem object for each part, and
                        % recombine them into a single one
                        rgem = gem(rPart);
                        igem = gem(iPart);

                        % We copy the sum into the current object
                        this = rgem + 1i*igem;
                        return;
                    end

                    % We check that the string is of an acceptable form
                    if ~this.isValidString(varargin{1})
                        error('The format of this string is not supported');
                    end

                    % Now we set the precision and construct the c++ object
                    if isempty(forceDefaultPrecision) || (forceDefaultPrecision ~= 1)
                        precision = max(this.getWorkingPrecision, this.nbDigitsFromString(varargin{1}));
                    else
                        % We force the precision to be the one requested by
                        % the user
                        precision = this.getWorkingPrecision;
                    end
                    this.objectIdentifier = gem_mex('newFromMatlab', {varargin{1}}, precision);
                elseif iscell(varargin{1})
                    % First, we check that the cell only contains a real
                    % part
                    if ~this.isRealCell(varargin{1})
                        [rPart iPart] = this.separateRICell(varargin{1});
                        if ~this.isRealCell(rPart) || ~this.isRealCell(iPart) % To prevent an infinite loop
                            error('The complex string was not correctly separated.');
                        end

                        % We create a gem object for each part, and
                        % recombine them into a single one
                        rgem = gem(rPart);
                        igem = gem(iPart);

                        % We copy the sum into the current object
                        this = rgem + 1i*igem;
                        return;
                    end

                    % We check that the cell contains numbers in an
                    % acceptable form
                    if ~this.isValidCell(varargin{1})
                        error('The format of this cell is not supported');
                    end

                    % We verify that the precision high enough to make sure all
                    % numbers are well translated
                    precision = this.getWorkingPrecision;
                    for i = 1:numel(varargin{1})
                        if isnumeric(varargin{1}{i})
                            precision = max(precision, 15);
                        elseif ischar(varargin{1}{i})
                            precision = max(precision, this.nbDigitsFromString(varargin{1}{i}));
                        else
                            error('Only doubles and strings are supported in a cell array.');
                        end
                    end
                    % Now we can construct the gem object from these
                    % strings
                    this.objectIdentifier = gem_mex('newFromMatlab', varargin{1}, precision);
                else
                    error('Wrong instruction upon creation of a new gem object.');
                end
            elseif nargin == 2
                if isequal(lower(varargin{1}),'encapsulate') && isequal(class(varargin{2}), 'uint64')
                    % If the second argument is of type 'uint64', then we interpret
                    % it as pointing to an existing instance of a C++ class, so we
                    % encapsulate it into the current gem instance.

                    % But since this should be a private constructor, we
                    % first we check that the caller is the current file
                    % (i.e. gem.m)
                    [ST I] = dbstack('-completenames');
                    if (length(ST) < 2) || (isempty(strfind(ST(2).file,'/@gem/')) && isempty(strfind(ST(2).file,'\@gem\')) && isempty(strfind(ST(2).file,'/@sgem/'))&& isempty(strfind(ST(2).file,'\@sgem\')) && isempty(strfind(ST(2).name,'gemRand')))
                        error('Only gem.m and sgem.m is allowed to encapsulate an integer into a new gem object.');
                    end

                    % This creates a matlab object which points to the C++
                    % object referred to by this number.
                    this.objectIdentifier = varargin{2};

                    % We check that the created object is valid (that it really points
                    % to a proper c++ object). Otherwise we produce an error.
                    if ~(this.checkIdentifierValidity)
                        this.objectIdentifier = 0; % We reset the handle...
                        error('Invalid reference given upon construction of a new gem object.');
                    end
                elseif isnumeric(varargin{2}) && isequal(size(varargin{2}), [1 1])
                    % Then we interpret this call as a call for the library to
                    % create an instance of this class from some numerical
                    % matlab data (e.g. a numerical number), together with
                    % a specific precision. We thus set the precision
                    % accordingly and call the default constructor.
                    
                    if (varargin{2} < 1)
                        error('At least one digit of precision is required');
                    end
                    
                    % Save default precision
                    previousPrecision = this.getWorkingPrecision;
                    
                    % Assigned desired precision
                    this.setWorkingPrecision(varargin{2});

                    % We tell the library that the precision must be
                    % enforced
                    forceDefaultPrecision = 1;
                    
                    % Create object
                    this = gem(varargin{1});
                    
                    % We don't need precision enforcement anymore
                    forceDefaultPrecision = 0;
                    
                    % We restore the default precision
                    this.setWorkingPrecision(previousPrecision);
                else
                    error('Wrong instruction upon creation of a new gem object.');
                end
            else
                error('Too many parameters in the creation of a new gem object.');
            end
        end

        % Destructor - Destroy the C++ class instance
        function delete(this)
            % Latest versions of matlab might try to delete object which
            % were not fully constructed yet...
            if ~isempty(this.objectIdentifier)
                gem_mex('delete', this.objectIdentifier);
            end
        end    
    end

    methods(Access = protected)
        % We offer an alternative to the default copy operation 'y=x;' to 
        % performs a deep copy rather than a shallow one
        function cp = copy(el)
            cp = gem(el);
        end
    end
    
    % Since the load function does not depend on a class instance (but creates one),
    % it needs to be a static method, and so we need to define it here...
    methods (Static)
        function result = loadobj(structure)
            fprintf('loading...');
            % The result should be an instance of a gem object with the data contained in the provided structure
            if structure.dataVersion > 1
                error('The object was saved with a newer version of the library. Please upgrade the library to load it again.');
            else
                % loading instruction...
                newObjectIdentifier = gem_mex('loadobj', structure);

                % since this is a static function, we still need to
                % encapsulate the object created into a matlab instance of
                % this class
                result = gem('encapsulate', newObjectIdentifier);
            end
            disp('done');
        end
    end


    %% Here come the methods which define the class-wide variables
    % 'workingPrecision' and 'displayPrecision'
    methods (Static)%, Access = private)
        % This method implements a static variable for the whole class that
        % defines the default precision of new numbers
        function value = workingPrecision(newValue)
            persistent precision;
            if isempty(precision)
                precision = 50;
                % We call the mex file to set the default working precision
                gem_mex('setWorkingPrecision', precision);
            end
            if nargin >= 1
                if precision < 1
                    error('The precision need to be larger or equal to 1    ');
                end
                precision = double(newValue);
                % We call the mex interface to make this the default
                % working precision
                gem_mex('setWorkingPrecision', precision);
            end
            value = precision;
        end
        % This method implements a static variable for the whole class that
        % defines the default precision at which numbers are displayed
        function value = displayPrecision(newValue)
            persistent precision;
            if isempty(precision)
                precision = 20;
            end
            if nargin >= 1
                if newValue >= 1
                    precision = double(newValue);
                else
                    % This means plotting as many nonzero digits as
                    % available
                    precision = -1;
                end
            end
            value = precision;
        end
    end
    
    
    %% The following methods are needed in the class, they are implemented as
    % methods rather than functions, because octave doesn't support declaration 
    % of functions into a class file
    methods (Static, Access = private)
        
        % This function verifies that a string doesn't contain any imaginary part
        function bool = isRealString(str)
            bool = 1 - (sum(str == 'i') > 0);
        end


        % This function verifies that a cell array doesn't contain any imaginary
        % part (but the cell array could contain either doubles, or strings)
        function bool = isRealCell(x)
            for i = 1:numel(x)
                if ischar(x{i})
                    if sum(x{i} == 'i') > 0
                        bool = 0;
                        return;
                    end
                elseif isnumeric(x{i}) && (numel(x{i}) == 1)
                    if ~isreal(x{i})
                        bool = 0;
                        return;
                    end
                else
                    error('Cell arrays can only contain doubles or strings.');
                end
            end
            bool = 1;
        end


        % This function separates the real and imaginary parts in a string
        % Here, we expect the imaginary part to appear after the real part.
        function [rPart iPart] = separateRIString(str)
            if this.isRealString(str)
                rPart = str;
                iPart = '';
            else
                % First, we remove spaces and force lower case
                str = lower(str(find(str~=' ')));

                % Now, let us identify the imaginary part...
                % maybe it is written in scientific notation...
                imagStart = regexp(str,'(+|-|)[0-9]+(\.|)[0-9]*e(+|-|)[0-9]+i');

                if isempty(imagStart)
                    % If not, we check if the imaginary part is written as a simple
                    % number
                    imagStart = regexp(str,'(+|-|)[0-9]+(\.|)[0-9]*i');
                end

                if isempty(imagStart)
                    % There is no imaginary part
                    rPart = str;
                    iPart = '';
                else
                    rPart = str(1:imagStart-1);
                    iPart = str(imagStart:end);
                    % We remove the 'i' from the imaginary part
                    iPart = iPart(find(iPart~='i'));
                end
            end
        end


        % This function separates the real and imaginary parts in a cell array
        function [rPart iPart] = separateRICell(x)
            rPart = cell(size(x));
            iPart = cell(size(x));
            for i = 1:numel(x)
                if ischar(x{i})
                    [rPart{i} iPart{i}] = this.separateRIString(x{i});
                elseif isnumeric(x{i}) && (numel(x{i}) == 1)
                    rPart{i} = real(x{i});
                    iPart{i} = imag(x{i});
                else
                    error('Cell arrays can only contain doubles or strings.');
                end
            end
        end


        % This function verifies that the string str defines a valid number
        % NOTE : We should also make sure that the spaces are located at the right
        %        place, but for the moment, we only check things up to the spaces.
        %        If mpfr doesn't like some spaces, we could also remove them here.
        function bool = isValidString(str)
            % We only support the following characters : 0-9,+,-,e,E,.
            if ~isempty(regexp(str,'[^ 0-9eE+-\.]'))
                bool = 0;
                return;
            end

            % We remove spaces, and lower the e
            str = lower(str);
            str = str(find(str~=' '));

            % Now we check if the string contains a number in scientific notation
            str2 = regexprep(str,'(+|-|)[0-9]+(\.|)[0-9]*e(+|-|)[0-9]+','');
            if length(str2) ~= length(str)
                if length(str2) == 0
                    % The string denotes a number in scientific notation
                    bool = 1;
                    return;
                else
                    % Something remains in the string after we removed the number
                    % in scientific notation. This is not an acceptable string
                    bool = 0;
                    return;
                end
            end

            % nothing was found, so we check if there is just a standard number
            str3 = regexprep(str,'(+|-|)[0-9]+(\.|)[0-9]*','');
            bool = (length(str3) == 0);
        end


        % This function verifies that a cell is in a suitable form to be passed to
        % the c++ library
        function bool = isValidCell(x)
            for i = 1:numel(x)
                if ischar(x{i})
                    if ~this.isValidString(x{i})
                        bool = 0;
                        return;
                    end
                elseif ~(isnumeric(x{i}) && (numel(x{i}) == 1))
                    error('Cell arrays can only contain doubles or strings.');
                end
            end
            bool = 1;
        end


        % This function counts how many digits significant are present in the
        % number described in the provided string
        function n = nbDigitsFromString(str)

            % We don't want to consider a possible exponent
            exponent = find(lower(str) == 'e', 1);
            if isempty(exponent)
                exponent = length(str)+1;
            end

            % We count spaces
            spaces = find(str(1:exponent-1) == ' ');
            % dots
            dots = find(str(1:exponent-1) == '.');
            % pluses
            pluses = find(str(1:exponent-1) == '+');
            % minuses
            minuses = find(str(1:exponent-1) == '-');

            % And return the number of digits
            n = exponent-1 - length(spaces) - length(dots) - length(pluses) - length(minuses);
        end

        % This function checks whether the library binaries are in the path.
        function value = checkForBinaries()
            persistent binariesOk;
            if isempty(binariesOk) %|| (binariesOk ~= 1)
                % In the source file version of this library, we start by
                % checking whether the c++ library was compiled. If not, we 
                % suggest to download the binaries.
                tmp = mfilename('fullpath');
                if (exist([tmp(1:end-8), 'gem_mex.', mexext], 'file') ~= 3)
                    warning('The library binaries were not found. You may wish to download them online at https://www.github.com/jdbancal/gem/releases .');
                    binariesOk = 0;
                else
                    binariesOk = 1;
                end
            end
            value = binariesOk;
        end
        
    end

end


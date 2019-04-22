% sGEM : an implementation of sparse GMP Eigen Mantrices for MATLAB
%
% sgem is a MATLAB class wrapper to the sparse GmpEigenMatrix C++ class
%
% It can be used as follows:
%   a = sgem(3.5)
%   b = sgem(2.5+1i)
%   c = a+b
%
%   sgem(pi)        % uses matlab 15 digit pi, embedded with default precision
%   sgem(pi, 2)     % gives 3.1 (only two digits taken from matlab's 15 digit pi)
%   sgem(gem('pi')) % uses gem constant pi with full default precision
%
%   S = rand(3)
%   [i,j,s] = find(S);
%   [m,n] = size(S);
%   sgem(i,j,s,m,n)            % Constructs a sparse matrix from nonzero elements
%   sgem(i,j,s)                % same if m=max(i) and n=max(j);
%   sgem(i,j,s,m,n, precision) % Same but with given precision
%   sgem(i,j,s,precision)      % same if m=max(i) and n=max(j);

% The precision of gem objects is set at construction. Either explicitely,
% or by using the global function gemSetWorkingPrecision. The default value
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

% Note: An sgem object is a handle-type object. This garantees that
% matlab manages its memory cleanly (calling the destructor whenever
% needed). This also means that "y=x" performes a soft copy of the object:
% modifying x afterwards _also_ modifies y. For this reason, any procedure
% that modifies an sgem object should always work on a copy! (or at least
% produce the result in a new object)

% Usage : Best usage is to always start by setting the working precision to
%   whatever precision we want the library to work. This default precision
%   enters at several stages. For instance, when computing log10(x), even
%   if x has 100 digits, the result is only precise up to the precision
%   specified throught gemWorkingPrecision.

classdef sgem < handle
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
        function this = sgem(varargin)

            if nargin == 0
                % Without further argument we construct a new empty instance
                this.objectIdentifier = sgem_mex('new');
            elseif nargin == 1
                % Then we check that the argument is of the same class
                % before copying it
                if isequal(class(varargin{1}), 'sgem')
                    % If it is an object of type 'sgem', we interpret 
                    % this call to a constructor as a call for copying the object
                    % into a new one
                    this.objectIdentifier = sgem_mex('new', varargin{1}.objectIdentifier);
                elseif isequal(class(varargin{1}), 'gem')
                    % The we create a sparse version of the provided dense 
                    % gem object with default threshold
                    threshold = min(size(varargin{1}))*10^gem(-this.getWorkingPrecision);
                    this.objectIdentifier = sgem_mex('sparse', varargin{1}.objectIdentifier, threshold.objectIdentifier);
                elseif isnumeric(varargin{1}) || islogical(varargin{1})
                    % Then we interpret this call as a call for the library to
                    % create an instance of this class from some numerical
                    % matlab data (e.g. a numerical number). We thus transfer
                    % this data to the C++ library. Since we are creating a
                    % sparse object, we transfer the information in an
                    % optimized way
                    [i j s] = find(varargin{1});
                    if size(i,1) < size(i,2) % for line vectors, matlab will return i,j,s with a different formatting... so we correct that 
                        i = i.';
                        j = j.';
                        s = s.';
                    end
                    if numel(s) == 0 % Size 1x0 is problematic for the c++ library, so we remove this possibility
                        i = [];
                        j = [];
                        s = [];
                    end
                    [m n] = size(varargin{1});
                    if isequal(class(s),'gem') % in principle this case should not occur
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, s.objectIdentifier, m, n, this.getWorkingPrecision);
                    elseif isa(varargin{1}, 'uint8') || isa(varargin{1}, 'uint16') || isa(varargin{1}, 'uint32') || isa(varargin{1}, 'uint64') ...
                        || isa(varargin{1}, 'int8') || isa(varargin{1}, 'int16') || isa(varargin{1}, 'int32') || isa(varargin{1}, 'int64')
                        % For integers, we set the values manually to make
                        % sure we don't loose any precision
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, zeros(size(s)), m, n, this.getWorkingPrecision);
                        s_values = gem(s);
                        this(:) = s_values;
                    else
                        if ~isa(s, 'double')
                            s = double(s);
                        end
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, s, m, n, this.getWorkingPrecision);
                    end
                elseif ischar(varargin{1}) || iscell(varargin{1})
                    % We first create a dense gem object, then return its
                    % sparse version
                    this = sgem(gem(varargin{1}));
                end
            elseif nargin == 2
                if isequal(lower(varargin{1}),'encapsulate') && isequal(class(varargin{2}), 'uint64')
                    % If the second argument is of type 'uint64', then we interpret
                    % it as pointing to an existing instance of a C++ class, so we
                    % encapsulate it into the current sgem instance.

                    % But since this should be a private constructor, we
                    % first we check that the caller is the current file
                    % (i.e. gem.m)
                    [ST I] = dbstack('-completenames');
                    if (length(ST) < 2) || (isempty(strfind(ST(2).file,'/@gem/')) && isempty(strfind(ST(2).file,'\@gem\')) && isempty(strfind(ST(2).file,'/@sgem/')) && isempty(strfind(ST(2).file,'\@sgem\')))
                        error('Only sgem.m is allowed to encapsulate an integer into a new sgem object.');
                    end

                    % This creates a matlab object which points to the C++
                    % object referred to by this number.
                    this.objectIdentifier = varargin{2};

                    % We check that the created object is valid (that it really points
                    % to a proper c++ object). Otherwise we produce an error.
                    if ~(this.checkIdentifierValidity)
                        this.objectIdentifier = 0; % We reset the handle...
                        error('Invalid reference given upon construction of a new sgem object.');
                    end
                elseif isequal(class(varargin{1}), 'gem')
                    % The we are asked to create a sparse version of
                    % a dense gem object, with custom zero threshold
                    threshold = varargin{2};
                    if ~isequal(class(threshold), 'gem')
                        threshold = gem(threshold);
                    end
                    this.objectIdentifier = sgem_mex('sparsify', varargin{1}.objectIdentifier, threshold);
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

                    % Create object
                    this = sgem(varargin{1});
                    
                    % We restore the default precision
                    this.setWorkingPrecision(previousPrecision);
                else
                    error('Wrong instruction upon creation of a new sgem object.');
                end
            elseif nargin == 3
                if isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3}) ...
                    && numel(varargin{1}) == numel(varargin{2}) && ((numel(varargin{1}) == numel(varargin{3})) || (numel(varargin{3}) == 1)) ...
                    % Then we interpret this call as a call for the library to
                    % create an instance of this class from some numerical
                    % matlab data with given row and column indices. We thus transfer
                    % this data to the C++ library. Since we are creating a
                    % sparse object, we transfer the information in an
                    % optimized way
                    i = double(varargin{1}); % We make sure i and j are double arrays
                    j = double(varargin{2});
                    s = full(varargin{3}); % We make sure s is not sparse
                    if numel(i) ~= numel(s)
                        s = ones(1,numel(i))*s;
                    end
                    if size(i,1) < size(i,2) % We make sure all vectors are in column form
                        i = i.';
                    end
                    if size(j,1) < size(j,2) % We make sure all vectors are in column form
                        j = j.';
                    end
                    if size(s,1) < size(s,2) % We make sure all vectors are in column form
                        s = s.';
                    end
                    % We verify that s contains no zero element
                    sel = find(s);
                    if ~isequal(sel',1:numel(s))
                        i = i(sel);
                        j = j(sel);
                        s = s(sel);
                    end
                    if numel(s) == 0 % Size 1x0 is problematic for the c++ library, so we remove this possibility
                        i = [];
                        j = [];
                        s = [];
                    end
                    m = max(i);
                    n = max(j);
                    if isequal(class(s),'gem') % in principle this case should not occur
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, s.objectIdentifier, m, n, this.getWorkingPrecision);
                    elseif isa(varargin{1}, 'uint8') || isa(varargin{1}, 'uint16') || isa(varargin{1}, 'uint32') || isa(varargin{1}, 'uint64') ...
                        || isa(varargin{1}, 'int8') || isa(varargin{1}, 'int16') || isa(varargin{1}, 'int32') || isa(varargin{1}, 'int64')
                        % For integers, we set the values manually to make
                        % sure we don't loose any precision
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, zeros(size(s)), m, n, this.getWorkingPrecision);
                        s_values = gem(s);
                        this(:) = s_values;
                    else
                        if ~isa(s, 'double')
                            s = double(s);
                        end
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, s, m, n, this.getWorkingPrecision);
                    end
                else
                    error('Wrong instruction upon creation of a new sgem object.');
                end
            elseif nargin == 4
                if isnumeric(varargin{4}) && isequal(size(varargin{4}), [1 1])
                    % Then we interpret this call as a call for the library to
                    % create an instance of this class from some numerical
                    % matlab data with given row and column indices, together with
                    % a specific precision. We thus set the precision
                    % accordingly and call the default constructor.
                    
                    if (varargin{4} < 1)
                        error('At least one digit of precision is required');
                    end
                    
                    % Save default precision
                    previousPrecision = this.getWorkingPrecision;
                    
                    % Assigned desired precision
                    this.setWorkingPrecision(varargin{4});

                    % Create object
                    this = sgem(varargin{1},varargin{2},varargin{3});
                    
                    % We restore the default precision
                    this.setWorkingPrecision(previousPrecision);
                else
                    error('Wrong instruction upon creation of a new sgem object.');
                end
            elseif nargin == 5
                if isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3}) ...
                    && numel(varargin{1}) == numel(varargin{2}) && ((numel(varargin{1}) == numel(varargin{3})) || (numel(varargin{3}) == 1)) ...
                    && isnumeric(varargin{4}) && isnumeric(varargin{5}) ...
                    && isequal(size(varargin{4}), [1 1]) && isequal(size(varargin{5}), [1 1])
                    % Then we interpret this call as a call for the library to
                    % create an instance of this class from some numerical
                    % matlab data with given row, column indices, and size. We thus transfer
                    % this data to the C++ library. Since we are creating a
                    % sparse object, we transfer the information in an
                    % optimized way
                    i = double(varargin{1}); % We make sure i and j are double arrays
                    j = double(varargin{2});
                    s = full(varargin{3}); % We make sure s is not sparse
                    if numel(i) ~= numel(s)
                        s = ones(1,numel(i))*s;
                    end
                    if size(i,1) < size(i,2) % We make sure all vectors are in column form
                        i = i.';
                    end
                    if size(j,1) < size(j,2) % We make sure all vectors are in column form
                        j = j.';
                    end
                    if size(s,1) < size(s,2) % We make sure all vectors are in column form
                        s = s.';
                    end
                    % We verify that s contains no zero element
                    sel = find(s);
                    if ~isequal(sel,1:numel(s))
                        i = i(sel);
                        j = j(sel);
                        s = s(sel);
                    end
                    if numel(s) == 0 % Size 1x0 is problematic for the c++ library, so we remove this possibility
                        i = [];
                        j = [];
                        s = [];
                    end
                    m = double(varargin{4}); % We make sure m and n are doubles
                    n = double(varargin{5});
                    if (numel(s) > 0) && ((min(i) < 1) || (max(i) > m) || (min(j) < 1) || (max(j) > n))
                        error('Incompatible sizes in construction of an sgem object.');
                    end
                    if isequal(class(s),'gem')
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, s.objectIdentifier, m, n, this.getWorkingPrecision);
                    else
                        this.objectIdentifier = sgem_mex('newFromMatlab', i, j, double(s), m, n, this.getWorkingPrecision);
                    end
                else
                    error('Wrong instruction upon creation of a new sgem object.');
                end
            elseif nargin == 6
                if isnumeric(varargin{6}) && isequal(size(varargin{6}), [1 1])
                    % Then we interpret this call as a call for the library to
                    % create an instance of this class from some numerical
                    % matlab data with given row, column indices, and size, together with
                    % a specific precision. We thus set the precision
                    % accordingly and call the default constructor.
                    
                    if (varargin{6} < 1)
                        error('At least one digit of precision is required');
                    end
                    
                    % Save default precision
                    previousPrecision = this.getWorkingPrecision;
                    
                    % Assigned desired precision
                    this.setWorkingPrecision(varargin{6});

                    % Create object
                    this = sgem(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
                    
                    % We restore the default precision
                    this.setWorkingPrecision(previousPrecision);
                else
                    error('Wrong instruction upon creation of a new sgem object.');
                end
            else
                error('Too many parameters in the creation of a new sgem object.');
            end
        end

        % Destructor - Destroy the C++ class instance
        function delete(this)
            % Latest versions of matlab might try to delete object which
            % were not fully constructed yet...
            if ~isempty(this.objectIdentifier)
                sgem_mex('delete', this.objectIdentifier);
            end
        end    
    end

    methods(Access = protected)
        % We offer an alternative to the default copy operation 'y=x;' to 
        % performs a deep copy rather than a shallow one
        function cp = copy(el)
            cp = sgem(el);
        end
    end
    
    % Since the load function does not depend on a class instance (but creates one),
    % it needs to be a static method, and so we need to define it here...
    methods (Static)
        function result = loadobj(structure)
            fprintf('loading...');
            % The result should be an instance of an sgem object with the data contained in the provided structure
            if structure.dataVersion > 1
                error('The object was saved with a newer version of the library. Please upgrade the library to load it again.');
            else
                % loading instruction...
                newObjectIdentifier = sgem_mex('loadobj', structure);

                % since this is a static function, we still need to
                % encapsulate the object created into a matlab instance of
                % this class
                result = sgem('encapsulate', newObjectIdentifier);
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
                sgem_mex('setWorkingPrecision', precision);
            end
            if nargin >= 1
                if precision < 1
                    error('The precision need to be larger or equal to 1    ');
                end
                precision = double(newValue);
                % We call the mex interface to make this the default
                % working precision
                sgem_mex('setWorkingPrecision', precision);
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
        % This method implements a static variable for the whole class that
        % defines whether we want the library to behave like matlab
        % regarding sparse operations.
        function value = sparseLikeMatlab(newValue)
            persistent likeMatlab;
            if isempty(likeMatlab)
                likeMatlab = 0;
            end
            if nargin >= 1
                if (newValue == 1) || (newValue == 0)
                    likeMatlab = double(newValue);
                end
            end
            value = likeMatlab;
        end
    end
end


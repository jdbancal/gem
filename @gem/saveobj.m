% Save the c++ object with matlab
function structure = saveobj(this)
    % This function should return a matlab structure that fully characterizes the current object
    structure = gem_mex('saveobj', this.objectIdentifier);
    % We add the format version
    structure.dataVersion = 1.0;
end

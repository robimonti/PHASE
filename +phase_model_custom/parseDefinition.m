function definition = parseDefinition(config)
%PARSEDEFINITION Validate a reproducible custom deterministic model.
definition = struct('expression',strtrim(config.customExpression), ...
    'functionFile',strtrim(config.customFunctionFile));
try
    parameters = jsondecode(config.customParameters);
catch ME
    error('PHASE:customParameters','Invalid parameter JSON: %s',ME.message);
end
if isempty(parameters), parameters = struct('name',{},'value',{},'estimate',{}); end
if ~isstruct(parameters)
    error('PHASE:customParameters','Parameters must be a JSON array of objects.');
end
reserved = {'x','y','s','t','pi','exp','expm1','log','log1p','sin','cos', ...
    'tan','tanh','sqrt','abs','erf','min','max'};
for k = 1:numel(parameters)
    p = parameters(k);
    if ~all(isfield(p,{'name','value','estimate'})) || ...
            ~isvarname(p.name) || strcmp(p.name,'phaseResidualMean') || ismember(p.name,reserved) || ...
            ~isnumeric(p.value) || ~isscalar(p.value) || ~isfinite(p.value) || ...
            ~islogical(p.estimate) || ~isscalar(p.estimate)
        error('PHASE:customParameters', ...
            'Each parameter needs a unique name, a finite value and estimate: true/false.');
    end
    if ~isfield(p,'lower') || isempty(p.lower), parameters(k).lower = -Inf; end
    if ~isfield(p,'upper') || isempty(p.upper), parameters(k).upper = Inf; end
    lo = parameters(k).lower; hi = parameters(k).upper;
    if ~isnumeric(lo) || ~isscalar(lo) || isnan(lo) || ...
            ~isnumeric(hi) || ~isscalar(hi) || isnan(hi) || lo >= hi || ...
            p.value < lo || p.value > hi || (p.estimate && (p.value == lo || p.value == hi))
        error('PHASE:customBounds','Parameter %s needs valid bounds and an interior starting value.',p.name);
    end
end
assigned=zeros(numel(parameters));
for k=1:numel(parameters)
    if isfield(parameters,'std') && ~isempty(parameters(k).std)
        sd=parameters(k).std;
        if ~isscalar(sd)||~isfinite(sd)||sd<0||parameters(k).estimate
            error('PHASE:customParameterUncertainty','Assigned standard deviations must be nonnegative and belong to fixed parameters.');
        end
        assigned(k,k)=sd^2;
    end
end
if isfield(config,'customParameterCovariance') && ~isempty(strtrim(config.customParameterCovariance))
    supplied=jsondecode(config.customParameterCovariance);
    if ~isempty(supplied)
        if ~isequal(size(supplied),size(assigned)) || any(~isfinite(supplied(:))) || ...
                norm(supplied-supplied','fro')>1e-10*max(1,norm(supplied,'fro')) || ...
                min(eig((supplied+supplied')/2)) < -1e-10*max(1,norm(supplied,'fro'))
            error('PHASE:customParameterUncertainty','Parameter covariance must be a finite symmetric positive semidefinite matrix in table order.');
        end
        assigned=(supplied+supplied')/2;
    end
end
free=find([parameters.estimate]);
if any(any(assigned(free,:)~=0))
    error('PHASE:customParameterUncertainty','Assigned covariance may refer only to fixed parameters.');
end
names = {parameters.name};
if numel(unique(names)) ~= numel(names)
    error('PHASE:customParameters','Parameter names must be unique.');
end
definition.parameters = parameters;
definition.assignedCovariance=assigned;
if ~isempty(definition.functionFile)
    [~,name,ext] = fileparts(definition.functionFile);
    if ~isfile(definition.functionFile) || ~strcmpi(ext,'.m') || ~isvarname(name)
        error('PHASE:customFunction','Select an existing MATLAB function .m file.');
    end
    definition.functionName = name;
else
    expression = definition.expression;
    % This is a numerical expression, not a MATLAB command/program.
    if isempty(expression) || ~isempty(regexp(expression,'[^\w\s.+*/^(),<>=~&|\-]','once')) || ...
            contains(expression,'..') || ~isempty(regexp(expression,'(?<![<>=~])=(?!=)','once'))
        error('PHASE:customExpression','Enter a numerical expression using x, y, s, t and the named parameters.');
    end
    tokens = regexp(regexprep(expression,'\d*\.?\d+[eE][+-]?\d+','0'), ...
        '[A-Za-z_]\w*','match');
    unknown = setdiff(tokens,[reserved,names]);
    if ~isempty(unknown)
        error('PHASE:customExpression','Unknown expression symbol: %s.',strjoin(unknown,', '));
    end
    % Function calls are allowed only for the small numerical whitelist.
    calls = regexp(expression,'([A-Za-z_]\w*)\s*\(','tokens');
    for k = 1:numel(calls)
        if ~ismember(calls{k}{1},reserved(6:end))
            error('PHASE:customExpression','Only numerical functions can be called.');
        end
    end
end
if isfield(config,'customResidualMean') && strcmp(config.customResidualMean,'constant')
    p=struct('name','phaseResidualMean','value',0,'estimate',true,'lower',-Inf,'upper',Inf);
    fields=fieldnames(parameters);
    for k=1:numel(fields), if ~isfield(p,fields{k}), p.(fields{k})=[]; end; end
    if isempty(parameters), definition.parameters=p;
    else, definition.parameters(end+1)=orderfields(p,parameters); end
    definition.assignedCovariance(end+1,end+1)=0;
end
end

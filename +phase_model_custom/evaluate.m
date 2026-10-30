function [value,jacobian] = evaluate(definition,values,coordinates,context)
%EVALUATE Relative LOS model. Columns: metric X, metric Y, along-track s, days.
% Reference each location to the same acquisition epoch as the observations.
value = relative(values);
if nargout > 1
    jacobian = zeros(size(coordinates,1),numel(values));
    mask=[definition.parameters.estimate];
    if isfield(definition,'assignedCovariance'), mask=mask | diag(definition.assignedCovariance)'>0; end
    for k = find(mask)
        h = eps^(1/3)*max(1,abs(values(k)));
        distance = min(values(k)-definition.parameters(k).lower, ...
            definition.parameters(k).upper-values(k));
        if definition.parameters(k).estimate && isfinite(distance), h=min(h,distance/4); end
        if h<=0, error('PHASE:customDerivative','Cannot differentiate a parameter on its bound.'); end
        plus = values; minus = values;
        plus(k) = plus(k)+h; minus(k) = minus(k)-h;
        jacobian(:,k) = (relative(plus)-relative(minus))/(2*h);
    end
end
    function result = relative(p)
        reference = coordinates; reference(:,4) = context.referenceDay;
        result = absolute(p,coordinates)-absolute(p,reference);
        nuisance=find(strcmp({definition.parameters.name},'phaseResidualMean'));
        if ~isempty(nuisance), result=result+p(nuisance); end
    end
    function result = absolute(p,c)
        coordinateXY=c(:,1:2);
        if isfield(context,'coordinateSystem') && strcmp(context.coordinateSystem,'lonlat')
            [latitude,longitude]=utm2deg(c(:,1),c(:,2),repmat(context.utmZone,size(c,1),1));
            coordinateXY=[longitude,latitude];
        end
        x = coordinateXY(:,1)-context.origin(1); y = coordinateXY(:,2)-context.origin(2);
        s = c(:,3); t = c(:,4);
        if isfield(context,'timeOffset'),t=t+context.timeOffset;end
        parameters = cell2struct(num2cell(p(:)),{definition.parameters.name},1);
        if ~isempty(definition.functionFile)
            % A .m plugin is executable user code, with an explicit contract.
            folder = fileparts(definition.functionFile);
            old = path; cleanup = onCleanup(@() path(old)); %#ok<NASGU>
            addpath(folder,'-begin');
            resolved = which(definition.functionName);
            if ~strcmp(char(java.io.File(resolved).getCanonicalPath()), ...
                    char(java.io.File(definition.functionFile).getCanonicalPath()))
                error('PHASE:customFunction','MATLAB resolved a different custom function. Rename the function.');
            end
            result = feval(definition.functionName,x,y,s,t,parameters);
        else
            names = {definition.parameters.name};
            argumentNames = strjoin([{'x','y','s','t'},names],',');
            evaluator = str2func(['@(' argumentNames ') ' definition.expression]);
            args = num2cell(p(:)');
            result = evaluator(x,y,s,t,args{:});
        end
        if isscalar(result), result = repmat(result,size(c,1),1); end
        if ~isnumeric(result) || ~isreal(result) || numel(result) ~= size(c,1) || any(~isfinite(result(:)))
            error('PHASE:customEvaluation','The model must return one finite real LOS displacement in mm per input row.');
        end
        result = result(:);
    end
end

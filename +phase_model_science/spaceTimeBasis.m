function B=spaceTimeBasis(c,state,referenceDay)
%SPACETIMEBASIS Native polynomial terms, scaled and referenced as displacements.
% A smooth spatial coefficient field multiplies each temporal power.
if strcmp(state.dimension,'1D'),z=c(:,3);else,z=c(:,1:2);end
z=(z-state.origin)./state.spatialScale;
if isfield(state,'spatialKind') && strcmp(state.spatialKind,'native cubic splines')
    n=state.spatialIntervals;
    native=struct('num_x',n,'num_y',n,'num_t',2, ...
        'knots_x',linspace(state.spatialBounds(1,1),state.spatialBounds(2,1),n+1), ...
        'knots_y',linspace(0,1,n+1),'knots_t',[0 .5 1]);
    if size(z,2)==1
        query=[z, .5*ones(size(z)),zeros(size(z)), .5*ones(size(z))];
    else
        native.knots_y=linspace(state.spatialBounds(1,2),state.spatialBounds(2,2),n+1);
        query=[z,zeros(size(z,1),1),.5*ones(size(z,1),1)];
    end
    design=full(phase_model_custom.nativeBasis2D(query,native));
    design=reshape(design,size(z,1),n+2,n+2,4);
    S=sum(design,4);
    if size(z,2)==1,S=reshape(sum(S,2),size(z,1),n+2);
    else,S=reshape(S,size(z,1),(n+2)^2);end
else
S=ones(size(c,1),1);
if state.spatialDegree>=1,S=[S z];end
if state.spatialDegree>=2
    if size(z,2)==1,S=[S z.^2];
    else,S=[S z(:,1).^2 z(:,2).^2 z(:,1).*z(:,2)];end
end
end
T=phase_model_science.trendBasis(c(:,4),referenceDay,state.timeScale,state.degree);
B=zeros(size(c,1),size(S,2)*size(T,2));
for k=1:size(T,2),B(:,(k-1)*size(S,2)+(1:size(S,2)))=S.*T(:,k);end
end

function scientificSelfTest()
%SCIENTIFICSELFTEST Small synthetic checks for model reporting and CRS logic.
% Run in MATLAB with the Mapping Toolbox before a scientific release.

% Most of this AOI is east of the 6-degree boundary (UTM zone 32).
aoi = [5.9 45; 6.4 45; 6.4 45.2; 5.9 45.2; 5.9 45];
[crs,~,epsgCode] = phase_model_beta.selectProjectUtm(aoi);
assert(epsgCode==32632,'The dominant UTM zone was not selected.');
samples = [5.95 45.05; 6.25 45.15];
[x,y] = phase_model_beta.projectWgs84(crs,samples);
[lat,lon] = projinv(crs,x,y);
assert(max(abs([lon,lat]-samples),[],'all') < 1e-7, ...
    'Coordinate round-trip exceeds tolerance.');

field = zeros(11,10);
field(11,4:5) = 10;
sigma = ones(size(field));
anomaly = phase_model_beta.modelAnomalies(field,sigma);
assert(all(strcmp(anomaly.class(1:10),'typical')));
assert(strcmp(anomaly.class{11},'persistent'));
assert(abs(anomaly.peakDeviationMm(11)-10)<1e-12);
assert(abs(anomaly.persistence(11)-0.2)<1e-12);

field(11,5) = 0;
anomaly = phase_model_beta.modelAnomalies(field,sigma);
assert(strcmp(anomaly.class{11},'typical'), ...
    'A single-epoch contrast must not create an anomaly label.');

fprintf('PHASE Model scientific self-test passed.\n');
end

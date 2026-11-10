function validatePreparedGeolocation(workDir)
%VALIDATEPREPAREDGEOLOCATION Never resume candidates at SNAP's (0,0) NoData.
folders=dir(fullfile(workDir,'PATCH_*'));folders=folders([folders.isdir]);
paths=[{workDir},arrayfun(@(v)fullfile(v.folder,v.name),folders,'UniformOutput',false)];
for k=1:numel(paths)
    for version={'ps1.mat','ps2.mat'}
        file=fullfile(paths{k},version{1});if ~isfile(file),continue;end
        vars=whos('-file',file);if ~ismember('lonlat',{vars.name}),continue;end
        data=load(file,'lonlat');xy=data.lonlat;
        if size(xy,2)~=2,continue;end
        invalid=~all(isfinite(xy),2)|abs(xy(:,1))>180|abs(xy(:,2))>=90|all(xy==0,2);
        if any(invalid)
            error('PHASE:InvalidPreparedGeolocation', ...
                ['%s contains %d candidates with invalid/NoData geolocation. ', ...
                 'Enable Prepare data and restart StaMPS from step 1 to apply the geolocation mask.'],file,nnz(invalid));
        end
    end
end
end

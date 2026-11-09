function [target,cfg] = psiAreaWorkDir(current,cfg)
%PSIAREAWORKDIR Keep each selection separate from the full-area PSI dataset.
parent=fileparts(current);[~,name]=fileparts(current);
if ~isempty(cfg.psi_source_dataset),name=cfg.psi_source_dataset;end
if isempty(regexp(name,'^[^/\\]+$','once')) || any(strcmp(name,{'.','..'}))
    error('PHASE:InvalidPsiDataset','Invalid PSI source dataset name.');
end
if strcmp(cfg.psi_area_mode,'full')
    target=fullfile(parent,name);cfg.psi_source_dataset='';return
end
polygon=phase_stamps_beta.validatePsiPolygon(cfg.psi_polygon);
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(uint8(jsonencode(struct('polygon',polygon,'master',cfg.master_date))));
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
target=fullfile(parent,[name '_AOI_' hash(1:12)]);
cfg.psi_source_dataset=name;cfg.psi_polygon=polygon;
% A new area must begin with fresh candidate preparation, never inherit PATCHes.
if ~isfolder(fullfile(target,'PATCH_1')) && ...
        (cfg.stamps_preparation~=0 || ~strcmp(cfg.stamps_first_step,'1'))
    error('PHASE:PsiAreaNeedsPreparation','A new PSI area requires Prepare data enabled and StaMPS first step 1.');
end
end

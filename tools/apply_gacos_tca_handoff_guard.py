"""Make the TRAIN/GACOS to StaMPS step-7 handoff explicit and recoverable.

PHASE_StaMPS.mlapp is an OOXML archive.  This one-shot migration edits only
matlab/document.xml via the shared round-trip helper, keeping the App Designer
archive layout intact.
"""

from pathlib import Path

from mlapp_roundtrip import edit_mlapp


REPO_ROOT = Path(__file__).resolve().parents[1]
MLAPP = REPO_ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"


FRESH_OLD = """\
                                setparm_aps('lambda', lambda, 1);
                                setparm_aps('heading', heading, 1);
                                aps_weather_model('gacos', 3, 3);
                                
                            end
                            
                            stamps(7,7); % execute StaMPS step 7
"""


FRESH_NEW = """\
                                setparm_aps('lambda', lambda, 1);
                                setparm_aps('heading', heading, 1);
                                % TRAIN derives its InSAR output path from pwd.  Force the
                                % processing root so tca&lt;psver&gt;.mat cannot be written elsewhere.
                                cd(cd_fullpath);
                                aps_weather_model('gacos', 3, 3);
                                cd(cd_fullpath);
                                
                            end

                            % begin TRAIN/TCA handoff guard
                            % ps_calc_scla (StaMPS step 7) loads ./tca&lt;psver&gt; as soon as
                            % subtr_tropo='y'.  Validate TRAIN's output here so a missing or
                            % incompatible correction never reaches that opaque load() call.
                            psver_tca = load(fullfile(cd_fullpath, 'psver.mat'), 'psver');
                            if strcmpi(strtrim(small_baseline_flag), 'y')
                                tca_name = sprintf('tca_sb%d.mat', psver_tca.psver);
                            else
                                tca_name = sprintf('tca%d.mat', psver_tca.psver);
                            end
                            tca_path = fullfile(cd_fullpath, tca_name);
                            tca_field = ['ph_tropo_' regexprep(strtrim(tropo_method), '^a_', '')];
                            if exist(tca_path, 'file') ~= 2
                                error('PHASE_StaMPS:gacosTcaMissing', ...
                                    ['TRAIN finished without creating %s. StaMPS step 7 cannot ', ...
                                     'apply tropo_method=%s without this file. aps_weather_model_InSAR ', ...
                                     'resolved to: %s'], tca_path, tropo_method, ...
                                     which('aps_weather_model_InSAR'));
                            end
                            tca_vars = who('-file', tca_path);
                            if ~ismember(tca_field, tca_vars)
                                error('PHASE_StaMPS:gacosTcaInvalid', ...
                                    ['%s exists, but does not contain %s for tropo_method=%s. ', ...
                                     'Restart PHASE at StaMPS step 6 to rebuild the correction.'], ...
                                     tca_path, tca_field, tropo_method);
                            end
                            updateOutput(app, ['TRAIN correction ready for StaMPS step 7: ' tca_path]);
                            clear psver_tca tca_name tca_path tca_field tca_vars
                            % end TRAIN/TCA handoff guard
                            
                            stamps(7,7); % execute StaMPS step 7
"""


RESUME_OLD = """\
                        else
                            stamps(stamps_first_step, stamps_last_step);
                            stamps(6,stamps_last_step); % subtract the computed corrections before the phase unwrapping
"""


RESUME_NEW = """\
                        else
                            % begin resumed TRAIN/TCA handoff guard
                            % A run resumed at step 7 used to assume that tca&lt;psver&gt;.mat was
                            % already present.  Rebuild a missing/incompatible GACOS correction
                            % from the staged .ztd maps before entering StaMPS.
                            cd(cd_fullpath);
                            psver_tca = load(fullfile(cd_fullpath, 'psver.mat'), 'psver');
                            if strcmpi(strtrim(small_baseline_flag), 'y')
                                tca_name = sprintf('tca_sb%d.mat', psver_tca.psver);
                            else
                                tca_name = sprintf('tca%d.mat', psver_tca.psver);
                            end
                            tca_path = fullfile(cd_fullpath, tca_name);
                            tca_field = ['ph_tropo_' regexprep(strtrim(tropo_method), '^a_', '')];
                            tca_vars = {};
                            if exist(tca_path, 'file') == 2
                                tca_vars = who('-file', tca_path);
                            end

                            if (exist(tca_path, 'file') ~= 2 || ~ismember(tca_field, tca_vars)) &amp;&amp; ...
                                    contains(tropo_method, 'a_gacos')
                                gacosFolder = fullfile(cd_fullpath, 'GACOS');
                                gacosZtd = dir(fullfile(gacosFolder, '**', '*.ztd'));
                                if isempty(gacosZtd)
                                    error('PHASE_StaMPS:gacosMapsMissing', ...
                                        ['Cannot rebuild %s: no GACOS .ztd maps were found below %s. ', ...
                                         'Restart PHASE at StaMPS step 6 to download/import them.'], ...
                                         tca_path, gacosFolder);
                                end
                                updateOutput(app, ['Rebuilding missing GACOS correction before StaMPS step 7: ' tca_path]);
                                setparm_aps('gacos_datapath', gacosFolder);
                                setparm_aps('UTC_sat', utc_time);
                                setparm_aps('lambda', lambda, 1);
                                setparm_aps('heading', heading, 1);
                                cd(cd_fullpath);
                                aps_weather_model('gacos', 3, 3);
                                cd(cd_fullpath);
                            end

                            if exist(tca_path, 'file') ~= 2
                                error('PHASE_StaMPS:gacosTcaMissing', ...
                                    ['StaMPS step 7 needs %s for tropo_method=%s, but TRAIN did ', ...
                                     'not create it. aps_weather_model_InSAR resolved to: %s'], ...
                                     tca_path, tropo_method, which('aps_weather_model_InSAR'));
                            end
                            tca_vars = who('-file', tca_path);
                            if ~ismember(tca_field, tca_vars)
                                error('PHASE_StaMPS:gacosTcaInvalid', ...
                                    ['%s does not contain %s. Restart PHASE at StaMPS step 6 ', ...
                                     'to rebuild the selected atmospheric correction.'], ...
                                     tca_path, tca_field);
                            end
                            updateOutput(app, ['TRAIN correction ready for resumed StaMPS processing: ' tca_path]);
                            clear psver_tca tca_name tca_path tca_field tca_vars gacosFolder gacosZtd
                            % end resumed TRAIN/TCA handoff guard

                            stamps(stamps_first_step, stamps_last_step);
                            stamps(6,stamps_last_step); % subtract the computed corrections before the phase unwrapping
"""


def patch(xml: str) -> str:
    marker = "% begin TRAIN/TCA handoff guard"
    if marker in xml:
        print("PHASE_StaMPS.mlapp already contains the TRAIN/TCA handoff guard.")
        return xml

    if xml.count(FRESH_OLD) != 1:
        raise RuntimeError("Fresh GACOS-to-step-7 anchor was not found exactly once")
    if xml.count(RESUME_OLD) != 1:
        raise RuntimeError("Resumed step-7 anchor was not found exactly once")

    xml = xml.replace(FRESH_OLD, FRESH_NEW, 1)
    return xml.replace(RESUME_OLD, RESUME_NEW, 1)


if __name__ == "__main__":
    if not MLAPP.is_file():
        raise SystemExit(f"Expected mlapp at {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"Patched {MLAPP.relative_to(REPO_ROOT)}")

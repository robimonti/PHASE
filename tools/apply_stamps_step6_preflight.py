"""Add a guarded Windows StaMPS step-6 preflight to PHASE_StaMPS.mlapp.

The historical Windows ps_unwrap fallback calls uw_nosnaphu, whose sparse
graph solver can run forever without throwing.  The maintained Windows fork
uses uw_3d plus a bundled snaphu.exe.  Detect a stale/shadowed StaMPS install
and a missing native binary before the synchronous App Designer callback
enters step 6, and persist full MATLAB exception reports for field debugging.
"""

from pathlib import Path

from mlapp_roundtrip import edit_mlapp


REPO_ROOT = Path(__file__).resolve().parents[1]
MLAPP = REPO_ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"


PREFLIGHT_ANCHOR = """\
                updateOutput(app, '----------------------- STEP 2: StaMPS processing started -----------------------');
                
                if contains(ph_output, 'unwrapped')
"""


PREFLIGHT_REPLACEMENT = """\
                updateOutput(app, '----------------------- STEP 2: StaMPS processing started -----------------------');

                % begin Windows StaMPS step-6 preflight
                % The old Windows ps_unwrap branch calls uw_nosnaphu.  Its
                % sparse-graph solver can loop forever without raising an error,
                % which looks exactly like a frozen PHASE app.  Refuse that stale
                % implementation and verify the maintained snaphu path up front.
                if contains(ph_output, 'unwrapped') && stamps_last_step >= 6
                    ps_unwrap_impl = which('ps_unwrap');
                    if isempty(ps_unwrap_impl)
                        error('PHASE_StaMPS:psUnwrapMissing', ...
                            ['ps_unwrap.m was not found. Check the StaMPS installation ', ...
                             'folder selected in PHASE: %s'], installation_folder);
                    end

                    if ispc
                        ps_unwrap_src = fileread(ps_unwrap_impl);
                        old_unwrap_call = regexp(ps_unwrap_src, ...
                            '(?m)^\s*\[ph_uw_some\]\s*=\s*uw_nosnaphu\s*\(', 'once');
                        if ~isempty(old_unwrap_call)
                            error('PHASE_StaMPS:staleWindowsPsUnwrap', ...
                                ['PHASE resolved an obsolete Windows ps_unwrap.m that is known ', ...
                                 'to loop indefinitely during step 6: %s. Update/reinstall the ', ...
                                 'pyccino/StaMPS master fork and make that folder the StaMPS ', ...
                                 'installation folder in PHASE.'], ps_unwrap_impl);
                        end

                        uw_stat_impl = which('uw_stat_costs');
                        if isempty(uw_stat_impl)
                            error('PHASE_StaMPS:uwStatCostsMissing', ...
                                'uw_stat_costs.m was not found below the selected StaMPS installation.');
                        end
                        uw_stat_src = fileread(uw_stat_impl);
                        unix_redirect_call = regexp(uw_stat_src, ...
                            '(?m)^\s*cmdstr\s*=.*>&\s*snaphu\.log', 'once');
                        if ~isempty(unix_redirect_call)
                            error('PHASE_StaMPS:staleWindowsUwStatCosts', ...
                                ['PHASE resolved an obsolete uw_stat_costs.m containing a ', ...
                                 'Unix-only snaphu command: %s. Update/reinstall the ', ...
                                 'pyccino/StaMPS master fork.'], uw_stat_impl);
                        end

                        [snaphu_status, snaphu_where] = dos('where snaphu');
                        if snaphu_status ~= 0 || isempty(strtrim(snaphu_where))
                            expected_snaphu = fullfile(installation_folder, ...
                                'external', 'snaphu', 'bin', 'snaphu.exe');
                            error('PHASE_StaMPS:snaphuMissing', ...
                                ['snaphu.exe is not available on PATH, so StaMPS step 6 ', ...
                                 'cannot unwrap safely on Windows. Expected installer output: ', ...
                                 '%s. Re-run the PHASE installer or restore the native StaMPS ', ...
                                 'binaries before retrying.'], expected_snaphu);
                        end

                        resolved_snaphu = regexprep(strtrim(snaphu_where), '[\\r\\n]+', '; ');
                        updateOutput(app, ['Step 6 preflight: ps_unwrap = ' ps_unwrap_impl]);
                        updateOutput(app, ['Step 6 preflight: snaphu = ' resolved_snaphu]);
                    end

                    updateOutput(app, ['Step 6 preflight passed. Phase unwrapping can take ', ...
                        'a long time; progress is printed in the MATLAB Command Window as ', ...
                        '''Processing IFG x of y''.']);
                    drawnow;
                    clear ps_unwrap_impl ps_unwrap_src old_unwrap_call uw_stat_impl
                    clear uw_stat_src unix_redirect_call snaphu_status snaphu_where
                    clear expected_snaphu resolved_snaphu
                end
                % end Windows StaMPS step-6 preflight
                
                if contains(ph_output, 'unwrapped')
"""


CATCH_ANCHOR = """\
                % If an error occurs
                updateOutput(app, ['An error occurred during script execution: ' ME.message]);
                updateOutput(app, 'Please check each step log file!');
"""


CATCH_REPLACEMENT = """\
                % Persist the complete stack trace.  The UI previously showed only
                % ME.message, which made a Windows native-tool failure impossible to
                % distinguish from a MATLAB/StaMPS source mismatch.
                try
                    phase_error_report = getReport(ME, 'extended', 'hyperlinks', 'off');
                    fprintf(2, '\\n%s\\n', phase_error_report);
                    phase_error_log = fullfile(pwd, 'PHASE_StaMPS_error.log');
                    phase_error_fid = fopen(phase_error_log, 'w');
                    if phase_error_fid ~= -1
                        fprintf(phase_error_fid, '%s\\n', phase_error_report);
                        fclose(phase_error_fid);
                    end
                catch
                    phase_error_log = '';
                end

                % If an error occurs
                if isempty(ME.identifier)
                    updateOutput(app, ['An error occurred during script execution: ' ME.message]);
                else
                    updateOutput(app, ['An error occurred during script execution [' ME.identifier ']: ' ME.message]);
                end
                if ~isempty(phase_error_log)
                    updateOutput(app, ['Full error report saved to: ' phase_error_log]);
                else
                    updateOutput(app, 'Please check each step log file!');
                end
"""


def patch(xml: str) -> str:
    marker = "% begin Windows StaMPS step-6 preflight"
    if marker in xml:
        # Repair archives produced by an early local draft of this migration:
        # Python interpreted MATLAB's \r/\n escape sequences while constructing
        # document.xml, splitting three MATLAB string literals across lines.
        repaired = xml.replace(
            "regexprep(strtrim(snaphu_where), '[\r\n]+', '; ')",
            r"regexprep(strtrim(snaphu_where), '[\r\n]+', '; ')",
        )
        repaired = repaired.replace(
            "fprintf(2, '\n%s\n', phase_error_report)",
            r"fprintf(2, '\n%s\n', phase_error_report)",
        )
        repaired = repaired.replace(
            "fprintf(phase_error_fid, '%s\n', phase_error_report)",
            r"fprintf(phase_error_fid, '%s\n', phase_error_report)",
        )
        repaired = repaired.replace(
            "[snaphu_status, snaphu_where] = system('where snaphu');",
            "[snaphu_status, snaphu_where] = dos('where snaphu');",
        )
        # App Designer stores MATLAB source in CDATA. XML entities inserted
        # into that CDATA remain literal MATLAB text and therefore break the
        # callback; operators must stay unescaped here.
        repaired = repaired.replace(
            "if stamps_first_step &lt;= 6 &amp;&amp; stamps_last_step >= 6",
            "if stamps_first_step <= 6 && stamps_last_step >= 6",
        )
        repaired = repaired.replace(
            "if stamps_first_step <= 6 && stamps_last_step >= 6",
            "if contains(ph_output, 'unwrapped') && stamps_last_step >= 6",
        )
        repaired = repaired.replace(
            r"(?m)^\s*cmdstr\s*=.*>&amp;\s*snaphu\.log",
            r"(?m)^\s*cmdstr\s*=.*>&\s*snaphu\.log",
        )
        if repaired != xml:
            print("Repaired MATLAB escape sequences in the step-6 preflight.")
        else:
            print("PHASE_StaMPS.mlapp already contains the step-6 preflight.")
        return repaired

    if xml.count(PREFLIGHT_ANCHOR) != 1:
        raise RuntimeError("StaMPS processing-start anchor was not found exactly once")
    if xml.count(CATCH_ANCHOR) != 1:
        raise RuntimeError("Start callback catch anchor was not found exactly once")

    xml = xml.replace(PREFLIGHT_ANCHOR, PREFLIGHT_REPLACEMENT, 1)
    return xml.replace(CATCH_ANCHOR, CATCH_REPLACEMENT, 1)


if __name__ == "__main__":
    if not MLAPP.is_file():
        raise SystemExit(f"Expected mlapp at {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"Patched {MLAPP.relative_to(REPO_ROOT)}")

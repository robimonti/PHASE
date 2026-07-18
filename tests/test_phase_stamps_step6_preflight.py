import zipfile


def _read_xml(path):
    with zipfile.ZipFile(path) as mlapp:
        return mlapp.read("matlab/document.xml").decode("utf-8")


def test_step6_preflight_precedes_every_stamps_branch(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")

    preflight = xml.index("% begin Windows StaMPS step-6 preflight")
    first_stamps_call = xml.index("stamps(stamps_first_step,6);")

    assert preflight < first_stamps_call
    guarded = xml[preflight:first_stamps_call]
    # Every unwrapped branch later invokes stamps(6, last), even when a run is
    # nominally resumed from step 7, so the preflight cannot key off first step.
    assert "contains(ph_output, 'unwrapped') && stamps_last_step >= 6" in guarded
    assert "which('ps_unwrap')" in guarded
    assert "which('uw_stat_costs')" in guarded
    assert "dos('where snaphu')" in guarded


def test_step6_preflight_rejects_both_known_windows_regressions(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    preflight = xml.index("% begin Windows StaMPS step-6 preflight")
    end = xml.index("% end Windows StaMPS step-6 preflight", preflight)
    guarded = xml[preflight:end]

    assert "uw_nosnaphu" in guarded
    assert "PHASE_StaMPS:staleWindowsPsUnwrap" in guarded
    assert r">&\s*snaphu\.log" in guarded
    assert "PHASE_StaMPS:staleWindowsUwStatCosts" in guarded
    assert "PHASE_StaMPS:snaphuMissing" in guarded


def test_start_callback_persists_extended_error_report(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")

    assert "getReport(ME, 'extended', 'hyperlinks', 'off')" in xml
    assert "PHASE_StaMPS_error.log" in xml
    assert "ME.identifier" in xml
    assert r"regexprep(strtrim(snaphu_where), '[\r\n]+', '; ')" in xml
    assert r"fprintf(2, '\n%s\n', phase_error_report)" in xml
    assert r"fprintf(phase_error_fid, '%s\n', phase_error_report)" in xml

    # Newlines inside these MATLAB character vectors would make the generated
    # App Designer source invalid even though document.xml remains valid XML.
    assert "fprintf(2, '\n%s\n', phase_error_report)" not in xml
    assert "fprintf(phase_error_fid, '%s\n', phase_error_report)" not in xml


def test_matlab_source_in_cdata_contains_no_literal_xml_entities(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp")
    code_start = xml.index("<![CDATA[")
    code_end = xml.index("]]>", code_start)
    matlab_source = xml[code_start:code_end]

    # App Designer stores source inside CDATA, where entity references are not
    # decoded. A literal &amp;lt; or &amp;amp; would therefore become invalid MATLAB.
    assert "&lt;" not in matlab_source
    assert "&gt;" not in matlab_source
    assert "&amp;" not in matlab_source

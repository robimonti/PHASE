"""Use one effective TRAIN/troposphere decision throughout PHASE_StaMPS.

The processing branch already required both the TRAIN checkbox and
subtr_tropo='y', but preparation and export still checked train_flag alone.
Consequently a no-correction run could later try to load tca2.mat. Derive one
flag after TRAIN recovery and use it for every TRAIN-dependent operation.
"""

from pathlib import Path

from mlapp_roundtrip import edit_mlapp


REPO_ROOT = Path(__file__).resolve().parents[1]
MLAPP = REPO_ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"


REQUEST_BEFORE = "                train_requested = (train_flag == 0);\n"
REQUEST_AFTER = (
    "                train_requested = (train_flag == 0) && "
    "strcmpi(strtrim(subtr_tropo), 'y');\n"
)

RECOVERY_END = """\
                % end TRAIN interactive recovery

                % Conversion of variables' format
"""

RECOVERY_END_WITH_FLAG = """\
                % end TRAIN interactive recovery

                % begin effective tropospheric-correction decision
                % TRAIN-dependent preparation, processing and export must all
                % agree. Checking train_flag alone made a subtr_tropo='n' run
                % attempt to export tca2.mat even though no correction existed.
                tropo_correction_enabled = (train_flag == 0) && ...
                    strcmpi(strtrim(subtr_tropo), 'y');
                if train_flag == 0 && ~tropo_correction_enabled
                    updateOutput(app, ['Tropospheric correction disabled by ', ...
                        'subtr_tropo=' char(string(subtr_tropo)) '.']);
                end
                % end effective tropospheric-correction decision

                % Conversion of variables' format
"""


REPLACEMENTS = (
    (
        "                    if train_flag == 0 % acquire the choice regarding the TRAIN software",
        "                    if tropo_correction_enabled % TRAIN correction is effectively enabled",
        2,
    ),
    (
        "                    if train_flag == 0 && strcmpi(strtrim(subtr_tropo), 'y')",
        "                    if tropo_correction_enabled",
        1,
    ),
    (
        "                    if train_flag == 0 % TRAIN tropo correction included",
        "                    if tropo_correction_enabled % TRAIN tropo correction included",
        1,
    ),
    (
        "                if stamps_last_step > 6 && train_flag == 0 && contains(ph_output, 'unwrapped')",
        "                if stamps_last_step > 6 && tropo_correction_enabled && contains(ph_output, 'unwrapped')",
        1,
    ),
)


def patch(xml: str) -> str:
    marker = "% begin effective tropospheric-correction decision"
    if marker in xml:
        print("PHASE_StaMPS.mlapp already uses the effective tropo flag.")
        return xml

    if xml.count(REQUEST_BEFORE) != 1:
        raise RuntimeError("TRAIN request-intent anchor was not found exactly once")
    if xml.count(RECOVERY_END) != 1:
        raise RuntimeError("TRAIN recovery-end anchor was not found exactly once")

    for old, _new, expected in REPLACEMENTS:
        if xml.count(old) != expected:
            raise RuntimeError(
                f"Expected {expected} occurrence(s) of effective-flag anchor: {old!r}"
            )

    xml = xml.replace(REQUEST_BEFORE, REQUEST_AFTER, 1)
    xml = xml.replace(RECOVERY_END, RECOVERY_END_WITH_FLAG, 1)
    for old, new, _expected in REPLACEMENTS:
        xml = xml.replace(old, new)
    return xml


if __name__ == "__main__":
    if not MLAPP.is_file():
        raise SystemExit(f"Expected mlapp at {MLAPP}")
    edit_mlapp(MLAPP, patch)
    print(f"Patched {MLAPP.relative_to(REPO_ROOT)}")

import zipfile


def _read_xml(path):
    with zipfile.ZipFile(path) as mlapp:
        return mlapp.read("matlab/document.xml").decode("utf-8")


def test_canonical_stamps_app_is_copied_not_moved_or_recycled(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing.mlapp")

    assert xml.count("The installed app is the single source of truth.") == 2
    assert xml.count("copyfile(canonical_stamps_app, dest_mlapp, 'f');") == 2
    assert "movefile(strcat(project_path_full, par, 'PHASE_StaMPS.mlapp')" not in xml
    assert "Copied PHASE_StaMPS.mlapp from previous folder:" not in xml
    assert "previous_folders = dir(fullfile(project_parent_path_full" not in xml
    assert xml.count("PHASE_Preprocessing:StaMPSAppMissing") == 2
    assert xml.count("copyfile(stamps_diagnostic, stamps_folder_full, 'f');") == 2


def test_existing_dataset_input_configuration_is_preserved(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing.mlapp")

    assert xml.count("if isfile(dst_input_mat)") == 2
    assert xml.count("Preserved the existing dataset input_StaMPS.mat.") == 2
    assert xml.count("copyfile(src_input_mat, dst_input_mat);") == 2
    assert "copyfile(src_input_mat, dst_input_mat, 'f');" not in xml
    assert "Refreshed input_StaMPS.mat in the StaMPS processing folder." not in xml


def test_handoff_still_launches_from_the_dataset_folder(phase_root):
    xml = _read_xml(phase_root / "PHASE_Preprocessing.mlapp")

    assert xml.count("stamps_app_full = fullfile(project_parent_path_full, stamps_folder);") == 2
    assert xml.count("cd(stamps_app_full);") == 2
    assert xml.count("run(stamps_app_file);") == 2

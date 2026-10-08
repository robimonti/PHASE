def test_hub_has_one_project_scoped_window_and_lazy_modules(phase_root):
    launcher = (phase_root / "PHASE_Hub.m").read_text(encoding="utf-8")
    hub = (phase_root / "+phase_hub" / "App.m").read_text(encoding="utf-8")

    assert "phase_hub.App(installation,projectRoot)" in launcher
    assert "phase_project.open(projectRoot)" in hub
    assert "uitabgroup(shell)" in hub
    assert "'Title','1 · Preprocessing'" in hub
    assert "'Title','2 · StaMPS'" in hub
    assert "'Title','3 · Model'" in hub
    assert "phase_preprocessing_beta.App( ..." in hub
    assert "phase_model_beta.App( ..." in hub
    assert "phase_stamps_beta.App(selected" in hub
    assert "obj.hasActiveWork()" in hub
    assert "obj.StampsApp.IsRunning" in hub


def test_hub_home_uses_english_and_phase_visual_style(phase_root):
    hub = (phase_root / "+phase_hub" / "App.m").read_text(encoding="utf-8")
    assert "'Your PHASE workspace'" in hub
    assert "'Open project'" in hub
    assert "'New project'" in hub
    assert "'Check for updates'" in hub
    assert "'Title','Project'" in hub
    assert "'ImageSource',logoPath" in hub
    assert "styleButton(newButton,true)" in hub
    assert "'Il tuo workspace PHASE'" not in hub
    assert "'Cerca update'" not in hub


def test_modules_can_embed_without_taking_ownership_of_hub_figure(phase_root):
    files = [
        phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta" / "App.m",
        phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta" / "App.m",
        phase_root / "+phase_model_beta" / "App.m",
    ]
    for path in files:
        controller = path.read_text(encoding="utf-8")
        assert "OwnsFigure = true" in controller
        assert "obj.OwnsFigure = false" in controller
        assert "obj.UIFigure = ancestor(parent,'figure')" in controller
        assert "uigridlayout(obj.HostContainer" in controller
        assert "if obj.OwnsFigure" in controller

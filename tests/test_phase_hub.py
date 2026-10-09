def test_hub_has_one_project_scoped_window_and_lazy_modules(phase_root):
    launcher = (phase_root / "PHASE_Hub.m").read_text(encoding="utf-8")
    hub = (phase_root / "+phase_hub" / "App.m").read_text(encoding="utf-8")

    assert "phase_hub.App(installation,projectRoot)" in launcher
    assert "phase_project.open(projectRoot)" in hub
    assert "obj.Navigation = uihtml(shell" in hub
    assert "obj.HomeView = uihtml(layout" in hub
    assert "obj.PreprocessingTab = uipanel(shell" in hub
    assert "obj.StampsTab = uipanel(shell" in hub
    assert "obj.ModelTab = uipanel(shell" in hub
    assert "phase_preprocessing_beta.App( ..." in hub
    assert "phase_model_beta.App( ..." in hub
    assert "phase_stamps_beta.App(selected" in hub
    assert "obj.hasActiveWork()" in hub
    assert "obj.StampsApp.IsRunning" in hub


def test_hub_home_uses_english_and_phase_visual_style(phase_root):
    hub = (phase_root / "+phase_hub" / "App.m").read_text(encoding="utf-8")
    view = (phase_root / "PHASE_Hub_UI.html").read_text(encoding="utf-8")
    assert "'PHASE_Hub_UI.html'" in hub
    assert "class=\"nav-pills\"" in view
    assert "data-action=\"updates\"" in view
    assert "data-action=\"open\"" in view
    assert "data-action=\"new\"" in view
    assert "PHASE_logo.png" in view
    assert "PHASE_mod1a.png" in view
    assert "PHASE_mod1b.png" in view
    assert "PHASE_mod2.png" in view
    assert "StaMPS PSI" in view
    assert "Displacement Modeling" in view
    assert "Open PHASE Preprocessing" in view
    assert "Open PHASE StaMPS" in view
    assert "Open Displacement Modeling" in view
    assert "border-radius:18px" in view
    assert "Your PHASE workspace" not in view
    assert "One project for Preprocessing" not in view
    assert "Il tuo workspace PHASE" not in view


def test_hub_progress_uses_published_outputs_and_highlights_next_stage(phase_root):
    hub = (phase_root / "+phase_hub" / "App.m").read_text(encoding="utf-8")
    view = (phase_root / "PHASE_Hub_UI.html").read_text(encoding="utf-8")
    progress = (phase_root / "+phase_project" / "workflowStatus.m").read_text(encoding="utf-8")
    assert "phase_project.workflowStatus(obj.ProjectRoot)" in hub
    assert "state.preprocessingComplete = progress.preprocessing" in hub
    assert "data-stage=\"preprocessing\"" in view
    assert "data-stage=\"stamps\"" in view
    assert "data-stage=\"model\"" in view
    assert "card.classList.toggle('done',complete)" in view
    assert "card.classList.toggle('next',recommended)" in view
    assert "<h2>Preprocessing</h2>" not in view
    assert "<h2>StaMPS PSI</h2>" not in view
    assert "<h2>Displacement Modeling</h2>" not in view
    assert "hasFiles(fullfile(folder,'diff0'))" in progress
    assert "isfile(fullfile(folder,'files','mat','PHASEresults.mat'))" in progress


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


def test_completed_modules_offer_next_section_inside_hub(phase_root):
    hub = (phase_root / "+phase_hub" / "App.m").read_text(encoding="utf-8")
    preprocessing = (phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta" / "App.m").read_text(encoding="utf-8")
    stamps = (phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta" / "App.m").read_text(encoding="utf-8")
    assert "function offerNextStage(obj, completed)" in hub
    assert "obj.refreshDatasets();" in hub
    assert "obj.showSection(next);" in hub
    assert "obj.UIFigure.UserData.offerNextStage('preprocessing')" in preprocessing
    assert "obj.UIFigure.UserData.offerNextStage('stamps')" in stamps
    assert "PHASE_Model_beta();" not in stamps

"""Custom routes share native residual algorithms and the original outlier filter."""

import subprocess


def test_custom_filter_and_automatic_filter_remain_verbatim(phase_root):
    for family in ("DET", "STC"):
        for dimension in ("1D", "2D"):
            path = f"MatlabFunctions/STmodel_{family}{dimension}.m"
            original = (phase_root / path).read_text()
            committed = subprocess.check_output(
                ["git", "show", f"HEAD:{path}"], cwd=phase_root, text=True
            )
            custom = (phase_root / path.replace(".m", "_custom.m")).read_text()
            auto_filter = original.split("%% 1) SPLINES-BASED OUTLIER REJECTION", 1)[1].split(
                "%% 2) PREPARE VARIABLES", 1
            )[0]
            custom_filter = custom.split("%% 1) SPLINES-BASED OUTLIER REJECTION", 1)[1].split(
                "%% 2) CUSTOM TWO-STAGE MODELLING", 1
            )[0]
            assert custom_filter == auto_filter
            assert auto_filter in committed


def test_automatic_mode_remains_default(phase_root):
    defaults = (phase_root / "+phase_model_beta/defaultConfig.m").read_text()
    assert "'deterministicModelMode','automatic'" in defaults
    engine = (phase_root / "+phase_model_beta/ProcessingEngine.m").read_text()
    for name in ("DET1D", "DET2D", "STC1D", "STC2D"):
        assert f"STmodel_{name}(" in engine
        assert f"STmodel_{name}_custom(" in engine


def test_custom_residual_methods_reuse_native_algorithms(phase_root):
    covariance = (phase_root / "+phase_model_custom/nativeCovariance.m").read_text()
    covariance_1d = (phase_root / "+phase_model_custom/nativeCovariance1D.m").read_text()
    assert "f2DEmpCovEst_accum(" in covariance
    for source in (covariance, covariance_1d):
        assert "f1DEmpCovEst_accum(" in source
        for fitter in ("fit_exp_auto", "fit_gaussian_auto", "fit_gaussian_cosine_auto"):
            assert fitter + "(" in source
        assert "r(ia)-r(ib)" not in source
    for dimension, adapter in (("1D", covariance_1d), ("2D", covariance)):
        original = (phase_root / f"MatlabFunctions/STmodel_STC{dimension}.m").read_text()
        expected = original.split("    % 3.7) Global Empirical Covariance Estimation", 1)[1].split(
            "    % Compute the covariance matrix of observations", 1
        )[0]
        expected = expected.replace("parfor i =", "for i =").replace(
            "0:round(lag_t(end)/100):round(lag_t(end))", "linspace(0,lag_t(end),101)"
        ).replace("0:round(lag_s(end)/100):round(lag_s(end))", "linspace(0,lag_s(end),101)")
        import re
        optional = r"    \[c_opt_[ts],.*?\]=phase_model_science.constrainCovariance.*?;\n"
        assert re.sub(optional, "", expected, flags=re.S) in re.sub(optional, "", adapter, flags=re.S)
    for dimension in ("1D", "2D"):
        original = (phase_root / f"MatlabFunctions/STmodel_DET{dimension}.m").read_text()
        expected = original.split("%% 4) SPLINES INTERPOLATION", 1)[1].split(
            "% 4.8) Splines solution estimation", 1
        )[0]
        adapter = (phase_root / f"+phase_model_custom/nativeDet{dimension}.m").read_text()
        assert expected in adapter

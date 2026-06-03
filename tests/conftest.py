def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "timelimit: marks tests that depend on timelimit and may behave unexpectedly on different machines (deselect with '-m \"not timelimit\"')"
    )
def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "timelimit: marks tests that depend on timelimit and may behave unexpectedly on different machines (deselect with '-m \"not timelimit\"')"
        "cpp: marks tests that depend on the c++ extension and will fail on the pure python package (deselect with '-m \"not cpp\"')"
    )

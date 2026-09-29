"""
Each module in this package builds exactly one tab of the app.

Every module exposes a single build_tab() function that creates its own
widgets, plot container, and callback as local variables and returns
(title, layout). Because nothing is created at module import time or
stored at module level, calling build_tab() twice (or building several
different tabs) never risks two tabs accidentally sharing the same
underlying widget or pane - each call gets entirely fresh objects.
"""

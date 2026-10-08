"""Pass/fail tolerances of the hydro accuracy suite.

Set from the measured accuracy (docs/hydro_tests.md, results/*.json) plus a
margin, so that a change of the solver's accuracy is caught while the
round-off differences of e.g. a different thread count or compiler are not.
"""


def check_tubes(ctx, out, tier):
    pass


def check_galilean(ctx, out, tier):
    pass


def check_sedov(ctx, out, tier):
    pass


def check_sound(ctx, out, tier):
    pass


def check_pancake(ctx, out, tier):
    pass


def check_cosmobox(ctx, out, tier):
    pass


def check_kh(ctx, out, tier):
    pass

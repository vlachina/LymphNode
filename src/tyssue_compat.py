"""Compatibility helpers for tyssue optional modules."""

import sys
import types


try:
    import tyssue._mesh_generation  # type: ignore # noqa: F401
except ModuleNotFoundError:
    # Some tyssue wheels miss the optional mesh extension but still import it.
    shim = types.ModuleType("tyssue._mesh_generation")

    def make_spherical(*args, **kwargs):
        raise NotImplementedError(
            "tyssue optional mesh-generation extension is unavailable in this "
            "installation; spherical mesh generation cannot be used."
        )

    shim.make_spherical = make_spherical
    sys.modules["tyssue._mesh_generation"] = shim


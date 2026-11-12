# SNAP geometry bridge

`PhaseCskGeometry.java` projects the requested ground AOI into the original
CSK radar grid using the installed SNAP `Orbit`, `SLCImage` and elevation model.
It reads metadata and DEM samples, never SAR amplitude or phase samples.
The DEM's native pixel centres and AOI boundaries are sampled; one 128-pixel
registration window is retained as a halo. PSI preparation masks that halo
using the DEM-corrected export coordinates.

The deliberately distributed `PhaseCskGeometry.class` contains portable Java 8
bytecode, not platform-native code. The original source and a SHA-256 manifest
are shipped alongside it. Rebuild after source edits with:

    python installer/build-csk-geometry.py

Only maintainers need a JDK. At runtime PHASE uses SNAP's bundled Java and
libraries, with a short manifest classpath to support Windows path limits.
External DEMs use the same default vertical-datum handling as SNAP's
`TopoPhaseRemoval`; normalization does not change the vertical datum.

The helper rejects missing elevations, invalid orbit results and geographic
DEMs with rotated grids. Requests partly outside a source scene are clipped
to the available intersection and reported explicitly. The original SAR
sampling and scientific phase operators remain unchanged.

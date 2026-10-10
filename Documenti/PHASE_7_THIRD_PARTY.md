# PHASE 7 third-party components

PHASE is distributed under GNU GPL v3 (see `LICENSE`). MATLAB and ESA SNAP
are separate prerequisites and are **not** included in the installers.

The macOS Apple Silicon DMG includes these separate runtime programs:

| Component | Source and license | Distribution note |
| --- | --- | --- |
| [StaMPS](https://github.com/pyccino/StaMPS/tree/7cabf05eddf8ebe8694e5346fe0f9d48aaef4962) | Pinned source included in the DMG; see its `LICENSE` | Native tools compiled from that source |
| [TRAIN](https://github.com/pyccino/TRAIN/tree/6d0273ae67d2a9f07a696b6a14298ef2c31607d8) | Pinned source included in the DMG; see its `LICENSE` and toolbox notices | MATLAB scripts bundled with PHASE |
| [SNAPHU 2.0.7](https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/) | `external/snaphu/LICENSE` included with StaMPS; matching `snaphu-v2.0.7.tar.gz` source supplied in the DMG and alongside the release downloads | Its CS2-derived parts are licensed for strictly noncommercial use; redistribution must remain free of charge |
| [Triangle 1.6](https://www.cs.cmu.edu/~quake/triangle.html) | Source and `external/triangle/LICENSE.txt` included with StaMPS | Free redistribution without compensation; commercial-product inclusion requires the author's arrangement |
| [GNU awk 5.4.0](https://ftp.gnu.org/gnu/gawk/gawk-5.4.0.tar.xz) | GNU GPL v3 or later; matching source `gawk-5.4.0.tar.xz` is supplied in the DMG and alongside the release downloads | The arm64 executable was built from the supplied source archive |

The GNU awk source archive has SHA-256
`3dd430f0cd3b4428c6c3f6afc021b9cd3c1f8c93f7a688dc268ca428a90b4ac1`.
The SNAPHU source archive has SHA-256
`c03ac126f9a964321bb5d6fb5b4004728368da268d7cb8407bb295a8abe5b262`.
Source links identify the exact upstream versions; the DMG also contains the
StaMPS, TRAIN and Triangle sources beside their binaries.

The Linux installer downloads the same pinned StaMPS and TRAIN forks and
builds the native tools locally. Its `snaphu`, `gawk` and `csh` prerequisites
are supplied separately by the user's operating system. The Windows wizard
downloads the pinned forks and runtime dependencies during installation.

PHASE's free release does not grant a commercial-use licence for Triangle or
the restricted portions of SNAPHU. Consult their original licence files
before redistributing or using those programs commercially.

#!/usr/bin/env python3
"""Rebuild the portable SNAP geometry bridge with Java 8 bytecode."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

root=Path(__file__).resolve().parents[1]
folder=root/'PHASE_Preprocessing/snap2stamps/java'
compiler=shutil.which('javac')
if not compiler:
    raise RuntimeError('A JDK is required to rebuild the geometry helper; end users use the bundled SNAP JRE.')
version=subprocess.run([compiler,'-version'],capture_output=True,text=True,check=True)
flags=['-source','8','-target','8'] if '1.8.' in version.stdout+version.stderr else ['--release','8']
subprocess.run([compiler,*flags,'-d',str(folder),str(folder/'PhaseCskGeometry.java')],check=True)
saved={'sourceSHA256':hashlib.sha256((folder/'PhaseCskGeometry.java').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
       'classSHA256':hashlib.sha256((folder/'PhaseCskGeometry.class').read_bytes()).hexdigest(),
       'classMajorVersion':52}
(folder/'PhaseCskGeometry.build.json').write_text(json.dumps(saved,indent=2)+'\n')
print('Rebuilt PHASE CSK geometry helper (Java 8; uses installed SNAP classes).')

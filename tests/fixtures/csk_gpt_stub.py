"""Test-only GPT substitute: validate substituted XML and create a DIM marker."""
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

root = ET.parse(sys.argv[1]).getroot()
source = root.find(".//node[@id='Read']/parameters/file").text
output = root.find(".//node[@id='Write']/parameters/file").text
assert Path(source).is_file(), source
print('TEST SNAP diagnostic: master reader started', flush=True)
if os.environ.get('PHASE_TEST_GPT_EXIT'):
    print('TEST SNAP error: simulated read/subset failure', flush=True)
    sys.exit(int(os.environ['PHASE_TEST_GPT_EXIT']))
if not os.environ.get('PHASE_TEST_GPT_NO_OUTPUT'):
    Path(output).write_text('<Dimap_Document/>')
    Path(output).with_suffix('.data').mkdir(exist_ok=True)
print('TEST SNAP diagnostic: subset complete', flush=True)

"""Required failures must be explicit, nonzero exits rather than silent fallback."""
import pathlib
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]


def fails(arguments,message):
    result = subprocess.run(arguments,cwd=ROOT,text=True,capture_output=True)
    assert result.returncode != 0, arguments
    assert message in result.stdout+result.stderr, result.stdout+result.stderr


demo = str(ROOT/'build/scatter_demo')
test = str(ROOT/'build/test_collision')
for energy in ('149999','19000001','NaN'):
    fails([demo,energy], 'nonfinite energy' if energy == 'NaN' else 'outside supported interval')
fails([demo,'1000000','0'],'N must be')
fails([demo,'not-a-number'],'Invalid incident energy')
fails([test,'bad_direction'],'input direction must be unit length')
fails([test,'bad_uniform'],'uniforms must be in [0,1)')
with tempfile.TemporaryDirectory() as directory:
    fails([demo,'1000000','1','1',directory],'cannot open')
    target = pathlib.Path(directory)
    (target/'total.endf').write_bytes((ROOT/'data/u238/total.endf').read_bytes()[:160])
    fails([demo,'1000000','1','1',directory],'truncated or unreadable')
print('PASS: invalid energies, direction, variates, CLI arguments, missing and truncated data')

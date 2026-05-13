
#!/usr/bin/env python3
import re
import glob
import os
import sys

nf = 'autophylo.nf'
if not os.path.exists(nf):
    print(f"{nf} not found in {os.getcwd()}")
    sys.exit(1)

text = open(nf).read()
m = re.search(r"params\.species\s*=\s*\[([^\]]+)\]", text, re.S)
if not m:
    print("Could not parse species from autophylo.nf")
    sys.exit(1)

items = []
for line in m.group(1).splitlines():
    s = line.strip()
    if not s:
        continue
    s = s.rstrip(',').strip()
    s = s.strip('\"\' ')
    if s:
        items.append(s)

expected = [s.replace(' ', '_') for s in items]

# locate most recent work data directory
paths = glob.glob('work/*/*/data') + glob.glob('work/*/data')
paths = [p for p in paths if os.path.isdir(p)]
if not paths:
    print('No work data dirs found')
    sys.exit(1)
paths = sorted(paths, key=lambda p: os.path.getmtime(p), reverse=True)
work_data = paths[0]
print('Using work data dir:', work_data)

report = []
for f in sorted(glob.glob(os.path.join(work_data, 'concatenated_*.fasta'))):
    gene = re.sub(r'.*concatenated_(.+)\.fasta$', r'\1', f)
    present = []
    try:
        with open(f) as fh:
            for line in fh:
                if line.startswith('>'):
                    mo = re.search(r"\[organism=([^\]]+)\]", line)
                    if mo:
                        sp = mo.group(1)
                    else:
                        sp = line[1:].strip().split()[0]
                        sp = sp.replace('_', ' ')
                    present.append(sp.replace(' ', '_'))
    except Exception as e:
        print('Error reading', f, e)
        continue
    present_set = set(present)
    missing = [s for s in expected if s not in present_set]
    report.append((gene, sorted(present), missing))

for gene, present, missing in report:
    print('GENE:', gene)
    print('  present_count:', len(present))
    print('  present:', ', '.join(present) if present else '(none)')
    print('  missing_count:', len(missing))
    print('  missing:', ', '.join(missing) if missing else '(none)')
    print()
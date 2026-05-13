#!/usr/bin/env python3
import sys, re

in_path, out_path = sys.argv[1], sys.argv[2]

best_len    = {}
best_header = {}
best_seq    = {}

org_pat = re.compile(r'\[organism=([^\]]+)\]')

current_org    = None
current_header = None
current_seq    = []

def flush():
    if current_org is None:
        return
    seq = "".join(current_seq)
    L   = len(seq)
    if L > best_len.get(current_org, -1):
        best_len[current_org]    = L
        best_header[current_org] = current_header
        best_seq[current_org]    = seq

with open(in_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            m = org_pat.search(line)
            current_org    = m.group(1) if m else line
            current_header = line
            current_seq    = []
        else:
            current_seq.append(line.replace(" ", ""))

flush()

with open(out_path, "w") as fh:
    for org, header in best_header.items():
        fh.write(header + "\n")
        seq = best_seq[org]
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i+70] + "\n")

print(f"Deduplicated {in_path}: kept {len(best_header)} species")
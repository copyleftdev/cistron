#!/usr/bin/env python3
# Differential corpus: vrs-python's fully-justified normalize + ga4gh_identify
# over random (sequence, variant). Validates cistron::fully_justified end-to-end.
import base64, random
from ga4gh.vrs import models
from ga4gh.vrs.dataproxy import _DataProxy
from ga4gh.vrs.normalize import normalize
from ga4gh.core import ga4gh_identify

random.seed(20240607)
ACC = "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"

class MemProxy(_DataProxy):
    def __init__(self, seq): self.seq = seq
    def get_sequence(self, identifier, start=None, end=None):
        return self.seq[start:end] if (start is not None or end is not None) else self.seq
    def get_metadata(self, identifier):
        return {"length": len(self.seq), "aliases": ["ga4gh:" + ACC], "alphabet": "ACGT"}

def rand_seq():
    alpha = "".join(random.sample("ACGT", random.choice([2, 2, 3, 4])))  # bias to repeats
    return "".join(random.choice(alpha) for _ in range(random.randrange(8, 21)))

rows, null_skipped = [], 0
target = 2000
while len(rows) < target:
    seq = rand_seq()
    start = random.randrange(0, len(seq))
    dellen = random.randrange(0, min(5, len(seq) - start + 1))
    end = start + dellen
    altlen = random.randrange(0, 5)
    alt = "".join(random.choice("ACGT") for _ in range(altlen))
    if seq[start:end] == alt:
        null_skipped += 1
        continue  # null: no change
    loc = models.SequenceLocation(sequenceReference=models.SequenceReference(refgetAccession=ACC), start=start, end=end)
    allele = models.Allele(location=loc, state=models.LiteralSequenceExpression(sequence=alt))
    try:
        n = normalize(allele, MemProxy(seq))
    except Exception:
        continue
    rows.append((seq, start, end, alt, ga4gh_identify(n)))

with open("identity/tests/vrs_norm_vectors.tsv", "w") as fh:
    for seq, s, e, alt, vid in rows:
        fh.write(f"{seq}\t{s}\t{e}\t{alt}\t{vid}\n")
print(f"wrote {len(rows)} vectors ({null_skipped} null skipped, RLE included)")

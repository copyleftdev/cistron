#!/usr/bin/env python3
# Generate a VRS conformance corpus from the ga4gh.vrs reference implementation.
# Run in a venv with `ga4gh.vrs` installed; commit identity/tests/vrs_vectors.tsv.
import base64, random
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify

random.seed(1729)  # deterministic corpus
BASES = "ACGT"

def rand_sq():
    raw = bytes(random.randrange(256) for _ in range(24))
    return "SQ." + base64.urlsafe_b64encode(raw).decode()

rows = [("SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 44908821, 44908822, "T")]  # official example
for _ in range(2000):
    start = random.randrange(0, 1_000_000)
    dellen = random.randrange(0, 6)
    altlen = random.randrange(0, 6)
    if dellen == 0 and altlen == 0:
        altlen = 1  # avoid a null (identity) allele
    alt = "".join(random.choice(BASES) for _ in range(altlen))
    rows.append((rand_sq(), start, start + dellen, alt))

with open("identity/tests/vrs_vectors.tsv", "w") as fh:
    for acc, start, end, alt in rows:
        loc = models.SequenceLocation(
            sequenceReference=models.SequenceReference(refgetAccession=acc),
            start=start, end=end)
        allele = models.Allele(location=loc, state=models.LiteralSequenceExpression(sequence=alt))
        fh.write(f"{acc}\t{start}\t{end}\t{alt}\t{ga4gh_identify(allele)}\n")
print(f"wrote {len(rows)} vectors")

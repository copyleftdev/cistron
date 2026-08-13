#!/usr/bin/env python3
# Differential corpus: biocommons hgvs 3'-normalization of genomic variants,
# validating cistron-hgvs::to_hgvs. Run with `hgvs` installed (see hgvsenv).
import random
import hgvs.parser, hgvs.normalizer
from hgvs.dataproviders.interface import Interface

ACC = "NC_TEST"
random.seed(4711)

class Mock(Interface):
    required_version = None
    def __init__(self, seq): self.seq = seq
    def get_seq(self, ac, start=None, end=None):
        return self.seq[start:end] if (start is not None or end is not None) else self.seq
    def data_version(self): return "mock"
    def schema_version(self): return "mock"
    def get_assembly_map(self, a): return {}
    def get_gene_info(self, g): return None
    def get_pro_ac_for_tx_ac(self, t): return None
    def get_acs_for_protein_seq(self, s): return []
    def get_similar_transcripts(self, t): return []
    def get_tx_exons(self, a, b, m): return []
    def get_tx_for_gene(self, g): return []
    def get_tx_for_region(self, a, m, s, e): return []
    def get_tx_identity_info(self, t): return None
    def get_tx_info(self, a, b, m): return None
    def get_tx_mapping_options(self, t): return []

hp = hgvs.parser.Parser()

def raw_hgvs(seq, start, end, alt):
    p1, dell, ref = start + 1, end - start, seq[start:end]
    if dell == 1 and len(alt) == 1:
        return f"{ACC}:g.{p1}{ref}>{alt}"
    if len(alt) == 0:
        return f"{ACC}:g.{p1}del" if dell == 1 else f"{ACC}:g.{p1}_{end}del"
    if dell == 0:
        return f"{ACC}:g.{start}_{start+1}ins{alt}"
    return f"{ACC}:g.{p1}delins{alt}" if dell == 1 else f"{ACC}:g.{p1}_{end}delins{alt}"

def rand_seq():
    alpha = "".join(random.sample("ACGT", random.choice([2, 2, 3, 4])))
    return "".join(random.choice(alpha) for _ in range(random.randrange(8, 21)))

rows, skipped = [], 0
while len(rows) < 2000:
    seq = rand_seq()
    start = random.randrange(0, len(seq))
    dellen = random.randrange(0, min(5, len(seq) - start + 1))
    end = start + dellen
    alt = "".join(random.choice("ACGT") for _ in range(random.randrange(0, 5)))
    if seq[start:end] == alt:            # null
        continue
    if dellen == 0 and start == 0:       # insertion at seq start: cistron can't 5'-flank
        continue
    hn = hgvs.normalizer.Normalizer(Mock(seq), shuffle_direction=3, cross_boundaries=False)
    try:
        vn = hn.normalize(hp.parse_hgvs_variant(raw_hgvs(seq, start, end, alt)))
    except Exception:
        skipped += 1
        continue
    bc = str(vn).split(":", 1)[1]        # strip accession -> bare "g.…"
    rows.append((seq, start, end, alt, bc))

with open("identity/../hgvs/tests/hgvs_vectors.tsv", "w") as fh:
    for seq, s, e, alt, bc in rows:
        fh.write(f"{seq}\t{s}\t{e}\t{alt}\t{bc}\n")
print(f"wrote {len(rows)} vectors ({skipped} skipped)")

------------------------------ MODULE Coords ------------------------------
(***************************************************************************)
(* The coordinate rules `cistron` is built on, checked by TLC.              *)
(*                                                                          *)
(* Three conventions name the same locus three ways and every off-by-one    *)
(* bug in genomics lives in the gaps between them:                          *)
(*                                                                          *)
(*   interbase        [lo, hi)   lo = bases before the span   (internal)    *)
(*   0-based half-open [lo, hi)   identical to interbase       (BED/UCSC)    *)
(*   1-based inclusive [a,  b ]   a = first base, b = last base (VCF/HGVS)   *)
(*                                                                          *)
(* Properties: the interbase<->1-based conversion round-trips, and overlap  *)
(* computed in 1-based-inclusive space agrees with interbase half-open      *)
(* overlap. The canary shows why half-open overlap must use strict `<`.      *)
(***************************************************************************)
EXTENDS Integers

N == 4   \* contig length; small enough for TLC to enumerate every interval

\* Interbase intervals on the contig: 0 <= lo <= hi <= N.
IV == { r \in [lo: 0..N, hi: 0..N] : r.lo <= r.hi }

\* Conversions (defined for NONEMPTY intervals, where a first/last base exists).
OneBasedOf(i) == [first |-> i.lo + 1, last |-> i.hi]
InterbaseOf(o) == [lo |-> o.first - 1, hi |-> o.last]

\* Half-open overlap: touching intervals (one ends where the next begins) do
\* NOT overlap. This is why the comparison is strict.
Overlaps(i, j) == (i.lo < j.hi) /\ (j.lo < i.hi)

\* Inclusive overlap: endpoints may coincide, so the comparison is <=.
OverlapsOB(o, p) == (o.first <= p.last) /\ (p.first <= o.last)

Nonempty(i) == i.lo < i.hi

RoundTrip ==
    \A i \in IV : Nonempty(i) => InterbaseOf(OneBasedOf(i)) = i

OverlapSym ==
    \A i \in IV, j \in IV : Overlaps(i, j) <=> Overlaps(j, i)

\* The convention-agreement law: the two overlap tests must give the same answer
\* on the same pair of loci, each expressed in its own convention.
Agreement ==
    \A i \in IV, j \in IV :
        (Nonempty(i) /\ Nonempty(j)) =>
            (Overlaps(i, j) <=> OverlapsOB(OneBasedOf(i), OneBasedOf(j)))

Inv == RoundTrip /\ OverlapSym /\ Agreement

\* Canary: half-open overlap written with <= (the classic bug). It counts
\* touching intervals as overlapping, so it must DISAGREE with inclusive space.
OverlapsBad(i, j) == (i.lo <= j.hi) /\ (j.lo <= i.hi)
AgreementBad ==
    \A i \in IV, j \in IV :
        (Nonempty(i) /\ Nonempty(j)) =>
            (OverlapsBad(i, j) <=> OverlapsOB(OneBasedOf(i), OneBasedOf(j)))
InvBad == AgreementBad   \* must be VIOLATED

VARIABLE tick
Init == tick = 0
Next == tick' = (IF tick = 0 THEN 1 ELSE 0)
===========================================================================

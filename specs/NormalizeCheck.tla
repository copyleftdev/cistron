-------------------------- MODULE NormalizeCheck --------------------------
(***************************************************************************)
(* An INDEPENDENT check of Normalize.tla's semantics, run by real TLC.      *)
(*                                                                          *)
(* The oracle (tlatools-rs) certifies that an implementation refines the    *)
(* spec, but it evaluates the spec with its own engine -- marking its own   *)
(* homework. TLC re-checks the two properties that make normalization a     *)
(* function you can content-address on:                                     *)
(*                                                                          *)
(*   Sound : the canonical form of an edit is UNIQUE                        *)
(*           (same denotation + both canonical  =>  identical)              *)
(*   Total : every variant HAS a canonical form                            *)
(*                                                                          *)
(* `Bases` is a model constant; the domain is bounded so TLC terminates.     *)
(***************************************************************************)
EXTENDS Integers, Sequences

CONSTANT Bases
MaxRefLen == 3
MaxAlleleLen == 2

\* All base sequences over `Bases` of length 0..n. In TLC a sequence is a
\* function on 1..k, so `[1..k -> Bases]` is exactly the length-k sequences.
BSeq(S, n) == UNION { [1..k -> S] : k \in 0..n }

Splice(s, pos, del, ins) ==
    SubSeq(s, 1, pos) \o ins \o SubSeq(s, pos + Len(del) + 1, Len(s))

Denotation(s, v) == Splice(s, v.pos, v.del, v.ins)

WellFormed(s, v) ==
    /\ v.pos \in 0..Len(s)
    /\ v.pos + Len(v.del) <= Len(s)
    /\ v.del = SubSeq(s, v.pos + 1, v.pos + Len(v.del))

Parsimonious(v) ==
    IF Len(v.del) = 0 \/ Len(v.ins) = 0
    THEN TRUE
    ELSE /\ Head(v.del) # Head(v.ins)
         /\ v.del[Len(v.del)] # v.ins[Len(v.ins)]

CanShiftLeft(s, v) ==
    IF v.pos = 0
    THEN FALSE
    ELSE IF Len(v.del) = 0 /\ Len(v.ins) = 0
         THEN TRUE
         ELSE IF Len(v.ins) = 0
              THEN s[v.pos] = s[v.pos + Len(v.del)]
              ELSE v.ins[Len(v.ins)] = s[v.pos + Len(v.del)]

IsCanonical(s, v) ==
    /\ WellFormed(s, v)
    /\ Parsimonious(v)
    /\ ~CanShiftLeft(s, v)

\* The broken definition the canary uses: left-alignment dropped. In a repeat,
\* two different positions both look "canonical", so Sound must fail on it.
IsCanonicalBad(s, v) == WellFormed(s, v) /\ Parsimonious(v)

Refs == { r \in BSeq(Bases, MaxRefLen) : Len(r) >= 1 }

VariantsOn(s) ==
    { v \in [ pos: 0..Len(s),
              del: BSeq(Bases, MaxAlleleLen),
              ins: BSeq(Bases, MaxAlleleLen) ] : WellFormed(s, v) }

Sound ==
    \A s \in Refs :
        \A v1 \in VariantsOn(s), v2 \in VariantsOn(s) :
            ( /\ IsCanonical(s, v1)
              /\ IsCanonical(s, v2)
              /\ Denotation(s, v1) = Denotation(s, v2) ) => v1 = v2

Total ==
    \A s \in Refs :
        \A v \in VariantsOn(s) :
            \E c \in VariantsOn(s) :
                IsCanonical(s, c) /\ Denotation(s, c) = Denotation(s, v)

SoundBad ==
    \A s \in Refs :
        \A v1 \in VariantsOn(s), v2 \in VariantsOn(s) :
            ( /\ IsCanonicalBad(s, v1)
              /\ IsCanonicalBad(s, v2)
              /\ Denotation(s, v1) = Denotation(s, v2) ) => v1 = v2

Inv == Sound /\ Total       \* must hold
InvBad == SoundBad          \* canary: must be VIOLATED

\* A minimal two-state machine so TLC has a graph to walk; the invariants are
\* closed predicates evaluated in each state.
VARIABLE tick
Init == tick = 0
Next == tick' = (IF tick = 0 THEN 1 ELSE 0)
===========================================================================

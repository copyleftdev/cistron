---------------------------- MODULE Normalize ----------------------------
(***************************************************************************)
(* The reference semantics of variant normalization for `cistron`.         *)
(*                                                                          *)
(* This spec does not *compute* a normal form; it *recognises* a correct    *)
(* one.  The Rust implementation supplies (raw -> out) as a two-state walk  *)
(* and the oracle certifies the step against `Next`.  Recognising is        *)
(* smaller and more trustworthy than reimplementing left-alignment here.    *)
(*                                                                          *)
(* Convention: LEFT-aligned + blunt/parsimonious (bcftools norm / VRS).     *)
(* HGVS would flip the shift direction to the RIGHT; that is a different    *)
(* module, deliberately not this one.                                       *)
(*                                                                          *)
(* Coordinates are interbase: `pos` is the number of reference bases before *)
(* the affected region (0-based).  Sequences are 1-indexed, per TLA+.       *)
(* A variant is a record [pos, del, ins]; bases are integers.               *)
(***************************************************************************)
EXTENDS Integers, Sequences

VARIABLES ref,   \* the reference sequence, a Seq of base integers
          raw,   \* the input variant the implementation was handed
          out,   \* the variant the implementation claims is normalized
          done   \* FALSE at the input, TRUE once normalization is claimed
vars == <<ref, raw, out, done>>

(* The alternate haplotype `v` produces against reference `s`: splice `ins`  *)
(* in place of the `Len(del)` reference bases starting at interbase `pos`.    *)
Splice(s, pos, del, ins) ==
    SubSeq(s, 1, pos) \o ins \o SubSeq(s, pos + Len(del) + 1, Len(s))

Denotation(s, v) == Splice(s, v.pos, v.del, v.ins)

(* A variant is well-formed only if its `del` is *exactly* the reference at   *)
(* the locus it claims to delete.  This is reference-agreement, the check     *)
(* that catches the largest class of real-world variant bugs.                 *)
WellFormed(s, v) ==
    /\ v.pos \in 0..Len(s)
    /\ v.pos + Len(v.del) <= Len(s)
    /\ v.del = SubSeq(s, v.pos + 1, v.pos + Len(v.del))

(* Blunt/parsimonious: del and ins share no leading and no trailing base.     *)
(* (When either allele is empty the condition is vacuously satisfied.)        *)
Parsimonious(v) ==
    IF Len(v.del) = 0 \/ Len(v.ins) = 0
    THEN TRUE
    ELSE /\ Head(v.del) # Head(v.ins)
         /\ v.del[Len(v.del)] # v.ins[Len(v.ins)]

(* Can the variant be shifted one base to the left with the same denotation? *)
(* Derived locally (no search) for the trimmed form:                          *)
(*   pure deletion  : the base leaving on the left equals the one entering.   *)
(*   ins present    : the last inserted base equals the base past the locus.  *)
CanShiftLeft(s, v) ==
    IF v.pos = 0
    THEN FALSE
    ELSE IF Len(v.del) = 0 /\ Len(v.ins) = 0
         THEN TRUE                                   \* null edit: canonical only at pos 0
         ELSE IF Len(v.ins) = 0
              THEN s[v.pos] = s[v.pos + Len(v.del)]  \* pure deletion (Len(del) > 0 here)
              ELSE v.ins[Len(v.ins)] = s[v.pos + Len(v.del)]

(* The canonical form: well-formed, parsimonious, and cannot be shifted left. *)
IsCanonical(s, v) ==
    /\ WellFormed(s, v)
    /\ Parsimonious(v)
    /\ ~CanShiftLeft(s, v)

(*-------------------------------- machine --------------------------------*)

Init ==
    /\ WellFormed(ref, raw)
    /\ out = raw
    /\ done = FALSE

(* The implementation's one move: it hands back `out'` and asserts it is a    *)
(* normalization of `raw`.  `Next` judges that claim -- it does not produce   *)
(* `out'`.  Two obligations, and only two:                                    *)
(*   preservation: `out'` denotes exactly what `raw` denotes;                 *)
(*   canonicity  : `out'` is the left-aligned, blunt representative.          *)
Normalize ==
    /\ done = FALSE
    /\ done' = TRUE
    /\ UNCHANGED <<ref, raw>>
    /\ Denotation(ref, out') = Denotation(ref, raw)
    /\ IsCanonical(ref, out')

Next == Normalize

Spec == Init /\ [][Next]_vars

(*-------------------------------- theorems -------------------------------*)
(* Recorded as the properties the harness enumerates the domain to check.     *)
(*   Soundness   : same denotation => same normal form  (identity is sound)   *)
(*   Idempotence : normalizing a canonical variant is a stutter (free here)   *)
=============================================================================

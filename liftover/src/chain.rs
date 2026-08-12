//! Parsing UCSC `.chain` files into [`Chain`]s.
//!
//! A chain file is a sequence of chains. Each begins with a header line
//! (`chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart
//! qEnd id`) followed by alignment lines: `size dt dq` for every block but the
//! last, then a lone `size`. `dt`/`dq` are the gap sizes in the target (t) and
//! query (q) before the next block.
//!
//! We lift **from t to q**: source = t (always `+` strand in liftOver chains),
//! target = q (which may be `-`). A `-` query keeps its coordinates in the
//! reverse-strand frame, which is exactly what [`Chain`]'s `Minus` handling
//! reflects. All arithmetic is checked so malformed input errors, never panics.

use crate::{Block, BuildError, Chain, Strand};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ChainParseError {
    /// A `chain` header line without the expected fields.
    BadHeader(String),
    /// An alignment line that is neither `size dt dq` nor a lone `size`.
    BadDataLine(String),
    /// A field that should have been a non-negative integer.
    BadNumber(String),
    /// Data appeared before any `chain` header.
    DataBeforeHeader(String),
    /// The target (t) strand was not `+` (liftOver chains require it).
    UnsupportedTargetStrand,
    /// An unrecognized query (q) strand (not `+` or `-`).
    BadQueryStrand(String),
    /// A coordinate computation overflowed.
    Overflow,
    /// The assembled chain was not well-formed.
    Build(BuildError),
}

fn uint(s: &str) -> Result<usize, ChainParseError> {
    s.parse::<usize>()
        .map_err(|_| ChainParseError::BadNumber(s.into()))
}

fn advance(pos: usize, size: usize, gap: usize) -> Result<usize, ChainParseError> {
    pos.checked_add(size)
        .and_then(|p| p.checked_add(gap))
        .ok_or(ChainParseError::Overflow)
}

/// Parse every chain in `text`. Blank lines and `#` comments are ignored.
pub fn parse_chains(text: &str) -> Result<Vec<Chain>, ChainParseError> {
    let mut chains = Vec::new();
    let mut lines = text.lines().peekable();

    while let Some(raw) = lines.next() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some(rest) = line.strip_prefix("chain") else {
            return Err(ChainParseError::DataBeforeHeader(line.into()));
        };

        let f: Vec<&str> = rest.split_whitespace().collect();
        // score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd [id]
        if f.len() < 11 {
            return Err(ChainParseError::BadHeader(line.into()));
        }
        let t_size = uint(f[2])?;
        if f[3] != "+" {
            return Err(ChainParseError::UnsupportedTargetStrand);
        }
        let mut t_pos = uint(f[4])?;
        let q_size = uint(f[7])?;
        let strand = match f[8] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            other => return Err(ChainParseError::BadQueryStrand(other.into())),
        };
        let mut q_pos = uint(f[9])?;

        let mut blocks = Vec::new();
        while let Some(peek) = lines.peek() {
            let dl = peek.trim();
            if dl.is_empty() || dl.starts_with('#') || dl.starts_with("chain") {
                break;
            }
            lines.next();
            let parts: Vec<&str> = dl.split_whitespace().collect();
            match parts.as_slice() {
                [size] => {
                    blocks.push(Block {
                        src: t_pos,
                        tgt: q_pos,
                        len: uint(size)?,
                    });
                    break; // a lone size ends the chain's block list
                }
                [size, dt, dq] => {
                    let (size, dt, dq) = (uint(size)?, uint(dt)?, uint(dq)?);
                    blocks.push(Block {
                        src: t_pos,
                        tgt: q_pos,
                        len: size,
                    });
                    t_pos = advance(t_pos, size, dt)?;
                    q_pos = advance(q_pos, size, dq)?;
                }
                _ => return Err(ChainParseError::BadDataLine(dl.into())),
            }
        }

        chains.push(Chain::new(t_size, q_size, strand, blocks).map_err(ChainParseError::Build)?);
    }

    Ok(chains)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Interval, Lifted};

    #[test]
    fn parses_a_plus_chain_and_lifts() {
        let text = "chain 100 chr1 1000 + 10 30 chr2 2000 + 500 520 1\n10 5 5\n5\n";
        let chains = parse_chains(text).unwrap();
        assert_eq!(chains.len(), 1);
        // block 1: src[10,20) -> tgt[500,510); block 2: src[25,30) -> tgt[515,520)
        assert_eq!(
            chains[0].lift(Interval { start: 12, end: 14 }),
            Ok(Lifted {
                start: 502,
                end: 504,
                strand: Strand::Plus
            })
        );
        assert_eq!(
            chains[0].lift(Interval { start: 26, end: 28 }),
            Ok(Lifted {
                start: 516,
                end: 518,
                strand: Strand::Plus
            })
        );
    }

    #[test]
    fn parses_a_minus_chain_and_reflects() {
        let text = "chain 100 chr1 1000 + 10 20 chr2 2000 - 500 510 1\n10\n";
        let chains = parse_chains(text).unwrap();
        let lifted = chains[0].lift(Interval { start: 12, end: 14 }).unwrap();
        // chain-frame tgt [502,504) reflects to [2000-504, 2000-502).
        assert_eq!(
            lifted,
            Lifted {
                start: 1496,
                end: 1498,
                strand: Strand::Minus
            }
        );
    }

    #[test]
    fn comments_and_blank_lines_are_ignored() {
        let text = "# a comment\n\nchain 1 c 100 + 0 5 d 100 + 0 5 1\n5\n\n";
        assert_eq!(parse_chains(text).unwrap().len(), 1);
    }

    #[test]
    fn terminates_blocks_on_blank_comment_or_next_header() {
        // Headers are 11-field (no id). Each chain's block list ends on a
        // `size dt dq` line whose *direct* successor is, in turn, a blank line,
        // a comment, and the next header — so every terminator branch must fire
        // on its own, not be masked by another.
        let text = concat!(
            "chain 1 c 100 + 0 30 d 100 + 0 30\n",
            "10 5 5\n",
            "\n", // blank terminates chain 1's blocks
            "chain 2 c 100 + 0 30 d 100 + 0 30\n",
            "10 5 5\n",
            "# comment\n", // comment terminates chain 2's blocks
            "chain 3 c 100 + 0 30 d 100 + 0 30\n",
            "10 5 5\n",
            "chain 4 c 100 + 50 60 d 100 + 50 60\n", // header terminates chain 3's blocks
            "10\n",
        );
        let chains = parse_chains(text).unwrap();
        assert_eq!(chains.len(), 4);
        assert_eq!(
            chains[0].lift(Interval { start: 3, end: 5 }),
            Ok(Lifted {
                start: 3,
                end: 5,
                strand: Strand::Plus
            })
        );
        assert_eq!(
            chains[3].lift(Interval { start: 52, end: 54 }),
            Ok(Lifted {
                start: 52,
                end: 54,
                strand: Strand::Plus
            })
        );
    }

    #[test]
    fn malformed_input_errors_never_panics() {
        assert!(matches!(
            parse_chains("chain 1 c 100 +\n5\n"),
            Err(ChainParseError::BadHeader(_))
        ));
        assert!(matches!(
            parse_chains("chain 1 c 100 - 0 5 d 100 + 0 5 1\n5\n"),
            Err(ChainParseError::UnsupportedTargetStrand)
        ));
        assert!(matches!(
            parse_chains("chain 1 c 100 + 0 5 d 100 ? 0 5 1\n5\n"),
            Err(ChainParseError::BadQueryStrand(_))
        ));
        assert!(matches!(
            parse_chains("5\n"),
            Err(ChainParseError::DataBeforeHeader(_))
        ));
        assert!(matches!(
            parse_chains("chain 1 c 100 + 0 5 d 100 + 0 5 1\n5 5\n"),
            Err(ChainParseError::BadDataLine(_))
        ));
        assert!(matches!(
            parse_chains("chain 1 c 100 + 0 5 d 100 + 0 5 1\nx\n"),
            Err(ChainParseError::BadNumber(_))
        ));
    }
}

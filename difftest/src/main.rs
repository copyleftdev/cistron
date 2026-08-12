//! Differential test: does `cistron::normalize` agree with `bcftools norm`?
//!
//! Our own validators (the TLA+ oracle, the brute-force reference, TLC) all
//! check cistron against cistron's self-authored spec. They cannot catch a
//! wrong *convention*. `bcftools norm` is the convention the community trusts,
//! so this shells out to it and compares byte-for-byte.
//!
//! For every (reference, variant) we emit a VCF record, let bcftools left-align
//! it, and compare its POS/REF/ALT against `cistron::normalize` rendered back to
//! VCF. Any disagreement is printed.

use std::collections::HashMap;
use std::fs;
use std::process::Command;

use cistron::{Base, Interbase, Variant};
use cistron_vcf::to_vcf;

const ALPHABET: [Base; 3] = [Base::A, Base::C, Base::G];
const MIN_REF: usize = 3;
const MAX_REF: usize = 6;
const MAX_ALLELE: usize = 2;

struct Expect {
    pos: u64,
    reference: String,
    alt: String,
    contig: String,
    input: String, // POS REF>ALT as fed to bcftools
}

fn render(bases: &[Base]) -> String {
    bases
        .iter()
        .map(|b| ['A', 'C', 'G', 'T'][b.index() as usize])
        .collect()
}

fn seqs(max: usize) -> Vec<Vec<Base>> {
    let mut out = vec![vec![]];
    let mut frontier = vec![vec![]];
    for _ in 0..max {
        let mut next = Vec::new();
        for s in &frontier {
            for &b in &ALPHABET {
                let mut e = s.clone();
                e.push(b);
                next.push(e);
            }
        }
        out.extend(next.iter().cloned());
        frontier = next;
    }
    out
}

fn variants_on(reference: &[Base]) -> Vec<Variant> {
    let mut out = Vec::new();
    for pos in 0..=reference.len() {
        let dmax = (reference.len() - pos).min(MAX_ALLELE);
        for dlen in 0..=dmax {
            let del = reference[pos..pos + dlen].to_vec();
            for ins in seqs(MAX_ALLELE) {
                out.push(Variant {
                    pos: Interbase::new(pos),
                    del: del.clone(),
                    ins,
                });
            }
        }
    }
    out
}

fn main() {
    let bcftools = std::env::var("BCFTOOLS").unwrap_or_else(|_| {
        format!(
            "{}/.local/src/bcftools-1.21/bcftools",
            std::env::var("HOME").unwrap()
        )
    });

    let dir = std::env::temp_dir().join("cistron-difftest");
    fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("ref.fa");
    let vcf_in = dir.join("in.vcf");
    let vcf_out = dir.join("out.vcf");

    let refs: Vec<Vec<Base>> = (MIN_REF..=MAX_REF)
        .flat_map(seqs)
        .filter(|r| r.len() >= MIN_REF)
        .collect();

    let mut fasta = String::new();
    let mut header = String::from("##fileformat=VCFv4.2\n");
    let mut body = String::new();
    let mut expect: HashMap<String, Expect> = HashMap::new();
    let mut id = 0u64;
    let mut skipped = 0u64;

    for (ci, reference) in refs.iter().enumerate() {
        let contig = format!("c{ci}");
        fasta.push_str(&format!(">{contig}\n{}\n", render(reference)));
        header.push_str(&format!(
            "##contig=<ID={contig},length={}>\n",
            reference.len()
        ));

        for raw in variants_on(reference) {
            if raw.denotation(reference) == *reference {
                continue; // no-op edit, not a variant
            }
            let vin = match to_vcf(reference, &raw) {
                Ok(a) => a,
                Err(_) => {
                    skipped += 1;
                    continue;
                }
            };
            let canonical = raw.normalize(reference).unwrap();
            let vexp = match to_vcf(reference, &canonical) {
                Ok(a) => a,
                Err(_) => {
                    skipped += 1;
                    continue;
                }
            };
            let name = format!("v{id}");
            id += 1;
            body.push_str(&format!(
                "{contig}\t{}\t{name}\t{}\t{}\t.\t.\t.\n",
                vin.pos,
                vin.ref_str(),
                vin.alt_str()
            ));
            expect.insert(
                name,
                Expect {
                    pos: vexp.pos,
                    reference: vexp.ref_str(),
                    alt: vexp.alt_str(),
                    contig: contig.clone(),
                    input: format!("{}:{}>{}", vin.pos, vin.ref_str(), vin.alt_str()),
                },
            );
        }
    }

    header.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    fs::write(&fa, &fasta).unwrap();
    fs::write(&vcf_in, format!("{header}{body}")).unwrap();
    let _ = fs::remove_file(dir.join("ref.fa.fai"));

    eprintln!(
        "prepared {} records ({skipped} skipped: null/unanchorable)",
        expect.len()
    );

    let out = Command::new(&bcftools)
        .args(["norm", "-c", "w", "-f"])
        .arg(&fa)
        .arg(&vcf_in)
        .arg("-o")
        .arg(&vcf_out)
        .output()
        .unwrap_or_else(|e| panic!("failed to run {bcftools}: {e}"));
    eprintln!("bcftools: {}", String::from_utf8_lossy(&out.stderr).trim());
    if !out.status.success() {
        eprintln!("bcftools exited nonzero");
        std::process::exit(2);
    }

    let output = fs::read_to_string(&vcf_out).unwrap();
    let mut compared = 0u64;
    let mut agree = 0u64;
    let mut mismatches: Vec<String> = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for line in output.lines() {
        if line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 5 {
            continue;
        }
        let (pos, name, refb, altb) = (c[1], c[2], c[3], c[4]);
        seen.insert(name.to_string());
        let Some(e) = expect.get(name) else { continue };
        compared += 1;
        if pos == e.pos.to_string() && refb == e.reference && altb == e.alt {
            agree += 1;
        } else if mismatches.len() < 20 {
            mismatches.push(format!(
                "  {} in={} | cistron={}:{}>{} bcftools={pos}:{refb}>{altb}",
                e.contig, e.input, e.pos, e.reference, e.alt
            ));
        }
    }

    let dropped: Vec<&String> = expect
        .keys()
        .filter(|k| !seen.contains(k.as_str()))
        .collect();

    println!("\n=== cistron vs bcftools norm ===");
    println!("compared : {compared}");
    println!("agree    : {agree}");
    println!("mismatch : {}", compared - agree);
    println!(
        "dropped  : {} (records bcftools did not emit)",
        dropped.len()
    );
    if !mismatches.is_empty() {
        println!("\nfirst mismatches (contig in=POS:REF>ALT):");
        for m in &mismatches {
            println!("{m}");
        }
    }
    if compared - agree == 0 && dropped.is_empty() {
        println!("\nAGREEMENT: cistron::normalize matches bcftools norm on all {agree} records.");
        std::process::exit(0);
    }
    std::process::exit(1);
}

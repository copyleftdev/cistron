#![no_main]
use libfuzzer_sys::fuzz_target;
use cistron_vcf::parse_line;

fuzz_target!(|data: &[u8]| {
    let s = String::from_utf8_lossy(data);
    let _ = parse_line(&s);
});

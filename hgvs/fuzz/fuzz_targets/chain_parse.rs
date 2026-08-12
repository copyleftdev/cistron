#![no_main]
use libfuzzer_sys::fuzz_target;
use cistron_liftover::chain::parse_chains;

// The chain parser takes untrusted files; it must never panic.
fuzz_target!(|data: &[u8]| {
    let s = String::from_utf8_lossy(data);
    let _ = parse_chains(&s);
});

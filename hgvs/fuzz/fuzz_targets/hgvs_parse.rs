#![no_main]
use libfuzzer_sys::fuzz_target;
use cistron::Base;
use cistron_hgvs::from_hgvs;

// The parser takes untrusted community strings; it must never panic.
fuzz_target!(|data: &[u8]| {
    let reference = [
        Base::A, Base::C, Base::G, Base::T, Base::A, Base::C, Base::G, Base::T,
    ];
    let s = String::from_utf8_lossy(data);
    let _ = from_hgvs(&reference, &s);
});

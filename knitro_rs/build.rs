extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    // Tell cargo to tell rustc to link the system bzip2
    // shared library.
    println!("cargo:rustc-link-lib=knitro");
    #[cfg(feature = "knitro_13")]
    println!("cargo:rustc-link-search=/opt/Knitro/knitro-13.0.1-Linux-64/lib/");
    #[cfg(feature = "knitro_12")]
    println!("cargo:rustc-link-search=/opt/Knitro/knitro-12.3.0-Linux-64/lib/");

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    #[cfg(feature = "knitro_13")]
    println!("cargo:rerun-if-changed=knitro13.h");
    #[cfg(feature = "knitro_12")]
    println!("cargo:rerun-if-changed=knitro12.h");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    #[cfg(feature = "knitro_13")]
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("knitro13.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");
    #[cfg(feature = "knitro_12")]
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("knitro12.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

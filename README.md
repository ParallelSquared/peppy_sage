Work in progress.

Minimal Python bindings for core Sage functionality. Intended as an internal library for JMod. Allows us to build a peptide index directly from a library and quickly compute hyperscore and other PSM characteristics in Rust.

To compile for testing install maturin (via `pip install maturin` or equivalent) and the official Rust language, then navigate to the peppy_sage directory containing the Rust cargo.toml (Cargo.toml) and run the command `maturin develop --release`. This will add a Python library to your current venv with the bindings laid out in peppy_sage/peppy_sage/__init__.py.

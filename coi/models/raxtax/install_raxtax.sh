#!/usr/bin/env bash

# first install rust (may not be needed if already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --no-modify-path --profile minimal
# this reported version of rustc 1.92.0 (ded5c06cf 2025-12-08)

# install raxtax
$HOME/.cargo/bin/cargo install raxtax --version 1.3.0 --profile=ultra
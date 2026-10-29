//! J1 consumes verified saved S1 evidence, never realigns or changes its report.

use super::*;
use crate::observed_transcript_structures as s1;
use serde_json::{Value, json};

#[cfg(test)]
#[path = "junction_assessment/tests.rs"]
mod tests;

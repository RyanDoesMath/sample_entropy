use serde::{Deserialize, Serialize};

/// Struct to store the name along with the entropy values.
#[derive(Debug, Serialize, Deserialize)]
pub struct VitalEntropies {
    pub name: String,
    pub sbp_sampen: f32,
    pub mbp_sampen: f32,
    pub dbp_sampen: f32,
}

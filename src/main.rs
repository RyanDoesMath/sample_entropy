use std::error::Error;
use csv;

use std::time::{Duration, Instant};

fn main() {
    const M: usize = 2;
    const R: f32 = 1.5;
    let mut data: Vec<f32> = Vec::new();
    for i in 1..3000 {
        data.push(i as f32);
    }
    let start = Instant::now();
    let sampen: f32 = sample_entropy(M, R, &data);
    let duration = start.elapsed();
    println!("Sample Entropy: {sampen}");
    println!("Ran in {:?} seconds.", duration);
}

/// Constructs the template vectors for a given time series.
///
/// # Arguments
///
/// * `m` - the window size for a single template.
/// * `data` - the time series data.
///
fn construct_templates(m: usize, data: &Vec<f32>) -> Vec<Vec<f32>> {
    let mut templates: Vec<Vec<f32>> = Vec::new();
    let mut new_template: Vec<f32>;
    for i in m..data.len()+1 {
        let mut new_template: Vec<f32> = Vec::new();
        for j in (i-m)..i {
            new_template.push(data[j]);
        }
        templates.push(new_template);
    }
    return templates;
}

/// Gets the number of matches for a vector of templates.
///
/// # Arguments
///
/// * `templates` - an immutable reference to the a vector containing all templates.
/// * `r` - the distance threshold over which a match does not occur.
///
fn get_matches(templates: &Vec<Vec<f32>>, r: &f32) -> u32 {
    let mut matches: u32 = 0;
    let mut distance: f32;
    
    for i in 0..templates.len() {
        for j in i+1..templates.len() {
            if is_match(&templates[i], &templates[j], &r) {
                matches += 1;
            }
        }
    }
    return matches*2;
}

/// Determines if two templates match.
///
/// The chebyshev distance is a distance metric between two vectors. It is
/// defined as the largest elementwise difference between the vectors.
/// A match occurs between two vectors when their chebyshev distance is
/// less than 'r'. Thus, if at any point the difference between two elements
/// is greater than 'r', we don't need to check any more of the vector.
///
/// # Arguments
///
/// * `vec_1` - an immutable reference to a template vector.
/// * `vec_2` - another immutable reference to a template vector.
/// * `r` - the distance threshold over which a match does not occur.
///
fn is_match(vec_1: &Vec<f32>, vec_2: &Vec<f32>, r: &f32) -> bool{
    for i in 0..vec_1.len() {
        if (vec_1[i] - vec_2[i]).abs() >= *r {
            return false;
        }
    }
    return true;
}

/// Computes sample entropy for a waveform.
///
/// # Arguments
/// * `m` - the smaller of the two template sizes.
/// * `r` - the distance threshold over which a match does not occur.
/// * `data` - a vector containing the waveform data.
///
fn sample_entropy(m: usize, r: f32, data: &Vec<f32>) -> f32 {
    let templates_size_m: Vec<Vec<f32>> = construct_templates(m, &data);
    let m_plus_one = m + 1;
    let templates_size_m_plus_1: Vec<Vec<f32>> = construct_templates(m_plus_one, &data);
    let length_m_template_matches: f32 = get_matches(&templates_size_m, &r) as f32;
    let length_m_plus_1_template_matches: f32 = get_matches(&templates_size_m_plus_1, &r) as f32;
    let ratio: f32= length_m_plus_1_template_matches/length_m_template_matches;
    let sampen: f32 = -(ratio).ln();
    return sampen;
}

/// Reads waveform data from a file in to a vector.
///
/// This is unused for now. Going to rewrite in a much more general way.

fn read_csv(path: &str) -> Result<Vec<i32>, Box<dyn Error>> {
    // Read data from path.
    let mut reader = csv::Reader::from_path(path)?;

    // Initialize vectors.
    let mut record_names: Vec<String> = vec![];
    let mut mean_blood_pressures: Vec<i32> = vec![];
    let mut systolic_blood_pressures: Vec<i32> = vec![];
    let mut diastolic_blood_pressures: Vec<i32> = vec![];
    let mut this_name: &str = "";
    let mut last_name: &str = "";
    // Read the values into the arrays.
    for result in reader.records() {
        let record = result?;
        this_name = &record[0];

        let name = &record[0];
        let mbp = record[1].parse::<f32>()? as i32;
        let sbp = record[2].parse::<f32>()? as i32;
        let dbp = record[3].parse::<f32>()? as i32;
        
        record_names.push(name.to_owned());
        mean_blood_pressures.push(mbp.to_owned());
        systolic_blood_pressures.push(sbp.to_owned());
        diastolic_blood_pressures.push(dbp.to_owned());

        last_name = &record[0];
    }
    Ok(mean_blood_pressures)
}

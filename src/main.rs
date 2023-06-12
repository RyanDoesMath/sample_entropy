use std::error::Error;
use csv;

use std::time::{Instant};

fn main() {
    const M: usize = 2;
    const R: f32 = 0.2;

    let vital_file_result = read_csv("D:/datasets/vitaldb_individual_csvs/0001.csv");

    let start = Instant::now();
    let vital_file_test = match vital_file_result {
        Ok(result) => result,
        Err(error) => panic!("Problem opening the csv file: {:?}", error),
    };
    let duration = start.elapsed();

    println!("{:?}", sample_entropy(M, R, &vital_file_test.sbp));
    println!("{:?}", duration);
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
        new_template = Vec::new();
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

/// Vital file struct for holding the data.
pub struct VitalFile {
    name: String,
    sbp: Vec<f32>,
    mbp: Vec<f32>,
    dbp: Vec<f32>,
}

/// Reads waveform data from a file into a vector.
/// 
/// Due to waves being different length, they cannot be put into a single csv
/// file without doing awkward things. For convenience, csv files for each
/// vital filename was made. The vital_file struct holds this data.
/// 
/// # Arguments
/// * `path` - a reference to a string filepath to a csv file.
/// 

fn read_csv(path: &str) -> Result<VitalFile, Box<dyn Error>> {
    // Read data from path.
    let mut reader = csv::Reader::from_path(path)?;

    // Initialize vectors.
    let mut record_names: Vec<String> = vec![];
    let mut mean_blood_pressures: Vec<f32> = vec![];
    let mut systolic_blood_pressures: Vec<f32> = vec![];
    let mut diastolic_blood_pressures: Vec<f32> = vec![];
    // Read the values into the arrays.
    for result in reader.records() {
        let record = result?;

        let name = &record[0];
        let mbp = record[1].parse::<f32>()?;
        let sbp = record[2].parse::<f32>()?;
        let dbp = record[3].parse::<f32>()?;
        
        record_names.push(name.to_string());
        mean_blood_pressures.push(mbp);
        systolic_blood_pressures.push(sbp);
        diastolic_blood_pressures.push(dbp);
    }

    let new_vital_file = VitalFile{
        name: record_names[0].clone(),
        sbp: systolic_blood_pressures,
        mbp: mean_blood_pressures,
        dbp: diastolic_blood_pressures,
    };

    Ok(new_vital_file)
}

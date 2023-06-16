use std::error::Error;
use csv;
use glob::glob;
use tqdm::tqdm;
use rayon::prelude::*;
use indicatif::ParallelProgressIterator;
use std::time::Instant;
use std::fs::File;
use std::io::Write;

fn main() -> std::io::Result<()> {
    let glob_pattern: String = String::from("D:/datasets/vitaldb_individual_csvs/*.csv");
    
    println!("Reading vital files...");
    let mut vital_files: Vec<VitalFile> = Vec::new();
    for file in tqdm(glob(&glob_pattern).expect("Failed to read glob pattern.")) {
        let path: String = match file {
            Ok(path) => path.into_os_string().into_string().unwrap(),
            Err(error) => panic!("{:?}", error),
        };
        let vital_file = read_csv(&path);
        let vital_file = match vital_file {
            Ok(result) => result,
            Err(error) => panic!("Problem opening the csv file: {:?}", error),
        };
        vital_files.push(vital_file);
    }
    const M: usize = 2;

    println!("Computing sample entropy...");
    let start = Instant::now();
    let sample_entropies: Vec<VitalEntropies> = {
        vital_files.par_iter().progress()
                   .map(|vf| compute_sampen_for_vital_file(M, &vf))
                   .collect::<Vec<VitalEntropies>>()
                   
    };
    let duration = start.elapsed();
    println!("{:?}", duration);

    println!("Saving to csv...");
    let entropy_csv: String = {
        sample_entropies.iter()
                        .map(|ve| vital_entropy_to_csv_line(&ve))
                        .collect::<Vec<String>>()
                        .join("")
    };
    let mut file = File::create("vitaldb_entropies_rust.csv")?;
    write!(file, "{}", entropy_csv);
    Ok(())
}

/// Vital file struct for holding the data.
pub struct VitalFile {
    name: String,
    sbp: Vec<f32>,
    mbp: Vec<f32>,
    dbp: Vec<f32>,
}

/// Struct to store the name along with the entropy values.
pub struct VitalEntropies {
    name: String,
    sbp_sampen: f32,
    mbp_sampen: f32,
    dbp_sampen: f32,
}

/// Computes sample entropy for a single VitalFile struct.
fn compute_sampen_for_vital_file(m: usize, vitalf: &VitalFile) -> VitalEntropies {
    let sbp_sampen: f32 = compute_sampen_for_wave(m, detrend_data(vitalf.sbp.clone()));
    let mbp_sampen: f32 = compute_sampen_for_wave(m, detrend_data(vitalf.mbp.clone()));
    let dbp_sampen: f32 = compute_sampen_for_wave(m, detrend_data(vitalf.dbp.clone()));

    return VitalEntropies {
        name: vitalf.name.clone(),
        sbp_sampen: sbp_sampen,
        mbp_sampen: mbp_sampen,
        dbp_sampen: dbp_sampen,
    }
}

fn compute_sampen_for_wave(m: usize, data: Vec<f32>) -> f32 {
    let stdev: f32 = standard_deviation(&data);
    let r: f32 = stdev*0.2;
    return sample_entropy(m, r, &data);
}

fn vital_entropy_to_csv_line(ve: &VitalEntropies) -> String {
    let comma: String = String::from(", ");
    let newline: String = String::from("\n");

    let mut line: String = ve.name.clone();
    line = line + &comma;
    line = line + &ve.sbp_sampen.to_string() + &comma;
    line = line + &ve.mbp_sampen.to_string() + &comma;
    line = line + &ve.dbp_sampen.to_string() + &newline;

    return line;
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
    let ratio: f32 = length_m_plus_1_template_matches/length_m_template_matches;
    let sampen: f32 = -(ratio).ln();
    return sampen;
}

/// Vectorized one liner for computing the mean of a vector.
fn mean(data: &Vec<f32>) -> f32 {
    data.iter().sum::<f32>() as f32 / data.len() as f32
}

/// Vectorized read-only code that computes standard deviation.
fn standard_deviation(data: &Vec<f32>) -> f32 {
    let xbar: f32 = mean(data);
    let squared_err: Vec<f32> = data.iter().map(|x| (x - xbar).powf(2.0)).collect();
    return ((squared_err.iter().sum::<f32>())/((data.len() as f32))).sqrt();
}

/// Detrends the data via a linear detrending.
/// 
/// Fits an ordinary least squares regression line to the data, then subtracts
/// the estimation from the model to detrend the data. This is done at the
/// suggestion of the 1994 paper by Pincus, S.M.; Goldberger, A.L. titled:
/// "Physiological time-series analysis: what does regularity quantify?"
///
/// In theory there is a nice closed form expression for denominator. It might
/// be useful to speed the program up, but honestly it is already fairly fast.
///
/// # Arguments
/// `data` - a mutable reference to the waveform data.
///
fn detrend_data(data: Vec<f32>) -> Vec<f32> {
    let xbar: f32 = ((data.len() as f32)+1.0)/2.0;
    let ybar: f32 = mean(&data);
    // beta hat is the estimate of the slope parameter.
    let beta_hat: f32 = {
        let data_enum = &data.iter().enumerate().collect::<Vec<_>>();
        let numerator: f32 = data_enum
            .iter()
            .map(|(x, y)| ((*x as f32)+1.0-xbar)*(*y-ybar))
            .collect::<Vec<_>>()
            .into_iter()
            .sum::<f32>();
        let denominator: f32 = data_enum
            .iter()
            .map(|(x, _y)| ((*x as f32)+1.0-xbar).powf(2.0))
            .collect::<Vec<_>>()
            .into_iter()
            .sum::<f32>();
        numerator/denominator 
    };
    // alpha hat is the estimate of the intercept parameter.
    let alpha_hat: f32 = &ybar - &beta_hat*&xbar;
    
    let detrended_data = {
        let data_enum = &data.iter().enumerate().collect::<Vec<_>>();
        data_enum
            .iter()
            .map(|(ix, val)| *val - &alpha_hat - (&beta_hat*((*ix as f32)+1.0)))
            .collect::<Vec<f32>>()
    };
    return detrended_data
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

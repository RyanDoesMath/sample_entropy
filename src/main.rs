use glob::glob;
use indicatif::{ParallelProgressIterator, ProgressBar};
use rayon::prelude::*;
use std::error::Error;
use std::time::Instant;
mod stats;
mod vital_entropies;
use csv::Writer;

use vital_entropies::VitalEntropies;

fn main() -> std::io::Result<()> {
    let glob_pattern: String = String::from("D:/datasets/vitaldb_individual_csvs/*.csv");
    println!("Reading vital files...");
    let vital_files = read_glob_into_vitalfiles(&glob_pattern);
    const M: usize = 2;

    println!("Computing sample entropy...");
    let start = Instant::now();
    let sample_entropies: Vec<VitalEntropies> = {
        vital_files
            .par_iter()
            .progress()
            .map(|vf| compute_sampen_for_vital_file(M, vf))
            .collect::<Vec<VitalEntropies>>()
    };
    let duration = start.elapsed();
    println!("Sample entropy computation finished in: {:?}", duration);

    println!("Saving to csv...");
    let mut writer = Writer::from_path("vitaldb_entropies_rust.csv")?;
    for element in sample_entropies.iter() {
        writer.serialize(element)?;
    }
    writer.flush()?;

    Ok(())
}

/// Vital file struct for holding the data.
pub struct VitalFile {
    name: String,
    sbp: Vec<f32>,
    mbp: Vec<f32>,
    dbp: Vec<f32>,
}

/// Computes sample entropy for a single VitalFile struct.
fn compute_sampen_for_vital_file(m: usize, vitalf: &VitalFile) -> VitalEntropies {
    let sbp_sampen: f32 = compute_sampen_for_wave(m, stats::detrend_data(&vitalf.sbp));
    let mbp_sampen: f32 = compute_sampen_for_wave(m, stats::detrend_data(&vitalf.mbp));
    let dbp_sampen: f32 = compute_sampen_for_wave(m, stats::detrend_data(&vitalf.dbp));

    VitalEntropies {
        name: vitalf.name.clone(),
        sbp_sampen,
        mbp_sampen,
        dbp_sampen,
    }
}

fn compute_sampen_for_wave(m: usize, data: Vec<f32>) -> f32 {
    let stdev: f32 = stats::standard_deviation(&data);
    let r: f32 = stdev * 0.2;
    stats::sample_entropy(m, r, &data)
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

    let new_vital_file = VitalFile {
        name: record_names[0].clone(),
        sbp: systolic_blood_pressures,
        mbp: mean_blood_pressures,
        dbp: diastolic_blood_pressures,
    };

    Ok(new_vital_file)
}

/// Reads all the files from the glob pattern into a vector of VitalFiles.
///
/// # Arguments
/// * `glob_pattern` - a String pattern for glob to use.
///

fn read_glob_into_vitalfiles(glob_pattern: &str) -> Vec<VitalFile> {
    let bar = {
        let glob_files = glob(glob_pattern).expect("Failed to read glob pattern.");
        ProgressBar::new(glob_files.count() as u64)
    };

    let glob_files = glob(glob_pattern).expect("Failed to read glob pattern.");
    let mut vital_files: Vec<VitalFile> = Vec::new();
    for file in glob_files {
        let path: String = match file {
            Ok(path) => path.into_os_string().into_string().unwrap(),
            Err(error) => panic!("{:?}", error),
        };
        let vital_file = match read_csv(&path) {
            Ok(result) => result,
            Err(error) => panic!("Problem opening the csv file: {:?}", error),
        };
        vital_files.push(vital_file);
        bar.inc(1);
    }

    vital_files
}

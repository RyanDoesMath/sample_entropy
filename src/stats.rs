use itertools::Itertools;

/// Constructs the template vectors for a given time series.
///
/// # Arguments
///
/// * `window_size` - the window size for a single template.
/// * `ts_data` - the time series data.
///
fn construct_templates(window_size: usize, ts_data: &Vec<f32>) -> Vec<Vec<f32>> {
    let num_windows = ts_data.len() - window_size + 1;
    (0..num_windows)
        .map(|x| ts_data[x..x + window_size].to_vec())
        .collect::<Vec<Vec<f32>>>()
}

/// Returns 2 times the number of unique pairs of template vectors where the
/// chebyshev distance between each pair of vectors is less than the given
/// threshold.
///
/// # Arguments
///
/// * `templates` - an immutable reference to the a vector containing all templates.
/// * `threshold` - the distance threshold over which a match does not occur.
///
fn get_matches(templates: &[Vec<f32>], threshold: &f32) -> usize {
    let mut matches: u32 = 0;

    for i in 0..templates.len() {
        for j in i + 1..templates.len() {
            if is_match(&templates[i], &templates[j], &threshold) {
                matches += 1;
            }
        }
    }
    (matches * 2).try_into().unwrap()
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
fn is_match(vec_1: &[f32], vec_2: &Vec<f32>, r: &f32) -> bool {
    let threshold = *r;
    return vec_1
        .iter()
        .zip(vec_2)
        .all(|x: (&f32, &f32)| (x.0 - x.1).abs() < threshold);
}

/// Computes sample entropy for a waveform.
///
/// # Arguments
/// * `m` - the smaller of the two template sizes.
/// * `r` - the distance threshold over which a match does not occur.
/// * `data` - a vector containing the waveform data.
///
pub fn sample_entropy(m: usize, r: f32, data: &Vec<f32>) -> f32 {
    let templates_size_m: Vec<Vec<f32>> = construct_templates(m, data);
    let m_plus_one = m + 1;
    let templates_size_m_plus_1: Vec<Vec<f32>> = construct_templates(m_plus_one, data);
    let length_m_template_matches: f32 = get_matches(&templates_size_m, &r) as f32;
    let length_m_plus_1_template_matches: f32 = get_matches(&templates_size_m_plus_1, &r) as f32;
    let ratio: f32 = length_m_plus_1_template_matches / length_m_template_matches;
    let sampen: f32 = -(ratio).ln();
    sampen
}

/// Vectorized one liner for computing the mean of a vector.
pub fn mean(data: &[f32]) -> f32 {
    data.iter().sum::<f32>() / data.len() as f32
}

/// Vectorized read-only code that computes standard deviation.
pub fn standard_deviation(data: &[f32]) -> f32 {
    let xbar: f32 = mean(data);
    let squared_err_sum: f32 = data
        .iter()
        .fold(0_f32, |acc, x| acc + ((x - xbar).powf(2.0)));
    (squared_err_sum / (data.len() as f32)).sqrt()
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
/// `data` - an immutable vector slice of waveform data.
///
pub fn detrend_data(data: &[f32]) -> Vec<f32> {
    let xbar: f32 = (data.len() + 1) as f32 / 2.0;
    let ybar: f32 = mean(data);
    // beta hat is the estimate of the slope parameter.
    let beta_hat: f32 = {
        let (numerator, denominator): (f32, f32) =
            data.iter()
                .enumerate()
                .fold((0_f32, 0_f32), |acc, (index, value)| {
                    let temp = (index + 1) as f32 - xbar;
                    let num_acc = acc.0 + (temp * (value - ybar));
                    let den_acc = acc.1 + (temp.powf(2.0));
                    (num_acc, den_acc)
                });
        numerator / denominator
    };
    // alpha hat is the estimate of the intercept parameter.
    let alpha_hat: f32 = ybar - beta_hat * xbar;

    data.iter()
        .enumerate()
        .map(|(ix, val)| val - alpha_hat - (beta_hat * ((ix as f32) + 1.0)))
        .collect::<Vec<f32>>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constuct_templates_1() {
        let expected: Vec<Vec<f32>> = vec![vec![1_f32], vec![2f32], vec![3_f32]];
        assert_eq!(expected, construct_templates(1, &vec![1_f32, 2_f32, 3_f32]));
    }

    #[test]
    fn test_constuct_templates_2() {
        let expected: Vec<Vec<f32>> = vec![
            vec![1_f32, 2_f32],
            vec![2f32, 3_f32],
            vec![3_f32, 4f32],
            vec![4_f32, 5_f32],
        ];
        assert_eq!(
            expected,
            construct_templates(2, &vec![1_f32, 2_f32, 3_f32, 4_f32, 5_f32])
        );
    }
}

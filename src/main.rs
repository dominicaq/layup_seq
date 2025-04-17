use std::time::Instant;
use plotters::prelude::*;

/*
Layup Sequence

Runtime & Space Complexity Summary:
Time Complexity: O(n)
Space Complexity: O(1)

Report:
A naive solution would use recursion, resulting in exponential runtime O(2^n)
and linear space complexity O(n) due to call stack overhead. Instead, an
iterative approach avoids this overhead and reduces the runtime to O(n).

We can optimize further by observing that only the last two computed values
(S(n-1) and S(n-2)) are needed at any step. By caching these in just two variables
and using a third to hold the current result (S(n)), we achieve constant space complexity, O(1).

Overhead Reduction:
Due to the potential size of the result being too big for standard types, we should return a
`BigInt`, however that would be costly (unwanted memory allocs, complex calculations, etc).
To minimize this overhead, we can instead return a `u64` value instead, using a
large prime modulus (`1,000,000,007`) to represent the result within bounds. This
changes our negation arithmetic however, since we need to handle potential negative
values for correctness. We simply add the modulus before applying the
final modulo to ensure the result stays non-negative:

    Recurrence: 2 * S(n-1) - S(n-2)
    (2 * prev - prev_prev) % MOD

    -prev_prev turns into -> MOD - prev_prev

    (2 * prev + (MOD - prev_prev)) % MOD

This guarantees correctness while keeping all computations within native integer range.

Results:
The analysis of the file `runtime_analysis.png`` shows that our code observations
and optimizations align with the expected O(n) reference.
*/

fn main() {
    // Parameters
    let start = 100;
    let n = 10_000;

    let step = 100;
    let num_trials = 100;

    // Run and plot the algorithm
    create_runtime_plot(&measure_execution_times(start, n, step, num_trials));
}

pub fn layup_sequence(n: usize) -> u64 {
    if n == 1 {
        return 1;
    } else if n == 2 {
        return 2;
    }

    // Using a large prime to prevent int overflow (max: 2^64 - 1)
    const MOD: u64 = 1_000_000_007;

    // Cache the 3 variables
    let mut prev_prev: u64 = 1; // S(n-2)
    let mut prev: u64 = 2;      // S(n-1)
    let mut current: u64 = 0;   // S(n)

    for i in 3..=n {
        current = if i % 2 == 0 {
            // Even: S(n-1) + S(n-2)
            (prev + prev_prev) % MOD
        } else {
            // Odd: 2 * S(n-1) - S(n-2)
            // Use MOD addition to handle possible negative result
            // (2 * prev - prev_prev) % MOD ≡ (2 * prev + MOD - prev_prev) % MOD
            (2 * prev % MOD + (MOD - prev_prev)) % MOD
        };

        prev_prev = prev;
        prev = current;
    }

    current
}

/*
 * Plot Deliverable
 */
fn measure_execution_times(start: usize, end: usize, step: usize, trials: usize) -> Vec<(usize, f64)> {
    let n_values: Vec<usize> = (start..=end).step_by(step).collect();
    let mut times = Vec::new();

    for &n in &n_values {
        let mut total = 0.0;
        // Multiple trials per n value to reduce variability
        for _ in 0..trials {
            let t0 = Instant::now();
            let _ = layup_sequence(n);
            let elapsed = t0.elapsed().as_secs_f64() * 1000.0;
            total += elapsed;
        }
        let avg_time = total / trials as f64;
        times.push((n, avg_time));
        println!("n = {}, avg time: {:.3} ms", n, avg_time);

        if n == end {
            println!("Reached maximum value n = {}", end);
        }
    }

    times
}

fn create_runtime_plot(times: &[(usize, f64)]) {
    // Create a plot
    let root = BitMapBackend::new("runtime_analysis.png", (800, 600))
        .into_drawing_area();
    root.fill(&WHITE).unwrap();

    let max_n = times.iter().map(|(n, _)| *n).fold(0, usize::max);
    let max_t = times.iter().map(|(_, t)| *t).fold(0.0, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Runtime Analysis of layup_sequence (up to n={})", max_n), ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0..max_n, 0.0..max_t * 1.1)
        .unwrap();

    chart.configure_mesh()
        .x_labels(10)
        .y_labels(10)
        .x_desc("Input Size (n)")
        .y_desc("Runtime (milliseconds)")
        .draw()
        .unwrap();

    // Draw data points
    chart.draw_series(LineSeries::new(
        times.iter().map(|(n, t)| (*n, *t)),
        &RED,
    ))
    .unwrap()
    .label("Runtime")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Draw a line representing perfect O(n) growth for comparison
    let scale_factor = max_t / max_n as f64;
    let reference_points: Vec<(usize, f64)> = (0..=10).map(|i| {
        let n = i * 1000;
        (n, n as f64 * scale_factor)
    }).collect();

    chart.draw_series(LineSeries::new(
        reference_points,
        &BLUE.mix(0.5),
    ))
    .unwrap()
    .label("O(n) reference")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE.mix(0.5)));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    // Save the result
    root.present().unwrap();
    println!("Plot has been saved as runtime_analysis.png");
}

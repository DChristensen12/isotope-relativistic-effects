use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// Bins a mass spectrometer scan's Intensity column into peak areas. Used by
// The_Method.calculate_areas_and_isotopes to turn a CSV's raw Mass/Intensity
// columns into the areas (and optionally masses) for each peak in a scan.

/// Sums `intensity` into consecutive bins covering the row range `[row_start,
/// row_end]`, and optionally reports which `mass` sits at each bin's peak
/// Intensity.
///
/// row_start / row_end use the 1-based "spreadsheet row number" convention rather
/// than 0-based array indices -- these are the row numbers as they'd appear in the
/// original spreadsheet. `pandas_difference_index` and `header_affect_on_index` are
/// the two offsets (normally 1 and 1) needed to convert a spreadsheet row number
/// into a 0-based array position: one for the 1-based-vs-0-based indexing
/// difference, one for the header row.
///
/// Two behaviors worth knowing about:
///
/// 1. The validity check (`row_start < 0 || row_end >= len || row_start >
///    row_end`) compares the raw, un-offset row_end against the array length
///    rather than the offset-adjusted position, so it's slightly more conservative
///    than the true usable bound. A range that fails this check returns
///    (empty, empty) rather than an error.
///
/// 2. The last bin in a row range that doesn't divide evenly by bin_size is not
///    clamped to row_end -- it can extend past the intended boundary, pulling in
///    whatever rows come next, all the way up to wherever the underlying array
///    actually ends (a bin that runs past the end of the array is truncated rather
///    than treated as an error).
///
/// Ties for the max Intensity within a bin resolve to the first occurrence.
pub fn bin_areas_core(
    intensity: &[f64],
    mass: Option<&[f64]>,
    row_start: i64,
    row_end: i64,
    bin_size: i64,
    pandas_difference_index: i64,
    header_affect_on_index: i64,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    if let Some(m) = mass {
        if m.len() != intensity.len() {
            return Err(format!(
                "mass length ({}) must match intensity length ({}) -- they should always \
                 come from the same dataframe",
                m.len(),
                intensity.len()
            ));
        }
    }
    // A non-positive bin_size would leave current_bin_start stuck in place forever
    // below, spinning while holding the GIL -- reject it up front instead.
    if bin_size <= 0 {
        return Err(format!("bin_size must be positive, got {bin_size}"));
    }

    let len = intensity.len() as i64;
    let offset = pandas_difference_index + header_affect_on_index;

    if row_start < 0 || row_end >= len || row_start > row_end {
        return Ok((Vec::new(), Vec::new()));
    }

    let mut areas = Vec::new();
    let mut masses_at_max = Vec::new();

    // Walk bins from row_start until current_bin_start passes row_end.
    let mut current_bin_start = row_start;
    while current_bin_start < row_end {
        // Not clamped to row_end -- see point 2 in the doc comment above.
        let current_bin_end = current_bin_start + bin_size;

        // Convert this bin's spreadsheet-row bounds into a 0-based array slice
        // (`+1` to make current_bin_end inclusive). A slice that runs past the
        // array's actual length is truncated rather than treated as an error, so
        // clamp both bounds into [0, len] instead of indexing out of range.
        let slice_start = (current_bin_start - offset).clamp(0, len);
        let slice_end_exclusive = (current_bin_end - offset + 1).clamp(0, len);

        if slice_start >= slice_end_exclusive {
            break;
        }
        let s = slice_start as usize;
        let e = slice_end_exclusive as usize; // exclusive

        let bin = &intensity[s..e];
        areas.push(bin.iter().sum());

        // Track which row in this bin had the peak Intensity so its Mass can be
        // recorded.
        if let Some(m) = mass {
            let (local_max_idx, _) = bin.iter().enumerate().fold(
                (0usize, f64::NEG_INFINITY),
                |(best_i, best_v), (i, &v)| if v > best_v { (i, v) } else { (best_i, best_v) },
            );
            masses_at_max.push(m[s + local_max_idx]);
        }

        current_bin_start += bin_size + 1;
    }

    Ok((areas, masses_at_max))
}

/// PyO3-facing wrapper, imported from Python as `spectral_kernels.bin_areas(...)`.
/// Pulls raw f64 slices out of the numpy arrays, hands them to bin_areas_core, and
/// returns the results as numpy arrays. Pass mass=None to skip the argmax work when
/// only the summed areas are needed.
#[pyfunction]
#[pyo3(signature = (
    intensity, row_start, row_end, bin_size, mass=None,
    pandas_difference_index=1, header_affect_on_index=1
))]
#[allow(clippy::too_many_arguments)]
fn bin_areas<'py>(
    py: Python<'py>,
    intensity: PyReadonlyArray1<f64>,
    row_start: i64,
    row_end: i64,
    bin_size: i64,
    mass: Option<PyReadonlyArray1<f64>>,
    pandas_difference_index: i64,
    header_affect_on_index: i64,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let intensity_slice = intensity
        .as_slice()
        .map_err(|e| PyValueError::new_err(format!("intensity array: {e}")))?;
    let mass_slice = mass
        .as_ref()
        .map(|m| {
            m.as_slice()
                .map_err(|e| PyValueError::new_err(format!("mass array: {e}")))
        })
        .transpose()?;

    let (areas, masses_at_max) = bin_areas_core(
        intensity_slice,
        mass_slice,
        row_start,
        row_end,
        bin_size,
        pandas_difference_index,
        header_affect_on_index,
    )
    .map_err(PyValueError::new_err)?;

    Ok((
        areas.into_pyarray_bound(py),
        masses_at_max.into_pyarray_bound(py),
    ))
}

// Registers the module so `import spectral_kernels` works after building with
// maturin. Add more kernels here later if others earn their way out of Python too
// (e.g. the mass-tolerance row search from the rows_are_masses branch).
#[pymodule]
fn spectral_kernels(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bin_areas, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_evenly_behaved_bins() {
        // Mass = 100..119, Intensity = 1..20
        let mass: Vec<f64> = (0..20).map(|i| 100.0 + i as f64).collect();
        let intensity: Vec<f64> = (0..20).map(|i| (i + 1) as f64).collect();
        let (areas, masses) =
            bin_areas_core(&intensity, Some(&mass), 3, 13, 2, 1, 1).unwrap();
        assert_eq!(areas, vec![9.0, 18.0, 27.0, 36.0]);
        assert_eq!(masses, vec![103.0, 106.0, 109.0, 112.0]);
    }

    #[test]
    fn last_bin_can_extend_past_row_end() {
        // row_end=10 with bin_size=3 doesn't divide evenly, so the last bin isn't
        // clamped and pulls in rows past row_end that still exist further in the
        // array.
        let mass: Vec<f64> = (0..20).map(|i| 100.0 + i as f64).collect();
        let intensity: Vec<f64> = (0..20).map(|i| (i + 1) as f64).collect();
        let (areas, masses) =
            bin_areas_core(&intensity, Some(&mass), 3, 10, 3, 1, 1).unwrap();
        assert_eq!(areas, vec![14.0, 30.0]);
        assert_eq!(masses, vec![104.0, 108.0]);
    }

    #[test]
    fn overrunning_bin_is_truncated_at_the_array_end() {
        // Mass = 200..209, Intensity = 1..10 (only 10 rows this time). The
        // unclamped bin_end would extend past the array's actual length, so the
        // out-of-range slice is truncated rather than erroring.
        let mass: Vec<f64> = (0..10).map(|i| 200.0 + i as f64).collect();
        let intensity: Vec<f64> = (0..10).map(|i| (i + 1) as f64).collect();
        let (areas, masses) =
            bin_areas_core(&intensity, Some(&mass), 3, 9, 4, 1, 1).unwrap();
        assert_eq!(areas, vec![20.0, 34.0]);
        assert_eq!(masses, vec![205.0, 209.0]);
    }

    #[test]
    fn invalid_ranges_return_empty_instead_of_erroring() {
        // A bad row range signals nothing to compute, not a programming mistake,
        // so it returns empty results rather than an Err.
        let mass: Vec<f64> = (0..20).map(|i| 100.0 + i as f64).collect();
        let intensity: Vec<f64> = (0..20).map(|i| (i + 1) as f64).collect();
        let (areas, masses) = bin_areas_core(&intensity, Some(&mass), -1, 10, 2, 1, 1).unwrap();
        assert!(areas.is_empty() && masses.is_empty());
        let (areas, masses) = bin_areas_core(&intensity, Some(&mass), 5, 100, 2, 1, 1).unwrap();
        assert!(areas.is_empty() && masses.is_empty());
    }

    #[test]
    fn ties_in_intensity_resolve_to_first_occurrence() {
        let intensity = [5.0, 5.0, 1.0, 0.0, 0.0, 0.0];
        let mass = [300.0, 301.0, 302.0, 303.0, 304.0, 305.0];
        let (areas, masses) = bin_areas_core(&intensity, Some(&mass), 2, 3, 1, 1, 1).unwrap();
        assert_eq!(areas, vec![10.0]);
        assert_eq!(masses, vec![300.0]); // not 301.0, despite the tie
    }

    #[test]
    fn mass_is_optional_and_skips_argmax_work_when_omitted() {
        let intensity = [5.0, 5.0, 1.0, 0.0, 0.0, 0.0];
        let (areas, masses) = bin_areas_core(&intensity, None, 2, 3, 1, 1, 1).unwrap();
        assert_eq!(areas, vec![10.0]);
        assert!(masses.is_empty());
    }

    #[test]
    fn rejects_mismatched_mass_length() {
        let intensity = [1.0, 2.0, 3.0];
        let mass = [1.0, 2.0];
        assert!(bin_areas_core(&intensity, Some(&mass), 0, 2, 1, 1, 1).is_err());
    }

    #[test]
    fn rejects_non_positive_bin_size_instead_of_hanging() {
        let intensity = [1.0, 2.0, 3.0];
        assert!(bin_areas_core(&intensity, None, 0, 2, 0, 1, 1).is_err());
    }
}

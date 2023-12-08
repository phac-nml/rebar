use crate::dataset::{Dataset, SearchResult};
use crate::recombination::Recombination;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::{debug, warn};
use std::collections::BTreeMap;
use std::fmt;
use std::str::FromStr;

// ----------------------------------------------------------------------------
// Validate

#[derive(Clone, Debug)]
pub struct Validate {
    pub status: Status,
    pub details: Vec<Details>,
}

// ----------------------------------------------------------------------------
// Status

#[derive(Clone, Debug)]
pub enum Status {
    Pass,
    Fail,
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let result = match self {
            Status::Pass => "pass",
            Status::Fail => "fail",
        };
        write!(f, "{}", result)
    }
}

impl FromStr for Status {
    type Err = Report;
    fn from_str(input: &str) -> Result<Self, Report> {
        let result = match input {
            "pass" => Status::Pass,
            "fail" => Status::Fail,
            _ => return Err(eyre!("Unknown status: {input}")),
        };
        Ok(result)
    }
}

// ----------------------------------------------------------------------------
// Details

#[derive(Clone, Debug, PartialEq)]
pub enum Details {
    IncorrectRecombinant,
    IncorrectParent,
    IncorrectPopulation,
    NoRecombinationDetected,
}

impl fmt::Display for Details {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let result = match self {
            Details::IncorrectRecombinant => "incorrect_recombinant",
            Details::IncorrectParent => "incorrect_parent",
            Details::IncorrectPopulation => "incorrect_population",
            Details::NoRecombinationDetected => "no_recombination_detected",
        };
        write!(f, "{}", result)
    }
}

impl FromStr for Details {
    type Err = Report;
    fn from_str(input: &str) -> Result<Self, Report> {
        let result = match input {
            "incorrect_recombinant" => Details::IncorrectRecombinant,
            "incorrect_parent" => Details::IncorrectParent,
            "incorrect_population" => Details::IncorrectPopulation,
            "no_recombination_detected" => Details::NoRecombinationDetected,
            _ => return Err(eyre!("Unknown details: {input}")),
        };
        Ok(result)
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

/// Validate a best match and recombination to expected dataset values.
pub fn validate(
    dataset: &Dataset,
    best_match: &SearchResult,
    recombination: &Recombination,
) -> Result<Option<Validate>, Report> {
    // Use the sequence ID as the expected population
    let expected_population = best_match.sequence_id.replace("population_", "");

    // If the sequence ID was not in the dataset, return no validation results
    if !dataset.populations.contains_key(&expected_population) {
        let validate: Option<Validate> = None;
        Ok(validate)
    } else {
        // ----------------------------------------------------------------
        // Population
        //
        // Is the expected population (sequence id) exactly the same as
        // the best match (consensus population)?

        let observed_population = &best_match.consensus_population;
        let validate_population = expected_population == *observed_population;

        debug!(
            "Validate population: observed={observed_population}, expected={expected_population}"
        );

        // ----------------------------------------------------------------
        // Recombinant
        //
        // Is this a recombinant, and if so, is the recombinant ancestor
        // exactly as expected?

        let observed_recombinant = &recombination.recombinant;
        let expected_recombinant =
            dataset.phylogeny.get_recombinant_ancestor(&expected_population)?;

        let validate_recombinant = expected_recombinant == *observed_recombinant;

        debug!("Validate recombinant: observed={observed_recombinant:?}, expected={expected_recombinant:?}");

        // ----------------------------------------------------------------
        // Parent Validation

        let observed_parents = &recombination.parents;
        let expected_parents = if let Some(expected_recombinant) = expected_recombinant {
            dataset.phylogeny.get_parents(&expected_recombinant)?
        } else {
            Vec::new()
        };

        debug!("Validate parents: observed={observed_parents:?}, expected={expected_parents:?}");
        // parent validation is already done in recombination::detect_recombination
        // to decide whether it's a novel variant or not. It seems wasteful to run it again...

        let validate_parent =
            compare_parents(observed_parents, &expected_parents, dataset)?;

        // ----------------------------------------------------------------
        // Recombination Validation
        // Were parents and breakpoints detected at all?

        // ----------------------------------------------------------------
        // Summary

        let validate = if validate_population && validate_recombinant && validate_parent {
            Validate {
                status: Status::Pass,
                details: Vec::new(),
            }
        } else {
            let mut validate = Validate {
                status: Status::Fail,
                details: Vec::new(),
            };
            if !validate_population {
                validate.details.push(Details::IncorrectPopulation);
            }
            if !validate_parent {
                // Were parents and breakpoints detected at all?
                if !expected_parents.is_empty() && observed_parents.is_empty() {
                    validate.details.push(Details::NoRecombinationDetected);
                } else {
                    validate.details.push(Details::IncorrectParent);
                }
            }
            if !validate_recombinant
                && !validate.details.contains(&Details::NoRecombinationDetected)
            {
                validate.details.push(Details::IncorrectRecombinant);
            }
            warn!(
                "{expected_population} failed validation: {details}",
                details = validate.details.iter().join(", ")
            );
            validate
        };

        Ok(Some(validate))
    }
}

/// Compare expected and observed recombination parents.
pub fn compare_parents(
    observed: &Vec<String>,
    expected: &[String],
    dataset: &Dataset,
) -> Result<bool, Report> {
    debug!("Observed parents: {observed:?}");
    debug!("Expected parents: {expected:?}");

    // check if the expected parents are actually in the dataset populations
    // ie. we actually have sequence data for them
    let expected_filter = expected
        .iter()
        .map(|p| dataset.get_ancestor_with_sequence(p).unwrap_or(p.clone()))
        .collect_vec();

    let expected = &expected_filter;

    // one is empty, other is not
    if (observed.is_empty() && !expected.is_empty())
        || (!observed.is_empty() && expected.is_empty())
    {
        return Ok(false);
    }
    // exact match
    let mut observed_exact = observed.clone();
    observed_exact.retain(|p| expected.contains(p));
    if observed_exact.len() == expected.len() {
        return Ok(true);
    }

    // ------------------------------------------------------------------------
    // allow some flexibility with regards to match parents/child

    let mut expected_fulfilled: BTreeMap<&str, Option<String>> = BTreeMap::new();
    expected.iter().for_each(|e| {
        expected_fulfilled.insert(e, None);
    });

    let recombination = false;

    // ------------------------------------------------------------------------
    // expected is ancestor of observed, but no recombinant in between
    // Ex. XE: expected = BA.1,BA.2; observed=BA.1.17.2.1,BA.2
    let observed_ancestors = observed
        .iter()
        .map(|pop| {
            let paths = dataset.phylogeny.get_ancestors(pop, recombination)?;
            let ancestors = paths.into_iter().flatten().unique().collect_vec();
            Ok((pop, ancestors))
        })
        .collect::<Result<Vec<_>, Report>>()?;

    observed_ancestors.into_iter().for_each(|(obs, ancestors)| {
        debug!("Observed {obs} ancestors: {ancestors:?}");
        expected.iter().for_each(|exp| {
            if ancestors.contains(exp) {
                expected_fulfilled.insert(exp, Some(obs.clone()));
            }
        });
    });
    debug!("expected_fulfilled: {expected_fulfilled:?}");

    // if no expected parents remain to be fullfilled
    if !expected_fulfilled.values().contains(&None) {
        return Ok(true);
    }

    // ------------------------------------------------------------------------
    // observed is ancestor of expected, but no recombinant in between
    let expected_ancestors = expected
        .iter()
        .map(|pop| {
            let paths = dataset.phylogeny.get_ancestors(pop, recombination)?;
            let ancestors = paths.into_iter().flatten().unique().collect_vec();
            Ok((pop, ancestors))
        })
        .collect::<Result<Vec<_>, Report>>()?;

    expected_ancestors.into_iter().for_each(|(exp, ancestors)| {
        debug!("Expected {exp} ancestors: {ancestors:?}");
        observed.iter().for_each(|obs| {
            if ancestors.contains(obs) {
                expected_fulfilled.insert(exp, Some(obs.clone()));
            }
        });
    });
    debug!("expected_fulfilled: {expected_fulfilled:?}");

    // if no expected parents remain to be fullfilled
    if !expected_fulfilled.values().contains(&None) {
        return Ok(true);
    }

    Ok(false)
}

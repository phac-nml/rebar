use crate::dataset::{Dataset, SearchResult};
use crate::recombination::Recombination;
use color_eyre::eyre::{eyre, Report, Result};
use itertools::Itertools;
use log::warn;
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

#[derive(Clone, Debug)]
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

        // ----------------------------------------------------------------
        // Recombinant
        //
        // Is this a recombinant, and if so, is the recombinant ancestor
        // exactly as expected?

        let observed_recombinant = &best_match.recombinant;
        let expected_recombinant = dataset
            .phylogeny
            .get_recombinant_ancestor(&expected_population)?;

        let validate_recombinant = expected_recombinant == *observed_recombinant;

        // ----------------------------------------------------------------
        // Parent Validation

        let observed_parents = &recombination.parents;
        let expected_parents = if let Some(expected_recombinant) = expected_recombinant {
            dataset.phylogeny.get_parents(&expected_recombinant)?
        } else {
            Vec::new()
        };

        let validate_parent = if expected_parents.is_empty() {
            observed_parents.is_empty()
        } else {
            let mut expected_parents_descendants = expected_parents
                .iter()
                .flat_map(|parent| dataset.phylogeny.get_descendants(parent).unwrap())
                .unique()
                .collect_vec();

            expected_parents_descendants.retain(|p| observed_parents.contains(p));
            expected_parents_descendants.len() == observed_parents.len()
        };

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
            if !validate_recombinant {
                validate.details.push(Details::IncorrectRecombinant);
            }
            if !validate_parent {
                // Were parents and breakpoints detected at all?
                if !expected_parents.is_empty() && observed_parents.is_empty() {
                    validate.details.push(Details::NoRecombinationDetected);
                } else {
                    validate.details.push(Details::IncorrectParent);
                }
            }
            if !validate_population {
                validate.details.push(Details::IncorrectPopulation);
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

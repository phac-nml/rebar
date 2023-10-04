use crate::cli::run;
//use crate::phylogeny::Phylogeny;
use color_eyre::eyre::{Report, Result};
//use itertools::Itertools;
use log::debug;

/// Create default SARS-CoV-2 recombinant edge cases.
pub fn default() -> Result<Vec<run::Args>, Report> {
    let mut edge_cases: Vec<run::Args> = Vec::new();

    // // --------------------------------------------------------------------
    // // Designated Recombinants

    // let manual = vec!["XCF", "XCG"].into_iter().map(String::from).collect_vec();

    // for recombinant in &phylogeny.recombinants {
    //     if manual.contains(recombinant) {
    //         continue
    //     }
    //     debug!("Creating auto edge case: {recombinant}");
    //     let mut edge_case = run::Args::default();
    //     edge_case.population = Some(recombinant.to_string());
    //     edge_case.parents = Some(phylogeny.get_parents(recombinant)?);
    //     edge_cases.push(edge_case);
    // }

    // --------------------------------------------------------------------
    // Manual

    // XR: BA.2 and BA.1 with no unique subs from BA.1
    let recombinant = "XR".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let mut edge_case = run::Args::default();
    edge_case.min_subs = 0;
    edge_case.population = Some(recombinant);
    edge_cases.push(edge_case);

    // XCF: XBB and FE.1 (XBB.1.18.1) with no unique subs from XBB
    let recombinant = "XCF".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let mut edge_case = run::Args::default();
    edge_case.min_subs = 0;
    edge_case.min_consecutive = 2;
    edge_case.population = Some(recombinant);
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // XCG: BA.5.2 and XBB.1 with only 2 consecutive bases

    let recombinant = "XCG".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let mut edge_case = run::Args::default();
    edge_case.min_subs = 1;
    edge_case.min_consecutive = 2;
    edge_case.population = Some(recombinant);
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // Finish

    Ok(edge_cases)
}

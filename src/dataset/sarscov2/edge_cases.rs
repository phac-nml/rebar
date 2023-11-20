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
    let edge_case = run::Args {
        min_subs: 0,
        population: Some(recombinant),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // XP: BA.1 and BA.2 with a single sub from BA.2
    let recombinant = "XP".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_consecutive: 1,
        min_length: 1,
        parents: Some(vec!["BA.1.1".to_string(), "BA.2".to_string()]),
        population: Some(recombinant),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // XAZ: BA.5 and BA.2.5 with a single sub from BA.2.5
    let recombinant = "XAZ".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_consecutive: 1,
        min_length: 1,
        min_subs: 0,
        parents: Some(vec!["BA.5".to_string(), "BA.2.5".to_string()]),
        population: Some(recombinant),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // XBK: BA.5.2 and CJ.2 with only 2 consecutive alleles from BA.5.2
    let recombinant = "XBK".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_consecutive: 2,
        population: Some(recombinant),
        parents: Some(vec!["BA.5.2".to_string(), "CJ.1".to_string()]),
        ..Default::default()
    };
    edge_cases.push(edge_case);
    // XBQ: BA.5.2 and CJ.2 with only 2 consecutive alleles from BA.5.2
    let recombinant = "XBQ".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_consecutive: 2,
        population: Some(recombinant),
        parents: Some(vec!["BA.5.2".to_string(), "CJ.1".to_string()]),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // XCF: XBB and FE.1 (XBB.1.18.1) with no unique subs from XBB
    let recombinant = "XCF".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_subs: 0,
        min_consecutive: 2,
        parents: Some(vec!["XBB".to_string(), "FE.1".to_string()]),
        population: Some(recombinant),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // XCG: BA.5.2 and XBB.1 with only 2 consecutive bases

    let recombinant = "XCG".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_subs: 1,
        min_consecutive: 2,
        population: Some(recombinant),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // Finish

    Ok(edge_cases)
}

use crate::cli::run;
use color_eyre::eyre::{Report, Result};
use log::debug;

/// Create default SARS-CoV-2 recombinant edge cases.
pub fn default() -> Result<Vec<run::Args>, Report> {
    let mut edge_cases: Vec<run::Args> = Vec::new();

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

    // --------------------------------------------------------------------
    // XAV: BA.5.1.24 and BA.2 with only two consecutive BA.2 subs

    let recombinant = "XAV".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_consecutive: 2,
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

    // --------------------------------------------------------------------
    // XBZ: EF.1.3 with very small 5' from BA.5.2
    let recombinant = "XBZ".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_length: 200,
        min_consecutive: 2,
        parents: Some(vec!["EF.1.3".to_string(), "BA.5.2".to_string()]),
        population: Some(recombinant),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // XCF: FE.1 (XBB.1.18.1) and XBB.1.5.44 with few consecutive
    let recombinant = "XCF".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_subs: 0,
        min_consecutive: 2,
        parents: Some(vec!["FE.1".to_string(), "XBB".to_string()]),
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
    // XCL: GK.2.1 and FY.4.2 (only pop with C7417T, 2023-11-20)

    let recombinant = "XCL".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        population: Some(recombinant),
        parents: Some(vec!["GK.2.1.1".to_string(), "FY.4.2".to_string()]),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // XCM: XBB.2.3.13 and DV.7.1.2

    let recombinant = "XCM".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        population: Some(recombinant),
        parents: Some(vec!["DV.7.1.2".to_string(), "XBB.2.3.13".to_string()]),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // XCQ: XBB.2.3 and XBB.1.5 with only 2 consecutive

    let recombinant = "XCQ".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_consecutive: 2,
        min_length: 400,
        population: Some(recombinant),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // XDB: XBB.1.16.19 with double reversion assumed from XBB

    let recombinant = "XDB".to_string();
    debug!("Creating manual edge case: {recombinant:?}");
    let edge_case = run::Args {
        min_subs: 0,
        min_consecutive: 2,
        population: Some(recombinant),
        parents: Some(vec!["XBB.1.16.19".to_string(), "XBB".to_string()]),
        ..Default::default()
    };
    edge_cases.push(edge_case);

    // --------------------------------------------------------------------
    // Finish

    Ok(edge_cases)
}

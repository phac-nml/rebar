use crate::cli;
use crate::dataset::attributes::Name;
use crate::utils::table::Table;
use color_eyre::eyre::{Report, Result};
use itertools::Itertools;
use strum::{EnumProperty, IntoEnumIterator};

/// List datasets
pub fn datasets(args: &cli::dataset::list::Args) -> Result<(), Report> {
    // table of name, tag, cli_version
    let mut table = Table::new();
    table.headers = vec![
        "Name",
        "CLI Version",
        "Minimum Tag Date",
        "Maximum Tag Date",
    ]
    .into_iter()
    .map(String::from)
    .collect_vec();

    for name in Name::iter() {
        // Check if this was not the name requested by CLI args
        if let Some(args_name) = &args.name {
            if &name != args_name {
                continue;
            }
        }

        // check if this datset name is actually implemented currently
        if name.get_str("implemented").unwrap_or("false") != "true" {
            continue;
        }

        // Extract compatibility attributes
        let compatibility = name.compatibility()?;

        let cli_version = compatibility.cli.version.unwrap_or(String::new());
        let min_date = if let Some(min_date) = compatibility.dataset.min_date {
            min_date.format("%Y-%m-%d").to_string()
        } else {
            "latest".to_string()
        };
        let max_date = if let Some(max_date) = compatibility.dataset.max_date {
            max_date.format("%Y-%m-%d").to_string()
        } else {
            "latest".to_string()
        };

        // Add to row
        let row = vec![
            name.to_string(),
            cli_version.to_string(),
            min_date.to_string(),
            max_date.to_string(),
        ];
        table.rows.push(row);
    }

    println!("\n{}", table.to_markdown()?);

    Ok(())
}

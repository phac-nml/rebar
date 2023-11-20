use clap::{Args as ClapArgs, Parser};
use color_eyre::eyre::{Report, Result, WrapErr};
use either::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Run on input alignment or dataset population.
#[derive(Clone, Debug, Deserialize, Parser, Serialize)]
#[clap(verbatim_doc_comment)]
pub struct Args {
    /// Dataset directory.
    #[clap(short = 'd', long, required = true)]
    #[serde(
        skip_serializing_if = "Args::is_default_dataset_dir",
        skip_deserializing
    )]
    pub dataset_dir: PathBuf,

    #[command(flatten)]
    #[serde(skip_serializing_if = "Args::is_default_input", skip_deserializing)]
    pub input: Input,

    // Hidden attribute, will be used for edge cases.
    #[arg(hide = true)]
    pub population: Option<String>,

    /// Restrict parent search to just these candidate parents.
    #[arg(long, value_delimiter = ',')]
    pub parents: Option<Vec<String>>,

    /// Remove these populations from the dataset.
    ///
    /// Regardless of whether you use '*' or not, all descendants of the
    /// specified populations will be removed.
    #[arg(short = 'k', long, value_delimiter = ',')]
    pub knockout: Option<Vec<String>>,

    /// Number of bases to mask at the 5' and 3' ends.
    ///
    /// Comma separated.
    #[arg(short = 'm', long)]
    #[arg(long, value_delimiter = ',')]
    pub mask: Vec<usize>,

    /// Maximum number of search iterations to find each parent.
    #[arg(short = 'i', long, default_value_t = Args::default().max_iter)]
    pub max_iter: usize,

    /// Maximum number of parents.
    #[arg(long, default_value_t = Args::default().max_parents)]
    pub max_parents: usize,

    /// Minimum number of parents.
    #[arg(long, default_value_t = Args::default().min_parents)]
    pub min_parents: usize,

    /// Minimum number of consecutive bases in a parental region.
    #[arg(short = 'c', long, default_value_t = Args::default().min_consecutive)]
    pub min_consecutive: usize,

    /// Minimum length of a parental region.
    #[arg(short = 'l', long, default_value_t = Args::default().min_length)]
    pub min_length: usize,

    /// Minimum number of substitutions in a parental region.
    #[arg(short = 's', long, default_value_t = Args::default().min_subs)]
    pub min_subs: usize,

    /// Run a naive search, which does not use information about designated recombinant parents.
    #[arg(short = 'u', long, default_value_t = Args::default().naive)]
    pub naive: bool,

    /// Output directory.
    ///
    /// If the directory does not exist, it will be created.
    #[clap(short = 'o', long, required = true)]
    #[serde(
        skip_serializing_if = "Args::is_default_output_dir",
        skip_deserializing
    )]
    pub output_dir: PathBuf,

    /// Number of CPU threads to use.
    #[clap(short = 't', long, default_value_t = Args::default().threads)]
    #[serde(skip)]
    pub threads: usize,
}

impl Default for Args {
    fn default() -> Self {
        Args {
            dataset_dir: PathBuf::new(),
            input: Input::default(),
            knockout: None,
            mask: vec![100, 200],
            max_iter: 3,
            min_parents: 2,
            max_parents: 2,
            min_consecutive: 3,
            min_length: 500,
            min_subs: 1,
            naive: false,
            output_dir: PathBuf::new(),
            parents: None,
            population: None,
            threads: 1,
        }
    }
}

impl Args {
    pub fn new() -> Self {
        Args {
            dataset_dir: PathBuf::new(),
            input: Input::default(),
            knockout: None,
            mask: vec![0, 0],
            max_iter: 0,
            min_parents: 0,
            max_parents: 0,
            min_consecutive: 0,
            min_length: 0,
            min_subs: 0,
            output_dir: PathBuf::new(),
            parents: None,
            population: None,
            threads: 0,
            naive: false,
        }
    }

    /// Check if input is default.
    pub fn is_default_dataset_dir(path: &Path) -> bool {
        path == Args::default().dataset_dir
    }

    /// Check if input is default.
    pub fn is_default_input(input: &Input) -> bool {
        input == &Args::default().input
    }
    /// Check if output directory is default.
    pub fn is_default_output_dir(path: &Path) -> bool {
        path == Args::default().output_dir
    }

    /// Override Args for edge case handling of particular recombinants.
    pub fn apply_edge_case(&self, new: &Args) -> Result<Args, Report> {
        let mut output = self.clone();

        output.max_iter = new.max_iter;
        output.max_parents = new.max_parents;
        output.min_consecutive = new.min_consecutive;
        output.min_length = new.min_length;
        output.min_subs = new.min_subs;
        output.parents = new.parents.clone();
        output.naive = new.naive;

        Ok(output)
    }

    /// Read args from file.
    ///
    /// Returns an Either with options:
    ///     Left: Args
    ///     Right: Vec<args>
    pub fn read(path: &Path, multiple: bool) -> Result<Either<Args, Vec<Args>>, Report> {
        let input = std::fs::read_to_string(path)
            .wrap_err_with(|| "Failed to read file: {path:?}.")?;
        let output: Vec<Args> = serde_json::from_str(&input)
            .wrap_err_with(|| format!("Failed to parse file: {path:?}"))?;

        if multiple {
            Ok(Right(output))
        } else {
            Ok(Left(output[0].clone()))
        }
    }

    /// Write args to file.
    pub fn write(args: &[Args], path: &Path) -> Result<(), Report> {
        // create file
        let mut file = File::create(path)
            .wrap_err_with(|| format!("Failed to create file: {path:?}"))?;

        // parse to string
        let args = serde_json::to_string_pretty(args)
            .wrap_err_with(|| format!("Failed to parse: {args:?}"))?;

        // write to file
        file.write_all(format!("{}\n", args).as_bytes())
            .wrap_err_with(|| format!("Failed to write file: {path:?}"))?;

        Ok(())
    }
}

#[derive(ClapArgs, Clone, Debug, Deserialize, Serialize, PartialEq)]
#[group(required = true, multiple = true)]
pub struct Input {
    /// Input fasta alignment.
    #[arg(long, value_delimiter = ',')]
    pub populations: Option<Vec<String>>,

    /// Input dataset population.
    #[arg(long)]
    pub alignment: Option<PathBuf>,
}

impl Default for Input {
    fn default() -> Self {
        Self::new()
    }
}

impl Input {
    pub fn new() -> Self {
        Input {
            populations: None,
            alignment: None,
        }
    }
}

use crate::dataset::attributes::Tag;
use crate::utils::remote_file::RemoteFile;
use chrono::Local;
use color_eyre::eyre::{Report, Result, WrapErr};
use indoc::formatdoc;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Create and write Toy1 reference sequence.
pub fn reference(_tag: &Tag, path: &Path) -> Result<RemoteFile, Report> {
    let sequences = formatdoc!(
        "
        >Reference
        AAAAAAAAAAAAAAAAAAAA
        "
    );

    let mut file = File::create(path)
        .wrap_err_with(|| format!("Unable to create file: {path:?}"))?;
    file.write_all(sequences.as_bytes())
        .wrap_err_with(|| format!("Unable to write file: {path:?}"))?;

    let remote_file = RemoteFile {
        local_path: path.to_owned(),
        date_created: Local::now().into(),
        ..Default::default()
    };

    Ok(remote_file)
}

/// Create and write Toy1 populations sequence.
pub fn populations(_tag: &Tag, path: &Path) -> Result<RemoteFile, Report> {
    let sequences = formatdoc!(
        "
        >A
        CCCCCCAACCCCCCCCCCCC
        >B
        TTTTTTTTTTTTTTTTTTAA
        >C
        AAGGGGGGGGGGGGGGGGGG
        >D
        CCCCCCAACCCTTTTTTTAA
        >E
        AAGCCCAACCCTTTTTTTAA
        "
    );

    let mut file = File::create(path)
        .wrap_err_with(|| format!("Unable to create file: {path:?}"))?;
    file.write_all(sequences.as_bytes())
        .wrap_err_with(|| format!("Unable to write file: {path:?}"))?;

    let remote_file = RemoteFile {
        local_path: path.to_owned(),
        date_created: Local::now().into(),
        ..Default::default()
    };

    Ok(remote_file)
}

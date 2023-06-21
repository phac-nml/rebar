use bio::io::fasta;
use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use log::debug;
use std::fs::{remove_file, write, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use tempfile::TempDir;
use zstd::stream::read::Decoder;

pub async fn download_file(
    url: &str,
    output_path: &PathBuf,
    decompress: bool,
) -> Result<(), Report> {
    let ext = Path::new(&url).extension().unwrap().to_str().unwrap();

    let response = reqwest::get(url).await?;
    if response.status() != 200 {
        return Err(eyre!(
            "Unable to download file: {url}\nStatus code {}.",
            response.status()
        ));
    }

    if decompress {
        // Write bytes to a tmp file
        let tmp_dir = TempDir::new()?;
        let tmp_path = PathBuf::from(tmp_dir.path()).join(format!("tmpfile.{ext}"));
        let content = response.bytes().await?;
        write(&tmp_path, content)
            .wrap_err(format!("Unable to write file: {:?}", tmp_path))?;
        decompress_file(&tmp_path, output_path, true)?;
    } else {
        let content = response.text().await?;
        write(output_path, content)
            .wrap_err(format!("Unable to write file: {:?}", output_path))?;
    }

    Ok(())
}

pub fn decompress_file(
    input: &PathBuf,
    output: &PathBuf,
    inplace: bool,
) -> Result<(), Report> {
    let ext = input.extension().unwrap();

    match ext.to_str().unwrap() {
        "zst" => {
            let reader = File::open(input)?;
            let mut decoder = Decoder::new(reader)?;
            let mut buffer = String::new();
            decoder.read_to_string(&mut buffer)?;
            write(output, buffer)
                .wrap_err(format!("Unable to write file: {:?}", output))?;

            if inplace {
                remove_file(input)?;
            }
        }
        _ => return Err(eyre!("Decompression for .{ext:?} is not implemented yet.")),
    };

    Ok(())
}

pub fn identify_diagnostic_mutations(alignment: &Path) -> Result<(), Report> {
    debug!("Loading alignment: {:?}", alignment);
    let alignment_reader = fasta::Reader::from_file(&alignment).expect("Unable to read alignment");
    Ok(())
}

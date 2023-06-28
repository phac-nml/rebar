pub mod constants;
pub mod data;
pub mod text;

use color_eyre::eyre::{eyre, Report, Result};
use color_eyre::Help;
use itertools::Itertools;
use log::debug;
use raqote::*;

pub fn draw_polygon(
    canvas: &mut raqote::DrawTarget,
    x_coords: &[f32],
    y_coords: &[f32],
    fill: &raqote::Source,
    stroke: &raqote::Source,
    stroke_style: &raqote::StrokeStyle,
) -> Result<(), Report> {
    // construct a list of path operations (drawing instructions)
    let mut ops: Vec<PathOp> = Vec::new();
    let winding = raqote::Winding::EvenOdd;

    for (i, it) in x_coords.iter().zip(y_coords.iter()).enumerate() {
        let (x, y) = it;
        let point = Point::new(*x, *y);
        let op = match i {
            0 => raqote::PathOp::MoveTo(point),
            _ => raqote::PathOp::LineTo(point),
        };
        ops.push(op);
    }

    ops.push(raqote::PathOp::Close);

    let path = raqote::Path { ops, winding };

    // fill polygon
    canvas.fill(&path, fill, &DrawOptions::new());
    // outline polygon
    canvas.stroke(&path, stroke, stroke_style, &DrawOptions::new());

    Ok(())
}

pub fn get_base_rgba(base: &String, ref_base: &String, pal_i: usize) -> [u8; 4] {
    // default WHITE
    let mut rgba = [255, 255, 255, 255];

    // convert to char for alphabet comparison
    let base_char = base.chars().next().unwrap();

    // same as reference, light palette
    if base == ref_base {
        rgba = constants::PALETTE_LIGHT[pal_i];
    }
    // mutation, dark palette
    else if constants::ALPHABET.contains(&base_char) {
        rgba = constants::PALETTE_DARK[pal_i];
    }

    rgba
}

pub fn create(
    barcodes_path: &std::path::Path,
    linelist_path: &std::path::Path,
    output_path: &std::path::Path,
) -> Result<(), Report> {
    // ------------------------------------------------------------------------
    // Import Data
    // ------------------------------------------------------------------------

    let mut linelist = data::Table::from_tsv(linelist_path)?;
    let barcodes = data::Table::from_tsv(barcodes_path)?;

    // check for mandatory columns and header pos
    let coord_i = barcodes.header_position("coord")?;
    let origin_i = barcodes.header_position("origin")?;
    let reference_i = barcodes.header_position("Reference")?;

    // get coords
    let coords = barcodes
        .rows
        .iter()
        .map(|row| &row[coord_i])
        .unique()
        .collect_vec();

    // get parents (origins column)
    let parents = barcodes
        .rows
        .iter()
        .map(|row| row[origin_i].to_string())
        .unique()
        .collect_vec();

    if parents.len() > constants::PALETTE_DARK.len() {
        return Err(
            eyre!("There are more parents than colors in the palette!")
            .suggestion(format!("Are you sure you want to plot recombination involving {} parents?", parents.len()))
            .suggestion("If so, please contact the developer to expand the color palette options."));
    }

    // get sequence ids (columns after mandatory cols and parents)
    let sequence_ids = barcodes
        .headers
        .iter()
        .skip(3 + parents.len())
        .cloned()
        .collect_vec();

    // search for sequence_ids in the linelist
    let strain_i = linelist.header_position("strain")?;

    // filter linelist down to just these samples
    linelist.rows = linelist
        .rows
        .into_iter()
        .filter(|row| sequence_ids.contains(&row[strain_i]))
        .collect_vec();

    // make a list of all the populations/sequences we will plot
    // Reference, Parents, Sequences
    let mut populations = vec!["Reference".to_string()];
    populations.extend(parents.clone());
    populations.extend(sequence_ids.clone());

    // ------------------------------------------------------------------------
    // Section
    // ------------------------------------------------------------------------

    // the size of the plot is going to be determined by:
    //   - number of coords in barcodes (width)
    //   - longest sequence_id name (width on left-hand side)
    //   - longest coordinate (height)
    //   - consants::X_INC, which will be width/height white-space for borders
    //   - number of parents and sequences (height)
    //   - ...? amount of white-space on top and bottom

    let num_coords = barcodes.rows.len();
    //let longest_seq_id = sequence_ids.iter().map(|id| id.len()).max()?;
    //println!("{longest_seq_id:?}");

    // longest sequence id (in pixels)
    let text_color = raqote::Color::new(255, 0, 0, 0);
    let longest_sequence_id = sequence_ids
        .iter()
        .map(|id| text::TextProps {
            text: id.to_owned(),
            size: constants::FONT_SIZE,
        })
        .map(|text_props| text::text(&text_props, text_color))
        .map(|text| text.width)
        .max()
        .ok_or_else(|| eyre!("Failed to calculated the maximum sequence ID length"))?;

    let section_gap = constants::X_INC / 2;

    // this is the x position each section will start at
    let section_x = constants::X_INC             // white-space left
        + longest_sequence_id                    // labels on left-hand-side
        + constants::X_INC; // white-space buffer labels-sections

    let canvas_width = section_x                 // x coord of section start
        + (num_coords as i32 * constants::X_INC) // coord boxes
        + constants::X_INC; // white-space right

    let mut section_y = constants::X_INC; // white-space top

    let canvas_height = section_y                        // white-space top
        + constants::X_INC                               // reference bases
        + (constants::X_INC * parents.len() as i32)      // parent bases
        + section_gap                                    // gap between parents and samples
        + (constants::X_INC * sequence_ids.len() as i32) // sequence/sample bases
        + constants::X_INC; // white-space bottom

    println!("{canvas_width} x {canvas_height}");

    // add white space between sub boxes by making them smaller than X_INC
    let sub_box_w = constants::X_INC as f32 * 0.8;

    // ------------------------------------------------------------------------
    // Draw
    // ------------------------------------------------------------------------

    let mut canvas = DrawTarget::new(canvas_width, canvas_height);

    // draw white background
    let mut background = PathBuilder::new();
    background.rect(0., 0., canvas_width as f32, canvas_height as f32);
    let background = background.finish();
    canvas.fill(&background, &constants::WHITE, &DrawOptions::new());

    // ------------------------------------------------------------------------
    // Regions
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // Sub Box
    // ------------------------------------------------------------------------

    debug!("Drawing sub base boxes.");

    section_y += section_gap;

    // iterate through sub coordinates
    for (coord_i, _coord) in coords.iter().enumerate() {
        // absolute x coord
        let x = (section_x + (constants::X_INC * (coord_i as i32))) as f32;

        // reference genome base
        let ref_base = barcodes.rows[coord_i][reference_i].to_string();

        // origins for samples
        let origin = barcodes.rows[coord_i][origin_i].to_string();

        // iterate through samples
        for (pop_i, population) in populations.iter().enumerate() {
            // absolute y coord
            let mut y = (section_y + (constants::X_INC * (pop_i as i32))) as f32;

            // shift samples down, to make gap between parents
            // 1 for reference + variable number of parents
            if pop_i > parents.len() {
                y += section_gap as f32;
            }

            // On the first coord, write pop label
            if coord_i == 0 {
                let text_properties = text::TextProps {
                    text: population.to_string(),
                    size: constants::FONT_SIZE,
                };
                let text_buffer = text::text(&text_properties, text_color);
                let text_x = (section_x - longest_sequence_id - constants::X_INC / 2
                    + (longest_sequence_id - text_buffer.width))
                    as f32;
                let text_y =
                    y + (constants::X_INC as f32 / 2.) - (text_buffer.height as f32 / 2.);
                text_buffer.render(&mut canvas, Point::new(text_x, text_y));
            }

            // adjust box coord based on width/height of sub box
            let box_x = x + (constants::X_INC as f32 / 2. - sub_box_w / 2.);
            let box_y = y + (constants::X_INC as f32 / 2. - sub_box_w / 2.);

            // draw sub box
            let draw_x = vec![box_x, box_x, box_x + sub_box_w, box_x + sub_box_w];
            let draw_y = vec![box_y, box_y + sub_box_w, box_y + sub_box_w, box_y];

            // box color, whether base is reference or not
            let pop_header_i = barcodes.header_position(population)?;
            let pop_base = barcodes.rows[coord_i][pop_header_i].to_string();
            let pop_color: Source;

            // is this the reference genome?
            if &**population == "Reference" {
                pop_color = constants::GREY;
            }
            // is this a parent?
            else if parents.contains(population) {
                let [r, g, b, a] = get_base_rgba(&pop_base, &ref_base, pop_i - 1);
                pop_color = Source::Solid(SolidSource { r, g, b, a });
            }
            // otherwise, it's a sequence
            else {
                // color by origin
                let parent_i = parents.iter().position(|p| *p == origin).unwrap();
                let [r, g, b, a] = get_base_rgba(&pop_base, &ref_base, parent_i);
                pop_color = Source::Solid(SolidSource { r, g, b, a });
            }

            draw_polygon(
                &mut canvas,
                &draw_x,
                &draw_y,
                &pop_color,
                &constants::TRANSPARENT,
                &constants::BASIC_STROKE_STYLE,
            )?;

            // draw sub text
            let pop_header_i = barcodes.header_position(population)?;
            let pop_base = barcodes.rows[coord_i][pop_header_i].to_string();

            let text_properties = text::TextProps {
                text: pop_base,
                size: constants::FONT_SIZE,
            };
            let text_buffer = text::text(&text_properties, text_color);
            let text_x = (box_x + (sub_box_w / 2.)) - ((text_buffer.width / 2) as f32);
            let text_y = (box_y + (sub_box_w / 2.)) - ((text_buffer.height / 2) as f32);
            text_buffer.render(&mut canvas, Point::new(text_x, text_y));
        }
    }

    // ------------------------------------------------------------------------
    // Export
    // ------------------------------------------------------------------------

    canvas.write_png(output_path).unwrap();

    Ok(())
}

pub mod constants;
pub mod text;

use crate::utils;
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

#[allow(unused_variables)]
pub fn create(
    barcodes_path: &std::path::Path,
    linelist_path: &std::path::Path,
    annotations_path: Option<&std::path::Path>,
    output_path: &std::path::Path,
) -> Result<(), Report> {
    // ------------------------------------------------------------------------
    // Import Data
    // ------------------------------------------------------------------------

    // mandatory import data
    let mut linelist = utils::Table::from_tsv(linelist_path)?;
    let barcodes = utils::Table::from_tsv(barcodes_path)?;

    // optional import data
    let mut annotations = utils::Table::new();
    if let Some(annotations_path) = annotations_path {
        annotations = utils::Table::from_tsv(annotations_path)?
    }

    // check for mandatory columns and header pos
    let genome_length_i = linelist.header_position("genome_length")?;
    let breakpoints_i = linelist.header_position("breakpoints")?;

    let coord_i = barcodes.header_position("coord")?;
    let origin_i = barcodes.header_position("origin")?;
    let reference_i = barcodes.header_position("Reference")?;

    // get genome_length, just use first
    let genome_length = linelist
        .rows
        .iter()
        .map(|row| &row[genome_length_i])
        .next()
        .unwrap()
        .parse::<usize>()?;

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

    let section_gap = constants::X_INC;
    let label_gap = constants::X_INC / 2.;

    // this is the x position each section will start at
    let section_x = constants::X_INC             // white-space left
        + longest_sequence_id as f32             // labels on left-hand-side
        + label_gap                              // gap between labels and sections
        + constants::X_INC; // white-space buffer labels-sections

    let canvas_width = section_x                 // x coord of section start
        + (num_coords as f32 * constants::X_INC) // coord boxes
        + constants::X_INC; // white-space right

    let mut section_y = constants::X_INC; // white-space top

    let mut canvas_height = constants::X_INC             // white-space top
        + (constants::X_INC * 2.) + section_gap           // parent regions and text labels
        + constants::X_INC  + section_gap                // guide section
        + constants::X_INC                               // reference bases
        + (constants::X_INC * parents.len() as f32)      // parent bases
        + section_gap                                    // gap between parents and samples
        + (constants::X_INC * sequence_ids.len() as f32) // sequence/sample bases
        + constants::X_INC; // white-space bottom

    // add extra height for optional annotations section
    if !annotations.rows.is_empty() {
        canvas_height += (constants::X_INC * 2.) + section_gap
    }

    debug!("Creating canvas: {canvas_width} x {canvas_height}");

    // add white space between sub boxes by making them smaller than X_INC
    let sub_box_w = constants::X_INC * 0.8;

    // convert genomic coordinates to pixels, based on coords
    let pixels_per_base = (num_coords as f32 * constants::X_INC) / genome_length as f32;

    // ------------------------------------------------------------------------
    // Background
    // ------------------------------------------------------------------------

    let mut canvas = DrawTarget::new(canvas_width as i32, canvas_height as i32);

    // draw white background
    let mut background = PathBuilder::new();
    background.rect(0., 0., canvas_width, canvas_height);
    let background = background.finish();
    canvas.fill(&background, &constants::WHITE, &DrawOptions::new());

    // ------------------------------------------------------------------------
    // Regions
    // ------------------------------------------------------------------------

    debug!("Drawing genomic coordinates guide.");

    // draw section label
    let (label, x) = (String::from("Regions"), section_x);
    let y = section_y + constants::X_INC;
    draw_section_label(&mut canvas, label, text_color, x, y, label_gap)?;

    // draw grey box as background (with cross hatching eventually)
    let box_x = section_x;
    let box_y = section_y + constants::X_INC;
    let box_w = num_coords as f32 * constants::X_INC;
    let box_h = constants::X_INC;

    let draw_x = vec![box_x, box_x, box_x + box_w, box_x + box_w];
    let draw_y = vec![box_y, box_y + box_h, box_y + box_h, box_y];
    draw_polygon(
        &mut canvas,
        &draw_x,
        &draw_y,
        &constants::GREY,
        &constants::TRANSPARENT,
        &constants::BASIC_STROKE_STYLE,
    )?;

    // iterate over parental regions in linelist
    let regions_i = linelist.header_position("regions")?;
    // they *should be all the same, just grab first
    let regions = &linelist
        .rows
        .iter()
        .map(|row| row[regions_i].to_string())
        .next()
        .unwrap();
    // 0-1000|parent1,1000-2000|parent2;
    let regions_split = regions.split(',').collect_vec();

    for (region_i, region) in regions_split.iter().enumerate() {
        // 0-1000|parent1
        let region_parts = region.split('|').collect_vec();
        let parent = region_parts[1];
        let regions_coords = region_parts[0].split('-').collect_vec();

        let mut region_start = regions_coords[0].parse::<usize>()?;
        if region_i == 0 {
            region_start = 0;
        }
        let mut region_end = regions_coords[1].parse::<usize>()?;
        if region_i == regions_split.len() - 1 {
            region_end = genome_length;
        }

        // draw the region box, leave X_INC gap at top for parent text labels
        // convert genomic coordinates to pixel coordinates
        let box_x = section_x + (region_start as f32 * pixels_per_base);
        let box_w = (region_end - region_start) as f32 * pixels_per_base;

        let box_y = section_y + constants::X_INC;
        let box_h = constants::X_INC;

        // region text label
        let text_properties = text::TextProps {
            text: parent.to_string(),
            size: constants::FONT_SIZE,
        };
        let text_buffer = text::text(&text_properties, text_color);
        let text_x = (box_x + (box_w / 2.)) - ((text_buffer.width / 2) as f32);
        let text_y = section_y;
        text_buffer.render(&mut canvas, Point::new(text_x, text_y));

        // region text line
        let draw_x = vec![box_x + (box_w / 2.), box_x + (box_w / 2.)];
        let draw_y = vec![box_y, text_y + text_buffer.height as f32];
        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &constants::TRANSPARENT,
            &constants::BLACK,
            &constants::BASIC_STROKE_STYLE,
        )?;

        // color
        let parent_i = parents.iter().position(|p| *p == parent).unwrap();
        let [r, g, b, a] = constants::PALETTE_DARK[parent_i];
        let color = Source::Solid(SolidSource { r, g, b, a });

        // draw region box
        let draw_x = vec![box_x, box_x, box_x + box_w, box_x + box_w];
        let draw_y = vec![box_y, box_y + box_h, box_y + box_h, box_y];
        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &color,
            &constants::TRANSPARENT,
            &constants::BASIC_STROKE_STYLE,
        )?;
    }

    section_y += constants::X_INC * 2.;

    // ------------------------------------------------------------------------
    // Guide
    // ------------------------------------------------------------------------

    debug!("Drawing genomic guide.");

    section_y += section_gap;

    // draw section label
    let (label, x) = (String::from("Genome"), section_x);
    let y = section_y + constants::X_INC;
    draw_section_label(&mut canvas, label, text_color, x, y, label_gap)?;

    // draw grey box
    let box_x = section_x;
    let box_y = section_y + constants::X_INC;
    let box_w = num_coords as f32 * constants::X_INC;
    let box_h = constants::X_INC;

    let draw_x = vec![box_x, box_x, box_x + box_w, box_x + box_w];
    let draw_y = vec![box_y, box_y + box_h, box_y + box_h, box_y];
    draw_polygon(
        &mut canvas,
        &draw_x,
        &draw_y,
        &constants::GREY,
        &constants::TRANSPARENT,
        &constants::BASIC_STROKE_STYLE,
    )?;

    // ------------------------------------------------------------------------
    // Annotations (Optional)

    // special pallete for annotations, that interleaves dark and light
    // skip colors reserverd for parents, x 2 for interleaved palette
    // todo!() raise error if no colors left, parents consumed it all
    let annot_palette = constants::PALETTE_DARK
        .iter()
        .zip(constants::PALETTE_LIGHT.iter())
        .skip(parents.len())
        .flat_map(|(dark, light)| vec![*dark, *light])
        .collect_vec();

    for (i, row) in annotations.rows.iter().enumerate() {
        let abbrev_i = annotations.header_position("abbreviation")?;
        let start_i = annotations.header_position("start")?;
        let end_i = annotations.header_position("end")?;

        let abbreviation = &annotations.rows[i][abbrev_i];
        let start = annotations.rows[i][start_i].parse::<usize>()?;
        let end = annotations.rows[i][end_i].parse::<usize>()?;

        // use colors from the color palette that are not reserved for pops
        let mut color_i = i;
        if color_i >= annot_palette.len() {
            color_i -= annot_palette.len();
        }
        let [r, g, b, a] = annot_palette[color_i];
        let color = Source::Solid(SolidSource { r, g, b, a });
        // draw the region box, leave X_INC gap at top for annotation labels
        // convert genomic coordinates to pixel coordinates
        let box_x = section_x + (start as f32 * pixels_per_base);
        let box_y = section_y + constants::X_INC;
        let box_w = (end - start) as f32 * pixels_per_base;
        let box_h = constants::X_INC;
        let draw_x = vec![box_x, box_x, box_x + box_w, box_x + box_w];
        let draw_y = vec![box_y, box_y + box_h, box_y + box_h, box_y];

        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &color,
            &constants::TRANSPARENT,
            &constants::BASIC_STROKE_STYLE,
        )?;

        // text label
        let text_properties = text::TextProps {
            text: abbreviation.to_string(),
            size: constants::FONT_SIZE - 5,
        };
        let text_buffer = text::text(&text_properties, text_color);
        let text_x = (box_x + (box_w / 2.)) - ((text_buffer.width / 2) as f32);
        // interleave height of adjacent labels
        let text_y = if let 0 = i % 2 {
            section_y
        } else {
            section_y - (constants::X_INC / 2.)
        };
        text_buffer.render(&mut canvas, Point::new(text_x, text_y));

        // text line
        let draw_x = vec![box_x + (box_w / 2.), box_x + (box_w / 2.)];
        let draw_y = vec![box_y, text_y + text_buffer.height as f32];
        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &constants::TRANSPARENT,
            &constants::BLACK,
            &constants::BASIC_STROKE_STYLE,
        )?;
    }

    // ------------------------------------------------------------------------
    // Sub markers

    // draw coord black lines
    for coord in coords.iter() {
        // convert genomic coord to numeric then to pixels
        let coord = coord.parse::<f32>().unwrap();

        let line_x = section_x + (coord * pixels_per_base);
        let line_y1 = section_y + constants::X_INC;
        let line_y2 = section_y + constants::X_INC * 2.;
        let draw_x = vec![line_x, line_x];
        let draw_y = vec![line_y1, line_y2];
        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &constants::TRANSPARENT,
            &constants::BLACK,
            &constants::BASIC_STROKE_STYLE,
        )?;
    }

    section_y += constants::X_INC * 2.;

    // ------------------------------------------------------------------------
    // Guide to Sub Polyons
    // ------------------------------------------------------------------------

    // draw coord black lines
    for (i, coord) in coords.iter().enumerate() {
        // convert genomic coord to numeric then to pixels
        let coord = coord.parse::<f32>().unwrap();

        // x coordinate for top of triangle, connects with guide
        let guide_x = section_x + (coord * pixels_per_base);
        let guide_y = section_y;
        // x coordinates for bottom of triangle, connects with sub bases
        // adjust based on sub_box_w
        let sub_x_buff = (constants::X_INC - sub_box_w) / 2.;
        let sub_x1 = section_x + (i as f32 * constants::X_INC) + sub_x_buff;
        let sub_x2 = section_x + ((i + 1) as f32 * (constants::X_INC)) - sub_x_buff;
        let sub_y1 = section_y + (constants::X_INC * 3.) + sub_x_buff;
        let sub_y2 = sub_y1;
        // Draw triangle from guide to top row of subs
        let draw_x = vec![guide_x, sub_x1, sub_x2];
        let draw_y = vec![guide_y, sub_y1, sub_y2];

        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &constants::LIGHT_GREY,
            &constants::TRANSPARENT,
            &constants::BASIC_STROKE_STYLE,
        )?;
    }

    // ------------------------------------------------------------------------
    // Axis coordinates (on top of sub triangle polygons)

    // check how long the longest coord label is (in pixels)
    let longest_coord = coords
        .iter()
        .map(|c| text::TextProps {
            text: (*c).clone(),
            size: constants::FONT_SIZE - 5,
        })
        .map(|text_props| text::text(&text_props, text_color))
        .map(|text| text.width)
        .max()
        .ok_or_else(|| eyre!("Failed to calculated the maximum coord length"))?;
    // how many x_inc are needed for it, add 1 extra x_inc for buffer
    let longest_coord_x_inc = (longest_coord as f32 / constants::X_INC).ceil() + 1.0;
    // maximum number of coord labels we can fit
    let max_num_coords = num_coords as f32 / longest_coord_x_inc;
    // calculate interval, round up to next pretty number (ex. 500)
    let coord_interval =
        (((genome_length as f32 / max_num_coords) / 500.).ceil() * 500.) as usize;

    let mut ax_coords = (0..genome_length).step_by(coord_interval).collect_vec();
    ax_coords.push(genome_length);

    for coord in ax_coords {
        // convert genomic units to pixel coords
        let coord_x = section_x + (coord as f32 * pixels_per_base);

        // draw line
        let line_y1 = section_y;
        let line_y2 = line_y1 + (constants::X_INC / 4.);
        let draw_x = vec![coord_x, coord_x];
        let draw_y = vec![line_y1, line_y2];
        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &constants::TRANSPARENT,
            &constants::BLACK,
            &constants::BASIC_STROKE_STYLE,
        )?;

        // draw text
        let text_properties = text::TextProps {
            text: coord.to_string(),
            size: constants::FONT_SIZE - 5,
        };
        let text_buffer = text::text(&text_properties, text_color);
        let text_x = coord_x - (text_buffer.width / 2) as f32;
        let text_y = line_y2;
        text_buffer.render(&mut canvas, Point::new(text_x, text_y));
    }

    // ------------------------------------------------------------------------
    // Breakpoints
    // ------------------------------------------------------------------------

    // draw section label
    let (label, x) = (String::from("Breakpoints"), section_x);
    let y = section_y + constants::X_INC;
    draw_section_label(&mut canvas, label, text_color, x, y, label_gap)?;

    let breakpoints = linelist
        .rows
        .iter()
        .flat_map(|row| row[breakpoints_i].split(','))
        .unique()
        .collect_vec();

    let dash_stroke_style: StrokeStyle = StrokeStyle {
        cap: LineCap::Square,
        join: LineJoin::Miter,
        width: constants::LINE_WIDTH,
        miter_limit: 2.,
        dash_array: vec![8., 20.],
        dash_offset: 0.,
    };

    for (i, breakpoint) in breakpoints.iter().enumerate() {
        // get the region start/end
        let breakpoint_parts = breakpoint.split('-').collect_vec();
        let prev_region_end = breakpoint_parts[0].parse::<f32>()? - 1.;
        let next_region_start = breakpoint_parts[1].parse::<f32>()? + 1.;
        //println!("{prev_region_end} {next_region_start}");

        // which subs does this fall between
        let coord_prev_i = coords
            .iter()
            .position(|c| **c == prev_region_end.to_string())
            .unwrap();
        let coord_next_i = coords
            .iter()
            .position(|c| **c == next_region_start.to_string())
            .unwrap();
        //println!("\t{coord_prev_i} {coord_next_i}");

        // middle will depend on breakpoints uncertainy
        let line_x = section_x
            + (coord_prev_i + 1 + (coord_next_i - coord_prev_i) / 2) as f32
                * constants::X_INC;
        // top of the line will be the sme
        let line_y1 = if let 0 = i % 2 {
            section_y + (constants::X_INC * 2.)
        } else {
            section_y + constants::X_INC
        };

        // option 1: draw line if coords right next to each other
        if coord_next_i - coord_prev_i == 1 {
            // the line bottom is alll the way at the bottom of ref, parents, sample
            let sub_y_buff = (constants::X_INC - sub_box_w) / 2.;
            let line_y2 = section_y
                + (constants::X_INC * 4.) // this height of the breakpoints section
                + (constants::X_INC * populations.len() as f32)
                - sub_y_buff;
            let draw_x = vec![line_x, line_x];
            let draw_y = vec![line_y1, line_y2];
            draw_polygon(
                &mut canvas,
                &draw_x,
                &draw_y,
                &constants::TRANSPARENT,
                &constants::BLACK,
                &dash_stroke_style,
            )?;
        }
        // option 2: draw box around coords if multiple
        else {
            todo!();
        }

        // breakpoint label, prep buffer, draw box background first before render
        let label = format!("Breakpoint #{}", i + 1);
        let text_properties = text::TextProps {
            text: label,
            size: constants::FONT_SIZE,
        };
        let text_buffer = text::text(&text_properties, text_color);
        let text_x = line_x - ((text_buffer.width / 2) as f32);
        let text_y = line_y1;

        // draw breakpoint box background, add several pixels for buffer
        let box_x = line_x - (text_buffer.width as f32 / 2.) - 5.;
        let box_y = line_y1 - 5.;
        let box_w = text_buffer.width as f32 + 10.;
        let box_h = text_buffer.height as f32 + 10.;
        let draw_x = vec![box_x, box_x, box_x + box_w, box_x + box_w];
        let draw_y = vec![box_y, box_y + box_h, box_y + box_h, box_y];

        draw_polygon(
            &mut canvas,
            &draw_x,
            &draw_y,
            &constants::WHITE,
            &constants::BLACK,
            &constants::BASIC_STROKE_STYLE,
        )?;

        // render breakpoint label
        text_buffer.render(&mut canvas, Point::new(text_x, text_y));
    }

    section_y += constants::X_INC * 2.;

    // ------------------------------------------------------------------------
    // Sub Box
    // ------------------------------------------------------------------------

    debug!("Drawing sub base boxes.");

    section_y += section_gap;

    // iterate through sub coordinates
    for (coord_i, _coord) in coords.iter().enumerate() {
        // absolute x coord
        let x = section_x + (constants::X_INC * coord_i as f32);

        // reference genome base
        let ref_base = barcodes.rows[coord_i][reference_i].to_string();

        // origins for samples
        let origin = barcodes.rows[coord_i][origin_i].to_string();

        // iterate through samples
        for (pop_i, population) in populations.iter().enumerate() {
            // absolute y coord
            let mut y = section_y + (constants::X_INC * pop_i as f32);

            // shift samples down, to make gap between parents
            // 1 for reference + variable number of parents
            if pop_i > parents.len() {
                y += section_gap;
            }

            // On the first coord, write pop label.
            // pop label has same pos as section_label, use that function
            if coord_i == 0 {
                let (label, x, y) = (String::from(population), section_x, y);
                draw_section_label(&mut canvas, label, text_color, x, y, label_gap)?;
            }

            // adjust box coord based on width/height of sub box
            let box_x = x + (constants::X_INC / 2.) - (sub_box_w / 2.);
            let box_y = y + (constants::X_INC / 2.) - (sub_box_w / 2.);

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

pub fn draw_section_label(
    canvas: &mut raqote::DrawTarget,
    label: String,
    text_color: raqote::Color,
    section_x: f32,
    section_y: f32,
    label_gap: f32,
) -> Result<(), Report> {
    // draw section label
    let text_properties = text::TextProps {
        text: label,
        size: constants::FONT_SIZE,
    };
    let text_buffer = text::text(&text_properties, text_color);
    let text_x = section_x - label_gap - text_buffer.width as f32;
    let text_y = section_y + (constants::X_INC / 2.) - (text_buffer.height as f32 / 2.);
    text_buffer.render(canvas, Point::new(text_x, text_y));

    Ok(())
}

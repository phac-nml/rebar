pub mod basic;
pub mod constants;
pub mod text;

use crate::utils;
use color_eyre::eyre::{eyre, Report, Result};
use color_eyre::Help;
use itertools::Itertools;
use log::debug;
use raqote::*;
use usvg::{NodeExt, TextToPath, TreeTextToPath, TreeWriting};

/// possible output extension
#[derive(Debug, Clone)]
pub enum Ext {
    Png,
    Svg,
}

impl Ext {
    pub fn from_str(ext: &str) -> Result<Self, Report> {
        match ext {
            "png" => Ok(Ext::Png),
            "svg" => Ok(Ext::Svg),
            _ => Err(eyre!("Plot extension .{ext:?} is not implemented yet.")),
        }
    }
}

pub fn draw_polygon(
    utree: &mut usvg::Tree,
    x_coords: &[f32],
    y_coords: &[f32],
    fill: Option<usvg::Fill>,
    stroke: Option<usvg::Stroke>,
) -> Result<(), Report> {
    let mut path = usvg::tiny_skia_path::PathBuilder::new();

    for (i, it) in x_coords.iter().zip(y_coords.iter()).enumerate() {
        let (x, y) = it;
        match i {
            0 => path.move_to(*x, *y),
            _ => path.line_to(*x, *y),
        }
    }
    path.close();
    let path = path.finish().unwrap();

    // convert from tiny_skia to usvg Path
    let mut path = usvg::Path::new(std::rc::Rc::new(path));

    // set aesthetic attributes
    path.fill = fill;
    path.stroke = stroke;

    // add this path to the svg node
    utree.root.append_kind(usvg::NodeKind::Path(path));

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
    // check if output_path extension
    let output_ext = output_path.extension().unwrap().to_str().unwrap();
    let output_ext = Ext::from_str(output_ext)?;

    let mut fontdb = usvg::fontdb::Database::new();
    fontdb.load_system_fonts();

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
    // Setup
    // ------------------------------------------------------------------------

    // the size of the plot is going to be determined by:
    //   - number of coords in barcodes (width)
    //   - longest sequence_id name (width on left-hand side)
    //   - longest coordinate (height)
    //   - consants::X_INC, which will be width/height white-space for borders
    //   - number of parents and sequences (height)
    //   - ...? amount of white-space on top and bottom

    let num_coords = barcodes.rows.len();

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
    let label_gap = constants::X_INC / 2.0;

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
        + (constants::X_INC * 2.0) + section_gap         // parent regions and text labels
        + constants::X_INC  + section_gap                // guide section
        + constants::X_INC                               // reference bases
        + (constants::X_INC * parents.len() as f32)      // parent bases
        + section_gap                                    // gap between parents and samples
        + (constants::X_INC * sequence_ids.len() as f32) // sequence/sample bases
        + constants::X_INC; // white-space bottom

    // add extra height for optional annotations section
    if !annotations.rows.is_empty() {
        canvas_height += (constants::X_INC * 2.0) + section_gap;
    }

    // add white space between sub boxes by making them smaller than X_INC
    let sub_box_w = constants::X_INC as f32 * 0.8;

    // convert genomic coordinates to pixels, based on coords
    let pixels_per_base = (num_coords as f32 * constants::X_INC) / genome_length as f32;

    // ------------------------------------------------------------------------
    // Canvas Background
    // ------------------------------------------------------------------------

    debug!("Creating canvas: {canvas_width} x {canvas_height}");

    let size = usvg::Size::from_wh(canvas_width, canvas_height).unwrap();
    // create a default usvg tree of the specified size
    let mut utree = usvg::Tree {
        size,
        view_box: usvg::ViewBox {
            rect: size.to_non_zero_rect(0.0, 0.0),
            aspect: usvg::AspectRatio::default(),
        },
        root: usvg::Node::new(usvg::NodeKind::Group(usvg::Group::default())),
    };

    // fill with opaque white background
    let fill = usvg::Fill {
        paint: usvg::Paint::Color(usvg::Color::new_rgb(255, 255, 255)),
        opacity: usvg::tiny_skia_path::NormalizedF32::new(1.).unwrap(),
        rule: usvg::FillRule::NonZero,
    };
    let draw_x = vec![0., 0., canvas_width, canvas_width];
    let draw_y = vec![0., canvas_height, canvas_height, 0.];
    draw_polygon(&mut utree, &draw_x, &draw_y, Some(fill), None)?;

    // ------------------------------------------------------------------------
    // Regions
    // ------------------------------------------------------------------------

    debug!("Drawing genomic coordinates guide.");

    // draw section label
    // let (label, x) = (String::from("Regions"), section_x);
    // let y = section_y + (constants::X_INC as f32 * 1.5) as i32;
    // draw_section_label(&mut canvas, label, text_color, x, y, label_gap)?;

    // draw grey box as background (with cross hatching eventually)
    let [r, g, b] = constants::GREY;
    let fill = usvg::Fill {
        paint: usvg::Paint::Color(usvg::Color::new_rgb(r, g, b)),
        opacity: usvg::tiny_skia_path::NormalizedF32::new(1.).unwrap(),
        rule: usvg::FillRule::NonZero,
    };
    let box_x = section_x;
    let box_y = section_y + constants::X_INC;
    let box_w = num_coords as f32 * constants::X_INC;
    let box_h = constants::X_INC;
    let draw_x = vec![box_x, box_x, box_x + box_w, box_x + box_w];
    let draw_y = vec![box_y, box_y + box_h, box_y + box_h, box_y];
    draw_polygon(&mut utree, &draw_x, &draw_y, Some(fill), None)?;

    // test text draw
    let label = "XBB.1.16.1".to_string();

    // absolute position of label
    let chunks = vec![usvg::TextChunk {
        x: Some(section_x),
        y: Some((box_y + box_h / 2.) + 10.),
        anchor: usvg::TextAnchor::End,
        spans: vec![basic::TextSpan::default().default],
        text_flow: usvg::TextFlow::Linear,
        text: label.clone(),
    }];

    // relative position of each char in string
    // hold x, y, and dy constant. just update  dx
    let [x, y, dy] = [Some(0.), Some(0.), Some(0.)];
    let positions = label
        .chars()
        .map(|c| usvg::CharacterPosition {
            x,
            y,
            dy,
            dx: Some(constants::LETTER_GAP),
        })
        .collect_vec();

    let text = usvg::Text {
        id: "test".to_string(),
        transform: usvg::Transform::default(),
        rendering_mode: usvg::TextRendering::OptimizeLegibility,
        rotate: vec![0.0],
        writing_mode: usvg::WritingMode::LeftToRight,
        // individual char positions
        positions,
        chunks,
    };
    // convert to path
    debug!("{text:?}");
    //utree.root.append_kind(usvg::NodeKind::Text(text));

    // ------------------------------------------------------------------------
    // Export
    // ------------------------------------------------------------------------

    //let mut text_nodes = Vec::new();
    // We have to update text nodes in clipPaths, masks and patterns as well.
    // for node in utree.root.descendants() {
    //     if let usvg::NodeKind::Text(_) = *node.borrow() {

    //         let mut new_node = None;
    //         if let usvg::NodeKind::Text(ref text) = *node.borrow() {
    //             debug!("text: {text:?}");
    //             let mut absolute_ts = node.parent().unwrap().abs_transform();
    //             debug!("absolute_ts: {absolute_ts:?}");
    //             absolute_ts = absolute_ts.pre_concat(text.transform);
    //             debug!("absolute_ts: {absolute_ts:?}");
    //             new_node = text.convert(&fontdb, absolute_ts);
    //             debug!("new_node: {new_node:?}");
    //         }

    //     }
    // }

    // for child in utree.root.children() {
    //     debug!("child: {child:?}");
    // }

    // convert text to paths
    utree.convert_text(&fontdb);

    match output_ext {
        // export to svg, simple write of tree to str
        Ext::Svg => {
            // convert utree to string, then write to file
            let xml_options = usvg::XmlOptions::default();
            let utree_data = utree.to_string(&xml_options);
            std::fs::write(output_path, utree_data);
        }
        // export to png, convert tree to img bytes
        Ext::Png => {
            // convert utree to a render tree (rtree)
            let rtree = resvg::Tree::from_usvg(&utree);

            // construct a pixmap to hold bytes
            let pixmap_size = rtree.size.to_int_size();
            let mut pixmap =
                resvg::tiny_skia::Pixmap::new(pixmap_size.width(), pixmap_size.height())
                    .unwrap();

            // render to file
            rtree.render(usvg::Transform::default(), &mut pixmap.as_mut());
            pixmap.save_png(output_path).unwrap();
        }
    }

    Ok(())
}

pub fn draw_section_label(
    canvas: &mut raqote::DrawTarget,
    label: String,
    text_color: raqote::Color,
    section_x: i32,
    section_y: i32,
    label_gap: i32,
) -> Result<(), Report> {
    // draw section label
    let text_properties = text::TextProps {
        text: label,
        size: constants::FONT_SIZE,
    };
    let text_buffer = text::text(&text_properties, text_color);
    let text_x = (section_x - label_gap - text_buffer.width) as f32;
    let text_y = section_y as f32 + (constants::X_INC as f32 / 2.)
        - (text_buffer.height as f32 / 2.);
    text_buffer.render(canvas, Point::new(text_x, text_y));

    Ok(())
}

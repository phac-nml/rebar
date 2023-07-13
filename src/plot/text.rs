use color_eyre::eyre::{eyre, Report, Result, WrapErr};
use image::{imageops, ImageBuffer, Rgba};
use itertools::Itertools;
//use log::debug;
use std::path::Path;

#[derive(Debug)]
pub enum HorizontalAlignment {
    Left,
    Center,
    Right,
}

#[derive(Debug)]
pub enum VerticalAlignment {
    Top,
    Center,
    Bottom,
}

/// Load font from file path
pub fn load_font(path: &Path) -> Result<rusttype::Font, Report> {
    let font_bytes = std::fs::read(path).wrap_err_with(|| {
        format!("Could not load font from file path: {}", path.display())
    })?;
    let font = rusttype::Font::try_from_vec(font_bytes)
        .ok_or_else(|| eyre!("Could not convert file to Font: {}", path.display()))?;

    Ok(font)
}

/// Convert text string to an image::ImageBuffer
pub fn to_image(
    text: &str,
    font_path: &Path,
    font_size: f32,
    color: &image::Rgba<u8>,
) -> Result<ImageBuffer<Rgba<u8>, Vec<u8>>, Report> {
    // get rgba channels of text color
    let [r, g, b, a] = [color.0[0], color.0[1], color.0[2], color.0[3]];

    // load font from file path
    let font = load_font(font_path)?;

    // ------------------------------------------------------------------------
    // Image Dimensions

    // Set font size and get metrics
    let scale = rusttype::Scale::uniform(font_size);
    let metrics = font.v_metrics(scale);

    // layout the glyphs in the text horizontally
    let glyphs: Vec<_> = font
        .layout(text, scale, rusttype::point(0., 0. + metrics.ascent))
        .collect();

    // get output image height from the font metrics, since height is only dependent on font
    let height = (metrics.ascent - metrics.descent).ceil();

    // Get output image widths from the pixel bounding boxes, since width is dependent
    // on font + the text to write (for horizontal layout)
    let min_x = glyphs
        .iter()
        .map(|g| {
            if let Some(bounding_box) = g.pixel_bounding_box() {
                bounding_box.min.x
            } else {
                0
            }
        })
        .min()
        .unwrap();

    let max_x = glyphs
        .iter()
        .map(|g| {
            if let Some(bounding_box) = g.pixel_bounding_box() {
                bounding_box.max.x
            } else {
                0
            }
        })
        .max()
        .unwrap();

    let width = if min_x >= 0 { max_x } else { max_x - min_x };

    //debug!("width: {width}, height: {height}");

    // ------------------------------------------------------------------------
    // Rasterize Glyphs

    // construct an image buffer to hold text pixels
    let mut image_buffer =
        ImageBuffer::<Rgba<u8>, Vec<_>>::new(width as u32, height as u32);
    let default_pixel: Rgba<u8> = Rgba([0, 0, 0, 0]);

    // iterate through each glyph ('letter')
    for glyph in glyphs {
        //debug!("glyph: {glyph:?}");

        if let Some(bounding_box) = glyph.pixel_bounding_box() {
            //debug!("\tbounding_box: {bounding_box:?}");

            // rasterize each glyph, by iterating through the pixels
            // x, y are relative to bounding box, v is 'coverage'
            glyph.draw(|x, y, v| {
                //debug!("\t\tx: {x}, y: {y}, v: {v}");
                let y = y as i32 + bounding_box.min.y;

                // sometimes x bounding box is negative, because kerning is applied
                // ex. the letter 'T' in isolation
                // in this case, force 0 to be the start point
                let x = if bounding_box.min.x >= 0 {
                    x as i32 + bounding_box.min.x
                } else {
                    x as i32
                };

                //debug!("\t\tx: {x}, y: {y}, v: {v}");

                // construct a pixel
                let pixel = Rgba([
                    (r as f32 * v) as u8,
                    (g as f32 * v) as u8,
                    (b as f32 * v) as u8,
                    (a as f32 * v) as u8,
                ]);
                // add pixel to image buffer, if that pixel is still the default
                if image_buffer.get_pixel(x as u32, y as u32) == &default_pixel {
                    image_buffer.put_pixel(x as u32, y as u32, pixel);
                }
            });
        }
    }

    Ok(image_buffer)
}

/// Convert image::ImageBuffer data text to u32 for raqote
pub fn to_raqote_data(
    image: &ImageBuffer<Rgba<u8>, Vec<u8>>,
) -> Result<Vec<u32>, Report> {
    // convert from u8 to u32
    let data = image
        .pixels()
        .map(|rgba| u32::from_be_bytes([rgba.0[3], rgba.0[0], rgba.0[1], rgba.0[2]]))
        .collect_vec();

    Ok(data)
}

//#[derive(Debug)]
pub struct DrawRaqoteArgs<'canvas> {
    pub canvas: &'canvas mut raqote::DrawTarget,
    pub text: String,
    pub font_path: std::path::PathBuf,
    pub font_size: f32,
    pub color: image::Rgba<u8>,
    pub x: f32,
    pub y: f32,
    pub vertical_alignment: VerticalAlignment,
    pub horizontal_alignment: HorizontalAlignment,
    pub rotate: u32,
}

impl<'canvas> DrawRaqoteArgs<'canvas> {
    pub fn from_canvas(canvas: &'canvas mut raqote::DrawTarget) -> Self {
        DrawRaqoteArgs {
            canvas,
            text: String::new(),
            font_path: std::path::PathBuf::from(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            ),
            font_size: 12.0,
            color: image::Rgba([0, 0, 0, 255]),
            x: 0.0,
            y: 0.0,
            vertical_alignment: VerticalAlignment::Top,
            horizontal_alignment: HorizontalAlignment::Left,
            rotate: 0,
        }
    }
}

/// Draw text string onto raqote canvas.
pub fn draw_raqote(
    args: &mut DrawRaqoteArgs,
) -> Result<ImageBuffer<Rgba<u8>, Vec<u8>>, Report> {
    let image = to_image(&args.text, &args.font_path, args.font_size, &args.color)?;

    // optional rotate
    let image = match args.rotate {
        90 => imageops::rotate90(&image),
        180 => imageops::rotate180(&image),
        270 => imageops::rotate270(&image),
        _ => image,
    };

    let data = to_raqote_data(&image)?;

    let x = match args.horizontal_alignment {
        HorizontalAlignment::Left => args.x,
        HorizontalAlignment::Center => args.x - (image.width() as f32 / 2.),
        HorizontalAlignment::Right => args.x - image.width() as f32,
    };

    let y = match args.vertical_alignment {
        VerticalAlignment::Top => args.y,
        VerticalAlignment::Center => args.y - (image.height() as f32 / 2.),
        VerticalAlignment::Bottom => args.y - image.height() as f32,
    };

    let point = raqote::Point::new(x, y);
    let draw_image = raqote::Image {
        width: image.width() as i32,
        height: image.height() as i32,
        data: &data,
    };
    args.canvas
        .draw_image_at(point.x, point.y, &draw_image, &raqote::DrawOptions::new());

    Ok(image)
}

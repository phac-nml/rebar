// Author: @J-Cake
// Source: https://gist.github.com/J-Cake/ddccf99d3f7d6fc947fc60204aa41e09
// Edited by: Katherine Eaton
// Last edited: 2023-06-28

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct TextRenderBuffer {
    pub width: i32,
    pub height: i32,
    pub baseline: i32,
    pub data: Vec<u32>,
    pub text: String,
}

impl TextRenderBuffer {
    pub fn to_image(&self) -> raqote::Image {
        raqote::Image {
            width: self.width,
            height: self.height,
            data: &self.data,
        }
    }

    pub fn render(&self, ctx: &mut raqote::DrawTarget, point: raqote::Point) -> &Self {
        ctx.draw_image_at(
            point.x,
            point.y,
            &self.to_image(),
            &raqote::DrawOptions {
                blend_mode: raqote::BlendMode::SrcAtop,
                alpha: 1.,
                antialias: raqote::AntialiasMode::Gray,
            },
        );

        self
    }

    // pub fn measure(&self) -> Size {
    //     return Size::new(self.width as f64, self.height as f64, self.baseline as f64);
    // }
}

impl std::fmt::Debug for TextRenderBuffer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "RenderBuffer {{ {}x{}_{} '{}' }}",
            self.width, self.height, self.baseline, self.text
        )
        .unwrap();
        Ok(())
    }
}

pub(crate) fn get_font<'a>() -> rusttype::Font<'a> {
    // TODO: Replace with dynamic font-lookup
    return rusttype::Font::try_from_vec(
        std::fs::read("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf").unwrap(),
    )
    .unwrap();
    //return rusttype::Font::try_from_vec(std::fs::read("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf").unwrap()).unwrap();
}

#[derive(Eq, PartialEq, Debug, Clone, Hash)]
pub struct TextProps {
    pub text: String,
    pub size: i32,
}

unsafe impl Send for TextProps {}

unsafe impl Send for TextRenderBuffer {}

fn render_text(props: TextProps, colour: raqote::Color) -> TextRenderBuffer {
    let font = get_font();

    //println!("props: {props:?}");
    // the scaling factor depends on the props.size?
    //let size = rusttype::Scale::uniform(props.size as f32 * 1.333f32);
    let size = rusttype::Scale::uniform(props.size as f32);
    let metrics = font.v_metrics(size);

    let glyphs: Vec<_> = font
        .layout(&props.text, size, rusttype::point(0., 0. + metrics.ascent))
        .collect();

    let mut image = {
        let width: i32 = {
            let min_x = glyphs
                .first()
                .map(|g| g.pixel_bounding_box().unwrap().min.x)
                .unwrap();
            let max_x = glyphs
                .last()
                .map(|g| g.pixel_bounding_box().unwrap().max.x)
                .unwrap();
            //max_x - min_x + 1
            max_x - min_x + 2
        };
        let height: i32 = (metrics.ascent - metrics.descent).ceil() as i32;

        TextRenderBuffer {
            width,
            height,
            text: props.text.clone(),
            baseline: metrics.ascent as i32,
            data: vec![0u32; (width * height) as usize],
        }
    };

    // Loop through the glyphs in the text, positing each one on a line
    for (_i, glyph) in glyphs.iter().enumerate() {
        //println!("\t{i}: {glyph:?}");

        if let Some(bounding_box) = glyph.pixel_bounding_box() {
            //println!("\t\tbounding_box: {bounding_box:?}");

            glyph.draw(|x, y, v| {
                // original, usize conversion causes overflow, when bounding_box.min.x is < 0
                //let index = ((y as usize + bounding_box.min.y as usize) * image.width as usize) + (x as usize + bounding_box.min.x as usize);
                // custom, do usize conversion after

                let index = (((y as isize + bounding_box.min.y as isize)
                    * image.width as isize)
                    + (x as isize + bounding_box.min.x as isize))
                    as usize;
                // let index = if bounding_box.min.x < 0 {
                //     (((y as isize + bounding_box.min.y as isize) * image.width as isize) + (x as isize + bounding_box.min.x as isize)) as usize
                // } else {
                //     ((y as usize + bounding_box.min.y as usize) * image.width as usize) + (x as usize + bounding_box.min.x as usize)
                // };
                //println!("\t\t\tx: {x}, y: {y}, v: {v}, index: {index}");

                image.data[index] = u32::from_be_bytes([
                    (colour.a() as f32 * v) as u8,
                    (colour.r() as f32 * v) as u8,
                    (colour.g() as f32 * v) as u8,
                    (colour.b() as f32 * v) as u8,
                ]);
            });
        }
    }

    image
}

lazy_static::lazy_static!(static ref CACHE: std::sync::Mutex<std::collections::HashMap<TextProps, std::sync::Arc<TextRenderBuffer>>> = std::sync::Mutex::new(std::collections::HashMap::new()););
pub(crate) fn text(
    props: &TextProps,
    colour: raqote::Color,
) -> std::sync::Arc<TextRenderBuffer> {
    let mut cache = CACHE.lock().unwrap();

    if let Some(texture) = cache.get(props) {
        return texture.clone();
    }

    let buffer = std::sync::Arc::new(render_text(props.clone(), colour));

    cache.insert(props.clone(), buffer);
    return cache.get(props).unwrap().clone();
}

// pub struct Size {
//     pub width: f64,
//     pub height: f64,
//     pub baseline: f64,
// }

// impl Size {
//     pub fn new(width: f64, height: f64, baseline: f64) -> Self {
//         Size { width, height, baseline }
//     }
// }

// fn measure_text(props: &TextProps) -> Size {
//     let font = get_font();
//     let size = rusttype::Scale::uniform(props.size as f32 * 1.333f32);
//     let metrics = font.v_metrics(size);

//     let glyphs: Vec<_> = font
//         .layout(&props.text, size, rusttype::point(0., 0. + metrics.ascent))
//         .collect();

//     return Size::new({
//                          let min_x = glyphs
//                              .first()
//                              .map(|g| g.pixel_bounding_box().unwrap().min.x)
//                              .unwrap();
//                          let max_x = glyphs
//                              .last()
//                              .map(|g| g.pixel_bounding_box().unwrap().max.x)
//                              .unwrap();
//                          (max_x - min_x + 1) as f64
//                      }, (metrics.ascent - metrics.descent).ceil() as f64, metrics.ascent as f64);
// }

// pub(crate) fn measure(props: &TextProps) -> Size {
//     let cache = CACHE.lock().unwrap();

//     if let Some(texture) = cache.get(props) {
//         return Size::new(texture.width as f64, texture.height as f64, texture.baseline as f64);
//     }

//     return measure_text(props)
// }

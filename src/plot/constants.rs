use raqote::*;

pub const ALPHABET: [char; 4] = ['A', 'C', 'G', 'T'];

// D3 Palettes: https://github.com/d3/d3-3.x-api-reference/blob/master/Ordinal-Scales.md#categorical-colors
pub const PALETTE_DARK: [[u8; 4]; 9] = [
    [31, 119, 180, 255],  // dark blue
    [255, 127, 14, 255],  // dark orange
    [44, 160, 44, 255],   // dark green
    [214, 39, 39, 255],   // dark red
    [148, 103, 189, 255], // dark purple
    [140, 86, 59, 255],   // dark brown
    [227, 119, 195, 255], // dark pink
    //[127, 127, 127, 255], // dark grey
    [188, 189, 34, 255], // dark green alt
    [23, 189, 207, 255], // dark teal
];

pub const PALETTE_LIGHT: [[u8; 4]; 9] = [
    [174, 199, 232, 255], // light blue
    [255, 187, 120, 255], // light orange
    [152, 223, 138, 255], // light green
    [255, 152, 150, 255], // light red
    [197, 176, 213, 255], // light purple
    [196, 156, 148, 255], // light brown
    [247, 182, 210, 255], // light pink
    //[199, 199, 199, 255], // light grey
    [219, 219, 141, 255], // light green alt
    [158, 218, 229, 255], // light teal
];

pub const WHITE: Source = Source::Solid(SolidSource {
    r: 255,
    g: 255,
    b: 255,
    a: 255,
});
pub const BLACK: Source = Source::Solid(SolidSource {
    r: 0,
    g: 0,
    b: 0,
    a: 255,
});

pub const BLACK_ALT: image::Rgba<u8> = image::Rgba([0, 0, 0, 255]);

pub const GREY: Source = Source::Solid(SolidSource {
    r: 225,
    g: 225,
    b: 225,
    a: 255,
});

pub const LIGHT_GREY: Source = Source::Solid(SolidSource {
    r: 240,
    g: 240,
    b: 240,
    a: 255,
});

pub const DARK_GREY: Source = Source::Solid(SolidSource {
    r: 175,
    g: 175,
    b: 175,
    a: 255,
});

pub const TRANSPARENT: Source = Source::Solid(SolidSource {
    r: 0,
    g: 0,
    b: 0,
    a: 0,
});

pub const TEXT_COLOR: image::Rgba<u8> = image::Rgba([0, 0, 0, 255]);
pub const FONT_SIZE: f32 = 25.0;
// want at least 2, a single pixel causes blurry artifacts
pub const LINE_WIDTH: f32 = 2.0;

// needs to be just large enough to fit single char in FONT_SIZE
pub const X_INC: f32 = 50.0;
// small buffer for legibility
pub const BUFFER: f32 = 5.0;

pub const BASIC_STROKE_STYLE: StrokeStyle = StrokeStyle {
    cap: LineCap::Square,
    join: LineJoin::Miter,
    width: LINE_WIDTH,
    miter_limit: 2.,
    dash_array: Vec::new(),
    dash_offset: 16.,
};

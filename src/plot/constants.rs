use raqote::*;

pub const ALPHABET: [char; 4] = ['A', 'C', 'G', 'T'];

// D3 Palettes: https://github.com/d3/d3-3.x-api-reference/blob/master/Ordinal-Scales.md#categorical-colors
pub const PALETTE_DARK: [[u8; 4]; 3] = [
    [31, 119, 180, 255], // dark blue
    [255, 127, 14, 255], // dark orange
    [44, 160, 44, 255],  // dark green
];

pub const PALETTE_LIGHT: [[u8; 4]; 3] = [
    [174, 199, 232, 255], // light blue
    [255, 187, 120, 255], // light orange
    [152, 223, 138, 255], // light green
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
pub const GREY: Source = Source::Solid(SolidSource {
    r: 196,
    g: 196,
    b: 196,
    a: 255,
});
pub const TRANSPARENT: Source = Source::Solid(SolidSource {
    r: 0,
    g: 0,
    b: 0,
    a: 0,
});
pub const FONT_SIZE: i32 = 12;
// want at least 2, a single pixel causes blurry artifacts
pub const LINE_WIDTH: f32 = 2.;

// needs to be just large enough to fit single char in FONT_SIZE
pub const X_INC: i32 = 25;

// pseudo-contants
pub const BASIC_STROKE_STYLE: StrokeStyle = StrokeStyle {
    cap: LineCap::Square,
    join: LineJoin::Miter,
    width: LINE_WIDTH,
    miter_limit: 2.,
    dash_array: Vec::new(),
    dash_offset: 16.,
};
//static text_color: raqote::Color = raqote::Color::new(255, 0, 0, 0);

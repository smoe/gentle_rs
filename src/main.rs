// use std::io::stdout;

use lazy_static::lazy_static;
use dna_sequence::DNAsequence;
use render_dna::RenderSVG;

use crate::facility::Facility;

pub mod error;
pub mod facility;
// pub mod location;
pub mod restriction_enzyme;
pub mod protease;
pub mod enzymes;
pub mod dna_sequence;
pub mod amino_acids;
pub mod render_dna;

lazy_static! {
    pub static ref FACILITY: Facility = Facility::new();
}

struct CircularDnaImageData {
    image: slint::Image,
}


fn main() {
    let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
    let dna = dna.get(0).unwrap();
    let r = RenderSVG::from_dna_sequence(dna);
    // let _ = svg::write(stdout(),&r.document);
    let image = CircularDnaImageData {
        image: slint::Image::load_from_svg_data(r.document.to_string().as_bytes()).expect("Can't create image from SVG"),
    };
    println!("{:?}",image.image.size());
    let main_window = MainWindow::new().unwrap();
    main_window.set_render(image.image);
    main_window.run().unwrap();
}

slint::slint! {
    component CircularDnaRender inherits Rectangle {
        callback clicked;

        in property <image> render;

        width: 512px;
        height: 512px;
    
        Image {
            source: render;
            width: parent.width;
            height: parent.height;
        }

        TouchArea {
            clicked => {
                root.clicked();
            }
        }
    }

    export component MainWindow inherits Window {
        in property <image> render;

        CircularDnaRender {
            render: render;
            clicked => {
                
            }
        }
    }
}
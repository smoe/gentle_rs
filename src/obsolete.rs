// pub mod render_dna_circular_svg;
// pub mod render_dna_sequence_svg;

// fn main() {
//     let config = Config::new().with_window(
//         WindowBuilder::default()
//             .with_title("GENtle 2")
//             .with_maximized(true)
//             // .with_inner_size(dioxus_desktop::LogicalSize::new(4000.0, 3000.0)),
//     );

//     dioxus_desktop::launch_cfg(app, config);
//     // dioxus_web::launch(app);
// }

// fn app(cx: Scope) -> Element {
//     let x = use_state(cx, || 0 as usize);
//     let window = use_window(cx);
//     println!("{}x{}",window.inner_size().width,window.inner_size().height);

//     let dna = DNAsequence::from_genbank_file("test_files/pGEX-3X.gb").unwrap();
//     let dna = dna.get(0).unwrap();
//     let plasmid_map = RenderCircularMapSVG::from_dna_sequence(dna);
//     let pasmid_map_svg_str = plasmid_map.document.to_string();

//     let sequence_map = RenderDnaSequenceSVG::from_dna_sequence(dna);
//     let sequence_svg_str = sequence_map.document.to_string();

//     cx.render(rsx! {

//         div {
//             display:"flex",
//             div {
//                 width: "10rem",
//                 height: "100%",
//                 format!("LEFT: {}",window.inner_size().width)
//             }
//             div {
//                 width: "100%",
//                 height: "100%",
//                 div {
//                     display: "flex",
//                     div {
//                         width: "30%",
//                         height: "50%",
//                         "DETAILS"
//                     }
//                     div {
//                         width: "70%",
//                         height: "300px",
//                         dangerous_inner_html: "{pasmid_map_svg_str}",
//                     }
//                 }
//                 div {
//                     width: "100%",
//                     height: "50%",
//                     position: "relative",
//                     overflow: "auto",
//                     dangerous_inner_html: "{sequence_svg_str}",
//                     onclick: move |event: Event<MouseData>| {
//                         let pos = event.data.coordinates().element();
//                         x.set((pos.x) as usize);
//                         let y = (pos.y) as usize;
//                         let n = sequence_map.get_click_location(*x.get(), y);
//                         println!("Clicked! Event: {event:#?}\n{pos:?}\n{n:?}");
//                     },
//                 }
//             }
//         }
//     })
// }

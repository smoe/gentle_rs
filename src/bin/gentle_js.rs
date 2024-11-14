use gentle::js_interface::JavaScriptInterface;
use rustyline::error::ReadlineError;
use rustyline::DefaultEditor;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // Create the JavaScript runtime
    let mut js = JavaScriptInterface::default();

    // Create line editor
    let mut rl = DefaultEditor::new()?;

    loop {
        // Read line from user
        let readline = rl.readline("GENtle> ");

        match readline {
            Ok(line) => {
                if line.trim() == "exit" {
                    break;
                }

                rl.add_history_entry(line.as_str())?;

                js.run(line);

                // let result = runtime.execute_script("<usage>", line);
                // match result {
                //     Ok(_state) => {
                //         // println!("!{:?}", &value);
                //     }
                //     Err(e) => println!("Error: {}", e),
                // }
            }
            Err(ReadlineError::Interrupted) => {
                println!("CTRL-C");
                break;
            }
            Err(ReadlineError::Eof) => {
                println!("CTRL-D");
                break;
            }
            Err(err) => {
                println!("Error: {}", err);
                // break;
            }
        }
    }

    Ok(())
}

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

pub trait LineProcessor {
    fn process_line(&mut self, line: String);
}

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>+Clone,
{
    let file = match File::open(filename.clone()) {
        Err(why) => {
            println!(
                "Looking for {:?} in {:?}",
                filename.as_ref().to_str(),
                std::env::current_dir().unwrap().display()
            );
            panic!("Cannot open file: {}", why)
        }
        Ok(file) => file,
    };
    Ok(io::BufReader::new(file).lines())
}

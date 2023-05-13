use std::io::{BufRead, BufReader};

use sscanf::scanf;
use winter_math::StarkField;

use crate::errors::*;
use crate::io::LineProcessor;

#[derive(Debug)]

/// Parses .in or .wires file to produce an input value vector ("z").
pub struct JsnarkWireReaderParser<E: StarkField> {
    pub wires: Vec<E>,
}

/// Parses .in or .wires file to produce an input value vector ("z").
pub struct JsnarkWireParser<'a, E: StarkField> {
    pub verbose: bool,
    pub wires: &'a mut Vec<E>,
}

impl<'a, E: StarkField> JsnarkWireParser<'a, E> {
    pub fn new(wires: &'a mut Vec<E>) -> Result<Self, InputWireError> {
        Ok(JsnarkWireParser {
            verbose: false,
            wires: wires,
        })
    }

    fn handle_assign(&mut self, wire_id: usize, wire_val: E) {
        if self.verbose {
            println!("ASSIGN: {} {:?}", wire_id, wire_val);
        }

        if wire_id >= self.wires.len() {
            self.wires.resize(wire_id + 1, E::ZERO);
        }
        self.wires[wire_id] = wire_val;
    }
}

impl<'a, E: StarkField> LineProcessor for JsnarkWireParser<'a, E> {
    #[cfg_attr(feature = "flame_it", flame)]
    fn process_line(&mut self, line: String) {
        if self.verbose {
            println!("{}", line);
        }
        if line.starts_with("#") {
            return;
        }

        // Remove comments and trim end-whitespace.
        let mut parts = line.split("#");
        let mut buf = parts.next().unwrap();
        buf = buf.trim();

        match scanf!(buf, "{} {x}", usize, u128) {
            Some((wire_id, wire_value)) => {
                self.handle_assign(wire_id, E::from(wire_value));
                return;
            }
            None => {}
        }

        println!("FAILED WIRE: {}", line);
    }
}

impl<'a, E: StarkField> JsnarkWireReaderParser<E> {
    pub fn new() -> Result<Self, InputWireError> {
        Ok(JsnarkWireReaderParser {
            wires: Vec::<E>::new(),
        })
    }

    fn pad_power_two(&mut self) {
        let num_wires = self.wires.len();
        if !num_wires.is_power_of_two() {
            let padding = num_wires.next_power_of_two() - num_wires;
            for _ in 0..padding {
                self.wires.push(E::ZERO);
            }
        }
    }

    pub fn parse_wire_file(&mut self, wire_file: &str, verbose: bool) {
        if verbose {
            println!("Parse wire file {}", wire_file);
        }

        let mut wire_parser = JsnarkWireParser::<E>::new(&mut self.wires).unwrap();
        wire_parser.verbose = verbose;

        // if let Ok(lines) = crate::io::read_lines(wire_file) {
        //     for line in lines {
        //         match line {
        //             Ok(ip) => {
        //                 wire_parser.process_line(ip);
        //             }
        //             Err(e) => println!("{:?}", e),
        //         }
        //     }
        // }
        // match crate::io::read_lines(wire_file) {
        //     Ok(lines) => {
        //         for line in lines {
        //             match line {
        //                 Ok(ip) => {
        //                     wire_parser.process_line(ip);
        //                 }
        //                 Err(e) => println!("{:?}", e),
        //             }
        //         }
        //     }
        //     Err(e) => println!("{:?}", e),
        // }

        match crate::io::open_file(wire_file) {
            Ok(file) => {
                let reader = BufReader::new(file);
                for line in reader.lines() {
                    match line {
                        Ok(ip) => {
                            wire_parser.process_line(ip);
                        }
                        Err(e) => println!("{:?}", e),
                    }
                }
            }
            Err(e) => println!("{:?}", e),
        }

        self.pad_power_two();

        // if verbose {
        //     println!("{:?}", self.wires);
        // }
    }
}

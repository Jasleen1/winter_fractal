#![allow(dead_code,unused_imports)]

pub mod accumulator;
pub mod errors;

#[cfg(feature = "flame_it")]
extern crate flame;
#[cfg(feature = "flame_it")]
#[macro_use]
extern crate flamer;

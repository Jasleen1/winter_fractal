
#[cfg(feature = "flame_it")]
extern crate flame;
#[cfg(feature = "flame_it")]
#[macro_use]
extern crate flamer;

pub mod r1cs_orchestrator;

#[cfg(test)]
mod tests;
#![allow(dead_code,unused_imports)]

use super::flame_local::{merge_repeated_spans, dump_html_custom};

use std::io::Write;

#[cfg_attr(feature = "flame_it", flame)]
pub fn generate_flame_report(report_dir_opt: Option<&str>, filename_prefix: &str, focus_method: Option<&str>) {
    let report_dir = report_dir_opt.unwrap_or("gen/reports");
    std::fs::create_dir_all(report_dir).unwrap();

    {
        // Raw report.
        let filespec = format!("{report_dir}/{filename_prefix}-flame.html");
        let file = &mut std::fs::File::create(filespec).unwrap();
        flame::dump_html(file).unwrap_or_else(|why| { println!("! {:?}", why.kind()); });
    }

    if focus_method.is_some() {
        // Report with repeated spans merged.
        let method_filter = focus_method.unwrap();
        let html_filespec = format!("{report_dir}/{filename_prefix}-flame-merge.html");
        let html_file = &mut std::fs::File::create(html_filespec).unwrap();

        let mut spans = flame::threads().into_iter().next().unwrap().spans;
        let method_metric = gather_time(method_filter.to_string(), &spans);

        merge_repeated_spans(&mut spans);

        // HTML report of metrics.
        dump_html_custom(html_file, &spans).unwrap_or_else(|why| { println!("! {:?}", why.kind()); });

        // Text report of metrics.
        let text_filespec = format!("{report_dir}/{filename_prefix}-flame-merge.txt");
        let text_file = std::fs::File::create(text_filespec).unwrap();
        let text_writer = &mut std::io::BufWriter::new(&text_file);
        writeln!(text_writer, "timemetric: {:?}", method_metric).unwrap();
        writeln!(text_writer, "timepercent: {}", method_metric.cumul_ns as f64 * 100.0 / method_metric.elapsed_ns as f64).unwrap();
    }
}

// Calculate cumulative time in particular spans
pub fn gather_time(tag: String, spans: &Vec<flame::Span>) -> TimeMetric {
    let mut gatherer = TimeMetric::new(tag);
    gatherer.elapsed_ns = spans.iter().map(|sp| sp.delta).sum();
    gatherer.visit_vec(spans);
    gatherer
}

// Debugging
fn list_spans(spans: &mut Vec<flame::Span>) {
    spans.iter().for_each(|s| println!("{:?}", s));
}

// Visitor pattern for collecting statistics.
trait Visitor<T> {
    fn visit(&mut self, visitable: &T);

    fn visit_vec(&mut self, many_visitable: &Vec<T>) {
        for visitable in many_visitable.iter() {
            self.visit(visitable);
        }
    }

    fn visit_arr(&mut self, many_visitable: &[T]) {
        for visitable in many_visitable.iter() {
            self.visit(visitable);
        }
    }
}

// TimeMetric collects cumulative span delta time for matching spans.
#[derive(Debug, Clone)]
pub struct TimeMetric {
    tag: String,
    cumul_ns: u64,
    elapsed_ns: u64,
}

impl TimeMetric {
    fn new(tag: String) -> TimeMetric {
        TimeMetric { tag: tag, cumul_ns: 0, elapsed_ns: 0 }
    }
}

// Collect the aggregate delta time in the matching spans.
// Because recursion can happen, we only count the topmost
// call to a tagged method.
impl Visitor<flame::Span> for TimeMetric {
    fn visit(&mut self, visitable: &flame::Span) {
        if self.tag == visitable.name {
            self.cumul_ns += visitable.delta;
        } else {
            // Remark: avoid over-counting descendants who also match,
            // hence we don't blithely always descend.
            self.visit_vec(&visitable.children);
        }

        // Avoid doing this. It double-counts. Really want just the top level elapsed.
        //self.elapsed += visitable.delta;
    }
}
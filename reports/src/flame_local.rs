
// Wrappers and local forks of private flame-crate methods.

// Merge same-tagged children.
// Lifted from suggested approach in https://github.com/llogiq/flame/issues/33 .
pub fn merge_repeated_spans(spans: &mut Vec<flame::Span>) {
    if spans.is_empty() {
        return;
    }

    // Sort so spans to be merged are adjacent and spans with the most children are
    // merged into to minimise allocations.
    spans.sort_unstable_by(|s1, s2| {
        let a = (&s1.name, s1.depth, usize::max_value() - s1.children.len());
        let b = (&s2.name, s2.depth, usize::max_value() - s2.children.len());
        a.cmp(&b)
    });

    // Copy children and sum delta from spans to be merged
    let mut merge_targets = vec![0];
    {
        let mut spans_iter = spans.iter_mut().enumerate();
        let (_, mut current) = spans_iter.next().unwrap();
        // println!("* {},{}", current.name, current.depth);
        for (i, span) in spans_iter {
            if current.name == span.name && current.depth == span.depth {
                current.delta += span.delta;
                let mut children = std::mem::replace(&mut span.children, Vec::new());
                current.children.extend(children.into_iter());
            } else {
                current = span;
                // println!("+ {},{}", current.name, current.depth);
                merge_targets.push(i);
            }
        }
    }

    // Move merged spans to the front of the spans vector
    for (target_i, &current_i) in merge_targets.iter().enumerate() {
        spans.swap(target_i, current_i);
    }

    // Remove duplicate spans
    spans.truncate(merge_targets.len());

    // Merge children of the newly collapsed spans
    for span in spans {
        merge_repeated_spans(&mut span.children);
    }
}

// Lifted directly from flame-0.2.2 source, because it is private.
// flame-0.2.1-pre opens it up but flamer requires precise 0.2.2. No good.
// Ty Overby, flame maintainer, claimed in 2017 to intend to help with folding
// by deprecating the dump_html and using a serde approach but nothing arrived.
// So expect that if anything, dump_html_custom is deprecated and won't be
// exposed.
// -- your faithful scribe Don

use std::io::Write;
use std::io::Result as IoResult;
use flame::{Span};

pub fn dump_html_custom<W: Write>(mut out: W, spans: &[Span]) -> IoResult<()> {
    fn dump_spans<W: Write>(out: &mut W, span: &Span) -> IoResult<()> {
        writeln!(out, "{{")?;
        writeln!(out, r#"name: "{}","#, span.name)?;
        writeln!(out, "value: {},", span.delta)?;
        writeln!(out, "start: {},", span.start_ns)?;
        writeln!(out, "end: {},", span.end_ns)?;
        writeln!(out, "children: [")?;
        for child in &span.children {
            dump_spans(out, child)?;
            writeln!(out, ",")?;
        }
        writeln!(out, "],")?;
        writeln!(out, "}}")?;
        Ok(())
    }

    write!(out, r#"
<!doctype html>
<html>
    <head>
        <style>
            html, body {{
                width: 100%;
                height: 100%;
                margin: 0;
                padding: 0;
            }}
            {}
        </style>
        <script>
            {}
            {}
            {}
        </script>
    </head>
    <body>
        <script>
            var width = document.body.offsetWidth;
            var height = document.body.offsetHeight - 100;
            var flamegraph =
                d3.flameGraph()
                  .width(width)
                  .height(height)
                  .tooltip(false)
                  .sort(function(a, b){{
                    if (a.start < b.start) {{
                        return -1;
                    }} else if (a.start > b.start) {{
                        return 1;
                    }} else {{
                        return 0;
                    }}
                  }});
            d3.select("body").datum({{ children: [
"#, include_str!("../resources/flameGraph.css"), include_str!("../resources/d3.js"), include_str!("../resources/d3-tip.js"), include_str!("../resources/flameGraph.js"))?;

    for span in spans {
        dump_spans(&mut out, &span)?;
        writeln!(out, ",")?;
    }

    write!(out, r#"]}}).call(flamegraph);
         </script>
    </body>
</html>"#)?;

    Ok(())
}
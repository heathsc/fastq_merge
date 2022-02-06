use std::{io, str::FromStr};

use clap::{App, Arg};
use stderrlog::Timestamp;

use crate::config::Config;

struct LogLevel {
    level: usize,
}

impl FromStr for LogLevel {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "error" => Ok(LogLevel { level: 0 }),
            "warn" => Ok(LogLevel { level: 1 }),
            "info" => Ok(LogLevel { level: 2 }),
            "debug" => Ok(LogLevel { level: 3 }),
            "trace" => Ok(LogLevel { level: 4 }),
            "none" => Ok(LogLevel { level: 5 }),
            _ => Err("no match"),
        }
    }
}

impl LogLevel {
    fn get_level(&self) -> usize {
        if self.level > 4 {
            0
        } else {
            self.level
        }
    }
}

pub fn handle_cli() -> io::Result<Config> {
    let m = App::new("fastq_merge").version(crate_version!()).author("Simon Heath").about("Merge multiple FASTQ files (and optionally summary sequence files), removing duplicate reads by only keeping the first record found for each read ID")
	.arg(Arg::new("report_file").short('r').long("report-file").takes_value(true).value_name("FILE").help("Name for report file. Default stderr"))
	.arg(Arg::new("loglevel").short('l').long("loglevel").takes_value(true).value_name("LOGLEVEL").help("Set log level")
			.possible_values(&["none", "error", "warn", "info", "debug", "trace"]).ignore_case(true).default_value("info"))
	.arg(Arg::new("summary_files").short('s').long("summary-files").takes_value(true).multiple_occurrences(true).require_delimiter(true).value_name("FILE").help("Name of summary sequencing file(s)"))
	.arg(Arg::new("summary_name").short('S').long("summary-name").takes_value(true).value_name("FILE").help("Name of merged summary sequencing file"))
	.arg(Arg::new("input").takes_value(true).multiple_occurrences(true).value_name("FILE").help("Input FASTQ file(s) or directories. Default stdin"))
	.get_matches();

    let verbose = m.value_of_t::<LogLevel>("loglevel").unwrap();
    stderrlog::new()
        .verbosity(verbose.get_level())
        .timestamp(Timestamp::Second)
        .init()
        .unwrap();

    debug!("Handling command line inputs");
    let inputs = match m.values_of("input") {
        Some(v) => v.map(|s| s.to_owned()).collect(),
        None => Vec::new(),
    };

    let summary_files = match m.values_of("summary_files") {
        Some(v) => v.map(|s| s.to_owned()).collect(),
        None => Vec::new(),
    };

    let summary_name = m
        .value_of("summary_name")
        .unwrap_or("sequencing_summary.txt")
        .to_owned();

    let report_file = m.value_of("report_file").map(|s| s.to_owned());
    Config::new(inputs, summary_files, summary_name, report_file)
}

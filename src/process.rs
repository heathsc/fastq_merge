use chrono::Utc;
use std::{
    collections::HashMap,
    fmt,
    fs::File,
    io::{self, stderr, stdout, BufRead, BufReader, BufWriter, Error, ErrorKind, Read, Write},
};

use crate::config::{Config, InputFile, SummaryFile};
use crate::stats::{FileStats, FileStatsType};

/// An input file that can be either a fastq or a sequencing summary file.
struct MFile {
    reader: BufReader<Box<dyn Read>>,
    line: usize, // Current input line
    buf: String, // Input buffer
}

impl MFile {
    /// Reads next line (up to linefeed) from MFile into buf
    /// Return io::Err if an error was produced on reading, Ok<None> at EOF and
    /// Ok<Some<()>> if line succesfully read
    fn next_line(&mut self) -> io::Result<Option<()>> {
        self.buf.clear();
        self.reader.read_line(&mut self.buf).map(|l| {
            self.line += 1;
            if l == 0 {
                None
            } else {
                Some(())
            }
        })
    }

    /// Reads next line from MFile into buf and returns a vector of tab separated columns
    /// Return io::Err if an error was produced on reading, Ok<None> at EOF and
    /// Ok<Some<Vec<&str>>> on success
    fn next_parsed_line(&mut self) -> io::Result<Option<Vec<&str>>> {
        Ok(self
            .next_line()?
            .map(move |_| self.buf.trim_end().split('\t').collect()))
    }
}

/// Input fastq file
struct InFile<'a> {
    file: &'a InputFile,
    mfile: MFile,
}

impl<'a> InFile<'a> {
    /// Create new InFile from InputFile
    /// Returns error if file can not be opened
    fn new(file: &'a InputFile) -> io::Result<Self> {
        Ok(Self {
            file,
            mfile: MFile {
                reader: file.get_bufreader()?,
                line: 0,
                buf: String::new(),
            },
        })
    }

    /// Generate file truncted error
    fn trunc_err(&self) -> io::Error {
        Error::new(
            ErrorKind::Other,
            format!(
                "{}:{} Unexpected EOF - file truncated?",
                &self.file, self.mfile.line
            ),
        )
    }

    /// Generate unexpected character error
    fn unexp_char(&self, c: char) -> io::Error {
        Error::new(
            ErrorKind::Other,
            format!(
                "{}:{} Unexpected char at start of line (expecting '{}'",
                &self.file, self.mfile.line, c
            ),
        )
    }

    /// Read and process next fastq record (4 lines consisting of Read ID, Sequence, Spare, Quality)
    /// If the read ID is not already present in read_id_seen then record will output to w and the entry in read_id_seen set accordingly
    /// The number of reads, bases and the base quality distribution are stored in fs for output and skipped records separately
    fn next_fastq_record<W: Write, W1: Write>(
        &mut self,
        read_id_seen: &mut HashMap<Box<str>, ReadIdState>,
        fs: &mut FileStats,
        summ_avail: bool,
        w: &mut W,
        report: &mut W1,
    ) -> io::Result<Option<()>> {
        // Read ID line
        // uniq will be used for the remainder of the record to indicate whether this read shuold be output or not
        // An EOF at this stage is not an error but indicates the end of the file with no errors found
        let state = match self.mfile.next_line()? {
            Some(_) => {
                // Check that line starts with '@'
                match self.mfile.buf.strip_prefix('@') {
                    // The Read ID goes from the '@' to the first whitespace character exclusive
                    Some(s) => match s.split_ascii_whitespace().next() {
                        Some(id) => {
                            if let Some(x) = read_id_seen.get_mut(id) {
                                if x.seen(DUPLEX_READ) {
                                    FileStatsType::Duplex
                                } else if x.seen(SEEN_READ) {
                                    FileStatsType::NonUnique
                                } else {
                                    if summ_avail && !x.seen(SEEN_SUMMARY) {
                                        writeln!(report, "Warning: Read {} has a FASTQ record but is not present in the summary file", id)?;
                                    }
                                    x.set(SEEN_READ);
                                    FileStatsType::Unique
                                }
                            } else {
                                // Read ID entry not found so set read_id_state to show that we have seen this ID
                                read_id_seen.insert(Box::from(id), ReadIdState::new(SEEN_READ));
                                if summ_avail {
                                    writeln!(report, "Warning: Read {} has a FASTQ record but is not present in the summary file", id)?;
                                }
                                FileStatsType::Unique
                            }
                        }
                        // Read ID line only has the '@ line with no ID following it
                        None => {
                            return Err(Error::new(
                                ErrorKind::Other,
                                format!(
                                    "{}:{} Unexpected EOL - short line",
                                    &self.file, self.mfile.line
                                ),
                            ))
                        }
                    },
                    None => return Err(self.unexp_char('@')), // Read ID line does not start with '@'
                }
            }
            None => return Ok(None), // EOF
        };
        let uniq = state == FileStatsType::Unique;
        if uniq {
            write!(w, "{}", self.mfile.buf)?;
        }

        // An EOF for any of the subsequent stages is an error (truncated record)

        // Sequence line (we just check the length to compare against the quality line)
        let n_bases = match self.mfile.next_line()? {
            Some(_) => {
                if uniq {
                    write!(w, "{}", self.mfile.buf)?;
                }
                self.mfile.buf.trim_end().len()
            }
            None => return Err(self.trunc_err()),
        };

        // + line (just check that it is present)
        match self.mfile.next_line()? {
            Some(_) => {
                if !self.mfile.buf.starts_with('+') {
                    return Err(self.unexp_char('+'));
                }
                if uniq {
                    write!(w, "{}", self.mfile.buf)?
                }
            }
            None => return Err(self.trunc_err()),
        }

        // quality line
        match self.mfile.next_line()? {
            Some(_) => {
                let q = self.mfile.buf.trim_end();
                if q.len() != n_bases {
                    Err(Error::new(
                        ErrorKind::Other,
                        format!(
                            "{}:{} Quality and Sequence lines have different lengths",
                            self.file, self.mfile.line
                        ),
                    ))
                } else {
                    // Update stats
                    fs.add_record(q, state);
                    if uniq {
                        write!(w, "{}", self.mfile.buf)?;
                    }
                    Ok(Some(()))
                }
            }
            None => Err(self.trunc_err()),
        }
    }
}

impl<'a> fmt::Display for InFile<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.file)
    }
}

/// Look for the column index corresponding to a particular header column.  Returns error if not found
fn chk_col(
    s: &str,
    columns: &HashMap<Box<str>, usize>,
    col_ix: &[Option<usize>],
    file: &SummaryFile,
) -> io::Result<usize> {
    let i = columns
        .get(s)
        .expect("required column not found in headers");
    match col_ix[*i].as_ref() {
        Some(x) => Ok(*x),
        None => Err(Error::new(
            ErrorKind::Other,
            format!(
                "Column '{}' not found in sequencing summary file {}",
                s, file
            ),
        )),
    }
}

/// Input sequence_summary file
struct SumFile<'a> {
    file: &'a SummaryFile,
    mfile: MFile,
    col_ix: Vec<Option<usize>>, // col_ix[i] gives the Some(j) where j is column index in the input for output column i, or None if this column does not exist in this input file
    n_cols: usize,              // number of columns in input
    read_id_ix: usize,          // Input column index of read_id column
    filter_ix: usize,           // Input column index of passes_filter column
    duplex_pair_ix: Option<usize>, // Input column of duplex_pair_read_id column (may not exist)
}

impl<'a> SumFile<'a> {
    /// Create new SumFile
    /// Generate output order using the HashMap columns as a template.  Afterwards we can iterate over
    /// col_ix to output the columns in the correct order.  If col_ix[i] is None then that columns does
    /// not exist in this datafile and so should be replaced by a dash ('-')
    fn new(file: &'a SummaryFile, columns: &HashMap<Box<str>, usize>) -> io::Result<Self> {
        let ix: Vec<_> = file
            .header_line()
            .trim_end()
            .split('\t')
            .map(|s| *columns.get(s).expect("Column not found"))
            .collect();
        let n_cols = ix.len();
        let mut col_ix = vec![None; columns.len()];
        for (i, j) in ix.into_iter().enumerate() {
            col_ix[j] = Some(i);
        }

        // Open file for reading
        let mut mfile = MFile {
            reader: file.get_reader()?,
            line: 0,
            buf: String::new(),
        };
        // Discard header line (we already have it in 'SummaryFile')
        let _ = mfile.next_line()?;
        let read_id_ix = chk_col("read_id", columns, &col_ix, file)?;
        let filter_ix = chk_col("passes_filtering", columns, &col_ix, file)?;
        let duplex_pair_ix = chk_col("duplex_pair_read_id", columns, &col_ix, file).ok();
        Ok(Self {
            file,
            mfile,
            col_ix,
            n_cols,
            read_id_ix,
            filter_ix,
            duplex_pair_ix,
        })
    }

    /// Read next line from input, checking that we have the expected number of columns.
    /// Retuens Ok(None) at EOF, OK(Some(_)) on a successful read and Err() otherwise.
    fn next_parsed_line(&mut self) -> io::Result<Option<Vec<&str>>> {
        let line = self.mfile.line + 1;
        let res = self.mfile.next_parsed_line();
        match res {
            Ok(Some(v)) if v.len() != self.n_cols => Err(Error::new(
                ErrorKind::Other,
                format!(
                    "{}: {} Unexpected number of columns found in sequencing summary file",
                    self.file, line
                ),
            )),
            _ => res,
        }
    }
}
impl<'a> fmt::Display for SumFile<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.file)
    }
}

/// Keeps track of whether we have seen the fastq record for a read ID and the summary record for a read ID
#[derive(Default, Copy, Clone)]
struct ReadIdState(u8);

const SEEN_READ: u8 = 1;
const SEEN_SUMMARY: u8 = 2;
const DUPLEX_READ: u8 = 4;
const PASSED: u8 = 8;
const SEEN_ALL: u8 = 15;

impl ReadIdState {
    fn new(mk: u8) -> Self {
        Self(mk & (SEEN_ALL))
    }

    fn seen(&self, mask: u8) -> bool {
        (self.0 & mask) != 0
    }

    fn check(&self, mask1: u8, mask2: u8) -> bool {
        (self.0 & mask1) == mask2
    }

    fn set(&mut self, mask: u8) -> &mut Self {
        self.0 |= mask & SEEN_ALL;
        self
    }
}

/// Initialize output of merged sequencing_summary output
/// As different files may not have the same columns (and maybe not in the same order), we first form a union of all
/// column headings that we find in all files.  Each new columns heading is given an index number that represents the output
/// column order.  These are assigned in the order that the headings are discovered, therefore the columns present in the
/// first input file will be in the same order in the output, and any subsequent columns are added at the end.
/// After this we check that the required headers 'read_id' and 'passes_filtering' are present and, if so, we write the column headers
/// in the final output order.
///
/// Returns Err() if required columns are not found or on write failure.
///
fn init_sequencing_output<W: Write>(
    sfiles: &[SummaryFile],
    wrt: &mut W,
) -> io::Result<HashMap<Box<str>, usize>> {
    // HashMap linking colulmn names to output column index
    let mut columns = HashMap::new();

    // Add column names in order to HashMap, using the current number of items
    // in columns as the hash value.
    for file in sfiles.iter() {
        for fd in file.header_line().trim_end().split('\t') {
            if !columns.contains_key(fd) {
                columns.insert(Box::from(fd), columns.len());
            }
        }
    }

    // Check for the required columns (hopefully the names are consistent across different Guppy versions!)
    if !columns.contains_key("read_id") {
        Err(Error::new(
            ErrorKind::Other,
            "Column 'read_id' not found in sequencing summary files".to_string(),
        ))
    } else if !columns.contains_key("passes_filtering") {
        Err(Error::new(
            ErrorKind::Other,
            "Column 'passes_filtering' not found in sequencing summary files".to_string(),
        ))
    } else {
        // HashMap does not store entries in insertion order, so we make a Vec with the columns heading in the required order
        // We use Option<_> as the vector element so that we can initialize v with None and then add the correct entries
        // as we iterate through the map.  At the end all of the elements should be Some(_), but if something catastrophic happens
        // we will detect it and panic! in the next section
        let mut v = vec![None; columns.len()];
        for (s, ix) in columns.iter() {
            v[*ix] = Some(s)
        }

        // Output the column headings in order as tab separated text
        // We use unwrap() on the element, because if all elements in v are not Some(&str) then the previous lines have a bug and we should panic!

        // First column
        write!(wrt, "{}", v[0].unwrap())?;

        // Subsequent columns
        for s in v[1..].iter().map(|x| x.unwrap()) {
            write!(wrt, "\t{}", s)?;
        }
        writeln!(wrt)?;

        // Return HashMap with column headings and indices
        Ok(columns)
    }
}

/// Ouput line from sequence_summary file in correct order, substituting a '-' if a column is missing
fn output_seq_sum_line<W: Write>(
    fields: &[&str],
    col_ix: &[Option<usize>],
    wrt: &mut W,
) -> io::Result<()> {
    // Map column index to appropriate column in input file or '-' if missing
    let map_fn = |x: &Option<usize>| match x {
        Some(ix) => fields[*ix],
        None => "-",
    };

    // First column
    write!(wrt, "{}", map_fn(&col_ix[0]))?;

    // Subsequent columns
    for fd in col_ix[1..].iter() {
        write!(wrt, "\t{}", map_fn(fd))?;
    }

    writeln!(wrt)
}

/// Process all input FASTQ files in turn, sending to stdout records that have not already been seem (unique)
/// and discarding those that have already been seen.
///
/// If present, we all merge the input summary sequence files in the same way, discarding entries wih the same
/// read ID as a previously output entry.  The merging of the summary files is slightly complicated by the fact that
/// the number (and potentially the order) of columns in the input files is not necessarily the same.  We therefore
/// take the union of the columns found in all of the input files, and use a common column ordering for the output.  Any
/// missing columns are indicated by a '-' character in the output.
///
/// The 'HashMap' read_id_seen is used to keep a track of the read ids that have already been processed.
pub fn process(cfg: &Config) -> io::Result<()> {
    // Get starting time
    let now = Utc::now();

    // Open output report file
    let mut report: Box<dyn Write> = if let Some(name) = cfg.report_file() {
        Box::new(File::create(name)?)
    } else {
        Box::new(stderr())
    };

    // We use read_id_seen to store the read ids that we have already processed
    // so we can skip reads with duplicate ids
    let mut read_id_seen: HashMap<Box<str>, ReadIdState> = HashMap::new();

    // Keep track of the total counts over all files
    let mut total_stats = FileStats::new();

    writeln!(report, "{}: fastq_merge starting", now)?;

    let summary_avail = !cfg.summary_files().is_empty();

    // Merge sequencing summary files if present
    if summary_avail {
        // Open sequencing summary output file
        let mut seq_sum_output = BufWriter::new(File::create(cfg.summary_file_name())?);

        // Setup column list and print output header line
        let columns = init_sequencing_output(cfg.summary_files(), &mut seq_sum_output)?;

        // Iterate over input files
        for file in cfg.summary_files().iter() {
            // Open file, check that required columns exist and generate output column ordering
            let mut sfile = SumFile::new(file, &columns)?;

            info!("Reading from sequencing summary file {}", sfile);

            let read_id_ix = sfile.read_id_ix;
            let filter_ix = sfile.filter_ix;
            let duplex_ix = sfile.duplex_pair_ix;
            let col_ix = sfile.col_ix.clone();

            // Iterate over input lines
            while let Some(fields) = sfile.next_parsed_line()? {
                let s = fields[read_id_ix];
                match read_id_seen.get(s) {
                    Some(x) => {
                        if !x.seen(SEEN_SUMMARY | DUPLEX_READ) {
                            output_seq_sum_line(&fields, &col_ix, &mut seq_sum_output)?;
                            read_id_seen.get_mut(s).unwrap().set(SEEN_SUMMARY);
                        }
                    }
                    None => {
                        let mut state = ReadIdState::new(SEEN_SUMMARY);
                        if fields[filter_ix] == "TRUE" {
                            state.set(PASSED);
                        }
                        read_id_seen.insert(Box::from(s), state);
                        output_seq_sum_line(&fields, &col_ix, &mut seq_sum_output)?;
                    }
                }
                // Check for duplex_pair_read_id field
                if let Some(ix) = duplex_ix {
                    let s = fields[ix];
                    match read_id_seen.get_mut(s) {
                        Some(x) => {
                            if !x.seen(DUPLEX_READ) {
                                writeln!(report, "Error: Read {} present as both as individual read and as a duplex pair", s)?;
                                x.set(DUPLEX_READ);
                            }
                        }
                        None => {
                            let mut state = ReadIdState::new(DUPLEX_READ);
                            if fields[filter_ix] == "TRUE" {
                                state.set(PASSED);
                            }
                            read_id_seen.insert(Box::from(s), state);
                        }
                    }
                }
            }
        }
    }

    // Open main output to stdout
    let mut output = BufWriter::new(stdout());

    // Iterate over inputs
    for file in cfg.inputs().iter().map(InFile::new) {
        let mut file = file?;
        info!("Reading from fastq file {}", file);

        // Counts for this file
        let mut stats = FileStats::new();

        // Iterate over inputs
        while file
            .next_fastq_record(
                &mut read_id_seen,
                &mut stats,
                summary_avail,
                &mut output,
                &mut report,
            )?
            .is_some()
        {}

        writeln!(report, "{}: {}", file, stats)?;

        // Accumulate total stats
        total_stats += stats;
    }

    if summary_avail {
        // Generate warnings if any output FASTQ records do not have a corresponding entry in the sequencing summary fileprint
        for (read_id, _) in read_id_seen
            .iter()
            .filter(|(_, state)| state.check(SEEN_READ | SEEN_SUMMARY | PASSED, SEEN_READ | PASSED))
        {
            writeln!(
                report,
                "Warning: Read {} has a FASTQ record but is not present in the summary file",
                read_id
            )?;
        }
    }

    // Output overall read statistics
    writeln!(report, "\nTotal Stats: {}", total_stats)
}
